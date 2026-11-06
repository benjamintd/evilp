#define _POSIX_C_SOURCE 200809L
#include "augment.h"
#include "exact.h"
#include "lp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void require(int yes, const char *message) {
    if (!yes) {
        fprintf(stderr, "augmentation contract: %s\n", message);
        exit(1);
    }
}
static LabAugmentLimits limits(void) {
    return (LabAugmentLimits){1024, 100000, 4};
}
static void read_model(LabModel *m) {
    char path[] = "/tmp/milplab-augmentation-XXXXXX";
    int fd = mkstemp(path);
    require(fd >= 0, "temporary MPS created");
    FILE *f = fdopen(fd, "w");
    require(f != NULL, "MPS stream created");
    fputs("NAME source\nROWS\n N obj\n E eq\n L cap\n G min\n"
          "COLUMNS\n x obj 0.1 eq 0.1\n x eq 0.2 cap 1\n"
          " y eq 0.3 min 1\n z cap 1\n"
          "RHS\n rhs eq 0.6 cap 2\n rhs min 0\n"
          "BOUNDS\n LO b x 0.1\n UP b x 2\n UP b y 2\n UP b z 2\nENDATA\n", f);
    require(!fclose(f), "MPS stream closed");
    require(!lab_model_read(m, path), "source model parsed");
    require(!unlink(path), "temporary MPS removed");
    require(m->exact && m->envelope && m->envelope->complete,
            "exact decimal and enclosure provenance retained");
}
static void dyadic_equal(const mpq_t q, double value, const char *message) {
    mpq_t binary;
    mpq_init(binary);
    mpq_set_d(binary, value);
    require(!mpq_cmp(q, binary), message);
    mpq_clear(binary);
}
static void source_preservation_and_lifetime(void) {
    LabModel source = {0}, owned = {0};
    read_model(&source);
    LabAugmentLimits cap = limits();
    double coefficient[] = {nextafter(.1, INFINITY), 0, 1.0 / 3.0};
    LabAddedRow row = {"cut_0", coefficient, nextafter(.2, -INFINITY), INFINITY};
    /* Reproduce presolve's borrowed row mask and bound overlays. */
    double lower[] = {.5, 0, 0}, upper[] = {2, 2, 2};
    double row_lower[] = {.6, -INFINITY, 0}, row_upper[] = {.6, INFINITY, INFINITY};
    LabModel view = source;
    view.col_lower = lower;
    view.col_upper = upper;
    view.row_lower = row_lower;
    view.row_upper = row_upper;
    require(lab_model_append_rows(&view, 1, &row, &cap, &owned) == LAB_AUGMENT_OK,
            "append to presolved borrowed view");
    require(owned.rows == 4 && owned.cols == 3 && owned.nnz == source.nnz + 2,
            "appended CSC dimensions");
    require(source.rows == 3 && source.col_lower[0] == .1 && source.row_upper[1] == 2,
            "original model immutable");
    require(owned.col_lower[0] == .5 && isinf(owned.row_upper[1]),
            "active LP bound overlays and row mask preserved");
    require(owned.exact != source.exact && owned.envelope != source.envelope &&
            owned.col_name[0] != source.col_name[0] && owned.row_name[0] != source.row_name[0],
            "all metadata and provenance independently owned");
    require(owned.exact->col_lower_source[0] == .1 &&
            !mpq_cmp(owned.exact->col_lower[0], source.exact->col_lower[0]),
            "original column source identity retained across overlay");
    require(!mpq_cmp(owned.exact->row_upper[1], source.exact->row_upper[1]) &&
            owned.envelope->row_upper[1].upper == source.envelope->row_upper[1].upper,
            "masked row's source sides remain original");
    for (size_t j = 0; j < source.cols; j++) {
        size_t next = owned.col_start[j];
        for (size_t k = source.col_start[j]; k < source.col_start[j + 1]; k++, next++) {
            require(owned.row_index[next] == source.row_index[k] &&
                    owned.value[next] == source.value[k] &&
                    !mpq_cmp(owned.exact->value[next], source.exact->value[k]) &&
                    !memcmp(&owned.envelope->value[next], &source.envelope->value[k],
                            sizeof(LabInterval)), "CSC source coefficients preserved exactly");
        }
        if (coefficient[j] != 0) {
            require(owned.row_index[next] == source.rows && owned.value[next] == coefficient[j],
                    "new entry's row index and binary coefficient");
            dyadic_equal(owned.exact->value[next], coefficient[j], "new coefficient exact dyadic");
            require(owned.envelope->value[next].lower == coefficient[j] &&
                    owned.envelope->value[next].upper == coefficient[j],
                    "new dyadic has a point enclosure");
        }
    }
    dyadic_equal(owned.exact->row_lower[3], row.lower, "new row's RHS exact dyadic");
    require(owned.envelope->complete && owned.parsed_decimal == source.parsed_decimal,
            "complete mixed source/dyadic provenance");
    coefficient[0] = 100;
    lab_model_free(&source);
    LabLpOptions opt = {.time_limit = 5, .iteration_limit = 10000,
                        .max_dimension = LAB_LP_SPARSE_DIMENSION_LIMIT};
    LabLpSession *session = lab_lp_session_create(&owned, &opt);
    require(session != NULL, "session borrows owned augmentation");
    LabLpResult result;
    lab_lp_session_solve(session, &result);
    require(!strcmp(result.status, "optimal") && result.has_primal,
            "owned model survives original and caller-row destruction");
    double objective = 0, violation = 0;
    require(!lab_exact_evaluate(&owned, result.primal, &objective, &violation) && violation <= 1e-7,
            "augmented point passes original-unit source evaluation");
    lab_lp_result_free(&result);
    lab_lp_session_free(session); /* borrower before owner */
    lab_model_free(&owned);
}
static void failure_atomicity(void) {
    LabModel source = {0}, owned = {0};
    read_model(&source);
    LabAugmentLimits cap = limits();
    require(lab_model_append_rows(&source, 0, NULL, &cap, &owned) == LAB_AUGMENT_OK,
            "zero-row clone supported");
    double a[] = {1, 0, 0};
    LabAddedRow rows[] = {{"new", a, 1, INFINITY}, {"new", a, 0, INFINITY}};
    double *saved = owned.value;
    cap.max_nonzeros = source.nnz;
    require(lab_model_append_rows(&source, 1, rows, &cap, &owned) == LAB_AUGMENT_RESOURCE &&
            owned.value == saved, "nonzero budget failure atomic");
    cap = limits();
    require(lab_model_append_rows(&source, 2, rows, &cap, &owned) == LAB_AUGMENT_INVALID &&
            owned.value == saved, "duplicate appended names rejected atomically");
    rows[0].name = source.row_name[0];
    require(lab_model_append_rows(&source, 1, rows, &cap, &owned) == LAB_AUGMENT_INVALID,
            "source-name collision rejected");
    rows[0].name = "new";
    a[0] = NAN;
    require(lab_model_append_rows(&source, 1, rows, &cap, &owned) == LAB_AUGMENT_INVALID,
            "NaN coefficient rejected");
    a[0] = INFINITY;
    require(lab_model_append_rows(&source, 1, rows, &cap, &owned) == LAB_AUGMENT_INVALID,
            "infinite coefficient rejected");
    a[0] = 1;
    rows[0].lower = INFINITY;
    require(lab_model_append_rows(&source, 1, rows, &cap, &owned) == LAB_AUGMENT_INVALID,
            "wrong infinity side rejected");
    rows[0].lower = 1;
    cap.max_rows = source.rows;
    require(lab_model_append_rows(&source, 1, rows, &cap, &owned) == LAB_AUGMENT_RESOURCE,
            "row budget rejected");
    cap = limits();
    require(lab_model_append_rows(&source, 0, NULL, &cap, &source) == LAB_AUGMENT_INVALID,
            "self replacement rejected");
    require(owned.value == saved && owned.rows == source.rows, "all failures preserve owner");
    lab_model_free(&owned);
    lab_model_free(&source);
}
static void constant_cut_and_binary_model(void) {
    LabModel source = {0}, owned = {0};
    read_model(&source);
    lab_exact_model_free(source.exact);
    source.exact = NULL;
    lab_source_envelope_free(source.envelope);
    source.envelope = NULL;
    source.parsed_decimal = 0;
    double a[] = {0, 0, 0};
    LabAddedRow cut = {"constant_cut", a, .5, INFINITY};
    LabAugmentLimits cap = limits();
    require(lab_model_append_rows(&source, 1, &cut, &cap, &owned) == LAB_AUGMENT_OK &&
            !owned.exact && !owned.envelope && !owned.parsed_decimal && owned.nnz == source.nnz,
            "binary-only model and zero-support contradiction retained");
    LabLpOptions opt = {.time_limit = 5, .iteration_limit = 10000,
                        .max_dimension = LAB_LP_SPARSE_DIMENSION_LIMIT};
    LabLpResult result;
    lab_lp_solve(&owned, &opt, &result);
    require(!strcmp(result.status, "infeasible"), "constant appended contradiction detected");
    lab_lp_result_free(&result);
    lab_model_free(&owned);
    lab_model_free(&source);
}
int main(void) {
    source_preservation_and_lifetime();
    failure_atomicity();
    constant_cut_and_binary_model();
    puts("owned augmentation contracts passed");
    return 0;
}

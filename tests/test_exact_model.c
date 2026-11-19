#define _POSIX_C_SOURCE 200809L
#include "model.h"
#include "augment.h"
#include "exact.h"
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void require(int condition, const char *text) {
    if (!condition) {
        fprintf(stderr, "exact model contract failed: %s\n", text);
        exit(1);
    }
}
static void read_text(LabModel *m, const char *text, int success) {
    char path[] = "/tmp/milp-exact-model-XXXXXX";
    int fd = mkstemp(path);
    require(fd >= 0, "temporary model");
    FILE *f = fdopen(fd, "w");
    require(f != NULL, "temporary stream");
    require(fputs(text, f) >= 0 && !fclose(f), "write model");
    int failed = lab_model_read(m, path);
    unlink(path);
    if (success && failed)
        fprintf(stderr, "%s\n", m->error);
    require(success ? !failed : failed != 0, "model read outcome");
    require(m->parsed_decimal, "parsed model provenance origin");
}
static void decimal_parser(void) {
    mpq_t q;
    mpq_init(q);
    require(!lab_exact_decimal(q, "+.100D+1") && mpq_cmp_ui(q, 1, 1) == 0,
            "D exponent parsed exactly");
    require(!lab_exact_decimal(q, "-12.5e-1") && mpq_cmp_si(q, -5, 4) == 0,
            "signed decimal parsed exactly");
    require(!lab_exact_decimal(q, "0.000") && mpq_sgn(q) == 0, "canonical zero");
    const char *invalid[] = {"", "+", ".", "1.2.3", "1e", "1e+", "1/3", "nan",
                             "0x1p0", "1e4097", "0e99999999999999999999999999"};
    for (size_t i = 0; i < sizeof invalid / sizeof *invalid; i++)
        require(lab_exact_decimal(q, invalid[i]) != 0, "invalid/bounded decimal rejected");
    mpq_clear(q);
}
static void source_arithmetic(void) {
    LabModel m = {0};
    read_text(&m,
        "NAME DECIMAL\nROWS\n N O\n L R\nCOLUMNS\n x O 1 R .1\n"
        " x R .2\n y O 0 R .7\nRHS\n b R .3\nRANGES\n g R .2\n"
        "BOUNDS\n LO b x .1\n UP b x .3\n FR b y\nENDATA\n", 1);
    require(m.exact != NULL, "small parsed model retains exact provenance");
    require(m.exact->rows == 1 && m.exact->cols == 2 && m.exact->nnz == 2,
            "exact dimensions match CSC");
    require(mpq_cmp_ui(m.exact->value[0], 3, 10) == 0,
            "duplicate decimal .1+.2 is exactly .3");
    require(m.value[0] == .3, "exact duplicate sum is rounded once into the numerical matrix");
    require(mpq_cmp_ui(m.exact->row_lower[0], 1, 10) == 0 &&
                mpq_cmp_ui(m.exact->row_upper[0], 3, 10) == 0,
            "ranged bounds retain exact subtraction");
    require(mpq_cmp_ui(m.exact->col_lower[0], 1, 10) == 0 &&
                mpq_cmp_ui(m.exact->col_upper[0], 3, 10) == 0,
            "column decimal sides exact");
    require(m.exact->col_lower_source[0] == m.col_lower[0] &&
                m.exact->col_upper_source[0] == m.col_upper[0],
            "source float bounds identify later overlays");
    require(!isfinite(m.col_lower[1]) && !isfinite(m.col_upper[1]) &&
                mpq_sgn(m.exact->col_lower[1]) == 0 && mpq_sgn(m.exact->col_upper[1]) == 0,
            "infinite sides use zero placeholders");
    lab_model_free(&m);
}
static void topology_disagreement(void) {
    const char *models[] = {
        "NAME LOST\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 10000000000000001\n"
        " x R -10000000000000000\nRHS\n b R 1\nENDATA\n",
        "NAME EXTRA\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R .1\n"
        " x R .2\n x R -.3\nRHS\n b R 1\nENDATA\n",
        "NAME RANGE\nROWS\n N O\n L R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R .3\nRANGES\n g R 1e-20\nENDATA\n",
        "NAME BOUND\nROWS\n N O\nCOLUMNS\n x O 0\n"
        "BOUNDS\n LO b x .1\n UP b x .100000000000000001\nENDATA\n"
    };
    for (size_t i = 0; i < sizeof models / sizeof *models; i++) {
        LabModel m = {0};
        read_text(&m, models[i], 1);
        require((m.exact != NULL) == (i < 2), "coalesced support is repaired; unresolved bound topology still disables certificates");
        if(i==0)require(m.nnz==1&&m.value[0]==1,"lost nonzero is restored");
        if(i==1)require(m.nnz==0,"exact cancellation removes the spurious binary nonzero");
        lab_model_free(&m);
    }
}
static void defaults_and_nonconstraints(void) {
    LabModel m = {0};
    read_text(&m,
        "NAME DEFAULTS\nOBJNAME O\nROWS\n N IGNORE\n N O\n E R\nCOLUMNS\n"
        " mark 'MARKER' 'INTORG'\n x O 0 R 1\n y O 0 R 1\n"
        " x IGNORE 10000000000000001\n x IGNORE -10000000000000000\n"
        " end 'MARKER' 'INTEND'\nRHS\n b R .3\nRANGES\n g R -.2\n"
        "BOUNDS\n LI b y .1\nENDATA\n", 1);
    require(m.exact != NULL, "ignored nonconstraint coefficients do not affect support");
    require(mpq_cmp_ui(m.exact->row_lower[0], 1, 10) == 0 &&
                mpq_cmp_ui(m.exact->row_upper[0], 3, 10) == 0,
            "negative equality range exact");
    require(mpq_sgn(m.exact->col_lower[0]) == 0 &&
                mpq_cmp_ui(m.exact->col_upper[0], 1, 1) == 0,
            "implicit marker binary upper bound retained");
    require(isinf(m.col_upper[1]) && mpq_cmp_ui(m.exact->col_lower[1], 1, 10) == 0,
            "explicit integer bound clears implicit upper bound");
    lab_model_free(&m);
    read_text(&m,
        "NAME NEGATIVE\nROWS\n N O\nCOLUMNS\n x O 0\n"
        "BOUNDS\n UP b x -.1\nENDATA\n", 1);
    require(m.exact != NULL && isinf(m.col_lower[0]) &&
                mpq_cmp_si(m.exact->col_upper[0], -1, 10) == 0,
            "negative upper-bound default lower side remains infinite");
    lab_model_free(&m);
}
static void bounded_storage_and_cleanup(void) {
    LabModel m = {0};
    read_text(&m,
        "NAME HUGEEXP\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 0e99999999\n"
        "RHS\n b R 0\nENDATA\n", 1);
    require(m.exact == NULL, "exact exponent resource limit leaves ordinary reader intact");
    lab_model_free(&m);
    read_text(&m,
        "NAME BAD\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R .1\n"
        "RHS\n b R .2\n b R .3\nENDATA\n", 0);
    lab_model_free(&m);
    size_t capacity = 1000000;
    char *text = calloc(capacity, 1);
    require(text != NULL, "storage test allocation");
    size_t used = (size_t)snprintf(text, capacity,
        "NAME LIMIT\nROWS\n N O\n G R\nCOLUMNS\n x O 0\n");
    for (size_t k = 0; k < 20001; k++)
        used += (size_t)snprintf(text + used, capacity - used, " x R .1\n");
    snprintf(text + used, capacity - used, "RHS\n b R 1\nENDATA\n");
    read_text(&m, text, 1);
    require(m.exact == NULL && m.nnz == 1, "entry budget disables only exact provenance");
    lab_model_free(&m);
    used = (size_t)snprintf(text, capacity, "NAME COLUMNS\nROWS\n N O\nCOLUMNS\n");
    for (size_t k = 0; k < 257; k++)
        used += (size_t)snprintf(text + used, capacity - used, " x%zu O 0\n", k);
    snprintf(text + used, capacity - used, "ENDATA\n");
    read_text(&m, text, 1);
    require(m.exact == NULL && m.cols == 257, "column budget leaves float model intact");
    lab_model_free(&m);
    used = (size_t)snprintf(text, capacity, "NAME ROWS\nROWS\n N O\n");
    for (size_t k = 0; k < 129; k++)
        used += (size_t)snprintf(text + used, capacity - used, " G R%zu\n", k);
    snprintf(text + used, capacity - used, "COLUMNS\n x O 0\nENDATA\n");
    read_text(&m, text, 1);
    require(m.exact == NULL && m.rows == 129, "row budget leaves float model intact");
    lab_model_free(&m);
    used = (size_t)snprintf(text, capacity,
        "NAME LONGBOUND\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R .1\n"
        "RHS\n b R 0\nBOUNDS\n LO b x .2\n UP b x 0.");
    memset(text + used, '0', 1100);
    used += 1100;
    snprintf(text + used, capacity - used, "\nENDATA\n");
    read_text(&m, text, 1);
    require(m.exact == NULL, "oversized bound lexeme permanently disables provenance");
    lab_model_free(&m);
    free(text);
}
static void original_source_evaluation(void) {
    LabModel m = {0};
    read_text(&m,
        "NAME FIXEDDECIMAL\nROWS\n N O\n G R\nCOLUMNS\n x O 1 R 1\n"
        " z O 0 R .1\nRHS\n b R 900719925474099.25\nBOUNDS\n"
        " UP b x 1\n FX b z 9007199254740992\nENDATA\n", 1);
    double point[] = {0, 0x1p53}, objective, violation;
    require(!lab_exact_evaluate(&m, point, &objective, &violation) && objective == 0 &&
                violation >= .05 && violation <= nextafter(.05, INFINITY),
            "decimal fixed contribution exposes the hidden one-twentieth violation");
    point[0] = .05;
    require(!lab_exact_evaluate(&m, point, &objective, &violation) && objective == .05 &&
                violation == 0,
            "one-twentieth point satisfies original decimal row");
    lab_model_free(&m);
    read_text(&m,
        "NAME COSTCANCEL\nOBJNAME O\nROWS\n N IGNORED\n N O\nCOLUMNS\n"
        " x O .1\n x O .2\n z O .1\n x IGNORED 99999999\n"
        "RHS\n b O 900719925474099.25\nBOUNDS\n LO b x .1\n UP b x 2\n"
        " FX b z 9007199254740992\nENDATA\n", 1);
    require(m.exact && mpq_cmp_ui(m.exact->objective[0], 3, 10) == 0,
            "selected objective duplicate coefficients retain decimal sums");
    point[0] = 2;
    require(!lab_exact_evaluate(&m, point, &objective, &violation) && objective == .55 &&
                violation == 0,
            "large fixed objective and offset cancellation retains the decimal remainder");
    point[0] = .1;
    require(!lab_exact_evaluate(&m, point, &objective, &violation) && violation == 0,
            "unchanged decimal lower side uses source provenance");
    m.col_lower[0] = .2;
    require(!lab_exact_evaluate(&m, point, &objective, &violation) && violation >= .1,
            "tightened lower overlay uses exact binary bound");
    m.col_lower[0] = .1;
    m.col_upper[0] = .2;
    point[0] = .3;
    require(!lab_exact_evaluate(&m, point, &objective, &violation) && violation > .09 &&
                violation < .11,
            "tightened upper overlay uses exact binary bound");
    point[0] = INFINITY;
    require(lab_exact_evaluate(&m, point, &objective, &violation) && isnan(objective) &&
                isinf(violation), "nonfinite source point fails explicitly");
    lab_model_free(&m);
}
static void compensated_binary_evaluation(void) {
    double lower[] = {0}, upper[] = {INFINITY};
    double column_lower[] = {-INFINITY, -INFINITY}, column_upper[] = {INFINITY, INFINITY};
    double cost[] = {1 + 0x1p-27, -1}, value[] = {1 + 0x1p-27, -1};
    size_t start[] = {0, 1, 2}, row[] = {0, 0};
    LabModel m = {.rows = 1, .cols = 2, .nnz = 2, .row_lower = lower, .row_upper = upper,
                  .col_lower = column_lower, .col_upper = column_upper, .objective = cost,
                  .col_start = start, .row_index = row, .value = value};
    double point[] = {1 - 0x1p-27, 1}, objective, violation;
    require(!lab_exact_evaluate(&m, point, &objective, &violation) && objective == -0x1p-54 &&
                violation >= 0x1p-54 && violation <= nextafter(0x1p-54, INFINITY),
            "FMA product errors retain cancellation in binary activity and cost");
    cost[0] = 0;
    cost[1] = 0;
    value[0] = DBL_MIN;
    value[1] = 0;
    point[0] = DBL_MIN;
    upper[0] = 0;
    lower[0] = -INFINITY;
    require(!lab_exact_evaluate(&m, point, &objective, &violation) && violation > 0 &&
                violation <= 8 * DBL_TRUE_MIN,
            "unrepresentable product tails get a conservative subnormal enclosure");
    cost[0] = DBL_MAX;
    point[0] = 2;
    require(lab_exact_evaluate(&m, point, &objective, &violation) && isnan(objective) &&
                isinf(violation), "binary product overflow fails explicitly");
    cost[0] = 0;
    value[0] = 0;
    point[0] = -DBL_MAX;
    column_lower[0] = DBL_MAX;
    require(lab_exact_evaluate(&m, point, &objective, &violation),
            "bound subtraction overflow cannot disappear through maximum");
    column_lower[0] = -INFINITY;
    point[0] = 0;
    int original_rounding = fegetround();
    require(fesetround(FE_DOWNWARD) == 0, "set test rounding mode");
    int rejected = lab_exact_evaluate(&m, point, &objective, &violation);
    require(fesetround(original_rounding) == 0, "restore rounding mode");
    require(rejected, "unsupported compensated rounding mode rejected");
    double many_cost[] = {1e100, 1e50, 1, -1e100, -1e50};
    double many_point[] = {1, 1, 1, 1, 1};
    double many_lower[] = {-INFINITY, -INFINITY, -INFINITY, -INFINITY, -INFINITY};
    double many_upper[] = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY};
    size_t no_entries[] = {0, 0, 0, 0, 0, 0};
    LabModel many = {.cols = 5, .objective = many_cost, .col_lower = many_lower,
                     .col_upper = many_upper, .col_start = no_entries};
    require(!lab_exact_evaluate(&many, many_point, &objective, &violation) && objective == 1 &&
                violation == 0,
            "discarded cost tails are recovered after cancellation across three scales");
}
static void original_objective_coalescing(void) {
    LabModel m = {0};
    read_text(&m,
        "NAME LOSTCOST\nROWS\n N O\nCOLUMNS\n x O 10000000000000000\n"
        " x O -10000000000000001\nBOUNDS\n UP b x 1\nENDATA\n", 1);
    require(m.exact && m.objective[0] == -1 &&
                mpq_cmp_si(m.exact->objective[0], -1, 1) == 0,
            "decision cost uses original decimal duplicate sum before binary rounding");
    lab_model_free(&m);
    read_text(&m,
        "NAME NEAREST\nROWS\n N O\nCOLUMNS\n x O .1\n"
        " x O .2\n x O -.2\nBOUNDS\n UP b x 1\nENDATA\n", 1);
    require(m.exact && m.objective[0] == .1,
            "source objective sum converts to nearest rather than truncating");
    lab_model_free(&m);
    read_text(&m,
        "NAME EXACTZERO\nROWS\n N O\nCOLUMNS\n x O .1\n"
        " x O .2\n x O -.3\nBOUNDS\n UP b x 1\nENDATA\n", 1);
    require(m.exact && m.objective[0] == 0 && mpq_sgn(m.exact->objective[0]) == 0,
            "exact zero objective is not a spurious cancellation residual");
    lab_model_free(&m);
    char text[2048];
    const char *prefix = "NAME UNDERFLOW\nROWS\n N O\nCOLUMNS\n x O 1.";
    size_t at = strlen(prefix);
    memcpy(text, prefix, at);
    memset(text + at, '0', 400);
    at += 400;
    snprintf(text + at, sizeof text - at, "1\n x O -1\nENDATA\n");
    read_text(&m, text, 0);
    require(strstr(m.error, "objective sum") != NULL,
            "nonzero exact objective underflow fails explicitly");
    lab_model_free(&m);
}
static void encloses(LabInterval interval, const char *expected, const char *message) {
    mpq_t q, endpoint;
    mpq_inits(q, endpoint, NULL);
    require(!mpq_set_str(q, expected, 10), "test rational syntax");
    mpq_canonicalize(q);
    int contained = !isnan(interval.lower) && !isnan(interval.upper) &&
                    interval.lower <= interval.upper;
    if (isfinite(interval.lower)) {
        mpq_set_d(endpoint, interval.lower);
        contained = contained && mpq_cmp(endpoint, q) <= 0;
    } else
        contained = contained && interval.lower == -INFINITY;
    if (isfinite(interval.upper)) {
        mpq_set_d(endpoint, interval.upper);
        contained = contained && mpq_cmp(endpoint, q) >= 0;
    } else
        contained = contained && interval.upper == INFINITY;
    mpq_clears(q, endpoint, NULL);
    require(contained, message);
}
static void source_envelopes(void) {
    LabInterval n = lab_source_number_interval("1", 1);
    require(n.lower == 1 && n.upper == 1, "integer literal envelope stays exact");
    n = lab_source_number_interval("9007199254740993", 9007199254740992.0);
    encloses(n, "9007199254740993", "rounded large integer is enclosed");
    n = lab_source_number_interval("18446744073709551615", 0x1p64);
    encloses(n, "18446744073709551615", "unsigned integer endpoint conversion is safe");
    n = lab_source_number_interval(".1", .1);
    encloses(n, "1/10", "decimal literal enclosed independently of binary rounding");
    n = lab_interval_multiply((LabInterval){1 + 0x1p-27, 1 + 0x1p-27},
                              (LabInterval){1 - 0x1p-27, 1 - 0x1p-27});
    encloses(n, "18014398509481983/18014398509481984", "FMA product tail enclosed");
    n = lab_interval_add((LabInterval){DBL_MAX, DBL_MAX},
                         (LabInterval){DBL_MAX, DBL_MAX});
    require(n.lower == DBL_MAX && n.upper == INFINITY, "sum overflow widens outward");
    n = lab_interval_multiply((LabInterval){0, 0}, (LabInterval){-INFINITY, INFINITY});
    require(n.lower == 0 && n.upper == 0, "zero coefficient times unbounded box is zero");
    n = lab_interval_divide((LabInterval){1, 1}, (LabInterval){-1, 1});
    require(n.lower == -INFINITY && n.upper == INFINITY,
            "uncertain denominator sign cannot tighten");
    LabModel m = {0};
    read_text(&m,
        "NAME ENVELOPE\nROWS\n N O\n L R\nCOLUMNS\n x O 0 R .1\n x R .2\n"
        "RHS\n b R .3\nRANGES\n g R -.2\nBOUNDS\n LO b x .1\n UP b x .3\nENDATA\n", 1);
    require(m.envelope && m.envelope->complete && m.envelope->nnz == m.nnz,
            "source envelope exists independently of exact provenance");
    encloses(m.envelope->value[0], "3/10", "coalesced source coefficient enclosed");
    encloses(m.envelope->row_lower[0], "1/10", "ranged lower source bound enclosed");
    encloses(m.envelope->row_upper[0], "3/10", "ranged upper source bound enclosed");
    encloses(m.envelope->col_lower[0], "1/10", "column lower source bound enclosed");
    require(m.envelope->col_lower_source[0] == m.col_lower[0],
            "envelope retains parsed bound identity for overlays");
    lab_model_free(&m);
    read_text(&m,
        "NAME HIDDENA\nROWS\n N O\n G R\nCOLUMNS\n"
        " x O 0 R 10000000000000001\n x R -10000000000000000\n"
        "RHS\n b R 1\nBOUNDS\n UP b x 1\nENDATA\n", 1);
    require(m.exact && m.nnz == 1 && m.value[0] == 1 && m.envelope->complete,
            "exact coalescing restores nonzero source support and its enclosure");
    encloses(m.envelope->value[0], "1", "lost nonzero source coefficient retained");
    lab_model_free(&m);
}
static void sparse_single_record_decimals(void) {
    /* A long, unrelated bound disables the complete snapshots. Three
     * different source values round to the same binary coefficient. */
    const char *prefix = "NAME EXCEPTIONS\nROWS\n N O\n L A\n L B\nCOLUMNS\n"
        " x A 1.1 B 1.10000000000000001\n y A 1.10000000000000002\n"
        " z B .1\n z B .2\nRHS\n r A 100 B 100\nBOUNDS\n UP b z 2.";
    char model[2048];
    size_t n = strlen(prefix);
    memcpy(model, prefix, n);
    memset(model + n, '0', 1100);
    strcpy(model + n + 1100, "\nENDATA\n");
    LabModel original = {0}, copy = {0};
    read_text(&original, model, 1);
    require(!original.exact && !original.native_exact &&
            original.envelope->exception_count == 3, "three bounded singleton tokens");
    double coefficients[] = {1, 0, 0};
    LabAddedRow added = {"extra", coefficients, 0, 1};
    LabAugmentLimits limits = {3, 5, 1};
    require(lab_model_append_rows(&original, 1, &added, &limits, &copy) == LAB_AUGMENT_OK,
            "append row between sparse tokens");
    lab_model_free(&original);
    mpq_t q, wanted;
    mpq_inits(q, wanted, NULL);
    const char *expected[] = {"1.1", "1.10000000000000001", "1", "1.10000000000000002"};
    for (size_t k = 0; k < 4; ++k)
        require(!lab_exact_decimal(wanted, expected[k]) &&
                !lab_source_exact_coefficient(&copy, k, q) && !mpq_cmp(wanted, q),
                "owned remapped original token or exact added coefficient");
    require(lab_source_exact_coefficient(&copy, 4, q) && !mpq_cmp(wanted, q),
            "duplicate group cannot inherit a constituent token");
    copy.value[3] = 2;
    require(lab_source_exact_coefficient(&copy, 3, q) && !mpq_cmp(wanted, q),
            "changed rounded value cannot reuse the source token");
    mpq_clears(q, wanted, NULL);
    lab_model_free(&copy);
}
int main(void) {
    decimal_parser();
    source_arithmetic();
    topology_disagreement();
    defaults_and_nonconstraints();
    bounded_storage_and_cleanup();
    original_source_evaluation();
    compensated_binary_evaluation();
    original_objective_coalescing();
    source_envelopes();
    sparse_single_record_decimals();
    puts("Original-decimal provenance, exact arithmetic, topology and resource contracts passed.");
    return 0;
}

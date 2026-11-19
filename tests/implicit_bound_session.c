#include "exact.h"
#include "lp.h"
#include "solution.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(c)                                                                                   \
    do {                                                                                           \
        if (!(c)) {                                                                                \
            fprintf(stderr, "implicit bound driver line%d: %s\n", __LINE__, #c);                   \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
static void emit(const LabModel *m, const LabLpResult *r, const char *directory, size_t step,
                 int cold) {
    if (strcmp(r->status, "optimal"))
        fprintf(stderr, "%s: %s\n", r->status, r->message);
    CHECK(!strcmp(r->status, "optimal") && r->has_primal);
    char path[4096];
    CHECK(snprintf(path, sizeof path, "%s/step%zu-cold%d.sol", directory, step, cold) <
          (int)sizeof path);
    FILE *f = fopen(path, "wx");
    CHECK(f);
    CHECK(!lab_solution_write_point(f, m, r->primal));
    CHECK(!fclose(f));
    printf("{\"step\":%zu,\"cold\":%d,\"warm\":%d,\"objective\":%.17g,\"rows\":%zu,\"columns\":%zu,"
           "\"implicit\":%zu,\"packing\":%zu,\"dual_iterations\":%zu}\n",
           step, cold, r->warm_start, r->objective, r->standard_rows, r->standard_columns,
           r->implicit_bound_rows, r->packing_bound_rows, r->dual_iterations);
}
int main(int argc, char **argv) {
    CHECK(argc == 8);
    LabModel m = {0};
    CHECK(!lab_model_read(&m, argv[1]));
    CHECK(m.cols >= 2);
    if (atoi(argv[7])) {
        lab_exact_model_free(m.exact);
        m.exact = NULL;
        lab_exact_model_free(m.native_exact);
        m.native_exact = NULL;
        CHECK(m.envelope);
        m.envelope->col_lower[0] = (LabInterval){-1, 1};
    }
    double *lo = malloc(m.cols * sizeof *lo), *hi = malloc(m.cols * sizeof *hi),
           *original_lo = malloc(m.cols * sizeof *lo), *original_hi = malloc(m.cols * sizeof *hi);
    CHECK(lo && hi && original_lo && original_hi);
    memcpy(lo, m.col_lower, m.cols * sizeof *lo);
    memcpy(hi, m.col_upper, m.cols * sizeof *hi);
    memcpy(original_lo, lo, m.cols * sizeof *lo);
    memcpy(original_hi, hi, m.cols * sizeof *hi);
    LabLpOptions options = {.time_limit = 60,
                            .iteration_limit = 100000,
                            .max_dimension = 8192,
                            .packing_bounds = 1,
                            .keep_bound_rows = atoi(argv[4]),
                            .dense = atoi(argv[5]),
                            .no_scaling = atoi(argv[6])};
    if (options.dense)
        options.max_dimension = 4096;
    LabLpSession *session = lab_lp_session_create(&m, &options);
    CHECK(session);
    FILE *steps = fopen(argv[2], "r");
    CHECK(steps);
    size_t step = 0;
    char token[4][128];
    while (fscanf(steps, "%127s %127s %127s %127s", token[0], token[1], token[2], token[3]) == 4) {
        for (size_t j = 0; j < 2; j++) {
            char *end;
            lo[j] = strtod(token[2 * j], &end);
            CHECK(!*end);
            hi[j] = strtod(token[2 * j + 1], &end);
            CHECK(!*end);
            CHECK(!lab_lp_session_set_bounds(session, j, lo[j], hi[j]));
        }
        LabModel view = m;
        view.col_lower = lo;
        view.col_upper = hi;
        LabLpResult result = {0};
        lab_lp_session_solve(session, &result);
        emit(&view, &result, argv[3], step, 0);
        lab_lp_result_free(&result);
        lab_lp_solve(&view, &options, &result);
        emit(&view, &result, argv[3], step, 1);
        lab_lp_result_free(&result);
        step++;
        CHECK(!memcmp(original_lo, m.col_lower, m.cols * sizeof *lo));
        CHECK(!memcmp(original_hi, m.col_upper, m.cols * sizeof *hi));
    }
    CHECK(step && feof(steps) && !ferror(steps));
    CHECK(!fclose(steps));
    lab_lp_session_free(session);
    free(lo);
    free(hi);
    free(original_lo);
    free(original_hi);
    lab_model_free(&m);
    return 0;
}

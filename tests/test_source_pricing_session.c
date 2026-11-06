#include "lp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define REQUIRE(test)                                                                              \
    do {                                                                                           \
        if (!(test)) {                                                                             \
            fprintf(stderr, "source session contract failed at line %d: %s\n", __LINE__, #test);   \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)

static void check(const LabModel *m, const LabLpResult *r, double target, const char *directory,
                  int dense, int scale, int step, int cold) {
    if (strcmp(r->status, "optimal"))
        fprintf(stderr, "unexpected %s: %s\n", r->status, r->message);
    REQUIRE(!strcmp(r->status, "optimal") && r->has_primal);
    REQUIRE(r->primal[0] == target && r->primal[1] == target);
    REQUIRE(fabs(r->objective - target * 1e-18) < 1e-15);
    char name[4096];
    REQUIRE(snprintf(name, sizeof name, "%s/d%d-s%d-step%d-cold%d.sol", directory, dense, scale,
                     step, cold) < (int)sizeof name);
    FILE *f = fopen(name, "w");
    REQUIRE(f != NULL);
    fprintf(f, "# objective %.17g\n", r->objective);
    for (size_t j = 0; j < m->cols; j++)
        fprintf(f, "%s %.17g\n", m->col_name[j], r->primal[j]);
    REQUIRE(!fclose(f));
    printf("{\"dense\":%d,\"scale\":%d,\"step\":%d,\"cold\":%d,"
           "\"target\":%.17g,\"objective\":%.17g,\"warm\":%d,"
           "\"source_attempts\":%zu,\"source_pivots\":%zu,\"source_scans\":%zu,"
           "\"binary_paths\":%zu}\n",
           dense, scale, step, cold, target, r->objective, r->warm_start,
           r->source_pricing_attempts, r->source_pricing_pivots, r->source_support_scans,
           r->source_support_binary_paths);
}
int main(int argc, char **argv) {
    REQUIRE(argc == 3);
    LabModel m = {0};
    REQUIRE(!lab_model_read(&m, argv[1]));
    REQUIRE(!m.exact && m.native_exact && m.objective_exact);
    REQUIRE(m.cols == 300 && m.rows == 200);
    LabExactModel *native = m.native_exact;
    double *lower = malloc(m.cols * sizeof *lower), *upper = malloc(m.cols * sizeof *upper);
    REQUIRE(lower && upper);
    for (int dense = 0; dense < 2; dense++)
        for (int scale = 0; scale < 2; scale++) {
            LabLpOptions opt = {.time_limit = 20,
                                .iteration_limit = 10000,
                                .max_dimension = 4096,
                                .dense = dense,
                                .no_scaling = !scale};
            LabLpSession *s = lab_lp_session_create(&m, &opt);
            REQUIRE(s != NULL);
            const double targets[] = {1e15, 5e14, 1e15, 2.5e14, 1e15};
            int warm = 0;
            for (int step = 0; step < 5; step++) {
                memcpy(lower, m.col_lower, m.cols * sizeof *lower);
                memcpy(upper, m.col_upper, m.cols * sizeof *upper);
                for (size_t j = 0; j < 2; j++) {
                    lower[j] = step == 3 ? targets[step] : 0;
                    upper[j] = targets[step];
                    REQUIRE(!lab_lp_session_set_bounds(s, j, lower[j], upper[j]));
                }
                LabLpResult r = {0};
                lab_lp_session_solve(s, &r);
                check(&m, &r, targets[step], argv[2], dense, scale, step, 0);
                warm += r.warm_start;
                lab_lp_result_free(&r);
                LabModel view = m;
                view.col_lower = lower;
                view.col_upper = upper;
                lab_lp_solve(&view, &opt, &r);
                check(&m, &r, targets[step], argv[2], dense, scale, step, 1);
                lab_lp_result_free(&r);
                REQUIRE(!m.exact && m.native_exact == native);
                REQUIRE(m.col_lower[0] == 0 && m.col_lower[1] == 0);
                REQUIRE(m.col_upper[0] == 1e15 && m.col_upper[1] == 1e15);
            }
            REQUIRE(warm >= 2);
            lab_lp_session_free(s);
        }
    free(lower);
    free(upper);
    lab_model_free(&m);
    return 0;
}

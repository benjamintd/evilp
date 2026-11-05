#define _POSIX_C_SOURCE 200809L
#include "lp.h"
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* These tests build original models independently of the optimizer reader.
 * Every warm result is checked against both a cold solve and a known optimum. */
static void require(int condition, const char *text) {
    if (!condition) {
        fprintf(stderr, "reoptimization contract failed: %s\n", text);
        exit(1);
    }
}
static LabLpOptions options(int dense) {
    LabLpOptions result = {0};
    result.time_limit = 10;
    result.iteration_limit = 10000;
    result.max_dimension = 512;
    result.model_sha256 = "0000000000000000000000000000000000000000000000000000000000000000";
    result.dense = dense;
    return result;
}
static void point(const LabModel *m, const LabLpResult *r, const double *lo, const double *hi) {
    require(r->has_primal, "optimal solve must return a primal");
    long double *activity = calloc(m->rows ? m->rows : 1, sizeof(long double));
    require(activity != NULL, "test allocation");
    long double objective = m->offset;
    for (size_t j = 0; j < m->cols; j++) {
        double x = r->primal[j];
        require(isfinite(x) && x >= lo[j] - 1e-7 && x <= hi[j] + 1e-7,
                "returned variable must satisfy current bounds");
        objective += (long double)m->objective[j] * x;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++)
            activity[m->row_index[k]] += (long double)m->value[k] * x;
    }
    for (size_t i = 0; i < m->rows; i++)
        require(isfinite(activity[i]) && activity[i] >= (long double)m->row_lower[i] - 1e-7L &&
                    activity[i] <= (long double)m->row_upper[i] + 1e-7L,
                "returned point must satisfy original rows");
    require(fabsl(objective - r->objective) < 1e-7L, "objective must match original point");
    free(activity);
}
static void set_bounds(LabLpSession *s, size_t j, double *lo, double *hi, double lower,
                       double upper) {
    require(lab_lp_session_set_bounds(s, j, lower, upper) == 0,
            "valid bound overlay must be accepted");
    lo[j] = lower;
    hi[j] = upper;
}
/* strict_warm=-1 leaves fallback choice open;0 is first cold session solve;
 *1 requires actual reuse on an unchanged finite-bound topology. */
static void compare(LabLpSession *s, const LabModel *m, const LabLpOptions *opt, const double *lo,
                    const double *hi, const char *status, double expected, int strict_warm,
                    int requires_dual) {
    LabModel view = *m;
    view.col_lower = (double *)lo;
    view.col_upper = (double *)hi;
    LabLpResult warm = {0}, cold = {0};
    lab_lp_session_solve(s, &warm);
    lab_lp_solve(&view, opt, &cold);
    if (strcmp(warm.status, status) || strcmp(cold.status, status)) {
        fprintf(stderr, "expected %s; warm=%s (%s),cold=%s (%s)\n", status, warm.status,
                warm.message, cold.status, cold.message);
        exit(1);
    }
    if (strict_warm >= 0)
        require(warm.warm_start == strict_warm, "warm-start diagnostic must reflect actual reuse");
    if (requires_dual)
        require(warm.dual_iterations > 0, "primal-infeasible reused basis must use dual pivots");
    if (!strcmp(status, "optimal")) {
        point(m, &warm, lo, hi);
        point(m, &cold, lo, hi);
        require(fabs(warm.objective - expected) < 1e-7, "warm objective must match known optimum");
        require(fabs(cold.objective - expected) < 1e-7, "cold objective must match known optimum");
        require(fabs(warm.objective - cold.objective) < 1e-7,
                "cold and warm objectives must agree");
    } else {
        require(!warm.has_primal && !cold.has_primal, "infeasible bound overlays have no primal");
    }
    lab_lp_result_free(&warm);
    lab_lp_result_free(&cold);
}
static void bounded_sequence(int dense) {
    char *row_names[] = {"capacity"}, *col_names[] = {"x", "y"};
    double row_lo[] = {-INFINITY}, row_hi[] = {7};
    double original_lo[] = {0, 0}, original_hi[] = {5, 5}, cost[] = {3, 2};
    size_t start[] = {0, 1, 2}, index[] = {0, 0};
    double value[] = {2, 1};
    unsigned char integer[] = {0, 0};
    LabModel m = {.rows = 1,
                  .cols = 2,
                  .nnz = 2,
                  .row_name = row_names,
                  .col_name = col_names,
                  .row_lower = row_lo,
                  .row_upper = row_hi,
                  .col_lower = original_lo,
                  .col_upper = original_hi,
                  .objective = cost,
                  .integer = integer,
                  .col_start = start,
                  .row_index = index,
                  .value = value,
                  .maximize = 1};
    LabLpOptions opt = options(dense);
    LabLpSession *s = lab_lp_session_create(&m, &opt);
    require(s != NULL, "bounded session creation");
    double lo[] = {0, 0}, hi[] = {5, 5};
    compare(s, &m, &opt, lo, hi, "optimal", 13, 0, 0);
    set_bounds(s, 0, lo, hi, 2, 5);
    compare(s, &m, &opt, lo, hi, "optimal", 12, 1, 1);
    set_bounds(s, 1, lo, hi, 0, 2);
    compare(s, &m, &opt, lo, hi, "optimal", 11.5, 1, 0);
    set_bounds(s, 0, lo, hi, 2, 2);
    compare(s, &m, &opt, lo, hi, "optimal", 10, 1, 0);
    set_bounds(s, 0, lo, hi, 0, 5);
    set_bounds(s, 1, lo, hi, 0, 5);
    compare(s, &m, &opt, lo, hi, "optimal", 13, 1, 0);
    set_bounds(s, 0, lo, hi, 0, 0);
    compare(s, &m, &opt, lo, hi, "optimal", 10, 1, 0);
    set_bounds(s, 0, lo, hi, 4, 5);
    compare(s, &m, &opt, lo, hi, "infeasible", 0, -1, 0);
    set_bounds(s, 0, lo, hi, 0, 5);
    compare(s, &m, &opt, lo, hi, "optimal", 13, -1, 0);
    require(lab_lp_session_set_bounds(s, 0, NAN, 5) != 0, "NaN overlay must be rejected");
    require(lab_lp_session_set_bounds(s, 0, 0, NAN) != 0, "NaN upper must be rejected");
    require(lab_lp_session_set_bounds(s, 0, -1, 5) != 0, "outside original lower must be rejected");
    require(lab_lp_session_set_bounds(s, 0, 0, 6) != 0, "outside original upper must be rejected");
    require(lab_lp_session_set_bounds(s, 2, 0, 1) != 0, "unknown column must be rejected");
    compare(s, &m, &opt, lo, hi, "optimal", 13, -1, 0);
    if (lab_lp_session_set_bounds(s, 0, 4, 3) == 0) {
        lo[0] = 4;
        hi[0] = 3;
        compare(s, &m, &opt, lo, hi, "infeasible", 0, -1, 0);
        set_bounds(s, 0, lo, hi, 0, 5);
        compare(s, &m, &opt, lo, hi, "optimal", 13, -1, 0);
    }
    require(original_lo[0] == 0 && original_lo[1] == 0 && original_hi[0] == 5 &&
                original_hi[1] == 5,
            "session overlays must not mutate original model bounds");
    lab_lp_session_free(s);
}
static void free_and_redundant(int dense) {
    char *row_names[] = {"equal", "twice", "zero"}, *col_names[] = {"x", "y"};
    double row_lo[] = {3, 6, 0}, row_hi[] = {3, 6, 0};
    double original_lo[] = {-INFINITY, 0}, original_hi[] = {INFINITY, 4}, cost[] = {1, -1};
    size_t start[] = {0, 2, 4}, index[] = {0, 1, 0, 1};
    double value[] = {1, 2, 2, 4};
    unsigned char integer[] = {0, 0};
    LabModel m = {.rows = 3,
                  .cols = 2,
                  .nnz = 4,
                  .row_name = row_names,
                  .col_name = col_names,
                  .row_lower = row_lo,
                  .row_upper = row_hi,
                  .col_lower = original_lo,
                  .col_upper = original_hi,
                  .objective = cost,
                  .integer = integer,
                  .col_start = start,
                  .row_index = index,
                  .value = value};
    LabLpOptions opt = options(dense);
    LabLpSession *s = lab_lp_session_create(&m, &opt);
    require(s != NULL, "free/redundant session creation");
    double lo[] = {-INFINITY, 0}, hi[] = {INFINITY, 4};
    compare(s, &m, &opt, lo, hi, "optimal", -9, 0, 0);
    set_bounds(s, 0, lo, hi, -2, INFINITY);
    compare(s, &m, &opt, lo, hi, "optimal", -4.5, -1, 0);
    set_bounds(s, 0, lo, hi, -2, 0);
    compare(s, &m, &opt, lo, hi, "optimal", -4.5, -1, 0);
    set_bounds(s, 1, lo, hi, 3, 4);
    compare(s, &m, &opt, lo, hi, "infeasible", 0, -1, 0);
    set_bounds(s, 0, lo, hi, -INFINITY, INFINITY);
    compare(s, &m, &opt, lo, hi, "optimal", -9, -1, 0);
    set_bounds(s, 1, lo, hi, 1, 1);
    compare(s, &m, &opt, lo, hi, "optimal", 0, -1, 0);
    set_bounds(s, 1, lo, hi, 0, 4);
    compare(s, &m, &opt, lo, hi, "optimal", -9, -1, 0);
    require(isinf(original_lo[0]) && original_lo[0] < 0 && isinf(original_hi[0]) &&
                original_hi[0] > 0 && original_lo[1] == 0 && original_hi[1] == 4,
            "free original domain remains immutable");
    lab_lp_session_free(s);
}
static void cumulative_deadline(int dense) {
    char *names[] = {"x"};
    double lo[] = {0}, hi[] = {2}, cost[] = {-1};
    size_t start[] = {0, 0};
    unsigned char integer[] = {0};
    LabModel m = {.cols = 1,
                  .col_name = names,
                  .col_lower = lo,
                  .col_upper = hi,
                  .objective = cost,
                  .integer = integer,
                  .col_start = start};
    LabLpOptions opt = options(dense);
    opt.time_limit = .5;
    LabLpSession *s = lab_lp_session_create(&m, &opt);
    require(s != NULL, "deadline session creation");
    LabLpResult first = {0}, expired = {0};
    lab_lp_session_solve(s, &first);
    require(!strcmp(first.status, "optimal"), "initial deadline solve");
    lab_lp_result_free(&first);
    struct timespec pause = {0, 550000000};
    while (nanosleep(&pause, &pause) != 0 && errno == EINTR) {
    }
    lab_lp_session_solve(s, &expired);
    require(!strcmp(expired.status, "time_limit"), "deadline must not reset between solves");
    require(!expired.has_primal, "expired session cannot claim a fresh optimum");
    lab_lp_result_free(&expired);
    lab_lp_session_free(s);
}
static void basis_identity_is_required_before_io(int dense) {
    char *names[] = {"x"};
    double lo[] = {0}, hi[] = {2}, cost[] = {-1};
    size_t start[] = {0, 0};
    unsigned char integer[] = {0};
    LabModel m = {.cols = 1,
                  .col_name = names,
                  .col_lower = lo,
                  .col_upper = hi,
                  .objective = cost,
                  .integer = integer,
                  .col_start = start};
    LabLpOptions opt = options(dense);
    opt.model_sha256 = NULL;
    opt.basis_in = "must-not-be-opened.bas";
    require(lab_lp_session_create(&m, &opt) == NULL,
            "persistent basis input requires a model identity before any IO");
    opt.model_sha256 = "not-a-64-character-identity";
    require(lab_lp_session_create(&m, &opt) == NULL, "malformed basis identity must be rejected");
    opt.model_sha256 = NULL;
    opt.basis_in = NULL;
    opt.basis_out = "must-not-be-created.bas";
    require(lab_lp_session_create(&m, &opt) == NULL,
            "persistent basis output requires a model identity before any IO");
    opt.basis_out = NULL;
    LabLpSession *s = lab_lp_session_create(&m, &opt);
    require(s != NULL, "in-memory sessions do not require a persisted-file identity");
    lab_lp_session_free(s);
}
static unsigned next_bound(uint32_t *seed) {
    *seed = *seed * UINT32_C(1664525) + UINT32_C(1013904223);
    return (unsigned)(*seed % 21);
}
static void generated_bound_sequences(int dense) {
    char *row_names[] = {"capacity"}, *col_names[] = {"x", "y"};
    double row_lo[] = {-INFINITY}, row_hi[] = {7};
    double original_lo[] = {0, 0}, original_hi[] = {5, 5}, cost[] = {3, 2};
    size_t start[] = {0, 1, 2}, index[] = {0, 0};
    double value[] = {2, 1};
    unsigned char integer[] = {0, 0};
    LabModel m = {.rows = 1,
                  .cols = 2,
                  .nnz = 2,
                  .row_name = row_names,
                  .col_name = col_names,
                  .row_lower = row_lo,
                  .row_upper = row_hi,
                  .col_lower = original_lo,
                  .col_upper = original_hi,
                  .objective = cost,
                  .integer = integer,
                  .col_start = start,
                  .row_index = index,
                  .value = value,
                  .maximize = 1};
    LabLpOptions opt = options(dense);
    LabLpSession *s = lab_lp_session_create(&m, &opt);
    require(s != NULL, "generated sequence session creation");
    uint32_t seed = UINT32_C(712931);
    double lo[] = {0, 0}, hi[] = {5, 5};
    compare(s, &m, &opt, lo, hi, "optimal", 13, 0, 0);
    for (size_t trial = 0; trial < 100; trial++) {
        for (size_t j = 0; j < 2; j++) {
            double a = (double)next_bound(&seed) / 4;
            double b = (double)next_bound(&seed) / 4;
            set_bounds(s, j, lo, hi, fmin(a, b), fmax(a, b));
        }
        if (2 * lo[0] + lo[1] > 7) {
            compare(s, &m, &opt, lo, hi, "infeasible", 0, -1, 0);
        } else {
            /* y contributes two objective units per capacity unit, x1.5.
             * This closed-form optimum is independent of both LP engines.
             * All arithmetic is on dyadic rationals and therefore exact. */
            double y = fmin(hi[1], 7 - 2 * lo[0]);
            double x = fmin(hi[0], (7 - y) / 2);
            compare(s, &m, &opt, lo, hi, "optimal", 3 * x + 2 * y, -1, 0);
        }
    }
    lab_lp_session_free(s);
}
int main(void) {
    for (int dense = 0; dense <= 1; dense++) {
        bounded_sequence(dense);
        free_and_redundant(dense);
        generated_bound_sequences(dense);
        cumulative_deadline(dense);
        basis_identity_is_required_before_io(dense);
    }
    puts("Cold/warm reoptimization, bound overlays, dual pivots and deadline contracts passed.");
    return 0;
}

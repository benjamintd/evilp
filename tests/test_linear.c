#include "linear.h"
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(condition)                                                                           \
    do {                                                                                           \
        if (!(condition)) {                                                                        \
            fprintf(stderr, "linear contract failed at %s:%d: %s\n", __FILE__, __LINE__,           \
                    #condition);                                                                   \
            return -1;                                                                             \
        }                                                                                          \
    } while (0)

static uint32_t state = 0x594fe381U;
static uint32_t next_random(void) {
    state = state * 1664525U + 1013904223U;
    return state;
}
static size_t csc(const double *a, size_t n, size_t *start, size_t *row, double *value) {
    size_t at = 0;
    for (size_t j = 0; j < n; j++) {
        start[j] = at;
        for (size_t i = 0; i < n; i++)
            if (a[i * n + j] != 0) {
                row[at] = i;
                value[at++] = a[i * n + j];
            }
    }
    start[n] = at;
    return at;
}
static int close_vector(const double *x, const double *expected, size_t n) {
    for (size_t i = 0; i < n; i++) {
        double scale = fabs(expected[i]);
        if (scale == 0)
            scale = 1;
        if (!isfinite(x[i]) || fabs(x[i] - expected[i]) > 1e-9 * scale)
            return 0;
    }
    return 1;
}
static int randomized_scaled_known_solutions(void) {
    for (size_t test = 0; test < 160; test++) {
        size_t n = 2 + next_random() % 18;
        double *base = calloc(n * n, sizeof(double));
        double *a = calloc(n * n, sizeof(double));
        double *value = malloc(n * n * sizeof(double));
        double *copy = malloc(n * n * sizeof(double));
        double *truth = malloc(n * sizeof(double));
        double *rhs = malloc(n * sizeof(double));
        double *x = malloc(n * sizeof(double));
        double *expected = malloc(n * sizeof(double));
        size_t *start = malloc((n + 1) * sizeof(size_t));
        size_t *row = malloc(n * n * sizeof(size_t));
        size_t *rp = malloc(n * sizeof(size_t)), *cp = malloc(n * sizeof(size_t));
        int *re = malloc(n * sizeof(int)), *ce = malloc(n * sizeof(int));
        CHECK(base && a && value && copy && truth && rhs && x && expected && start && row && rp &&
              cp && re && ce);
        for (size_t i = 0; i < n; i++) {
            rp[i] = cp[i] = i;
            double magnitude = 0;
            for (size_t j = 0; j < n; j++)
                if (j != i) {
                    double v = (double)(1 + next_random() % 4) / 8;
                    if (next_random() % 2)
                        v = -v;
                    base[i * n + j] = v;
                    magnitude += fabs(v);
                }
            base[i * n + i] = magnitude + 1;
            re[i] = (int)(next_random() % 801) - 400;
            ce[i] = (int)(next_random() % 801) - 400;
            truth[i] = ((double)(next_random() % 25) - 12) / 4;
        }
        for (size_t i = n; i-- > 1;) {
            size_t j = next_random() % (i + 1), t = rp[i];
            rp[i] = rp[j];
            rp[j] = t;
            j = next_random() % (i + 1);
            t = cp[i];
            cp[i] = cp[j];
            cp[j] = t;
        }
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                a[i * n + j] = scalbn(base[rp[i] * n + cp[j]], re[i] + ce[j]);
        size_t nnz = csc(a, n, start, row, value);
        memcpy(copy, value, nnz * sizeof(double));
        for (int dense = 0; dense < 2; dense++) {
            LabLinear s = {0};
            CHECK(lab_linear_factorize(&s, n, start, row, value, dense) == 0);
            CHECK(s.valid && s.scaling_changed && s.factor.sparse == !dense);
            CHECK(memcmp(copy, value, nnz * sizeof(double)) == 0);
            for (int transpose = 0; transpose < 2; transpose++) {
                for (size_t i = 0; i < n; i++) {
                    double sum = 0;
                    for (size_t j = 0; j < n; j++)
                        sum += (transpose ? base[rp[j] * n + cp[i]] : base[rp[i] * n + cp[j]]) *
                               truth[j];
                    rhs[i] = scalbn(sum, transpose ? ce[i] : re[i]);
                    expected[i] = scalbn(truth[i], -(transpose ? re[i] : ce[i]));
                }
                CHECK(lab_linear_solve(&s, rhs, x, transpose) == 0);
                /* Units can differ by 2^800. A true zero can have a large
                 * absolute original-unit forward error but tiny natural-unit
                 * error. Check the independently known unscaled variables. */
                for (size_t i = 0; i < n; i++) {
                    double natural = scalbn(x[i], transpose ? re[i] : ce[i]);
                    CHECK(fabs(natural - truth[i]) <= 1e-9 * (1 + fabs(truth[i])));
                    if (expected[i] != 0)
                        CHECK(fabs(x[i] - expected[i]) <= 1e-9 * fabs(expected[i]));
                }
                CHECK(s.max_backward_error <= 1e-10);
                CHECK(lab_linear_solve(&s, rhs, rhs, transpose) == 0);
                for (size_t i = 0; i < n; i++)
                    CHECK(fabs(scalbn(rhs[i], transpose ? re[i] : ce[i]) - truth[i]) <=
                          1e-9 * (1 + fabs(truth[i])));
            }
            CHECK(s.solve_calls == 4);
            lab_linear_free(&s);
        }
        free(base);
        free(a);
        free(value);
        free(copy);
        free(truth);
        free(rhs);
        free(x);
        free(expected);
        free(start);
        free(row);
        free(rp);
        free(cp);
        free(re);
        free(ce);
    }
    return 0;
}
static int equilibration_extremes(void) {
    const size_t start[] = {0, 1, 2}, row[] = {0, 1};
    const double value[] = {0x1p-900, 0x1p900};
    for (int dense = 0; dense < 2; dense++) {
        LabLinear s = {0};
        CHECK(lab_linear_factorize_ex(&s, 2, start, row, value, dense, 0, 1) != 0);
        CHECK(lab_linear_factorize(&s, 2, start, row, value, dense) == 0);
        CHECK(s.scaled_rows == 2);
        double rhs[] = {0x1p-900, -0x1p900}, x[2];
        CHECK(lab_linear_solve(&s, rhs, x, 0) == 0);
        CHECK(x[0] == 1 && x[1] == -1);
        CHECK(lab_linear_solve(&s, rhs, x, 1) == 0);
        CHECK(x[0] == 1 && x[1] == -1);
        CHECK(value[0] == 0x1p-900 && value[1] == 0x1p900);
        lab_linear_free(&s);
    }
    const size_t scalar_start[] = {0, 1}, scalar_row[] = {0};
    const double extremes[] = {DBL_TRUE_MIN, DBL_MAX};
    for (size_t i = 0; i < 2; i++) {
        LabLinear s = {0};
        double x;
        CHECK(lab_linear_factorize(&s, 1, scalar_start, scalar_row, &extremes[i], 0) == 0);
        CHECK(lab_linear_solve(&s, &extremes[i], &x, 0) == 0 && x == 1);
        CHECK(lab_linear_solve(&s, &extremes[i], &x, 1) == 0 && x == 1);
        lab_linear_free(&s);
    }
    return 0;
}
static int genuine_refinement_regression(void) {
    /* Own-generated LP-like bases: five dense integer rows, eight -I/+I bound
     * rows, and unit slack columns. No imported benchmark coefficients. The
     * unit rows expose cancellation hidden by the large coefficient rows. */
    enum { P = 5, K = 8, N = 21 };
    uint32_t seed = 0x78416U;
    size_t improved = 0, rejected_without = 0;
    for (size_t trial = 0; trial < 20; trial++) {
        double matrix[N * N] = {0}, truth[N] = {0}, rhs[N] = {0}, values[N * N], x[N];
        size_t start[N + 1], rows[N * N];
        for (size_t i = 0; i < P; i++)
            for (size_t j = 0; j < K; j++) {
                seed = seed * 1664525U + 1013904223U;
                matrix[i * N + j] = (double)(1 + seed % 10000);
            }
        for (size_t j = 0; j < K; j++) {
            matrix[(P + j) * N + j] = -1;
            matrix[(P + K + j) * N + j] = 1;
        }
        for (size_t j = K; j < N; j++)
            matrix[j * N + j] = 1;
        truth[1] = 1;
        truth[2] = -1;
        for (size_t j = 3; j < K; j++)
            truth[j] = (double)(j - 2) * 100000 * (j % 2 ? -1 : 1);
        for (size_t i = 0; i < N; i++)
            for (size_t j = 0; j < N; j++)
                rhs[i] += matrix[i * N + j] * truth[j];
        csc(matrix, N, start, rows, values);
        for (int dense = 0; dense < 2; dense++) {
            LabLinear raw = {0}, refined = {0};
            CHECK(lab_linear_factorize_ex(&raw, N, start, rows, values, dense, 0, 0) == 0);
            int raw_result = lab_linear_solve(&raw, rhs, x, 0);
            CHECK((raw_result == 0) == (raw.max_backward_error <= 1e-10));
            rejected_without += raw_result != 0;
            CHECK(lab_linear_factorize_ex(&refined, N, start, rows, values, dense, 0, 1) == 0);
            CHECK(lab_linear_solve(&refined, rhs, x, 0) == 0);
            CHECK(close_vector(x, truth, N));
            CHECK(refined.max_backward_error <= 1e-10);
            if (refined.initial_backward_error > 1e-14) {
                CHECK(refined.refinement_steps > 0);
                CHECK(refined.max_backward_error < refined.initial_backward_error * 0.001);
                improved++;
            }
            lab_linear_free(&raw);
            lab_linear_free(&refined);
        }
    }
    CHECK(improved > 0 && rejected_without > 0);
    return 0;
}
static int invalid_and_alias_cases(void) {
    const size_t empty[] = {0};
    LabLinear s = {0};
    CHECK(lab_linear_factorize(&s, 0, empty, NULL, NULL, 0) == 0);
    CHECK(lab_linear_solve(&s, NULL, NULL, 0) == 0);
    const size_t start[] = {0, 1, 2}, rows[] = {0, 1};
    const double values[] = {1, 1};
    CHECK(lab_linear_factorize(&s, 2, start, rows, values, 0) == 0);
    double rhs[] = {NAN, 1}, x[2];
    CHECK(lab_linear_solve(&s, rhs, x, 0) != 0);
    rhs[0] = 2;
    CHECK(lab_linear_solve(&s, rhs, x, 0) == 0 && x[0] == 2 && x[1] == 1);
    CHECK(lab_linear_factorize(&s, 1000001, empty, NULL, NULL, 0) != 0);
    CHECK(lab_linear_factorize(&s, 4097, empty, NULL, NULL, 1) != 0);
    const size_t bad_start[] = {0, 2, 1}, bad_rows[] = {0, 2};
    CHECK(lab_linear_factorize(&s, 2, bad_start, rows, values, 0) != 0);
    CHECK(lab_linear_factorize(&s, 2, start, bad_rows, values, 0) != 0);
    const double nonfinite[] = {1, INFINITY};
    CHECK(lab_linear_factorize(&s, 2, start, rows, nonfinite, 0) != 0);
    const size_t singular_row[] = {0, 0};
    CHECK(lab_linear_factorize(&s, 2, start, singular_row, values, 0) != 0);
    CHECK(lab_linear_solve(&s, rhs, x, 0) != 0);
    lab_linear_free(&s);
    return 0;
}
static int product_roundoff_is_not_an_exact_residual(void) {
    /* In binary arithmetic, rounded .1*.3 equals rounded .03, while the exact
     * product of those input doubles differs. The residual must retain it. */
    const size_t start[] = {0, 1}, row[] = {0};
    const double value[] = {0.1}, rhs[] = {0.03};
    for (int dense = 0; dense < 2; dense++) {
        LabLinear s = {0};
        double x;
        CHECK(lab_linear_factorize_ex(&s, 1, start, row, value, dense, 0, 0) == 0);
        CHECK(lab_linear_solve(&s, rhs, &x, 0) == 0);
        CHECK(value[0] * x == rhs[0]);
        CHECK(s.max_backward_error > 1e-20 && s.max_backward_error < 1e-16);
        lab_linear_free(&s);
    }
    return 0;
}
static int rounding_mode_contract(void) {
    const size_t start[] = {0, 1, 2}, row[] = {0, 1};
    const double value[] = {1, 1}, rhs[] = {1, 2};
    const int modes[] = {FE_UPWARD, FE_DOWNWARD, FE_TOWARDZERO};
    int saved = fegetround();
    CHECK(saved == FE_TONEAREST);
    for (size_t i = 0; i < sizeof modes / sizeof *modes; i++) {
        LabLinear s = {0};
        CHECK(fesetround(modes[i]) == 0);
        int code = lab_linear_factorize(&s, 2, start, row, value, 0);
        int restored = fesetround(saved);
        CHECK(restored == 0);
        CHECK(code != 0 && !s.valid && strstr(s.error, "round-to-nearest"));
        CHECK(lab_linear_factorize(&s, 2, start, row, value, 0) == 0);
        double x[] = {123, 456};
        CHECK(fesetround(modes[i]) == 0);
        code = lab_linear_solve(&s, rhs, x, 0);
        restored = fesetround(saved);
        CHECK(restored == 0);
        CHECK(code != 0 && s.valid && strstr(s.error, "round-to-nearest"));
        CHECK(x[0] == 123 && x[1] == 456);
        CHECK(lab_linear_solve(&s, rhs, x, 0) == 0 && x[0] == 1 && x[1] == 2);
        lab_linear_free(&s);
    }
    CHECK(fegetround() == saved);
    return 0;
}
static int check_updated_known_point(LabLinear *s, const double *matrix, size_t n) {
    double *truth = malloc(n * sizeof(double)), *rhs = malloc(n * sizeof(double));
    double *x = malloc(n * sizeof(double));
    CHECK(truth && rhs && x);
    for (size_t i = 0; i < n; i++)
        truth[i] = ((double)(i % 9) - 4) / 4;
    for (int transpose = 0; transpose < 2; transpose++) {
        for (size_t i = 0; i < n; i++) {
            rhs[i] = 0;
            for (size_t j = 0; j < n; j++)
                rhs[i] += (transpose ? matrix[j * n + i] : matrix[i * n + j]) * truth[j];
        }
        CHECK(lab_linear_solve(s, rhs, x, transpose) == 0);
        CHECK(close_vector(x, truth, n));
        CHECK(s->max_backward_error <= 1e-10);
        CHECK(lab_linear_solve(s, rhs, rhs, transpose) == 0);
        CHECK(close_vector(rhs, truth, n));
        CHECK(s->max_backward_error <= 1e-10);
    }
    double *retained = calloc(n * n, sizeof(double));
    CHECK(retained);
    for (size_t j = 0; j < n; j++)
        for (size_t k = s->col_start[j]; k < s->col_start[j + 1]; k++)
            retained[s->row_index[k] * n + j] += s->value[k];
    CHECK(memcmp(retained, matrix, n * n * sizeof(double)) == 0);
    free(retained);
    free(truth);
    free(rhs);
    free(x);
    return 0;
}
static int bounded_update_chains(void) {
    enum { N = 7, STEPS = 110 };
    for (int dense = 0; dense < 2; dense++)
        for (int scaling = 0; scaling < 2; scaling++) {
            double matrix[N * N] = {0}, values[N * N], initial[N * N], column[N], copy[N];
            double direction[N], saved_direction[N];
            size_t start[N + 1], row[N * N], update_rows[N];
            for (size_t i = 0; i < N; i++)
                matrix[i * N + i] = 8;
            size_t initial_nnz = csc(matrix, N, start, row, values);
            memcpy(initial, values, initial_nnz * sizeof(double));
            LabLinear s = {0};
            CHECK(lab_linear_factorize_ex(&s, N, start, row, values, dense, scaling, 1) == 0);
            size_t rebases = 0, accepted = 0;
            for (size_t step = 0; step < STEPS; step++) {
                size_t p = (step * 3 + 1) % N;
                for (size_t i = 0; i < N; i++) {
                    column[i] = i == p ? 8 + (double)(step % 4) / 4
                                      : ((double)((i + 2 * step) % 5) - 2) / 16;
                    update_rows[i] = N - 1 - i;
                }
                memcpy(copy, column, sizeof copy);
                int supplied = (step % 3) != 0;
                if (supplied) {
                    CHECK(lab_linear_solve(&s, column, direction, 0) == 0);
                    memcpy(saved_direction, direction, sizeof direction);
                }
                size_t before = s.update_count, calls = s.solve_calls;
                int code;
                if (step % 2) {
                    double reversed[N];
                    for (size_t i = 0; i < N; i++)
                        reversed[i] = column[update_rows[i]];
                    code = lab_linear_update(&s, p, N, update_rows, reversed,
                                             supplied ? direction : NULL);
                } else
                    code = lab_linear_update_dense(&s, p, column, supplied ? direction : NULL);
                CHECK(memcmp(copy, column, sizeof copy) == 0);
                if (supplied)
                    CHECK(memcmp(saved_direction, direction, sizeof direction) == 0);
                if (before == LAB_LINEAR_UPDATE_LIMIT) {
                    CHECK(code == 1 && s.update_count == before && s.valid);
                    CHECK(strstr(s.error, "limit"));
                    CHECK(check_updated_known_point(&s, matrix, N) == 0);
                    for (size_t i = 0; i < N; i++)
                        matrix[i * N + p] = column[i];
                    csc(matrix, N, start, row, values);
                    CHECK(lab_linear_factorize_ex(&s, N, start, row, values, dense, scaling, 1) == 0);
                    CHECK(s.update_count == 0);
                    rebases++;
                } else {
                    CHECK(code == 0 && s.update_count == before + 1 && s.valid);
                    CHECK(s.solve_calls == calls + (size_t)!supplied);
                    CHECK(s.accepted_updates == s.update_count);
                    for (size_t i = 0; i < N; i++)
                        matrix[i * N + p] = column[i];
                    accepted++;
                }
                CHECK(s.factor.sparse == !dense);
                CHECK(check_updated_known_point(&s, matrix, N) == 0);
                if (step < LAB_LINEAR_UPDATE_LIMIT)
                    CHECK(memcmp(initial, values, initial_nnz * sizeof(double)) == 0);
            }
            CHECK(rebases == STEPS / (LAB_LINEAR_UPDATE_LIMIT + 1));
            CHECK(accepted + rebases == STEPS);
            lab_linear_free(&s);
        }
    return 0;
}
static int scaled_update_units(void) {
    enum { N = 3 };
    const size_t permutation[] = {2, 0, 1};
    const int exponent[] = {-700, 500, 100};
    for (int dense = 0; dense < 2; dense++) {
        double matrix[N * N] = {0}, value[N * N], column[N], rhs[N], expected[N], x[N];
        size_t start[N + 1], row[N * N];
        for (size_t j = 0; j < N; j++)
            matrix[permutation[j] * N + j] = scalbn(1, exponent[j]);
        csc(matrix, N, start, row, value);
        LabLinear s = {0};
        CHECK(lab_linear_factorize(&s, N, start, row, value, dense) == 0);
        int original_row[N], original_col[N];
        memcpy(original_row, s.row_exponent, sizeof original_row);
        memcpy(original_col, s.col_exponent, sizeof original_col);
        for (size_t step = 0; step < 18; step++) {
            size_t p = step % N;
            for (size_t i = 0; i < N; i++)
                column[i] = matrix[i * N + p] * (step % 2 ? .5 : 2);
            CHECK(lab_linear_update_dense(&s, p, column, NULL) == 0);
            for (size_t i = 0; i < N; i++)
                matrix[i * N + p] = column[i];
            CHECK(memcmp(original_row, s.row_exponent, sizeof original_row) == 0);
            CHECK(memcmp(original_col, s.col_exponent, sizeof original_col) == 0);
            for (int transpose = 0; transpose < 2; transpose++) {
                const int point_exponent[] = {300, -200, -50};
                for (size_t i = 0; i < N; i++)
                    expected[i] = scalbn(i % 2 ? -1 : 1, point_exponent[i]);
                for (size_t i = 0; i < N; i++) {
                    rhs[i] = 0;
                    for (size_t j = 0; j < N; j++)
                        rhs[i] += (transpose ? matrix[j * N + i] : matrix[i * N + j]) * expected[j];
                }
                CHECK(lab_linear_solve(&s, rhs, x, transpose) == 0);
                CHECK(close_vector(x, expected, N));
                CHECK(s.max_backward_error <= 1e-10);
            }
        }
        CHECK(s.update_count == 18);
        lab_linear_free(&s);
    }
    return 0;
}
static int rejected_updates_preserve_basis(void) {
    const size_t start[] = {0, 1, 2, 3}, rows[] = {0, 1, 2};
    const double identity[] = {1, 1, 1};
    const double matrix[] = {1, 0, 0, 1, 1, 0, 0, 0, 1};
    for (int dense = 0; dense < 2; dense++) {
        LabLinear s = {0};
        CHECK(lab_linear_factorize(&s, 3, start, rows, identity, dense) == 0);
        const double column[] = {1, 1, 0};
        CHECK(lab_linear_update_dense(&s, 0, column, column) == 0);
        size_t saved_start[4], saved_rows[4];
        double saved_values[4];
        CHECK(s.nnz == 4);
        memcpy(saved_start, s.col_start, sizeof saved_start);
        memcpy(saved_rows, s.row_index, sizeof saved_rows);
        memcpy(saved_values, s.value, sizeof saved_values);
        const double singular[] = {0, 1, 0}, wrong_direction[] = {2, 0, 0};
        const double weak_column[] = {1e-12, 1 + 1e-12, 0}, weak_direction[] = {1e-12, 1, 0};
        const double nonfinite[] = {NAN, 1, 0}, nonfinite_direction[] = {1, INFINITY, 0};
        const size_t duplicate_rows[] = {0, 0}, outside_row[] = {3};
        CHECK(lab_linear_update_dense(&s, 0, singular, NULL) == 1);
        CHECK(lab_linear_update_dense(&s, 0, weak_column, weak_direction) == 1);
        CHECK(lab_linear_update_dense(&s, 0, column, wrong_direction) == 1);
        CHECK(lab_linear_update_dense(&s, 0, nonfinite, NULL) == -1);
        CHECK(lab_linear_update_dense(&s, 0, column, nonfinite_direction) == -1);
        CHECK(lab_linear_update_dense(&s, 3, column, NULL) == -1);
        CHECK(lab_linear_update(&s, 0, 2, duplicate_rows, identity, NULL) == -1);
        CHECK(lab_linear_update(&s, 0, 1, outside_row, identity, NULL) == -1);
        int saved_round = fegetround();
        CHECK(fesetround(FE_UPWARD) == 0);
        int code = lab_linear_update_dense(&s, 0, column, NULL);
        int restored = fesetround(saved_round);
        CHECK(restored == 0 && code == -1);
        CHECK(s.valid && s.update_count == 1 && s.rejected_updates == 9);
        CHECK(memcmp(saved_start, s.col_start, sizeof saved_start) == 0);
        CHECK(memcmp(saved_rows, s.row_index, sizeof saved_rows) == 0);
        CHECK(memcmp(saved_values, s.value, sizeof saved_values) == 0);
        CHECK(check_updated_known_point(&s, matrix, 3) == 0);
        /* Aliased sparse input is copied before the retained CSC is replaced. */
        size_t count = s.col_start[1] - s.col_start[0];
        CHECK(lab_linear_update(&s, 0, count, s.row_index, s.value, NULL) == 0);
        CHECK(check_updated_known_point(&s, matrix, 3) == 0);
        lab_linear_free(&s);
    }
    CHECK(lab_linear_update(NULL, 0, 0, NULL, NULL, NULL) == -1);
    CHECK(lab_linear_update_dense(NULL, 0, NULL, NULL) == -1);
    return 0;
}
static int large_sparse_updates(void) {
    const size_t n = 12000;
    size_t *start = malloc((n + 1) * sizeof(size_t)), *rows = malloc(n * sizeof(size_t));
    double *values = malloc(n * sizeof(double)), *column = calloc(n, sizeof(double));
    double *rhs = malloc(n * sizeof(double));
    CHECK(start && rows && values && column && rhs);
    for (size_t i = 0; i < n; i++) {
        start[i] = rows[i] = i;
        values[i] = 1;
    }
    start[n] = n;
    LabLinear s = {0};
    CHECK(lab_linear_factorize(&s, n, start, rows, values, 0) == 0);
    for (size_t step = 0; step < 3; step++) {
        size_t p = step * 4000;
        memset(column, 0, n * sizeof(double));
        column[p] = 2;
        column[p + 1] = .25;
        CHECK(lab_linear_update_dense(&s, p, column, NULL) == 0);
        CHECK(s.factor.sparse && s.factor.lu == NULL && s.factor.factor_nnz == n);
        CHECK(s.nnz == n + step + 1);
        for (int transpose = 0; transpose < 2; transpose++) {
            for (size_t i = 0; i < n; i++)
                rhs[i] = 1;
            for (size_t j = 0; j <= step; j++) {
                size_t q = j * 4000;
                rhs[q] = transpose ? 2.25 : 2;
                if (!transpose)
                    rhs[q + 1] = 1.25;
            }
            CHECK(lab_linear_solve(&s, rhs, rhs, transpose) == 0);
            for (size_t i = 0; i < n; i++)
                CHECK(fabs(rhs[i] - 1) <= 1e-12);
            CHECK(s.max_backward_error <= 1e-10);
        }
    }
    CHECK(s.update_count == 3);
    lab_linear_free(&s);
    free(start);
    free(rows);
    free(values);
    free(column);
    free(rhs);
    return 0;
}
int main(void) {
    if (randomized_scaled_known_solutions() || equilibration_extremes() ||
        genuine_refinement_regression() || invalid_and_alias_cases() ||
        product_roundoff_is_not_an_exact_residual() || rounding_mode_contract() ||
        bounded_update_chains() || scaled_update_units() || rejected_updates_preserve_basis() ||
        large_sparse_updates())
        return 1;
    puts("Linear equilibration/refinement/update contracts passed.");
    return 0;
}

#include "factor.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define CHECK(condition)                                                                           \
    do {                                                                                           \
        if (!(condition)) {                                                                        \
            fprintf(stderr, "factor contract failed at %s:%d: %s\n", __FILE__, __LINE__,           \
                    #condition);                                                                   \
            return -1;                                                                             \
        }                                                                                          \
    } while (0)

static uint32_t random_state = 0x213f49a7U;
static uint32_t next_random(void) {
    random_state = random_state * 1664525U + 1013904223U;
    return random_state;
}
static int close_vector(const double *actual, const double *expected, size_t n) {
    for (size_t i = 0; i < n; i++)
        if (!isfinite(actual[i]) || fabs(actual[i] - expected[i]) > 1e-10 * (1 + fabs(expected[i])))
            return 0;
    return 1;
}
static void multiply(const double *a, const double *x, double *rhs, size_t n, int transpose) {
    for (size_t i = 0; i < n; i++) {
        long double sum = 0;
        for (size_t j = 0; j < n; j++)
            sum += (long double)a[transpose ? j * n + i : i * n + j] * x[j];
        rhs[i] = (double)sum;
    }
}
static int matrix_case(const double *a, size_t n, int duplicates) {
    size_t *start = calloc(n + 1, sizeof *start);
    size_t *index = malloc((3 * n * n + 1) * sizeof *index);
    double *value = malloc((3 * n * n + 1) * sizeof *value);
    double *expected = malloc((n ? n : 1) * sizeof *expected);
    double *rhs = malloc((n ? n : 1) * sizeof *rhs);
    double *x = malloc((n ? n : 1) * sizeof *x);
    double *dense_x = malloc((n ? n : 1) * sizeof *dense_x);
    CHECK(start && index && value && expected && rhs && x && dense_x);
    size_t nnz = 0;
    for (size_t j = 0; j < n; j++) {
        start[j] = nnz;
        /* Reversed rows exercise unsorted CSC input. Split values test exact
         * coalescing rather than relying on the dense reference's parser. */
        for (size_t i = n; i-- > 0;)
            if (a[i * n + j] != 0) {
                index[nnz] = i;
                value[nnz++] = a[i * n + j] * (duplicates ? 0.5 : 1);
                if (duplicates) {
                    index[nnz] = i;
                    value[nnz++] = 0;
                    index[nnz] = i;
                    value[nnz++] = a[i * n + j] * 0.5;
                }
            }
        start[j + 1] = nnz;
        expected[j] = ((double)(next_random() % 41) - 20) / 4;
    }
    LabFactor sparse = {0}, dense = {0};
    CHECK(lab_factorize_sparse(&sparse, n, start, index, value) == 0);
    CHECK(sparse.sparse && sparse.valid && sparse.lu == NULL);
    CHECK(lab_factorize(&dense, a, n) == 0);
    CHECK(!dense.sparse && dense.valid);
    CHECK(sparse.factor_nnz <= sparse.peak_nnz);
    for (int transpose = 0; transpose < 2; transpose++) {
        multiply(a, expected, rhs, n, transpose);
        CHECK(lab_factor_solve(&sparse, rhs, x, transpose) == 0);
        CHECK(lab_factor_solve(&dense, rhs, dense_x, transpose) == 0);
        CHECK(close_vector(x, expected, n));
        CHECK(close_vector(dense_x, expected, n));
        CHECK(close_vector(x, dense_x, n));
        /* In-place RHS solves are part of the public contract. */
        CHECK(lab_factor_solve(&sparse, rhs, rhs, transpose) == 0);
        CHECK(close_vector(rhs, expected, n));
    }
    lab_factor_free(&sparse);
    lab_factor_free(&dense);
    free(start);
    free(index);
    free(value);
    free(expected);
    free(rhs);
    free(x);
    free(dense_x);
    return 0;
}
static int randomized_permuted_cases(void) {
    for (size_t test = 0; test < 200; test++) {
        size_t n = 1 + next_random() % 31;
        double *base = calloc(n * n, sizeof *base);
        double *a = calloc(n * n, sizeof *a);
        size_t *row = malloc(n * sizeof *row), *col = malloc(n * sizeof *col);
        CHECK(base && a && row && col);
        for (size_t i = 0; i < n; i++) {
            row[i] = col[i] = i;
            double magnitude = 0;
            for (size_t j = 0; j < n; j++)
                if (j != i && next_random() % 4 == 0) {
                    base[i * n + j] = ((double)(next_random() % 9) - 4) / 8;
                    magnitude += fabs(base[i * n + j]);
                }
            base[i * n + i] = 1 + magnitude;
        }
        for (size_t i = n; i-- > 1;) {
            size_t j = next_random() % (i + 1), t = row[i];
            row[i] = row[j];
            row[j] = t;
            j = next_random() % (i + 1);
            t = col[i];
            col[i] = col[j];
            col[j] = t;
        }
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                a[i * n + j] = base[row[i] * n + col[j]];
        CHECK(matrix_case(a, n, (int)(test % 2)) == 0);
        free(base);
        free(a);
        free(row);
        free(col);
    }
    return 0;
}
static int cancellation_and_pivot_cases(void) {
    const double a[] = {0, 2, 1, 0, 0, 0, 3, 1, 4, 1, 0, 2, 1, 0, 2, 6};
    CHECK(matrix_case(a, 4, 1) == 0);
    LabFactor f = {0};
    const size_t start[] = {0, 3, 5}, row[] = {0, 0, 0, 1, 1};
    const double value[] = {3, -3, 2, 1, 2}, rhs[] = {4, 9};
    double x[2];
    CHECK(lab_factorize_sparse(&f, 2, start, row, value) == 0);
    CHECK(f.factor_nnz == 2);
    CHECK(lab_factor_solve(&f, rhs, x, 0) == 0);
    CHECK(x[0] == 2 && x[1] == 3);
    /* Backend replacement must release the previous backend correctly. */
    const double identity[] = {1, 0, 0, 1};
    CHECK(lab_factorize(&f, identity, 2) == 0);
    CHECK(!f.sparse && !f.sparse_data && f.valid);
    CHECK(lab_factorize_sparse(&f, 2, start, row, value) == 0);
    lab_factor_free(&f);
    return 0;
}
static int fill_case(void) {
    const size_t n = 8;
    size_t start[9], row[22], count = 0;
    double value[22], rhs[8], x[8];
    for (size_t j = 0; j < n; j++) {
        start[j] = count;
        if (j == 0) {
            for (size_t i = 0; i < n; i++) {
                row[count] = i;
                value[count++] = i == 0 ? 8 : 1;
            }
        } else {
            row[count] = 0;
            value[count++] = 1;
            row[count] = j;
            value[count++] = 2;
        }
        rhs[j] = j == 0 ? 15 : 3;
    }
    start[n] = count;
    LabFactor f = {0};
    CHECK(count == 22);
    CHECK(lab_factorize_sparse(&f, n, start, row, value) == 0);
    CHECK(f.factor_nnz > count && f.factor_nnz <= n * n && f.peak_nnz >= f.factor_nnz);
    for (int transpose = 0; transpose < 2; transpose++) {
        CHECK(lab_factor_solve(&f, rhs, x, transpose) == 0);
        for (size_t i = 0; i < n; i++)
            CHECK(fabs(x[i] - 1) < 1e-12);
    }
    lab_factor_free(&f);
    return 0;
}
static int invalid_and_singular_cases(void) {
    LabFactor f = {0};
    const size_t empty[] = {0};
    CHECK(lab_factorize_sparse(&f, 0, empty, NULL, NULL) == 0);
    CHECK(lab_factor_solve(&f, NULL, NULL, 0) == 0);
    const size_t start[] = {0, 1, 2}, rows[] = {0, 0};
    const double values[] = {1, 1};
    double rhs[] = {1, 1}, x[2];
    CHECK(lab_factorize_sparse(&f, 2, start, rows, values) != 0);
    CHECK(!f.valid && lab_factor_solve(&f, rhs, x, 0) != 0);
    const size_t invalid_start[] = {0, 2, 1};
    const size_t one_row[] = {0};
    const double one_value[] = {1};
    CHECK(lab_factorize_sparse(&f, 2, invalid_start, one_row, one_value) != 0);
    const size_t scalar_start[] = {0, 1}, invalid_row[] = {1};
    CHECK(lab_factorize_sparse(&f, 1, scalar_start, invalid_row, one_value) != 0);
    const double nan_value[] = {NAN}, inf_value[] = {INFINITY};
    CHECK(lab_factorize_sparse(&f, 1, scalar_start, one_row, nan_value) != 0);
    CHECK(lab_factorize_sparse(&f, 1, scalar_start, one_row, inf_value) != 0);
    const size_t split_start[] = {0, 2}, split_row[] = {0, 0};
    const double overflow[] = {DBL_MAX, DBL_MAX};
    CHECK(lab_factorize_sparse(&f, 1, split_start, split_row, overflow) != 0);
    const size_t too_many[] = {0, 20000001};
    CHECK(lab_factorize_sparse(&f, 1, too_many, NULL, NULL) != 0);
    CHECK(lab_factorize(&f, NULL, SIZE_MAX) != 0);
    CHECK(lab_factorize_sparse(&f, 1, scalar_start, one_row, one_value) == 0);
    rhs[0] = NAN;
    CHECK(lab_factor_solve(&f, rhs, x, 0) != 0);
    lab_factor_free(&f);
    return 0;
}
static int large_banded_case(void) {
    const size_t n = 4096;
    size_t *start = malloc((n + 1) * sizeof *start);
    size_t *row = malloc(3 * n * sizeof *row);
    double *value = malloc(3 * n * sizeof *value);
    double *expected = malloc(n * sizeof *expected), *rhs = calloc(n, sizeof *rhs);
    double *transpose_rhs = calloc(n, sizeof *transpose_rhs), *x = malloc(n * sizeof *x);
    CHECK(start && row && value && expected && rhs && transpose_rhs && x);
    for (size_t i = 0; i < n; i++)
        expected[i] = ((double)(i % 17) - 8) / 4;
    size_t nnz = 0;
    for (size_t j = 0; j < n; j++) {
        start[j] = nnz;
        if (j) {
            row[nnz] = j - 1;
            value[nnz++] = -0.5;
        }
        row[nnz] = j;
        value[nnz++] = 4;
        if (j + 1 < n) {
            row[nnz] = j + 1;
            value[nnz++] = -1;
        }
        start[j + 1] = nnz;
        for (size_t e = start[j]; e < nnz; e++) {
            rhs[row[e]] += value[e] * expected[j];
            transpose_rhs[j] += value[e] * expected[row[e]];
        }
    }
    LabFactor f = {0};
    CHECK(lab_factorize_sparse(&f, n, start, row, value) == 0);
    CHECK(f.lu == NULL && f.sparse && f.valid);
    CHECK(f.factor_nnz == 3 * n - 2 && f.peak_nnz == f.factor_nnz);
    CHECK(lab_factor_solve(&f, rhs, x, 0) == 0 && close_vector(x, expected, n));
    CHECK(lab_factor_solve(&f, transpose_rhs, x, 1) == 0 && close_vector(x, expected, n));
    printf("Sparse banded factor: dimension=%zu, LU entries=%zu, dense capacity=%zu\n", n,
           f.factor_nnz, n * n);
    lab_factor_free(&f);
    free(start);
    free(row);
    free(value);
    free(expected);
    free(rhs);
    free(transpose_rhs);
    free(x);
    return 0;
}
int main(void) {
    if (randomized_permuted_cases() || cancellation_and_pivot_cases() || fill_case() ||
        invalid_and_singular_cases() || large_banded_case())
        return 1;
    puts("Sparse/dense LU contracts passed: 200 random permuted matrices, transpose/in-place "
         "solves, "
         "invalid/singular inputs, and sparse memory scaling.");
    return 0;
}

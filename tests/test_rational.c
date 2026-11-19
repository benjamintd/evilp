#include "rational.h"
#include "common.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

static void require(int condition, const char *message) {
    if (!condition) {
        fprintf(stderr, "exact reconstruction contract: %s\n", message);
        exit(1);
    }
}
static void dense_contract(void) {
    mpq_t a[9], b[3], sum, product;
    for (size_t i = 0; i < 9; i++)
        mpq_init(a[i]);
    for (size_t i = 0; i < 3; i++)
        mpq_init(b[i]);
    mpq_inits(sum, product, NULL);
    /* The first pivot requires a row exchange. Solutions are non-dyadic. */
    const long matrix[] = {0, 2, 1, 3, -1, 2, 1, 1, -4};
    for (size_t i = 0; i < 9; i++)
        mpq_set_si(a[i], matrix[i], 1);
    for (size_t i = 0; i < 3; i++)
        mpq_set_si(b[i], (long)i + 1, 1);
    require(!lab_rational_solve(3, a, b, lab_now() + 5), "exact nonsingular solve");
    for (size_t i = 0; i < 3; i++) {
        mpq_set_ui(sum, 0, 1);
        for (size_t j = 0; j < 3; j++) {
            mpq_set_si(product, matrix[i * 3 + j], 1);
            mpq_mul(product, product, b[j]);
            mpq_add(sum, sum, product);
        }
        require(mpq_cmp_si(sum, (long)i + 1, 1) == 0, "exact original equations");
    }
    for (size_t i = 0; i < 9; i++)
        mpq_set_ui(a[i], 0, 1);
    require(lab_rational_solve(3, a, b, lab_now() + 5), "singularity rejects proof");
    mpq_set_ui(a[0], 1, 1);
    mpz_mul_2exp(mpq_numref(a[0]), mpq_numref(a[0]), LAB_RATIONAL_MAX_BITS + 1);
    require(lab_rational_solve(1, a, b, lab_now() + 5), "bit budget is bounded");
    require(lab_rational_solve(0, NULL, NULL, lab_now() - 1), "expired budget");
    require(lab_rational_solve(LAB_RATIONAL_MAX_DIMENSION + 1, NULL, NULL, lab_now() + 1),
            "dimension cap");
    require(!lab_rational_solve(0, NULL, NULL, lab_now() + 1), "empty exact solve");
    for (size_t i = 0; i < 9; i++)
        mpq_clear(a[i]);
    for (size_t i = 0; i < 3; i++)
        mpq_clear(b[i]);
    mpq_clears(sum, product, NULL);
    puts("Bounded exact reconstruction contracts passed.");
}
static mpq_t *vector(size_t count) {
    mpq_t *q = calloc(count ? count : 1, sizeof *q);
    require(q != NULL, "test allocation");
    for (size_t i = 0; i < count; i++)
        mpq_init(q[i]);
    return q;
}
static void release(mpq_t *q, size_t count) {
    for (size_t i = 0; i < count; i++)
        mpq_clear(q[i]);
    free(q);
}
static void solve_known(size_t n, size_t nnz, const size_t *start, const size_t *row,
                         mpq_t *value) {
    mpq_t *expected = vector(n), *rhs = vector(n), *before = vector(n), *answer = vector(n);
    mpq_t *stored = vector(nnz), term;
    mpq_init(term);
    for (size_t j = 0; j < n; j++) {
        mpq_set_si(expected[j], (long)(j % 13) - 6, (unsigned long)(j % 5) + 1);
        mpq_set_si(answer[j], 42, 1);
        for (size_t k = start[j]; k < start[j + 1]; k++) {
            mpq_mul(term, value[k], expected[j]);
            mpq_add(rhs[row[k]], rhs[row[k]], term);
        }
    }
    for (size_t i = 0; i < n; i++)
        mpq_set(before[i], rhs[i]);
    for (size_t k = 0; k < nnz; k++)
        mpq_set(stored[k], value[k]);
    require(!lab_rational_sparse_solve(n, start, row, value, rhs, answer, lab_now() + 5),
            "sparse peeled solve succeeds");
    for (size_t i = 0; i < n; i++) {
        require(mpq_cmp(answer[i], expected[i]) == 0, "known exact sparse solution");
        require(mpq_cmp(rhs[i], before[i]) == 0, "sparse RHS stays immutable");
    }
    for (size_t k = 0; k < nnz; k++)
        require(mpq_cmp(value[k], stored[k]) == 0, "sparse coefficients stay immutable");
    mpq_clear(term);
    release(expected, n);
    release(rhs, n);
    release(before, n);
    release(answer, n);
    release(stored, nnz);
}
static void singleton_rows_and_columns(void) {
    /* First row is singleton; every column initially has degree at least two.
     * This exercises direct substitution before a nontrivial dense core. */
    size_t start[] = {0, 2, 4, 6}, row[] = {0, 1, 1, 2, 1, 2};
    const long entries[] = {2, 1, 2, 1, 3, 2};
    mpq_t *value = vector(6);
    for (size_t k = 0; k < 6; k++)
        mpq_set_si(value[k], entries[k], 1);
    solve_known(3, 6, start, row, value);
    release(value, 6);
    /* A deferred singleton column is removed before a singleton-row solve.
     * Back-substitution must not subtract already-substituted terms twice. */
    size_t mixed_start[] = {0, 1, 3, 6, 9};
    size_t mixed_row[] = {0, 0, 1, 0, 2, 3, 0, 2, 3};
    const long mixed_entries[] = {1, 2, 3, -1, 2, 1, 4, 3, 2};
    value = vector(9);
    for (size_t k = 0; k < 9; k++)
        mpq_set_si(value[k], mixed_entries[k], 7);
    solve_known(4, 9, mixed_start, mixed_row, value);
    release(value, 9);
}
static void long_triangular_basis(void) {
    const size_t n = 8191, nnz = 2 * n - 1;
    size_t *start = calloc(n + 1, sizeof *start), *row = calloc(nnz, sizeof *row);
    require(start && row, "triangular basis allocation");
    mpq_t *value = vector(nnz), *rhs = vector(n), *answer = vector(n);
    for (size_t j = 0; j < n; ++j) {
        start[j] = 2 * j;
        row[2 * j] = j;
        mpq_set_ui(value[2 * j], 1, 1);
        if (j + 1 < n) {
            row[2 * j + 1] = j + 1;
            mpq_set_si(value[2 * j + 1], -1, 1);
        }
        mpq_set_ui(rhs[j], 1, 1);
    }
    start[n] = nnz;
    require(!lab_rational_sparse_solve(n, start, row, value, rhs, answer, lab_now() + 5),
            "long singleton chain");
    for (size_t j = 0; j < n; ++j)
        require(mpq_cmp_ui(answer[j], j + 1, 1) == 0, "reconstructed chain values");
    /* The two columns of the final row compete for the same pivot. */
    mpq_set_ui(value[nnz - 3], 0, 1);
    require(lab_rational_sparse_solve(n, start, row, value, rhs, answer, lab_now() + 5),
            "singular chain rejected");
    for (size_t j = 0; j < n; ++j)
        require(mpq_cmp_ui(answer[j], j + 1, 1) == 0, "singular solve preserves output");
    release(answer, n);
    release(rhs, n);
    release(value, nnz);
    free(row);
    free(start);
}
static void large_peeled_basis(void) {
    size_t n = LAB_RATIONAL_SPARSE_MAX_DIMENSION;
    size_t *start = calloc(n + 1, sizeof *start), *row = calloc(n, sizeof *row);
    require(start && row, "large identity allocation");
    mpq_t *value = vector(n);
    for (size_t j = 0; j < n; j++) {
        start[j] = j;
        row[j] = n - 1 - j;
        mpq_set_si(value[j], (j % 2) ? -3 : 2, 7);
    }
    start[n] = n;
    solve_known(n, n, start, row, value);
    release(value, n);
    free(start);
    free(row);
    /* Two hundred equations, but only a three-dimensional nontrivial core.
     * The early unit columns depend on core variables through their rows. */
    n = 200;
    size_t nnz = 3 * (n - 3) + 9;
    start = calloc(n + 1, sizeof *start);
    row = calloc(nnz, sizeof *row);
    require(start && row, "large mixed basis allocation");
    value = vector(nnz);
    size_t at = 0;
    for (size_t j = 0; j < n; j++) {
        start[j] = at;
        if (j < n - 3) {
            row[at] = j;
            mpq_set_ui(value[at++], 1, 1);
        } else {
            if (j < n - 1)
                for (size_t i = 0; i < n - 3; i++) {
                    row[at] = i;
                    mpq_set_si(value[at++], j == n - 3 ? -2 : 3, 7);
                }
            for (size_t i = n - 3; i < n; i++) {
                row[at] = i;
                mpq_set_si(value[at++], i == j ? 5 : 1, 3);
            }
        }
    }
    start[n] = at;
    require(at == nnz, "mixed sparse entry count");
    solve_known(n, nnz, start, row, value);
    release(value, nnz);
    free(start);
    free(row);
}
static void generated_sparse_cases(void) {
    uint32_t state = UINT32_C(1877723);
    for (size_t sample = 0; sample < 80; sample++) {
        size_t n = 2 + sample % 11, count = 0;
        size_t start[14], row[169];
        mpq_t *value = vector(n * n);
        for (size_t j = 0; j < n; j++) {
            start[j] = count;
            for (size_t i = 0; i < n; i++) {
                state = UINT32_C(1664525) * state + UINT32_C(1013904223);
                long coefficient = i == j ? (long)(5 * n) : (long)(state % 5) - 2;
                if (!coefficient)
                    continue;
                row[count] = (i + sample) % n;
                mpq_set_si(value[count++], coefficient, 11);
            }
        }
        start[n] = count;
        solve_known(n, count, start, row, value);
        release(value, n * n);
    }
}
static void sparse_rejections(void) {
    size_t n = LAB_RATIONAL_SPARSE_MAX_CORE + 1;
    size_t *start = calloc(n + 1, sizeof *start), *row = calloc(2 * n, sizeof *row);
    mpq_t *value = vector(2 * n), *rhs = vector(n), *answer = vector(n);
    require(start && row, "rejection test allocation");
    for (size_t j = 0; j < n; j++) {
        start[j] = 2 * j;
        row[2 * j] = j;
        row[2 * j + 1] = (j + 1) % n;
        mpq_set_ui(value[2 * j], 2, 1);
        mpq_set_ui(value[2 * j + 1], 1, 1);
        mpq_set_ui(answer[j], 42, 1);
    }
    start[n] = 2 * n;
    require(lab_rational_sparse_solve(n, start, row, value, rhs, answer, lab_now() + 5),
            "large nonpeelable core rejected");
    for (size_t j = 0; j < n; j++)
        require(mpq_cmp_ui(answer[j], 42, 1) == 0, "failure preserves output");
    /* The old boundary now reconstructs this coupled cycle exactly. */
    size_t small = LAB_RATIONAL_MAX_DIMENSION + 1;
    row[2 * small - 1] = 0;
    solve_known(small, 2 * small, start, row, value);
    size_t one_start[] = {0, 1}, one_row[] = {0};
    require(lab_rational_sparse_solve(1, one_start, one_row, value, rhs, answer,
                                     lab_now() - 1), "expired sparse budget");
    require(lab_rational_sparse_solve(LAB_RATIONAL_SPARSE_MAX_DIMENSION + 1,
                                     NULL, NULL, NULL, NULL, NULL, lab_now() + 5),
            "sparse dimension cap");
    size_t many[] = {0, LAB_RATIONAL_SPARSE_MAX_NONZEROS + 1};
    require(lab_rational_sparse_solve(1, many, NULL, NULL, rhs, answer, lab_now() + 5),
            "sparse nonzero cap");
    one_start[0] = 1;
    require(lab_rational_sparse_solve(1, one_start, one_row, value, rhs, answer,
                                     lab_now() + 5), "nonzero CSC origin rejected");
    one_start[0] = 0;
    one_row[0] = 1;
    require(lab_rational_sparse_solve(1, one_start, one_row, value, rhs, answer,
                                     lab_now() + 5), "out-of-range sparse row rejected");
    one_row[0] = 0;
    size_t duplicate_start[] = {0, 2}, duplicate_row[] = {0, 0};
    require(lab_rational_sparse_solve(1, duplicate_start, duplicate_row, value, rhs, answer,
                                     lab_now() + 5), "duplicate CSC row rejected");
    mpz_mul_2exp(mpq_numref(value[0]), mpq_numref(value[0]), LAB_RATIONAL_MAX_BITS + 1);
    require(lab_rational_sparse_solve(1, one_start, one_row, value, rhs, answer,
                                     lab_now() + 5), "sparse bit budget");
    require(!lab_rational_sparse_solve(0, NULL, NULL, NULL, NULL, NULL, lab_now() + 5),
            "empty sparse solve");
    release(value, 2 * n);
    release(rhs, n);
    release(answer, n);
    free(start);
    free(row);
}
int main(void) {
    dense_contract();
    singleton_rows_and_columns();
    long_triangular_basis();
    large_peeled_basis();
    generated_sparse_cases();
    sparse_rejections();
    puts("Sparse exact singleton elimination, bounded dense cores and immutable inputs pass.");
    return 0;
}

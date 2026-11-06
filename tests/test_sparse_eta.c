/* Exercise the private eta kernel against its frozen dense calculation. The
 * public solve additionally establishes the whole-vector finite invariant. */
#include "../src/linear.c"
#include <gmp.h>
#include <stdio.h>

#define CHECK(c)                                                                                   \
    do {                                                                                           \
        if (!(c)) {                                                                                \
            fprintf(stderr, "sparse eta line%d: %s\n", __LINE__, #c);                              \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)
static uint64_t random_state = UINT64_C(0x826f234b97d18402);
static uint64_t random_word(void) {
    random_state ^= random_state << 13;
    random_state ^= random_state >> 7;
    random_state ^= random_state << 17;
    return random_state;
}
static int dense_reference(size_t n, size_t p, const double *d, double *x, int transpose) {
    double numerator = x[p], tail = 0;
    if (transpose) {
        for (size_t i = 0; i < n; i++) {
            if (i == p)
                continue;
            double product = d[i] * x[i], error = fma(d[i], x[i], -product);
            if (!isfinite(product) || !isfinite(error) ||
                (d[i] != 0 && x[i] != 0 && product == 0) || add(&numerator, &tail, -product) ||
                add(&numerator, &tail, -error))
                return -1;
        }
        numerator += tail;
    }
    double pivot = numerator / d[p];
    if (!isfinite(pivot) || (numerator != 0 && pivot == 0))
        return -1;
    if (!transpose)
        for (size_t i = 0; i < n; i++) {
            if (i == p)
                continue;
            double value = fma(-d[i], pivot, x[i]);
            if (!isfinite(value) || (x[i] == 0 && d[i] != 0 && pivot != 0 && value == 0))
                return -1;
            x[i] = value;
        }
    x[p] = pivot;
    return 0;
}
static int differential(void) {
    enum { N = 257, CASES = 30000 };
    double d[N], left[N], right[N];
    size_t row[N];
    double value[N];
    for (size_t sample = 0; sample < CASES; sample++) {
        size_t n = 1 + random_word() % N, p = random_word() % n, count = 0;
        int span = sample % 3 == 0 ? 1000 : sample % 3 == 1 ? 100 : 8;
        for (size_t i = 0; i < n; i++) {
            int exponent = (int)(random_word() % (uint64_t)(2 * span + 1)) - span;
            left[i] = scalbn((double)(1 + random_word() % 7) / 8, exponent);
            if (random_word() % 2)
                left[i] = -left[i];
            if (random_word() % 8 == 0)
                left[i] = (random_word() % 2) ? 0.0 : -0.0;
            right[i] = left[i];
            d[i] = (random_word() % 2) ? 0.0 : -0.0;
            if (i == p || random_word() % 16 == 0) {
                exponent = (int)(random_word() % (uint64_t)(2 * span + 1)) - span;
                d[i] = scalbn((double)(1 + random_word() % 7) / 8, exponent);
                if (random_word() % 2)
                    d[i] = -d[i];
                if (i != p) {
                    row[count] = i;
                    value[count++] = d[i];
                }
            }
        }
        LabLinearEta eta = {.column = p,
                            .count = count,
                            .row_index = row,
                            .value = value,
                            .pivot = d[p],
                            .sparse = 1};
        LabLinear state = {.n = n, .update_count = 1, .eta = &eta};
        int transpose = (int)(sample % 2);
        int a = dense_reference(n, p, d, left, transpose),
            b = eta_solve(&state, 0, right, transpose);
        CHECK(a == b);
        if (!a)
            for (size_t i = 0; i < n; i++)
                CHECK(left[i] == right[i]);
    }
    return 0;
}
static int exact_dyadic_oracle(void) {
    enum { N = 16 };
    double d[N] = {0}, x[N], actual[N];
    size_t row[] = {0, 12};
    double value[] = {.25, -.5};
    d[0] = .25;
    d[7] = 2;
    d[12] = -.5;
    LabLinearEta eta = {
        .column = 7, .count = 2, .row_index = row, .value = value, .pivot = 2, .sparse = 1};
    LabLinear state = {.n = N, .update_count = 1, .eta = &eta};
    mpq_t expected, pivot, term, coefficient, original;
    mpq_inits(expected, pivot, term, coefficient, original, NULL);
    for (int transpose = 0; transpose < 2; transpose++) {
        for (size_t i = 0; i < N; i++)
            actual[i] = x[i] = ((double)i - 8) / 8;
        CHECK(!eta_solve(&state, 0, actual, transpose));
        mpq_set_d(pivot, x[7]);
        if (transpose)
            for (size_t i = 0; i < N; i++)
                if (i != 7) {
                    mpq_set_d(coefficient, d[i]);
                    mpq_set_d(term, x[i]);
                    mpq_mul(term, term, coefficient);
                    mpq_sub(pivot, pivot, term);
                }
        mpq_set_ui(coefficient, 2, 1);
        mpq_div(pivot, pivot, coefficient);
        for (size_t i = 0; i < N; i++) {
            mpq_set_d(expected, x[i]);
            if (i == 7)
                mpq_set(expected, pivot);
            else if (!transpose) {
                mpq_set_d(coefficient, d[i]);
                mpq_mul(term, coefficient, pivot);
                mpq_sub(expected, expected, term);
            }
            mpq_set_d(original, actual[i]);
            CHECK(!mpq_cmp(expected, original));
        }
    }
    mpq_clears(expected, pivot, term, coefficient, original, NULL);
    return 0;
}
static int public_storage_and_guards(void) {
    enum { N = 16 };
    size_t start[N + 1], row[N];
    double identity[N], column[N] = {0}, rhs[N] = {0}, x[N];
    for (size_t i = 0; i < N; i++) {
        start[i] = row[i] = i;
        identity[i] = 1;
    }
    start[N] = N;
    LabLinear state = {0};
    CHECK(!lab_linear_factorize(&state, N, start, row, identity, 0));
    column[3] = 1;
    column[8] = DBL_TRUE_MIN;
    CHECK(!lab_linear_update_dense(&state, 3, column, column));
    CHECK(state.eta[0].sparse && state.eta[0].count == 1 && state.eta[0].row_index[0] == 8);
    CHECK(state.eta[0].value[0] == DBL_TRUE_MIN);
    for (int transpose = 0; transpose < 2; transpose++) {
        for (size_t i = 0; i < N; i++)
            x[i] = 17;
        rhs[15] = NAN;
        CHECK(lab_linear_solve(&state, rhs, x, transpose) == -1);
        CHECK(x[0] == 17);
        rhs[15] = INFINITY;
        CHECK(lab_linear_solve(&state, rhs, x, transpose) == -1);
        rhs[15] = 0;
        rhs[transpose ? 8 : 3] = .25;
        /* The nonzero subnormal must retain the old underflow rejection. */
        CHECK(lab_linear_solve(&state, rhs, x, transpose) == -1);
        rhs[transpose ? 8 : 3] = 0;
    }
    lab_linear_free(&state);
    CHECK(!lab_linear_factorize(&state, N, start, row, identity, 0));
    for (size_t i = 0; i < N; i++)
        column[i] = i == 3 ? 1 : .125;
    CHECK(!lab_linear_update_dense(&state, 3, column, column));
    CHECK(!state.eta[0].sparse && state.eta[0].count == N);
    lab_linear_free(&state);
    /* The old dense-equivalent cap still applies to an all-sparse chain. */
    size_t n = 125001;
    size_t *starts = lab_alloc(n + 1, sizeof(size_t)), *rows = lab_alloc(n, sizeof(size_t));
    double *ones = lab_alloc(n, sizeof(double)), *unit = lab_alloc(n, sizeof(double));
    for (size_t i = 0; i < n; i++) {
        starts[i] = rows[i] = i;
        ones[i] = 1;
    }
    starts[n] = n;
    unit[0] = 1;
    CHECK(!lab_linear_factorize(&state, n, starts, rows, ones, 0));
    for (size_t k = 0; k < 31; k++) {
        CHECK(!lab_linear_update(&state, 0, 1, rows, unit, unit));
        CHECK(state.eta[k].sparse && !state.eta[k].count);
    }
    CHECK(lab_linear_update(&state, 0, 1, rows, unit, unit) == 1 && state.update_count == 31);
    lab_linear_free(&state);
    free(starts);
    free(rows);
    free(ones);
    free(unit);
    return 0;
}
int main(void) {
    if (differential() || exact_dyadic_oracle() || public_storage_and_guards())
        return 1;
    puts("Sparse eta:30000dense comparisons,exact dyadic oracle,finite/underflow/storage guards "
         "passed");
    return 0;
}

#include "linear.h"
#include "common.h"
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Error-free transforms require binary64 operations without reassociation or
 * finite-only assumptions. Refuse compiler modes that invalidate that contract. */
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__)
#error "Linear residual refinement requires strict floating-point semantics"
#endif
#if FLT_RADIX != 2 || DBL_MANT_DIG != 53 || DBL_MIN_EXP != -1021 || DBL_MAX_EXP != 1024
#error "Linear residual refinement requires binary64 double precision"
#endif
#if FLT_EVAL_METHOD != 0
#error "Linear residual refinement requires double-precision expression evaluation"
#endif
#if defined(__clang__)
#pragma STDC FENV_ACCESS ON
#pragma STDC FP_CONTRACT OFF
#endif

enum {
    LINEAR_DIMENSION_LIMIT = 1000000,
    LINEAR_ENTRY_LIMIT = 20000000,
    LINEAR_DENSE_LIMIT = 4096,
    REFINEMENT_LIMIT = 8
};

static int failure(LabLinear *s, const char *message) {
    snprintf(s->error, sizeof s->error, "%s", message);
    return -1;
}
void lab_linear_free(LabLinear *s) {
    if (!s)
        return;
    lab_factor_free(&s->factor);
    free(s->col_start);
    free(s->row_index);
    free(s->value);
    free(s->row_exponent);
    free(s->col_exponent);
    memset(s, 0, sizeof *s);
}

/* Keep a normal input coefficient normal after each scaling stage. This clips
 * the usual max-norm power-of-two scale when a row/column has a very wide range.
 * An original subnormal can be scaled upward without losing any source bits. */
static int exponent_for(double smallest, double largest) {
    if (largest == 0)
        return 0;
    int maximum = ilogb(largest), minimum = ilogb(smallest);
    int low = (DBL_MIN_EXP - 1) - minimum;
    int high = (DBL_MAX_EXP - 1) - maximum;
    if (low > high)
        return 0;
    int exponent = -maximum;
    if (exponent < low)
        exponent = low;
    if (exponent > high)
        exponent = high;
    if (exponent < -1023)
        exponent = -1023;
    if (exponent > 1023)
        exponent = 1023;
    return exponent;
}
static int exact_scale(double input, int exponent, double *output) {
    double value = scalbn(input, exponent);
    if (!isfinite(value) || (input != 0 && value == 0) || scalbn(value, -exponent) != input)
        return -1;
    *output = value;
    return 0;
}
int lab_linear_factorize_ex(LabLinear *s, size_t n, const size_t *start, const size_t *row,
                            const double *value, int dense, int scaling, int refinement) {
    lab_linear_free(s);
    s->scaling_enabled = scaling != 0;
    s->refinement_enabled = refinement != 0;
    if (fegetround() != FE_TONEAREST)
        return failure(s, "Linear factorization requires round-to-nearest");
    if (n > LINEAR_DIMENSION_LIMIT || (dense && n > LINEAR_DENSE_LIMIT) || !start || start[0] ||
        start[n] > LINEAR_ENTRY_LIMIT || (start[n] && (!row || !value)))
        return failure(s, "Invalid basis dimensions or storage");
    for (size_t j = 0; j < n; j++)
        if (start[j] > start[j + 1])
            return failure(s, "Nonmonotone basis column offsets");
    size_t nnz = start[n];
    for (size_t k = 0; k < nnz; k++)
        if (row[k] >= n || !isfinite(value[k]))
            return failure(s, "Invalid original basis entry");
    s->n = n;
    s->nnz = nnz;
    s->col_start = lab_alloc(n + 1, sizeof(size_t));
    s->row_index = lab_alloc(nnz, sizeof(size_t));
    s->value = lab_alloc(nnz, sizeof(double));
    s->row_exponent = lab_alloc(n, sizeof(int));
    s->col_exponent = lab_alloc(n, sizeof(int));
    memcpy(s->col_start, start, (n + 1) * sizeof(size_t));
    if (nnz) {
        memcpy(s->row_index, row, nnz * sizeof(size_t));
        memcpy(s->value, value, nnz * sizeof(double));
    }
    double *scaled = lab_alloc(nnz, sizeof(double));
    double *smallest = lab_alloc(n, sizeof(double));
    double *largest = lab_alloc(n, sizeof(double));
    for (size_t i = 0; i < n; i++)
        smallest[i] = INFINITY;
    if (scaling) {
        for (size_t k = 0; k < nnz; k++) {
            double a = fabs(value[k]);
            if (a == 0)
                continue;
            smallest[row[k]] = fmin(smallest[row[k]], a);
            largest[row[k]] = fmax(largest[row[k]], a);
        }
        for (size_t i = 0; i < n; i++) {
            s->row_exponent[i] = exponent_for(smallest[i], largest[i]);
            s->scaled_rows += s->row_exponent[i] != 0;
        }
    }
    int bad = 0;
    for (size_t k = 0; k < nnz; k++)
        if (exact_scale(value[k], s->row_exponent[row[k]], &scaled[k])) {
            bad = 1;
            break;
        }
    if (!bad && scaling)
        for (size_t j = 0; j < n; j++) {
            double lo = INFINITY, hi = 0;
            for (size_t k = start[j]; k < start[j + 1]; k++) {
                double a = fabs(scaled[k]);
                if (a != 0) {
                    lo = fmin(lo, a);
                    hi = fmax(hi, a);
                }
            }
            s->col_exponent[j] = exponent_for(lo, hi);
            s->scaled_columns += s->col_exponent[j] != 0;
            for (size_t k = start[j]; k < start[j + 1]; k++)
                if (exact_scale(scaled[k], s->col_exponent[j], &scaled[k]))
                    bad = 1;
        }
    free(smallest);
    free(largest);
    if (bad) {
        free(scaled);
        return failure(s, "Basis equilibration overflow or information loss");
    }
    s->scaling_changed = s->scaled_rows || s->scaled_columns;
    int result;
    if (dense) {
        double *matrix = lab_alloc(n * n, sizeof(double));
        for (size_t j = 0; j < n; j++)
            for (size_t k = start[j]; k < start[j + 1]; k++)
                matrix[row[k] * n + j] += scaled[k];
        result = lab_factorize(&s->factor, matrix, n);
        free(matrix);
    } else {
        result = lab_factorize_sparse(&s->factor, n, start, row, scaled);
    }
    free(scaled);
    if (result)
        return failure(s, "Scaled basis is singular or ill-conditioned");
    s->valid = 1;
    return 0;
}
int lab_linear_factorize(LabLinear *s, size_t n, const size_t *start, const size_t *row,
                         const double *value, int dense) {
    return lab_linear_factorize_ex(s, n, start, row, value, dense, 1, 1);
}

/* Error-free sum of two finite doubles, subject to finite intermediate sums.
 * Together with fma(a,b,-a*b), this retains product and addition roundoff even
 * on platforms where long double has no more precision than double. */
static int add(double *hi, double *lo, double value) {
    double sum = *hi + value;
    double moved = sum - *hi;
    double error = (*hi - (sum - moved)) + (value - moved);
    double tail = *lo + error;
    double next = sum + tail;
    double displacement = next - sum;
    double remainder = (sum - (next - displacement)) + (tail - displacement);
    if (!isfinite(sum) || !isfinite(tail) || !isfinite(next) || !isfinite(remainder))
        return -1;
    *hi = next;
    *lo = remainder;
    return 0;
}
static void magnitude_add(double *scale, double *sum, double value) {
    value = fabs(value);
    if (value > *scale) {
        *sum = *sum * (*scale / value) + 1;
        *scale = value;
    } else if (value != 0) {
        *sum += value / *scale;
    }
}
static int residual(const LabLinear *s, const double *rhs, const double *x, double *r,
                    int transpose, double *error) {
    size_t n = s->n;
    double *lo = lab_alloc(n, sizeof(double));
    double *scale = lab_alloc(n, sizeof(double));
    double *magnitude = lab_alloc(n, sizeof(double));
    for (size_t i = 0; i < n; i++) {
        r[i] = rhs[i];
        scale[i] = magnitude[i] = 1;
        magnitude_add(&scale[i], &magnitude[i], rhs[i]);
    }
    int bad = 0;
    for (size_t j = 0; j < n && !bad; j++)
        for (size_t k = s->col_start[j]; k < s->col_start[j + 1]; k++) {
            size_t i = s->row_index[k], destination = transpose ? j : i;
            double a = s->value[k], b = x[transpose ? i : j];
            double product = a * b;
            double product_error = fma(a, b, -product);
            if (!isfinite(product) || !isfinite(product_error) ||
                add(&r[destination], &lo[destination], -product) ||
                add(&r[destination], &lo[destination], -product_error)) {
                bad = 1;
                break;
            }
            magnitude_add(&scale[destination], &magnitude[destination], product);
        }
    double worst = 0;
    for (size_t i = 0; i < n && !bad; i++) {
        r[i] += lo[i];
        double relative = fabs(r[i]) / scale[i] / magnitude[i];
        if (!isfinite(r[i]) || !isfinite(relative))
            bad = 1;
        worst = fmax(worst, relative);
    }
    free(lo);
    free(scale);
    free(magnitude);
    *error = bad ? INFINITY : worst;
    return bad ? -1 : 0;
}
static int raw_solve(LabLinear *s, const double *rhs, double *x, int transpose) {
    size_t n = s->n;
    double *scaled_rhs = lab_alloc(n, sizeof(double));
    double *y = lab_alloc(n, sizeof(double));
    int bad = 0;
    for (size_t i = 0; i < n; i++) {
        int exponent = transpose ? s->col_exponent[i] : s->row_exponent[i];
        if (exact_scale(rhs[i], exponent, &scaled_rhs[i])) {
            bad = 1;
            break;
        }
    }
    if (!bad && lab_factor_solve(&s->factor, scaled_rhs, y, transpose))
        bad = 1;
    for (size_t i = 0; i < n && !bad; i++) {
        int exponent = transpose ? s->row_exponent[i] : s->col_exponent[i];
        x[i] = scalbn(y[i], exponent);
        if (!isfinite(x[i]) || (y[i] != 0 && x[i] == 0))
            bad = 1;
    }
    free(scaled_rhs);
    free(y);
    return bad ? -1 : 0;
}
int lab_linear_solve(LabLinear *s, const double *rhs, double *x, int transpose) {
    s->solve_calls++;
    s->initial_backward_error = s->max_backward_error = INFINITY;
    s->error[0] = 0;
    if (fegetround() != FE_TONEAREST)
        return failure(s, "Linear solve requires round-to-nearest");
    if (!s->valid || (s->n && (!rhs || !x)))
        return failure(s, "Invalid linear solve state or vectors");
    size_t n = s->n;
    double *original_rhs = lab_alloc(n, sizeof(double));
    double *best = lab_alloc(n, sizeof(double));
    double *r = lab_alloc(n, sizeof(double));
    double *correction = lab_alloc(n, sizeof(double));
    double *candidate = lab_alloc(n, sizeof(double));
    if (n)
        memcpy(original_rhs, rhs, n * sizeof(double));
    int failed = 0;
    for (size_t i = 0; i < n; i++)
        if (!isfinite(original_rhs[i]))
            failed = 1;
    double best_error = INFINITY;
    if (failed || raw_solve(s, original_rhs, best, transpose) ||
        residual(s, original_rhs, best, r, transpose, &best_error)) {
        failure(s, "Initial original-unit basis solve failed");
        failed = 1;
        goto done;
    }
    s->initial_backward_error = best_error;
    for (int step = 0; s->refinement_enabled && step < REFINEMENT_LIMIT && best_error > 1e-14;
         step++) {
        if (raw_solve(s, r, correction, transpose))
            break;
        int changed = 0, finite = 1;
        for (size_t i = 0; i < n; i++) {
            candidate[i] = best[i] + correction[i];
            changed |= candidate[i] != best[i];
            finite &= isfinite(candidate[i]);
        }
        if (!changed || !finite)
            break;
        s->refinement_steps++;
        double next_error;
        if (residual(s, original_rhs, candidate, r, transpose, &next_error) ||
            next_error >= best_error)
            break;
        if (n)
            memcpy(best, candidate, n * sizeof(double));
        best_error = next_error;
    }
    if (n)
        memcpy(x, best, n * sizeof(double));
    s->max_backward_error = best_error;
    if (best_error > 1e-10) {
        failure(s, "Original-unit basis backward error exceeds 1e-10");
        failed = 1;
    }
done:
    free(original_rhs);
    free(best);
    free(r);
    free(correction);
    free(candidate);
    return failed ? -1 : 0;
}

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
    LINEAR_DIMENSION_LIMIT = LAB_LINEAR_DIMENSION_LIMIT,
    LINEAR_ENTRY_LIMIT = 128000000,
    LINEAR_DENSE_LIMIT = 4096,
    REFINEMENT_LIMIT = 8
};
#define LINEAR_ETA_ENTRY_LIMIT 4000000

static int failure(LabLinear *s, const char *message) {
    snprintf(s->error, sizeof s->error, "%s", message);
    return -1;
}
void lab_linear_csc_free(LabLinearCSC *copy) {
    if (!copy)
        return;
    free(copy->col_start);
    free(copy->row_index);
    free(copy->value);
    memset(copy, 0, sizeof *copy);
}
void lab_linear_free(LabLinear *s) {
    if (!s)
        return;
    lab_factor_free(&s->factor);
    free(s->base_col_start);
    free(s->base_row_index);
    free(s->base_value);
    for (size_t i = 0; i < s->overlay_count; i++) {
        free(s->overlay[i].row_index);
        free(s->overlay[i].value);
    }
    free(s->overlay);
    free(s->overlay_slot);
    free(s->row_exponent);
    free(s->col_exponent);
    for (size_t i = 0; i < s->update_count; i++)
        free(s->eta_value[i]);
    free(s->eta_column);
    free(s->eta_value);
    memset(s, 0, sizeof *s);
}
static void current_column(const LabLinear *s, size_t column, size_t *count,
                            const size_t **row, const double **value) {
    size_t slot = s->overlay_slot ? s->overlay_slot[column] : 0;
    if (slot) {
        const LabLinearColumn *c = &s->overlay[slot - 1];
        *count = c->nnz;
        *row = c->row_index;
        *value = c->value;
    } else {
        size_t start = s->base_col_start[column];
        *count = s->base_col_start[column + 1] - start;
        *row = s->base_row_index + start;
        *value = s->base_value + start;
    }
}
int lab_linear_get_column(const LabLinear *s, size_t column, size_t *count,
                          const size_t **row, const double **value) {
    if (count)
        *count = 0;
    if (row)
        *row = NULL;
    if (value)
        *value = NULL;
    if (!s || !s->valid || column >= s->n || !count || !row || !value)
        return -1;
    current_column(s, column, count, row, value);
    return 0;
}
static size_t counter_add(size_t value, size_t increment) {
    return value > SIZE_MAX - increment ? SIZE_MAX : value + increment;
}
int lab_linear_copy_csc(LabLinear *s, LabLinearCSC *copy) {
    if (!s || !s->valid || !copy || copy->n || copy->nnz || copy->col_start ||
        copy->row_index || copy->value)
        return -1;
    LabLinearCSC result = {.n = s->n, .nnz = s->nnz};
    result.col_start = lab_alloc(s->n + 1, sizeof(size_t));
    result.row_index = lab_alloc(s->nnz, sizeof(size_t));
    result.value = lab_alloc(s->nnz, sizeof(double));
    size_t at = 0;
    for (size_t j = 0; j < s->n; j++) {
        size_t count;
        const size_t *row;
        const double *value;
        current_column(s, j, &count, &row, &value);
        if (at > s->nnz || count > s->nnz - at) {
            lab_linear_csc_free(&result);
            return -1;
        }
        result.col_start[j] = at;
        if (count) {
            memcpy(result.row_index + at, row, count * sizeof(size_t));
            memcpy(result.value + at, value, count * sizeof(double));
        }
        at += count;
    }
    if (at != s->nnz) {
        lab_linear_csc_free(&result);
        return -1;
    }
    result.col_start[s->n] = at;
    *copy = result;
    s->csc_materializations = counter_add(s->csc_materializations, 1);
    s->csc_materialized_entries = counter_add(s->csc_materialized_entries, at);
    return 0;
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
    s->nnz = s->base_nnz = nnz;
    s->base_col_start = lab_alloc(n + 1, sizeof(size_t));
    s->base_row_index = lab_alloc(nnz, sizeof(size_t));
    s->base_value = lab_alloc(nnz, sizeof(double));
    s->row_exponent = lab_alloc(n, sizeof(int));
    s->col_exponent = lab_alloc(n, sizeof(int));
    memcpy(s->base_col_start, start, (n + 1) * sizeof(size_t));
    if (nnz) {
        memcpy(s->base_row_index, row, nnz * sizeof(size_t));
        memcpy(s->base_value, value, nnz * sizeof(double));
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
    s->factor.core_refinement = refinement != 0;
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
    for (size_t j = 0; j < n && !bad; j++) {
        size_t count;
        const size_t *row;
        const double *value;
        current_column(s, j, &count, &row, &value);
        for (size_t k = 0; k < count; k++) {
            size_t i = row[k], destination = transpose ? j : i;
            double a = value[k], b = x[transpose ? i : j];
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
static int base_solve(LabLinear *s, const double *rhs, double *x, int transpose) {
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
/* B_new = B_old E, with E equal to the identity except for column p=d.
 * E^-1 v uses x_p=v_p/d_p and x_i=v_i-d_i*x_p. E^-T only changes
 * component p: x_p=(v_p-sum_(i!=p) d_i*v_i)/d_p. The normal chain is
 * chronological after the base solve; the transpose chain is reversed before
 * the base transpose solve. All eta vectors use original basis coordinates. */
static int eta_solve(const LabLinear *s, size_t update, double *x, int transpose) {
    size_t p = s->eta_column[update];
    const double *d = s->eta_value[update];
    double numerator = x[p], tail = 0;
    if (transpose) {
        for (size_t i = 0; i < s->n; i++) {
            if (i == p)
                continue;
            double product = d[i] * x[i];
            double error = fma(d[i], x[i], -product);
            if (!isfinite(product) || !isfinite(error) ||
                (d[i] != 0 && x[i] != 0 && product == 0) ||
                add(&numerator, &tail, -product) || add(&numerator, &tail, -error))
                return -1;
        }
        numerator += tail;
    }
    double pivot_value = numerator / d[p];
    if (!isfinite(pivot_value) || (numerator != 0 && pivot_value == 0))
        return -1;
    if (!transpose)
        for (size_t i = 0; i < s->n; i++) {
            if (i == p)
                continue;
            double value = fma(-d[i], pivot_value, x[i]);
            if (!isfinite(value) ||
                (x[i] == 0 && d[i] != 0 && pivot_value != 0 && value == 0))
                return -1;
            x[i] = value;
        }
    x[p] = pivot_value;
    return 0;
}
static int raw_solve(LabLinear *s, const double *rhs, double *x, int transpose) {
    if (!transpose) {
        if (base_solve(s, rhs, x, 0))
            return -1;
        for (size_t i = 0; i < s->update_count; i++)
            if (eta_solve(s, i, x, 0))
                return -1;
        return 0;
    }
    double *work = lab_alloc(s->n, sizeof(double));
    if (s->n)
        memcpy(work, rhs, s->n * sizeof(double));
    int failed = 0;
    for (size_t i = s->update_count; i-- > 0;)
        if (eta_solve(s, i, work, 1)) {
            failed = 1;
            break;
        }
    if (!failed)
        failed = base_solve(s, work, x, 1);
    free(work);
    return failed ? -1 : 0;
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

static int update_failure(LabLinear *s, const char *message, int status) {
    if (s) {
        s->rejected_updates++;
        failure(s, message);
    }
    return status;
}
int lab_linear_update(LabLinear *s, size_t column, size_t count, const size_t *row,
                      const double *value, const double *direction) {
    if (!s)
        return -1;
    s->error[0] = 0;
    if (fegetround() != FE_TONEAREST)
        return update_failure(s, "Linear update requires round-to-nearest", -1);
    if (!s->valid || column >= s->n || count > s->n || (count && (!row || !value)))
        return update_failure(s, "Invalid basis update dimensions or storage", -1);
    if (s->update_count >= LAB_LINEAR_UPDATE_LIMIT ||
        s->n > LINEAR_ETA_ENTRY_LIMIT / (s->update_count + 1))
        return update_failure(s, "Basis update limit requires refactorization", 1);
    size_t n = s->n;
    double *new_column = lab_alloc(n, sizeof(double));
    unsigned char *seen = lab_alloc(n, sizeof(unsigned char));
    double *d = lab_alloc(n, sizeof(double));
    double *work = lab_alloc(n, sizeof(double));
    int status = -1;
    const char *message = "Invalid or duplicate basis update entry";
    for (size_t i = 0; i < count; i++) {
        if (row[i] >= n || seen[row[i]] || !isfinite(value[i]))
            goto done;
        seen[row[i]] = 1;
        new_column[row[i]] = value[i];
    }
    if (direction) {
        for (size_t i = 0; i < n; i++) {
            if (!isfinite(direction[i])) {
                message = "Nonfinite basis update direction";
                goto done;
            }
            d[i] = direction[i];
        }
    } else if (lab_linear_solve(s, new_column, d, 0)) {
        status = 1;
        message = "Basis update direction requires refactorization";
        goto done;
    }
    double error;
    status = 1;
    message = "Inaccurate basis update direction requires refactorization";
    /* A stricter local target leaves room for accumulation across the bounded
     * chain; every subsequent solve still checks the unchanged 1e-10 guard. */
    if (residual(s, new_column, d, work, 0, &error) || error > 1e-12)
        goto done;
    double largest = 0;
    for (size_t i = 0; i < n; i++)
        largest = fmax(largest, fabs(d[i]));
    message = "Weak basis update pivot requires refactorization";
    if (largest == 0 || d[column] == 0 || fabs(d[column]) < 1e-8 * largest ||
        !isfinite(1 / d[column]))
        goto done;
    for (size_t i = 0; i < n; i++)
        if (!isfinite(d[i] / d[column]))
            goto done;
    size_t old_count;
    const size_t *old_row;
    const double *old_value;
    current_column(s, column, &old_count, &old_row, &old_value);
    size_t nnz = s->nnz - old_count + count;
    message = "Updated basis storage limit requires refactorization";
    if (nnz > LINEAR_ENTRY_LIMIT)
        goto done;
    /* Retain only the replacement column. The base CSC and other overlays
     * stay untouched. Copy before releasing an old overlay because the input
     * row/value arrays may be borrowed from either current or base storage. */
    size_t *rows = lab_alloc(count, sizeof(size_t));
    double *values = lab_alloc(count, sizeof(double));
    if (count) {
        memcpy(rows, row, count * sizeof(size_t));
        memcpy(values, value, count * sizeof(double));
    }
    if (!s->overlay_slot) {
        s->overlay_slot = lab_alloc(n, sizeof(unsigned char));
        s->overlay = lab_alloc(LAB_LINEAR_UPDATE_LIMIT, sizeof(LabLinearColumn));
    }
    if (!s->eta_column) {
        s->eta_column = lab_alloc(LAB_LINEAR_UPDATE_LIMIT, sizeof(size_t));
        s->eta_value = lab_alloc(LAB_LINEAR_UPDATE_LIMIT, sizeof(double *));
    }
    s->eta_column[s->update_count] = column;
    s->eta_value[s->update_count] = d;
    s->update_count++;
    s->accepted_updates++;
    d = NULL;
    size_t slot = s->overlay_slot[column];
    if (!slot) {
        slot = ++s->overlay_count;
        s->overlay_slot[column] = (unsigned char)slot;
    }
    LabLinearColumn *replacement = &s->overlay[slot - 1];
    free(replacement->row_index);
    free(replacement->value);
    replacement->column = column;
    replacement->nnz = count;
    replacement->row_index = rows;
    replacement->value = values;
    s->nnz = nnz;
    s->column_entries_copied = counter_add(s->column_entries_copied, count);
    s->error[0] = 0;
    status = 0;
done:
    free(new_column);
    free(seen);
    free(d);
    free(work);
    return status ? update_failure(s, message, status) : 0;
}
int lab_linear_update_dense(LabLinear *s, size_t column, const double *new_column,
                            const double *direction) {
    if (!s)
        return -1;
    if (!s->valid || column >= s->n || !new_column)
        return update_failure(s, "Invalid dense basis update state or vector", -1);
    size_t *row = lab_alloc(s->n, sizeof(size_t));
    double *value = lab_alloc(s->n, sizeof(double));
    size_t count = 0;
    for (size_t i = 0; i < s->n; i++)
        if (new_column[i] != 0) {
            row[count] = i;
            value[count++] = new_column[i];
        }
    int status = lab_linear_update(s, column, count, row, value, direction);
    free(row);
    free(value);
    return status;
}

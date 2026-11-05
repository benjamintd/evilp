#include "lp.h"
#include "common.h"
#include "linear.h"
#include "exact.h"
#include "rational.h"
#include "pivot.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Sparse standard form and revised simplex. The dense factorization is retained
 * only as an explicitly selected reference backend. */
typedef struct {
    size_t count, capacity;
    size_t *row;
    double *value;
} Column;
typedef struct {
    size_t m, n, structural, total_rows;
    Column *a;
    double *b, *cost, *shift, *sign;
    size_t *basis, *original, *row_id, *art_row;
    size_t *source_row;
    double *row_scale;
    double *structural_norm;
    double *equation_scale;
    unsigned char *row_upper;
    unsigned char *artificial;
    unsigned char *basic;
    size_t *lower_row, *upper_row;
    unsigned char *transform;
    LabLinear linear;
    double *xb, *pi, *cb, *direction, *column;
    size_t scaled_rows, scaled_columns;
    const LabLpOptions *options;
    const LabModel *model;
    LabLpResult *result;
    double started;
} Work;
static int factor_basis(Work *w);
static void put(Work *w, size_t row, size_t col, double value) {
    if (value == 0)
        return;
    Column *c = &w->a[col];
    if (c->count == c->capacity) {
        size_t cap = c->capacity ? 2 * c->capacity : 4;
        size_t *ri = realloc(c->row, cap * sizeof(size_t));
        double *va = realloc(c->value, cap * sizeof(double));
        if (!ri || !va)
            exit(70);
        c->row = ri;
        c->value = va;
        c->capacity = cap;
    }
    c->row[c->count] = row;
    c->value[c->count++] = value;
}
static long double dot(const Work *w, size_t j, const double *v, long double *magnitude) {
    const Column *c = &w->a[j];
    long double sum = 0;
    for (size_t k = 0; k < c->count; k++) {
        long double term = (long double)c->value[k] * v[c->row[k]];
        sum += term;
        if (magnitude)
            *magnitude += fabsl(term);
    }
    return sum;
}
static void scatter(const Work *w, size_t j, double *v) {
    memset(v, 0, w->m * sizeof(double));
    const Column *c = &w->a[j];
    for (size_t k = 0; k < c->count; k++)
        v[c->row[k]] += c->value[k];
}
static long double column_norm(const Work *w, size_t j) {
    if (j < w->structural)
        return w->structural_norm[j];
    /* Slack and artificial columns have at most one unscaled unit entry. */
    return w->a[j].count ? fabs(w->a[j].value[0]) : 0;
}
static int linear_solve(Work *w, const double *rhs, double *x, int transpose) {
    double *saved = NULL;
    if (rhs == x) {
        saved = lab_alloc(w->m, sizeof(double));
        memcpy(saved, rhs, w->m * sizeof(double));
        rhs = saved;
    }
    size_t before = w->linear.refinement_steps;
    int failed = lab_linear_solve(&w->linear, rhs, x, transpose);
    w->result->refinement_steps += w->linear.refinement_steps - before;
    if (isfinite(w->linear.initial_backward_error))
        w->result->initial_backward_error =
            fmax(w->result->initial_backward_error, w->linear.initial_backward_error);
    if (isfinite(w->linear.max_backward_error))
        w->result->max_backward_error =
            fmax(w->result->max_backward_error, w->linear.max_backward_error);
    if (failed && w->linear.update_count && !factor_basis(w)) {
        w->result->basis_update_rejections++;
        failed = linear_solve(w, rhs, x, transpose);
    }
    free(saved);
    return failed;
}
static int transpose_solve(Work *w, const double *rhs, double *x) {
    if (linear_solve(w, rhs, x, 1))
        return -1;
    for (size_t j = 0; j < w->m; j++) {
        long double magnitude = fabs(rhs[j]);
        long double sum = dot(w, w->basis[j], x, &magnitude);
        if (!isfinite(sum) || !isfinite(magnitude) ||
            fabsl(sum - rhs[j]) > 1e-10L * (1 + magnitude))
            return w->linear.update_count && !factor_basis(w) ? transpose_solve(w, rhs, x) : -1;
    }
    return 0;
}
static void message(Work *w, const char *status, const char *text) {
    w->result->status = status;
    snprintf(w->result->message, sizeof w->result->message, "%s", text);
}
static int budget(Work *w) {
    if (lab_now() - w->started >= w->options->time_limit) {
        message(w, "time_limit", "LP time budget exhausted");
        return -1;
    }
    if (w->result->iterations >= w->options->iteration_limit) {
        message(w, "iteration_limit", "LP pivot budget exhausted");
        return -1;
    }
    return 0;
}
static double original_basic_tolerance(const Work *w, size_t position) {
    size_t j = w->basis[position];
    double tolerance = LAB_PRIMAL_TOLERANCE;
    if (j < w->structural)
        tolerance = fmin(tolerance, 1e-7 / fabs(w->sign[j]));
    else if (!w->artificial[j] && w->a[j].count == 1)
        tolerance = fmin(tolerance, 1e-7 * w->equation_scale[w->a[j].row[0]]);
    return tolerance;
}
static int factor_basis(Work *w) {
    memset(w->basic, 0, w->n);
    for (size_t i = 0; i < w->m; i++) {
        if (w->basis[i] >= w->n || w->basic[w->basis[i]]) {
            message(w, "numerical_failure", "Invalid or duplicated basis column");
            return -1;
        }
        w->basic[w->basis[i]] = 1;
    }
    w->result->refactorizations++;
    int failed;
    {
        size_t nnz = 0;
        for (size_t j = 0; j < w->m; j++)
            nnz += w->a[w->basis[j]].count;
        size_t *start = lab_alloc(w->m + 1, sizeof(size_t));
        size_t *row = lab_alloc(nnz, sizeof(size_t));
        double *value = lab_alloc(nnz, sizeof(double));
        size_t at = 0;
        for (size_t j = 0; j < w->m; j++) {
            start[j] = at;
            const Column *c = &w->a[w->basis[j]];
            for (size_t k = 0; k < c->count; k++) {
                row[at] = c->row[k];
                value[at++] = c->value[k];
            }
        }
        start[w->m] = at;
        failed = lab_linear_factorize_ex(&w->linear, w->m, start, row, value, w->options->dense,
                                         !w->options->no_scaling, !w->options->no_refinement);
        free(start);
        free(row);
        free(value);
    }
    if (w->linear.factor.peak_nnz > w->result->factor_nnz_peak)
        w->result->factor_nnz_peak = w->linear.factor.peak_nnz;
    if (w->linear.factor.core_dimension > w->result->factor_core_peak)
        w->result->factor_core_peak = w->linear.factor.core_dimension;
    if (w->linear.factor.eliminated_singletons > w->result->singleton_eliminations_peak)
        w->result->singleton_eliminations_peak = w->linear.factor.eliminated_singletons;
    if (failed) {
        message(w, "numerical_failure", "Singular or ill-conditioned basis");
        return -1;
    }
    return 0;
}
static int basis_values(Work *w, int primal_required) {
    if (linear_solve(w, w->b, w->xb, 0)) {
        message(w, "numerical_failure", "Basis solve failed residual check");
        return -1;
    }
    /* Recompute residuals rather than trusting a successful factorization. */
    long double *activity = lab_alloc(w->m, sizeof(long double));
    long double *magnitude = lab_alloc(w->m, sizeof(long double));
    for (size_t j = 0; j < w->m; j++) {
        const Column *c = &w->a[w->basis[j]];
        for (size_t k = 0; k < c->count; k++) {
            long double term = (long double)c->value[k] * w->xb[j];
            activity[c->row[k]] += term;
            magnitude[c->row[k]] += fabsl(term);
        }
    }
    int bad = 0;
    for (size_t i = 0; i < w->m; i++)
        if (!isfinite(activity[i]) || !isfinite(magnitude[i]) ||
            fabsl(activity[i] - w->b[i]) > 1e-10L * (1 + magnitude[i] + fabs(w->b[i])))
            bad = 1;
    free(activity);
    free(magnitude);
    if (bad) {
        message(w, "numerical_failure", "Basis solve failed residual check");
        return -1;
    }
    for (size_t i = 0; i < w->m; i++) {
        if (primal_required && w->xb[i] < -LAB_PRIMAL_TOLERANCE) {
            message(w, "numerical_failure", "Basis lost primal feasibility");
            return -1;
        }
        if (primal_required && w->xb[i] < 0)
            w->xb[i] = 0;
    }
    return 0;
}
static int refactor(Work *w, int primal_required) {
    return factor_basis(w) || basis_values(w, primal_required) ? -1 : 0;
}
static int is_basic(const Work *w, size_t j) {
    return w->basic[j] != 0;
}
static void change_basis(Work *w, size_t position, size_t column) {
    w->basic[w->basis[position]] = 0;
    w->basis[position] = column;
    w->basic[column] = 1;
}
/* The caller has changed its combinatorial basis. A rejected update leaves the
 * previous linear factor valid, so fallback builds the complete NEW basis. */
static int update_basis(Work *w, size_t position, const double *direction, int primal_required) {
    if (w->options->no_updates || !w->linear.valid || w->linear.n != w->m)
        return refactor(w, primal_required);
    const Column *c = &w->a[w->basis[position]];
    size_t before = w->linear.refinement_steps;
    int status = lab_linear_update(&w->linear, position, c->count, c->row, c->value, direction);
    w->result->refinement_steps += w->linear.refinement_steps - before;
    if (status == 0) {
        w->result->basis_updates++;
        if (!basis_values(w, primal_required))
            return 0;
    }
    w->result->basis_update_rejections++;
    return refactor(w, primal_required);
}
static void exact_coefficient(const LabModel *m, size_t k, mpq_t q) {
    if (m->exact)
        mpq_set(q, m->exact->value[k]);
    else
        mpq_set_d(q, m->value[k]);
}
static void exact_column_bound(const LabModel *m, size_t j, int upper, mpq_t q) {
    double current = upper ? m->col_upper[j] : m->col_lower[j];
    if (m->exact && current == (upper ? m->exact->col_upper_source[j]
                                      : m->exact->col_lower_source[j]))
        mpq_set(q, upper ? m->exact->col_upper[j] : m->exact->col_lower[j]);
    else
        mpq_set_d(q, current);
}
static int exact_source(const LabModel *m) {
    return m->exact || (!m->parsed_decimal && m->rows <= 128 && m->cols <= 256 &&
                       m->nnz <= 20000);
}
static void exact_shift(const Work *w, size_t j, mpq_t q) {
    const LabModel *m = w->model;
    int upper = w->transform[j] == 'U';
    if (w->transform[j] == 'F')
        mpq_set_ui(q, 0, 1);
    else if (m->exact && w->shift[j] == (upper ? m->exact->col_upper_source[j]
                                              : m->exact->col_lower_source[j]))
        mpq_set(q, upper ? m->exact->col_upper[j] : m->exact->col_lower[j]);
    else
        mpq_set_d(q, w->shift[j]);
}
static mpq_t *exact_shift_activity(const Work *w) {
    const LabModel *m = w->model;
    mpq_t *activity = lab_alloc(m->rows, sizeof(mpq_t));
    for (size_t i = 0; i < m->rows; i++)
        mpq_init(activity[i]);
    mpq_t shift, value;
    mpq_inits(shift, value, NULL);
    for (size_t j = 0; j < m->cols; j++) {
        exact_shift(w, j, shift);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            exact_coefficient(m, k, value);
            mpq_mul(value, value, shift);
            mpq_add(activity[m->row_index[k]], activity[m->row_index[k]], value);
        }
    }
    mpq_clears(shift, value, NULL);
    return activity;
}
static void exact_row_bound(const LabModel *m, size_t i, int upper, mpq_t q) {
    if (m->exact)
        mpq_set(q, upper ? m->exact->row_upper[i] : m->exact->row_lower[i]);
    else
        mpq_set_d(q, upper ? m->row_upper[i] : m->row_lower[i]);
}
static int exact_basis_value(const Work *w, size_t j, size_t entry, mpq_t value) {
    const LabModel *m = w->model;
    const Column *c = &w->a[j];
    size_t i = c->row[entry];
    if (j >= w->structural || w->source_row[i] == SIZE_MAX) {
        mpq_set_d(value, c->value[entry]);
        return 0;
    }
    size_t col = w->original[j], source = w->source_row[i];
    size_t lo = m->col_start[col], hi = m->col_start[col + 1];
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (m->row_index[mid] < source)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo == m->col_start[col + 1] || m->row_index[lo] != source)
        return -1;
    exact_coefficient(m, lo, value);
    mpq_t multiplier;
    mpq_init(multiplier);
    mpq_set_d(multiplier, w->sign[j]);
    mpq_mul(value, value, multiplier);
    mpq_set_d(multiplier, w->row_scale[i]);
    mpq_mul(value, value, multiplier);
    mpq_clear(multiplier);
    return 0;
}
static void exact_rhs(const Work *w, mpq_t *rhs) {
    const LabModel *m = w->model;
    mpq_t *activity = exact_shift_activity(w);
    mpq_t side, scale, shift;
    mpq_inits(side, scale, shift, NULL);
    for (size_t i = 0; i < w->m; i++)
        if (w->source_row[i] != SIZE_MAX) {
            size_t source = w->source_row[i];
            exact_row_bound(m, source, w->row_upper[i], side);
            mpq_sub(rhs[i], side, activity[source]);
            mpq_set_d(scale, w->row_scale[i]);
            mpq_mul(rhs[i], rhs[i], scale);
        }
    for (size_t j = 0; j < m->cols; j++) {
        exact_shift(w, j, shift);
        for (size_t i = 0; i < w->m; i++) {
            int upper = w->row_id[i] == w->upper_row[j];
            if (!upper && w->row_id[i] != w->lower_row[j])
                continue;
            exact_column_bound(m, j, upper, side);
            if (upper)
                mpq_sub(rhs[i], side, shift);
            else
                mpq_sub(rhs[i], shift, side);
            mpq_set_d(scale, w->equation_scale[i]);
            mpq_mul(rhs[i], rhs[i], scale);
        }
    }
    for (size_t i = 0; i < m->rows; i++)
        mpq_clear(activity[i]);
    free(activity);
    mpq_clears(side, scale, shift, NULL);
}
/* Reconstruct the original point before converting it to doubles. Singleton
 * elimination removes the many trivial bound/slack rows before the bounded
 * dense exact core. A genuinely negative basic value requires a basis repair. */
static int exact_primal(Work *w, double *point) {
    const LabModel *m = w->model;
    if (w->options->no_exact || !exact_source(m))
        return 0;
    w->result->exact_primal_attempts++;
    size_t n = w->m, nnz = 0;
    for (size_t t = 0; t < n; t++)
        nnz += w->a[w->basis[t]].count;
    if (n > 1024 || nnz > 100000)
        return -1;
    double deadline = fmin(w->started + w->options->time_limit, lab_now() + 1);
    size_t *start = lab_alloc(n + 1, sizeof(size_t));
    size_t *row = lab_alloc(nnz, sizeof(size_t));
    mpq_t *value = lab_alloc(nnz, sizeof(mpq_t));
    mpq_t *rhs = lab_alloc(n, sizeof(mpq_t)), *z = lab_alloc(n, sizeof(mpq_t));
    mpq_t *x = lab_alloc(m->cols, sizeof(mpq_t));
    mpq_t *activity = lab_alloc(m->rows, sizeof(mpq_t));
    for (size_t k = 0; k < nnz; k++)
        mpq_init(value[k]);
    for (size_t i = 0; i < n; i++) {
        mpq_init(rhs[i]);
        mpq_init(z[i]);
    }
    for (size_t j = 0; j < m->cols; j++)
        mpq_init(x[j]);
    for (size_t i = 0; i < m->rows; i++)
        mpq_init(activity[i]);
    size_t at = 0;
    int valid = 1;
    for (size_t t = 0; t < n; t++) {
        start[t] = at;
        size_t j = w->basis[t];
        for (size_t k = 0; k < w->a[j].count; k++) {
            row[at] = w->a[j].row[k];
            if (exact_basis_value(w, j, k, value[at]))
                valid = 0;
            at++;
        }
    }
    start[n] = at;
    exact_rhs(w, rhs);
    if (valid && lab_rational_sparse_solve(n, start, row, value, rhs, z, deadline))
        valid = 0;
    for (size_t i = 0; i < n && valid; i++)
        if (w->artificial[w->basis[i]] || mpq_sgn(z[i]) < 0)
            valid = 0;
    mpq_t term, bound;
    mpq_inits(term, bound, NULL);
    for (size_t j = 0; j < m->cols; j++)
        exact_shift(w, j, x[j]);
    for (size_t i = 0; i < n && valid; i++) {
        size_t j = w->basis[i];
        if (j < w->structural) {
            mpq_set_d(term, w->sign[j]);
            mpq_mul(term, term, z[i]);
            mpq_add(x[w->original[j]], x[w->original[j]], term);
        }
    }
    for (size_t j = 0; j < m->cols && valid; j++) {
        if (isfinite(m->col_lower[j])) {
            exact_column_bound(m, j, 0, bound);
            if (mpq_cmp(x[j], bound) < 0)
                valid = 0;
        }
        if (isfinite(m->col_upper[j])) {
            exact_column_bound(m, j, 1, bound);
            if (mpq_cmp(x[j], bound) > 0)
                valid = 0;
        }
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            exact_coefficient(m, k, term);
            mpq_mul(term, term, x[j]);
            mpq_add(activity[m->row_index[k]], activity[m->row_index[k]], term);
        }
    }
    for (size_t i = 0; i < m->rows && valid; i++) {
        if (isfinite(m->row_lower[i])) {
            exact_row_bound(m, i, 0, bound);
            if (mpq_cmp(activity[i], bound) < 0)
                valid = 0;
        }
        if (isfinite(m->row_upper[i])) {
            exact_row_bound(m, i, 1, bound);
            if (mpq_cmp(activity[i], bound) > 0)
                valid = 0;
        }
    }
    if (lab_now() >= deadline)
        valid = 0;
    double *candidate = lab_alloc(m->cols, sizeof(double));
    for (size_t j = 0; j < m->cols && valid; j++) {
        if (lab_exact_nearest(x[j], &candidate[j]) ||
            (mpq_sgn(x[j]) && candidate[j] == 0))
            valid = 0;
    }
    if (valid) {
        memcpy(point, candidate, m->cols * sizeof(double));
        w->result->exact_primal_recoveries++;
    }
    free(candidate);
    for (size_t k = 0; k < nnz; k++)
        mpq_clear(value[k]);
    for (size_t i = 0; i < n; i++) {
        mpq_clear(rhs[i]);
        mpq_clear(z[i]);
    }
    for (size_t j = 0; j < m->cols; j++)
        mpq_clear(x[j]);
    for (size_t i = 0; i < m->rows; i++)
        mpq_clear(activity[i]);
    free(start); free(row); free(value); free(rhs); free(z); free(x); free(activity);
    mpq_clears(term, bound, NULL);
    return valid ? 1 : -1;
}
/* Reconstruct a small Phase I basis over original decimal coefficients, then
 * validate the resulting row combination against the original model and box.
 * Stored standard-form RHS rounding is deliberately absent from this proof.
 * An unavailable/limited reconstruction cannot bypass an uncertain sign. */
static int phase_one_certificate(Work *w) {
    const LabModel *m = w->model;
    size_t n = w->m;
    if (w->options->no_exact || n > LAB_RATIONAL_MAX_DIMENSION ||
        (m->parsed_decimal && !m->exact))
        return 0;
    w->result->exact_attempts++;
    double deadline = fmin(w->started + w->options->time_limit, lab_now() + 1);
    mpq_t *a = lab_alloc(n * n, sizeof(mpq_t));
    mpq_t *rhs = lab_alloc(n, sizeof(mpq_t));
    mpq_t *y = lab_alloc(m->rows, sizeof(mpq_t));
    for (size_t k = 0; k < n * n; k++)
        mpq_init(a[k]);
    for (size_t i = 0; i < n; i++)
        mpq_init(rhs[i]);
    for (size_t i = 0; i < m->rows; i++)
        mpq_init(y[i]);
    mpq_t value, multiplier, margin, q, side;
    mpq_inits(value, multiplier, margin, q, side, NULL);
    int valid = 1;
    for (size_t t = 0; t < n; t++) {
        size_t j = w->basis[t];
        mpq_set_ui(rhs[t], w->artificial[j], 1);
        const Column *c = &w->a[j];
        for (size_t k = 0; k < c->count; k++) {
            size_t i = c->row[k];
            if (j < w->structural && w->source_row[i] != SIZE_MAX) {
                size_t col = w->original[j], source = w->source_row[i];
                size_t lo = m->col_start[col], hi = m->col_start[col + 1];
                while (lo < hi) {
                    size_t mid = lo + (hi - lo) / 2;
                    if (m->row_index[mid] < source)
                        lo = mid + 1;
                    else
                        hi = mid;
                }
                if (lo == m->col_start[col + 1] || m->row_index[lo] != source) {
                    valid = 0;
                    break;
                }
                exact_coefficient(m, lo, value);
                mpq_set_d(multiplier, w->sign[j]);
                mpq_mul(value, value, multiplier);
                mpq_set_d(multiplier, w->row_scale[i]);
                mpq_mul(value, value, multiplier);
            } else
                mpq_set_d(value, c->value[k]);
            mpq_add(a[t * n + i], a[t * n + i], value);
        }
        if (!valid)
            break;
    }
    if (valid && lab_rational_solve(n, a, rhs, deadline))
        valid = 0;
    for (size_t i = 0; i < n && valid; i++)
        if (w->source_row[i] != SIZE_MAX) {
            mpq_set_d(multiplier, w->row_scale[i]);
            mpq_mul(value, multiplier, rhs[i]);
            mpq_add(y[w->source_row[i]], y[w->source_row[i]], value);
        }
    for (size_t i = 0; i < m->rows && valid; i++) {
        int sign = mpq_sgn(y[i]);
        if (!sign)
            continue;
        double bound = sign > 0 ? m->row_lower[i] : m->row_upper[i];
        if (!isfinite(bound)) {
            valid = 0;
            break;
        }
        if (m->exact)
            mpq_set(side, sign > 0 ? m->exact->row_lower[i] : m->exact->row_upper[i]);
        else
            mpq_set_d(side, bound);
        mpq_mul(value, y[i], side);
        mpq_add(margin, margin, value);
    }
    for (size_t j = 0; j < m->cols && valid; j++) {
        mpq_set_ui(q, 0, 1);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            exact_coefficient(m, k, value);
            mpq_mul(value, value, y[m->row_index[k]]);
            mpq_add(q, q, value);
        }
        int sign = mpq_sgn(q);
        if (!sign)
            continue;
        double bound = sign > 0 ? m->col_upper[j] : m->col_lower[j];
        if (!isfinite(bound)) {
            valid = 0;
            break;
        }
        exact_column_bound(m, j, sign > 0, side);
        mpq_mul(value, q, side);
        mpq_sub(margin, margin, value);
        if (lab_now() >= deadline)
            valid = 0;
    }
    valid = valid && mpq_sgn(margin) > 0 && lab_now() < deadline;
    if (valid) {
        LabLpResult *r = w->result;
        r->farkas_rows = m->rows;
        r->farkas = lab_alloc(m->rows, sizeof(char *));
        for (size_t i = 0; i < m->rows; i++)
            r->farkas[i] = mpq_get_str(NULL, 10, y[i]);
        r->farkas_margin = mpq_get_str(NULL, 10, margin);
        r->exact_infeasibilities++;
    }
    for (size_t k = 0; k < n * n; k++)
        mpq_clear(a[k]);
    for (size_t i = 0; i < n; i++)
        mpq_clear(rhs[i]);
    for (size_t i = 0; i < m->rows; i++)
        mpq_clear(y[i]);
    free(a);
    free(rhs);
    free(y);
    mpq_clears(value, multiplier, margin, q, side, NULL);
    return valid;
}
/* 0 optimum, 1 unbounded, -1 explicit failure/limit. Bland's normal leaving
 * order is preserved. A rejected basis can retry a validated Harris pivot. */
static int simplex(Work *w, int phase_one) {
    int factor_ready = 0;
    for (;;) {
        if (budget(w) || (!factor_ready && refactor(w, 1)))
            return -1;
        factor_ready = 0;
        for (size_t i = 0; i < w->m; i++)
            w->cb[i] = phase_one ? (double)w->artificial[w->basis[i]] : w->cost[w->basis[i]];
        if (transpose_solve(w, w->cb, w->pi)) {
            message(w, "numerical_failure", "Transpose basis solve failed");
            return -1;
        }
        size_t entering = SIZE_MAX;
        double worst = 0;
        long double pi_norm = 0;
        for (size_t i = 0; i < w->m; i++)
            pi_norm += fabs(w->pi[i]);
        for (size_t j = 0; j < w->n; j++) {
            if ((!phase_one && w->artificial[j]) || is_basic(w, j))
                continue;
            long double reduced = phase_one ? (double)w->artificial[j] : w->cost[j];
            /* Scale solve uncertainty by this column, not a dimensionless
             * objective floor that can suppress a real tiny coefficient. */
            long double direct_cost = reduced, dot_magnitude = 0;
            reduced -= dot(w, j, w->pi, &dot_magnitude);
            long double scale = fabsl(direct_cost) + dot_magnitude;
            long double tolerance = fmaxl(1e-10L * fmaxl(scale, DBL_MIN),
                                          64 * DBL_EPSILON * pi_norm * column_norm(w, j));
            if (!isfinite((double)reduced) || !isfinite((double)scale) || !isfinite(tolerance)) {
                message(w, "numerical_failure", "Reduced cost overflow");
                return -1;
            }
            worst = fmax(worst, (double)-reduced);
            /* An explicit cost with no dual contribution has no cancellation
             * to discount. Ignoring it can hide an improving unbounded ray. */
            if (entering == SIZE_MAX &&
                (reduced < -tolerance || (direct_cost < 0 && dot_magnitude == 0)))
                entering = j;
        }
        if (!phase_one)
            w->result->dual_residual = worst;
        if (entering == SIZE_MAX) {
            if (phase_one) {
                double artificial_value = 0;
                for (size_t i = 0; i < w->m; i++)
                    if (w->artificial[w->basis[i]])
                        artificial_value += w->xb[i];
                if (artificial_value > 1e-8) {
                    if (phase_one_certificate(w)) {
                        message(w, "infeasible", "Exact original-model Farkas contradiction reconstructed");
                        return -1;
                    }
                    if (budget(w))
                        return -1;
                    if (worst > 0) {
                        message(w, "numerical_failure",
                                "Uncertain Phase I dual sign; infeasibility unproved");
                        return -1;
                    }
                }
            }
            return 0;
        }
        scatter(w, entering, w->column);
        if (linear_solve(w, w->column, w->direction, 0)) {
            message(w, "numerical_failure", "Pivot direction solve failed");
            return -1;
        }
        size_t leaving = SIZE_MAX;
        int ratio_status = lab_primal_ratio(w->m, w->xb, w->direction, w->basis,
                                            1, &leaving);
        if (ratio_status) {
            if (phase_one || ratio_status < 0) {
                message(w, "numerical_failure", "No stable leaving pivot");
                return -1;
            }
            return 1;
        }
        size_t alternative = SIZE_MAX;
        if (!w->options->bland_ratio &&
            lab_primal_ratio(w->m, w->xb, w->direction, w->basis, 0, &alternative))
            alternative = SIZE_MAX;
        size_t old_column = w->basis[leaving];
        const char *prior_status = w->result->status;
        char prior_message[sizeof w->result->message];
        memcpy(prior_message, w->result->message, sizeof prior_message);
        change_basis(w, leaving, entering);
        w->result->iterations++;
        if (!w->options->bland_ratio) {
            /* A ratio prediction cannot establish actual basis feasibility.
             * Check the newly factored basis before accepting it. Retain that
             * factor for the next iteration; the ordinary path factors once. */
            if (update_basis(w, leaving, w->direction, 1)) {
                w->result->rejected_pivots++;
                change_basis(w, leaving, old_column);
                if (alternative == SIZE_MAX || alternative == leaving || budget(w))
                    return -1;
                change_basis(w, alternative, entering);
                w->result->iterations++;
                if (refactor(w, 1)) {
                    w->result->rejected_pivots++;
                    return -1;
                }
                w->result->harris_rescues++;
                w->result->status = prior_status;
                memcpy(w->result->message, prior_message, sizeof prior_message);
            }
            factor_ready = 1;
        } else if (!w->options->no_updates) {
            if (update_basis(w, leaving, w->direction, 1))
                return -1;
            factor_ready = 1;
        }
    }
}
static void remove_row(Work *w, size_t row, size_t basis_position) {
    w->basic[w->basis[basis_position]] = 0;
    for (size_t j = 0; j < w->n; j++) {
        Column *c = &w->a[j];
        size_t keep = 0;
        double norm = 0;
        for (size_t k = 0; k < c->count; k++) {
            if (c->row[k] == row)
                continue;
            c->row[keep] = c->row[k] - (size_t)(c->row[k] > row);
            c->value[keep++] = c->value[k];
            norm = fmax(norm, fabs(c->value[keep - 1]));
        }
        c->count = keep;
        if (j < w->structural)
            w->structural_norm[j] = norm;
    }
    for (size_t i = row; i + 1 < w->m; i++) {
        w->b[i] = w->b[i + 1];
        w->row_id[i] = w->row_id[i + 1];
        w->source_row[i] = w->source_row[i + 1];
        w->row_scale[i] = w->row_scale[i + 1];
        w->equation_scale[i] = w->equation_scale[i + 1];
        w->row_upper[i] = w->row_upper[i + 1];
    }
    for (size_t i = basis_position; i + 1 < w->m; i++)
        w->basis[i] = w->basis[i + 1];
    for (size_t j = 0; j < w->n; j++)
        if (w->artificial[j] && w->art_row[j] > row)
            w->art_row[j]--;
    w->m--;
    w->result->redundant_rows++;
}
static int remove_artificials(Work *w) {
    for (size_t i = 0; i < w->m;) {
        if (!w->artificial[w->basis[i]]) {
            i++;
            continue;
        }
        if (budget(w) || refactor(w, 1))
            return -1;
        if (fabs(w->xb[i]) > 1e-9) {
            message(w, "numerical_failure", "Nonzero artificial variable after Phase I");
            return -1;
        }
        memset(w->cb, 0, w->m * sizeof(double));
        w->cb[i] = 1;
        if (transpose_solve(w, w->cb, w->pi))
            return -1;
        long double inverse_norm = 0;
        for (size_t k = 0; k < w->m; k++)
            inverse_norm += fabsl(w->pi[k]);
        size_t entering = SIZE_MAX;
        int ambiguous = 0;
        for (size_t j = 0; j < w->n; j++) {
            if (w->artificial[j] || is_basic(w, j))
                continue;
            long double v = 0, magnitude = 0, column_norm = 0;
            const Column *c = &w->a[j];
            for (size_t k = 0; k < c->count; k++) {
                long double term = (long double)w->pi[c->row[k]] * c->value[k];
                v += term;
                magnitude += fabsl(term);
                column_norm = fmaxl(column_norm, fabs(c->value[k]));
            }
            /* Include solve uncertainty, including components whose true
             * value is zero but whose computed multiplier is roundoff-sized. */
            long double uncertainty = fmaxl(magnitude, inverse_norm * column_norm);
            if (fabsl(v) > 1e-12L * fmaxl(uncertainty, DBL_MIN)) {
                entering = j;
                break;
            }
            if (v != 0)
                ambiguous = 1;
        }
        if (entering != SIZE_MAX) {
            change_basis(w, i, entering);
            w->result->iterations++;
            i++;
        } else if (ambiguous) {
            message(w, "numerical_failure", "Cannot establish redundant row with stable pivots");
            return -1;
        } else {
            /* The artificial column is a unit vector in its source equation.
             * Deleting that equation and this basic column preserves rank. */
            remove_row(w, w->art_row[w->basis[i]], i);
        }
    }
    return 0;
}
/* All chosen changes of units are exact powers of two. Exponents are clipped
 * toward zero to retain finite, reversible coefficients; scaling never drops
 * a nonzero or changes the caller's original model. */
static int limit_exponent(double value, int exponent) {
    if (value == 0)
        return exponent;
    int q = ilogb(fabs(value));
    if (exponent > 0 && q + exponent > DBL_MAX_EXP - 2)
        exponent = DBL_MAX_EXP - 2 - q > 0 ? DBL_MAX_EXP - 2 - q : 0;
    if (exponent < 0 && q + exponent < DBL_MIN_EXP - 1) {
        int minimum = DBL_MIN_EXP - 1 - q;
        exponent = minimum < 0 ? minimum : 0;
    }
    return exponent;
}
static int bounded_exponent(int exponent) {
    return exponent < -256 ? -256 : exponent > 256 ? 256 : exponent;
}
static int scale_form(Work *w, const LabModel *m) {
    if (w->options->no_scaling)
        return 0;
    double *norm = lab_alloc(w->m, sizeof(double));
    int *column_exp = lab_alloc(m->cols, sizeof(int));
    int *row_exp = lab_alloc(w->m, sizeof(int));
    for (size_t j = 0; j < w->structural; j++)
        for (size_t k = 0; k < w->a[j].count; k++) {
            size_t i = w->a[j].row[k];
            norm[i] = fmax(norm[i], fabs(w->a[j].value[k]));
        }
    for (size_t j = 0; j < m->cols; j++) {
        if (m->col_lower[j] == m->col_upper[j])
            continue;
        double width = m->col_upper[j] - m->col_lower[j];
        if (isfinite(width) && width > 0)
            column_exp[j] = bounded_exponent(ilogb(width));
    }
    for (size_t j = 0; j < w->structural; j++) {
        size_t original = w->original[j];
        if (!isfinite(m->col_lower[original]) || !isfinite(m->col_upper[original])) {
            double largest = 0;
            for (size_t k = 0; k < w->a[j].count; k++) {
                size_t i = w->a[j].row[k];
                if (w->source_row[i] != SIZE_MAX && norm[i] > 0)
                    largest = fmax(largest, fabs(w->a[j].value[k]) / norm[i]);
            }
            if (largest > 0)
                column_exp[original] = bounded_exponent(-ilogb(largest));
        }
    }
    /* One original free variable uses the same scale in both split columns. */
    for (size_t j = 0; j < w->structural; j++) {
        int e = column_exp[w->original[j]];
        e = limit_exponent(w->cost[j], e);
        for (size_t k = 0; k < w->a[j].count; k++)
            e = limit_exponent(w->a[j].value[k], e);
        column_exp[w->original[j]] = e;
    }
    int bad = 0;
    memset(norm, 0, w->m * sizeof(double));
    for (size_t j = 0; j < w->structural; j++) {
        int e = column_exp[w->original[j]];
        if (e)
            w->scaled_columns++;
        w->sign[j] = scalbn(w->sign[j], e);
        double cost = scalbn(w->cost[j], e);
        if (!isfinite(cost) || scalbn(cost, -e) != w->cost[j])
            bad = 1;
        w->cost[j] = cost;
        for (size_t k = 0; k < w->a[j].count; k++) {
            double old = w->a[j].value[k], value = scalbn(old, e);
            if (!isfinite(value) || scalbn(value, -e) != old)
                bad = 1;
            w->a[j].value[k] = value;
            norm[w->a[j].row[k]] = fmax(norm[w->a[j].row[k]], fabs(value));
        }
    }
    if (bad) {
        free(norm);
        free(column_exp);
        free(row_exp);
        message(w, "numerical_failure", "Irreversible or nonfinite column transformation");
        return -1;
    }
    for (size_t i = 0; i < w->m; i++) {
        int e = norm[i] > 0 ? bounded_exponent(-ilogb(norm[i])) : 0;
        row_exp[i] = limit_exponent(w->b[i], e);
    }
    for (size_t j = 0; j < w->structural; j++)
        for (size_t k = 0; k < w->a[j].count; k++) {
            size_t i = w->a[j].row[k];
            row_exp[i] = limit_exponent(w->a[j].value[k], row_exp[i]);
        }
    for (size_t j = 0; j < w->structural; j++)
        for (size_t k = 0; k < w->a[j].count; k++) {
            size_t i = w->a[j].row[k];
            double old = w->a[j].value[k], value = scalbn(old, row_exp[i]);
            if (!isfinite(value) || scalbn(value, -row_exp[i]) != old)
                bad = 1;
            w->a[j].value[k] = value;
        }
    for (size_t i = 0; i < w->m; i++) {
        int e = row_exp[i];
        if (e)
            w->scaled_rows++;
        double rhs = scalbn(w->b[i], e);
        if (!isfinite(rhs) || scalbn(rhs, -e) != w->b[i])
            bad = 1;
        w->b[i] = rhs;
        w->equation_scale[i] = scalbn(1, e);
        w->row_scale[i] *= w->equation_scale[i];
    }
    free(norm);
    free(column_exp);
    free(row_exp);
    if (bad) {
        message(w, "numerical_failure", "Irreversible or nonfinite scaling transformation");
        return -1;
    }
    return 0;
}
static int build(Work *w, const LabModel *m) {
    size_t dim = w->options->max_dimension, nz = 0, nr = 0;
    if (dim > LAB_LP_SPARSE_DIMENSION_LIMIT ||
        (w->options->dense && dim > LAB_LP_DENSE_DIMENSION_LIMIT)) {
        message(w, "unsupported", "Backend dimension ceiling exceeded");
        return -1;
    }
    for (size_t j = 0; j < m->cols; j++) {
        if (m->col_lower[j] > m->col_upper[j]) {
            message(w, "infeasible", "Conflicting column bounds");
            return -1;
        }
        if (m->col_lower[j] == m->col_upper[j])
            continue;
        nz += (!isfinite(m->col_lower[j]) && !isfinite(m->col_upper[j])) ? 2U : 1U;
        nr += (size_t)isfinite(m->col_lower[j]) + (size_t)isfinite(m->col_upper[j]);
    }
    for (size_t i = 0; i < m->rows; i++) {
        if (m->row_lower[i] > m->row_upper[i]) {
            message(w, "infeasible", "Conflicting row bounds");
            return -1;
        }
        nr += m->row_lower[i] == m->row_upper[i]
                  ? 1U
                  : (size_t)isfinite(m->row_lower[i]) + (size_t)isfinite(m->row_upper[i]);
    }
    if (nr > dim || nz > dim || nr > (dim - nz) / 2) {
        message(w, "unsupported", "Standard-form dimension limit exceeded");
        return -1;
    }
    if (m->nnz > (128000000U - 3 * nr) / 4) {
        message(w, "unsupported", "Sparse transformation entry budget exceeded");
        return -1;
    }
    w->structural = nz;
    w->m = w->total_rows = nr;
    w->n = nz + 2 * nr;
    size_t n = w->n, cap = nr ? nr : 1;
    w->a = lab_alloc(n, sizeof(Column));
    w->b = lab_alloc(cap, sizeof(double));
    w->cost = lab_alloc(n, sizeof(double));
    w->basis = lab_alloc(cap, sizeof(size_t));
    w->source_row = lab_alloc(cap, sizeof(size_t));
    w->row_scale = lab_alloc(cap, sizeof(double));
    w->equation_scale = lab_alloc(cap, sizeof(double));
    w->row_upper = lab_alloc(cap, 1);
    for (size_t i = 0; i < nr; i++)
        w->equation_scale[i] = 1;
    w->row_id = lab_alloc(cap, sizeof(size_t));
    w->artificial = lab_alloc(n, 1);
    w->basic = lab_alloc(n, 1);
    w->art_row = lab_alloc(n, sizeof(size_t));
    w->shift = lab_alloc(m->cols, sizeof(double));
    w->original = lab_alloc(nz, sizeof(size_t));
    w->sign = lab_alloc(nz, sizeof(double));
    w->lower_row = lab_alloc(m->cols, sizeof(size_t));
    w->upper_row = lab_alloc(m->cols, sizeof(size_t));
    w->transform = lab_alloc(m->cols, 1);
    size_t *first = lab_alloc(m->cols, sizeof(size_t));
    size_t *second = lab_alloc(m->cols, sizeof(size_t));
    size_t *lower = lab_alloc(m->rows, sizeof(size_t));
    size_t *upper = lab_alloc(m->rows, sizeof(size_t));
    long double *shifted = lab_alloc(m->rows, sizeof(long double));
    size_t z = 0;
    double sense = m->maximize ? -1 : 1;
    for (size_t j = 0; j < m->cols; j++) {
        first[j] = second[j] = w->lower_row[j] = w->upper_row[j] = SIZE_MAX;
        if (m->col_lower[j] == m->col_upper[j]) {
            w->shift[j] = m->col_lower[j];
            w->transform[j] = 'X';
        } else {
            int lf = isfinite(m->col_lower[j]), hf = isfinite(m->col_upper[j]);
            w->transform[j] = lf ? 'L' : hf ? 'U' : 'F';
            w->shift[j] = lf ? m->col_lower[j] : hf ? m->col_upper[j] : 0;
            first[j] = z;
            w->original[z] = j;
            w->sign[z] = lf || !hf ? 1 : -1;
            w->cost[z] = sense * m->objective[j] * w->sign[z];
            z++;
            if (!lf && !hf) {
                second[j] = z;
                w->original[z] = j;
                w->sign[z] = -1;
                w->cost[z++] = -sense * m->objective[j];
            }
        }
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++)
            shifted[m->row_index[k]] += (long double)m->value[k] * w->shift[j];
    }
    mpq_t *precise_shift = exact_source(m) ? exact_shift_activity(w) : NULL;
    mpq_t precise_rhs, precise_bound;
    mpq_inits(precise_rhs, precise_bound, NULL);
    size_t r = 0;
    for (size_t i = 0; i < m->rows; i++) {
        lower[i] = upper[i] = SIZE_MAX;
        int equality = m->row_lower[i] == m->row_upper[i];
        for (int side = 0; side < 2; side++) {
            double bound = side ? m->row_upper[i] : m->row_lower[i];
            if (!isfinite(bound) || (equality && side))
                continue;
            if (side)
                upper[i] = r;
            else
                lower[i] = r;
            long double rhs = bound - shifted[i];
            if (precise_shift) {
                exact_row_bound(m, i, side, precise_rhs);
                mpq_sub(precise_rhs, precise_rhs, precise_shift[i]);
                rhs = mpq_get_d(precise_rhs);
            }
            double flip = rhs < 0 ? -1 : 1;
            if (precise_shift && mpq_sgn(precise_rhs) < 0)
                flip = -1;
            w->b[r] = (double)(flip * rhs);
            w->source_row[r] = i;
            w->row_upper[r] = (unsigned char)side;
            w->row_scale[r] = flip;
            w->row_id[r] = r;
            double slack = equality ? 0 : flip * (side ? 1 : -1);
            put(w, r, nz + r, slack);
            put(w, r, nz + nr + r, 1);
            w->artificial[nz + nr + r] = 1;
            w->art_row[nz + nr + r] = r;
            w->basis[r] = slack == 1 ? nz + r : nz + nr + r;
            r++;
        }
    }
    for (size_t j = 0; j < m->cols; j++) {
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            size_t i = m->row_index[k];
            size_t rows[2] = {lower[i], upper[i]};
            for (int side = 0; side < 2; side++) {
                if (rows[side] == SIZE_MAX)
                    continue;
                size_t ri = rows[side];
                double value = w->row_scale[ri] * m->value[k];
                if (first[j] != SIZE_MAX)
                    put(w, ri, first[j], value * w->sign[first[j]]);
                if (second[j] != SIZE_MAX)
                    put(w, ri, second[j], -value);
            }
        }
        if (first[j] == SIZE_MAX)
            continue;
        /* Keep a slack equation for the initially implicit bound as well.
         * Tightening it changes only RHS; the root transformation stays fixed. */
        if (isfinite(m->col_lower[j])) {
            w->lower_row[j] = r;
            w->b[r] = 0;
            if (precise_shift) {
                exact_shift(w, j, precise_rhs);
                exact_column_bound(m, j, 0, precise_bound);
                mpq_sub(precise_rhs, precise_rhs, precise_bound);
                w->b[r] = mpq_get_d(precise_rhs);
            }
            put(w, r, first[j], -1);
            put(w, r, nz + r, 1);
            put(w, r, nz + nr + r, 1);
            w->artificial[nz + nr + r] = 1;
            w->art_row[nz + nr + r] = r;
            w->basis[r] = nz + r;
            w->source_row[r] = SIZE_MAX;
            w->row_id[r] = r;
            r++;
        }
        if (isfinite(m->col_upper[j])) {
            w->upper_row[j] = r;
            w->b[r] = w->transform[j] == 'U' ? 0 : m->col_upper[j] - w->shift[j];
            if (precise_shift) {
                exact_column_bound(m, j, 1, precise_rhs);
                exact_shift(w, j, precise_bound);
                mpq_sub(precise_rhs, precise_rhs, precise_bound);
                w->b[r] = mpq_get_d(precise_rhs);
            }
            put(w, r, first[j], w->transform[j] == 'U' ? -1 : 1);
            put(w, r, nz + r, 1);
            put(w, r, nz + nr + r, 1);
            w->artificial[nz + nr + r] = 1;
            w->art_row[nz + nr + r] = r;
            w->basis[r] = nz + r;
            w->source_row[r] = SIZE_MAX;
            w->row_id[r] = r;
            r++;
        }
    }
    free(first);
    free(second);
    free(lower);
    free(upper);
    free(shifted);
    if (precise_shift) {
        for (size_t i = 0; i < m->rows; i++)
            mpq_clear(precise_shift[i]);
        free(precise_shift);
    }
    mpq_clears(precise_rhs, precise_bound, NULL);
    unsigned char *has_structural = lab_alloc(nr, 1);
    for (size_t j = 0; j < nz; j++)
        for (size_t k = 0; k < w->a[j].count; k++)
            has_structural[w->a[j].row[k]] = 1;
    for (size_t i = 0; i < nr; i++) {
        if (!isfinite(w->b[i])) {
            free(has_structural);
            message(w, "numerical_failure", "Bound transformation overflow");
            return -1;
        }
        if (w->artificial[w->basis[i]] && w->b[i] > 0) {
            if (!has_structural[i]) {
                free(has_structural);
                message(w, "infeasible", "Contradictory constant row");
                return -1;
            }
        }
    }
    free(has_structural);
    for (size_t j = 0; j < n; j++) {
        w->result->standard_nnz += w->a[j].count;
        for (size_t k = 0; k < w->a[j].count; k++)
            if (!isfinite(w->a[j].value[k])) {
                message(w, "numerical_failure", "Coefficient transformation overflow");
                return -1;
            }
    }
    if (scale_form(w, m))
        return -1;
    /* Coefficients remain fixed between row deletions. Cache the same maximum
     * used by pricing's uncertainty guard instead of rescanning every column
     * after every pivot. Only structural columns need an allocated cache. */
    w->structural_norm = lab_alloc(nz, sizeof(double));
    for (size_t j = 0; j < nz; j++)
        for (size_t k = 0; k < w->a[j].count; k++)
            w->structural_norm[j] = fmax(w->structural_norm[j], fabs(w->a[j].value[k]));
    w->xb = lab_alloc(cap, sizeof(double));
    w->pi = lab_alloc(cap, sizeof(double));
    w->cb = lab_alloc(cap, sizeof(double));
    w->direction = lab_alloc(cap, sizeof(double));
    w->column = lab_alloc(cap, sizeof(double));
    w->result->standard_rows = nr;
    w->result->standard_columns = n;
    return 0;
}
static int load_basis(Work *w) {
    if (!w->options->basis_in)
        return 0;
    FILE *f = fopen(w->options->basis_in, "r");
    char magic[32], sha[128];
    size_t m = 0, n = 0, total = 0;
    if (!f || fscanf(f, "%31s %127s %zu %zu %zu", magic, sha, &m, &n, &total) != 5 ||
        strcmp(magic, "MILPLAB_BASIS_1") || strcmp(sha, w->options->model_sha256) || n != w->n ||
        total != w->total_rows || m != w->m) {
        if (f)
            fclose(f);
        message(w, "input_error", "Invalid or incompatible basis header");
        return -1;
    }
    size_t *ids = lab_alloc(m ? m : 1, sizeof(size_t));
    size_t *basis = lab_alloc(m ? m : 1, sizeof(size_t));
    int bad = 0;
    for (size_t i = 0; i < m; i++) {
        if (fscanf(f, "%zu %zu", &ids[i], &basis[i]) != 2 || ids[i] != w->row_id[i] ||
            (i && ids[i] <= ids[i - 1]) || basis[i] >= n || w->artificial[basis[i]]) {
            bad = 1;
            break;
        }
        for (size_t j = 0; j < i; j++)
            if (basis[i] == basis[j])
                bad = 1;
    }
    char extra[2];
    if (fscanf(f, "%1s", extra) == 1 || ferror(f))
        bad = 1;
    fclose(f);
    if (!bad)
        memcpy(w->basis, basis, m * sizeof(size_t));
    free(ids);
    free(basis);
    if (bad) {
        message(w, "input_error", "Invalid basis entries");
        return -1;
    }
    /* Phase I has independently established the active rows before loading.
     * The file cannot instruct the solver to drop any constraints. */
    return 1;
}
static int save_basis(Work *w) {
    if (!w->options->basis_out)
        return 0;
    FILE *f = fopen(w->options->basis_out, "wx");
    if (!f) {
        message(w, "input_error", "Cannot create fresh basis file");
        return -1;
    }
    fprintf(f, "MILPLAB_BASIS_1 %s %zu %zu %zu\n", w->options->model_sha256, w->m, w->n,
            w->total_rows);
    for (size_t i = 0; i < w->m; i++)
        fprintf(f, "%zu %zu\n", w->row_id[i], w->basis[i]);
    int bad = ferror(f);
    if (fclose(f))
        bad = 1;
    if (bad) {
        message(w, "input_error", "Cannot write basis file");
        return -1;
    }
    return 0;
}
/* Evaluate the existing support-gap safeguard in source arithmetic. This is
 * an optimizer rejection test, separate from the independent bound checker. */
static int exact_support_gap_bad(const Work *w, const double *primal, const double *dual,
                                 int *available) {
    const LabModel *m = w->model;
    mpq_t *activity = lab_alloc(m->rows, sizeof(mpq_t));
    for (size_t i = 0; i < m->rows; i++)
        mpq_init(activity[i]);
    mpq_t x, y, value, reduced, side, gap, variation, shift, threshold;
    mpq_inits(x, y, value, reduced, side, gap, variation, shift, threshold, NULL);
    for (size_t j = 0; j < m->cols; j++) {
        mpq_set_d(x, primal[j]);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            exact_coefficient(m, k, value);
            mpq_mul(value, value, x);
            mpq_add(activity[m->row_index[k]], activity[m->row_index[k]], value);
        }
        if (w->transform[j] != 'X') {
            exact_shift(w, j, shift);
            mpq_sub(value, x, shift);
            if (m->exact)
                mpq_set(reduced, m->exact->objective[j]);
            else
                mpq_set_d(reduced, m->objective[j]);
            mpq_mul(value, value, reduced);
            mpq_abs(value, value);
            mpq_add(variation, variation, value);
        }
    }
    *available = 1;
    for (size_t i = 0; i < m->rows && *available; i++) {
        int sign = (dual[i] > 0) - (dual[i] < 0);
        if (!sign)
            continue;
        if (!isfinite(sign > 0 ? m->row_lower[i] : m->row_upper[i])) {
            *available = 0;
            break;
        }
        exact_row_bound(m, i, sign < 0, side);
        mpq_sub(value, activity[i], side);
        mpq_set_d(y, dual[i]);
        mpq_mul(value, value, y);
        mpq_abs(value, value);
        mpq_add(gap, gap, value);
    }
    for (size_t j = 0; j < m->cols && *available; j++) {
        if (m->exact)
            mpq_set(reduced, m->exact->objective[j]);
        else
            mpq_set_d(reduced, m->objective[j]);
        if (m->maximize)
            mpq_neg(reduced, reduced);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            exact_coefficient(m, k, value);
            mpq_set_d(y, dual[m->row_index[k]]);
            mpq_mul(value, value, y);
            mpq_sub(reduced, reduced, value);
        }
        int sign = mpq_sgn(reduced);
        if (!sign)
            continue;
        if (!isfinite(sign > 0 ? m->col_lower[j] : m->col_upper[j])) {
            *available = 0;
            break;
        }
        exact_column_bound(m, j, sign < 0, side);
        mpq_set_d(x, primal[j]);
        mpq_sub(value, x, side);
        mpq_mul(value, value, reduced);
        mpq_abs(value, value);
        mpq_add(gap, gap, value);
    }
    mpq_set_ui(threshold, 1, 1);
    mpq_add(threshold, threshold, variation);
    mpq_set_ui(value, 1, 10000000);
    mpq_mul(threshold, threshold, value);
    int bad = *available && mpq_cmp(gap, threshold) > 0;
    for (size_t i = 0; i < m->rows; i++)
        mpq_clear(activity[i]);
    free(activity);
    mpq_clears(x, y, value, reduced, side, gap, variation, shift, threshold, NULL);
    return bad;
}
static int recover(Work *w, const LabModel *m) {
    LabLpResult *r = w->result;
    r->primal_repair_required = 0;
    r->primal = lab_alloc(m->cols, sizeof(double));
    r->dual = lab_alloc(m->rows, sizeof(double));
    memcpy(r->primal, w->shift, m->cols * sizeof(double));
    for (size_t i = 0; i < w->m; i++) {
        size_t j = w->basis[i];
        if (j < w->structural)
            r->primal[w->original[j]] += w->sign[j] * w->xb[i];
        if (w->source_row[i] != SIZE_MAX)
            r->dual[w->source_row[i]] += w->row_scale[i] * w->pi[i];
    }
    if (exact_primal(w, r->primal) < 0) {
        r->primal_repair_required = 1;
        if (lab_exact_evaluate(m, r->primal, &r->objective, &r->primal_residual))
            r->primal_residual = DBL_MAX;
        if (!budget(w))
            message(w, "numerical_failure", "Exact primal reconstruction unavailable or infeasible");
        return -1;
    }
    long double variation = 0;
    double residual = 0;
    long double *activity = lab_alloc(m->rows, sizeof(long double));
    for (size_t j = 0; j < m->cols; j++) {
        double x = r->primal[j];
        if (!isfinite(x)) {
            free(activity);
            message(w, "numerical_failure", "Primal reconstruction overflow");
            return -1;
        }
        variation += fabsl((long double)m->objective[j] * (x - w->shift[j]));
        residual = fmax(residual, m->col_lower[j] - x);
        residual = fmax(residual, x - m->col_upper[j]);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++)
            activity[m->row_index[k]] += (long double)m->value[k] * x;
    }
    for (size_t i = 0; i < m->rows; i++) {
        if (!isfinite(activity[i]) || !isfinite(r->dual[i])) {
            free(activity);
            message(w, "numerical_failure", "Row activity or multiplier overflow");
            return -1;
        }
        residual = fmax(residual, (double)(m->row_lower[i] - activity[i]));
        residual = fmax(residual, (double)(activity[i] - m->row_upper[i]));
    }
    double evaluated_objective, evaluated_residual;
    if (lab_exact_evaluate(m, r->primal, &evaluated_objective, &evaluated_residual)) {
        free(activity);
        r->primal_residual = DBL_MAX;
        message(w, "numerical_failure", "Original-model evaluation overflow or uncertainty");
        return -1;
    }
    r->objective = evaluated_objective;
    residual = evaluated_residual;
    r->primal_residual = residual;
    if (!isfinite(r->objective) || !isfinite(residual) || residual > 1e-7) {
        free(activity);
        message(w, "numerical_failure", "Original-model primal residual exceeds 1e-7");
        return -1;
    }
    /* A tolerance-small reduced cost can imply a large objective error over
     * wide bounds. This floating support calculation is a rejection safeguard,
     * not a certificate; milpcheck recomputes a rigorous bound independently. */
    long double sense = m->maximize ? -1 : 1;
    long double gap = 0;
    int available = 1;
    for (size_t i = 0; i < m->rows; i++) {
        double y = r->dual[i];
        if (y == 0)
            continue;
        double side = y > 0 ? m->row_lower[i] : m->row_upper[i];
        if (!isfinite(side)) {
            available = 0;
            break;
        }
        gap += fabsl((long double)y * (activity[i] - side));
    }
    for (size_t j = 0; j < m->cols && available; j++) {
        long double reduced = sense * m->objective[j];
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++)
            reduced -= (long double)m->value[k] * r->dual[m->row_index[k]];
        if (reduced == 0)
            continue;
        double side = reduced > 0 ? m->col_lower[j] : m->col_upper[j];
        if (!isfinite(side)) {
            available = 0;
            break;
        }
        gap += fabsl(reduced * (r->primal[j] - side));
    }
    free(activity);
    /* Complementarity differences avoid subtracting large nearly equal totals.
     * Fixed-variable contributions and objective offsets cannot hide a gap. */
    int gap_bad = exact_source(m) ? exact_support_gap_bad(w, r->primal, r->dual, &available)
                                 : available && (!isfinite(gap) || !isfinite(variation) ||
                                                 gap > 1e-7L * (1 + variation));
    if (gap_bad) {
        message(w, "numerical_failure", "Primal/support-bound gap exceeds stopping tolerance");
        return -1;
    }
    r->has_primal = 1;
    return 0;
}
void lab_lp_result_free(LabLpResult *result) {
    free(result->primal);
    free(result->dual);
    for (size_t i = 0; i < result->farkas_rows; i++)
        free(result->farkas[i]);
    free(result->farkas);
    free(result->farkas_margin);
    memset(result, 0, sizeof *result);
}
static void work_free(Work *w) {
    lab_linear_free(&w->linear);
    if (w->a)
        for (size_t j = 0; j < w->n; j++) {
            free(w->a[j].row);
            free(w->a[j].value);
        }
    free(w->a);
    free(w->b);
    free(w->cost);
    free(w->shift);
    free(w->sign);
    free(w->basis);
    free(w->original);
    free(w->row_id);
    free(w->art_row);
    free(w->source_row);
    free(w->row_scale);
    free(w->structural_norm);
    free(w->equation_scale);
    free(w->row_upper);
    free(w->artificial);
    free(w->basic);
    free(w->lower_row);
    free(w->upper_row);
    free(w->transform);
    free(w->xb);
    free(w->pi);
    free(w->cb);
    free(w->direction);
    free(w->column);
    memset(w, 0, sizeof *w);
}
static int dual_simplex(Work *w, int original_units);
static int restore_primal(Work *w, const LabModel *m) {
    if (recover(w, m) == 0)
        return 0;
    if (w->options->no_primal_repair ||
        (!w->result->primal_repair_required && !(w->result->primal_residual > 1e-7)))
        return -1;
    free(w->result->primal);
    free(w->result->dual);
    w->result->primal = w->result->dual = NULL;
    w->result->primal_repairs++;
    size_t before = w->result->dual_iterations;
    int status = dual_simplex(w, 1);
    w->result->primal_repair_iterations += w->result->dual_iterations - before;
    if (status)
        return -1;
    return recover(w, m);
}
static int cold_solve(Work *w, const LabModel *m) {
    w->model = m;
    if (build(w, m))
        return -1;
    if (simplex(w, 1))
        return -1;
    double violation = 0;
    for (size_t i = 0; i < w->m; i++)
        if (w->artificial[w->basis[i]])
            violation += w->xb[i];
    if (violation > 1e-8) {
        message(w, "infeasible", "Positive Phase I artificial objective (numerical claim)");
        return -1;
    }
    if (remove_artificials(w) || load_basis(w) < 0)
        return -1;
    int status = simplex(w, 0);
    if (status < 0 || restore_primal(w, m))
        return -1;
    if (status)
        message(w, "unbounded", "Improving LP ray found (numerical claim)");
    else {
        message(w, "optimal", "LP reduced costs satisfy numerical stopping tolerance");
        if (save_basis(w))
            return -1;
    }
    return 0;
}
/* Dual revised simplex repairs primal infeasibility caused by RHS changes while
 * preserving a dual-feasible basis. No new matrix or Phase I is needed. */
static int dual_simplex(Work *w, int original_units) {
    double *reduced = lab_alloc(w->n, sizeof(double));
    int outcome = -1;
    int factor_ready = !w->options->no_updates && w->linear.valid && w->linear.n == w->m &&
                       !basis_values(w, 0);
    for (;;) {
        if (budget(w) || (!factor_ready && refactor(w, 0)))
            break;
        factor_ready = 0;
        for (size_t i = 0; i < w->m; i++)
            w->cb[i] = w->cost[w->basis[i]];
        if (transpose_solve(w, w->cb, w->pi)) {
            message(w, "numerical_failure", "Warm transpose solve failed");
            break;
        }
        double worst = 0;
        int dual_bad = 0;
        long double pi_norm = 0;
        for (size_t i = 0; i < w->m; i++)
            pi_norm += fabs(w->pi[i]);
        for (size_t j = 0; j < w->n; j++) {
            if (w->artificial[j] || is_basic(w, j)) {
                reduced[j] = 0;
                continue;
            }
            long double scale = fabs(w->cost[j]);
            long double rc = w->cost[j] - dot(w, j, w->pi, &scale);
            long double tolerance = fmaxl(1e-10L * fmaxl(scale, DBL_MIN),
                                          64 * DBL_EPSILON * pi_norm * column_norm(w, j));
            if (!isfinite(rc) || !isfinite(scale) || !isfinite(tolerance) || rc < -tolerance)
                dual_bad = 1;
            reduced[j] = (double)fmaxl(rc, 0);
            worst = fmax(worst, (double)-rc);
        }
        w->result->dual_residual = worst;
        if (dual_bad) {
            message(w, "numerical_failure", "Warm basis lost dual feasibility");
            break;
        }
        size_t leaving = SIZE_MAX;
        for (size_t i = 0; i < w->m; i++)
            if (w->xb[i] < -(original_units ? original_basic_tolerance(w, i) : 1e-8) &&
                (leaving == SIZE_MAX || w->basis[i] < w->basis[leaving]))
                leaving = i;
        if (leaving == SIZE_MAX) {
            for (size_t i = 0; i < w->m; i++)
                if (w->xb[i] < 0)
                    w->xb[i] = 0;
            outcome = 0;
            break;
        }
        memset(w->cb, 0, w->m * sizeof(double));
        w->cb[leaving] = 1;
        if (transpose_solve(w, w->cb, w->column)) {
            message(w, "numerical_failure", "Warm pivot-row solve failed");
            break;
        }
        long double inverse_norm = 0;
        for (size_t i = 0; i < w->m; i++)
            inverse_norm += fabs(w->column[i]);
        if (!isfinite(inverse_norm)) {
            message(w, "numerical_failure", "Dual inverse-row norm overflow");
            break;
        }
        size_t entering = SIZE_MAX;
        long double ratio = INFINITY;
        int ambiguous = 0, pivot_bad = 0;
        for (size_t j = 0; j < w->n; j++) {
            if (w->artificial[j] || is_basic(w, j))
                continue;
            long double magnitude = 0, norm = 0;
            long double a = dot(w, j, w->column, &magnitude);
            for (size_t k = 0; k < w->a[j].count; k++)
                norm = fmaxl(norm, fabs(w->a[j].value[k]));
            long double tol = 1e-12L * fmaxl(fmaxl(magnitude, norm * inverse_norm), DBL_MIN);
            if (!isfinite(a) || !isfinite(magnitude) || !isfinite(norm) || !isfinite(tol)) {
                pivot_bad = 1;
                break;
            }
            if (a >= -tol) {
                if (a < 0)
                    ambiguous = 1;
                continue;
            }
            long double candidate = reduced[j] / -a;
            if (!isfinite(candidate)) {
                pivot_bad = 1;
                break;
            }
            long double tie = 1e-12L * (1 + fabsl(candidate));
            if (entering == SIZE_MAX || candidate < ratio - tie ||
                (fabsl(candidate - ratio) <= tie && j < entering)) {
                ratio = candidate;
                entering = j;
            }
        }
        if (pivot_bad) {
            message(w, "numerical_failure", "Dual ratio calculation overflow");
            break;
        }
        if (entering == SIZE_MAX) {
            if (ambiguous)
                message(w, "numerical_failure", "No stable dual ratio pivot");
            else {
                message(w, "infeasible", "Dual repair found no eligible column (numerical claim)");
                outcome = 1;
            }
            break;
        }
        change_basis(w, leaving, entering);
        w->result->iterations++;
        w->result->dual_iterations++;
        if (!w->options->no_updates) {
            if (update_basis(w, leaving, NULL, 0))
                break;
            factor_ready = 1;
        }
    }
    free(reduced);
    return outcome;
}
struct LabLpSession {
    const LabModel *original;
    LabModel model;
    LabLpOptions options, running;
    Work work;
    double deadline;
    int ready;
};
LabLpSession *lab_lp_session_create(const LabModel *m, const LabLpOptions *options) {
    if (!m || !options || !isfinite(options->time_limit) || !options->max_dimension ||
        !options->iteration_limit ||
        ((options->basis_in || options->basis_out) &&
         (!options->model_sha256 || strlen(options->model_sha256) != 64)))
        return NULL;
    LabLpSession *s = lab_alloc(1, sizeof *s);
    s->original = m;
    s->model = *m;
    s->options = *options;
    s->model.col_lower = lab_alloc(m->cols, sizeof(double));
    s->model.col_upper = lab_alloc(m->cols, sizeof(double));
    memcpy(s->model.col_lower, m->col_lower, m->cols * sizeof(double));
    memcpy(s->model.col_upper, m->col_upper, m->cols * sizeof(double));
    s->deadline = lab_now() + options->time_limit;
    return s;
}
int lab_lp_session_set_bounds(LabLpSession *s, size_t j, double lower, double upper) {
    if (!s || j >= s->model.cols || isnan(lower) || isnan(upper) || lower == INFINITY ||
        upper == -INFINITY || lower < s->original->col_lower[j] ||
        upper > s->original->col_upper[j])
        return -1;
    s->model.col_lower[j] = lower;
    s->model.col_upper[j] = upper;
    return 0;
}
static int update_rhs(LabLpSession *s) {
    Work *w = &s->work;
    const LabModel *m = &s->model;
    for (size_t j = 0; j < m->cols; j++) {
        double lo = m->col_lower[j], hi = m->col_upper[j];
        if (w->transform[j] == 'X' && (lo != w->shift[j] || hi != w->shift[j]))
            return -1;
        if (w->transform[j] == 'L' &&
            (!isfinite(lo) || lo < w->shift[j] || (isfinite(hi) != (w->upper_row[j] != SIZE_MAX))))
            return -1;
        if (w->transform[j] == 'U' && (!isfinite(hi) || hi > w->shift[j] || isfinite(lo)))
            return -1;
        if (w->transform[j] == 'F' && (isfinite(lo) || isfinite(hi)))
            return -1;
    }
    size_t *position = lab_alloc(w->total_rows, sizeof(size_t));
    for (size_t i = 0; i < w->total_rows; i++)
        position[i] = SIZE_MAX;
    for (size_t i = 0; i < w->m; i++)
        position[w->row_id[i]] = i;
    int bad = 0;
    for (size_t j = 0; j < m->cols; j++) {
        if (w->lower_row[j] != SIZE_MAX) {
            size_t i = position[w->lower_row[j]];
            if (i == SIZE_MAX) {
                bad = 1;
                break;
            }
            double raw = w->shift[j] - m->col_lower[j];
            int e = ilogb(w->equation_scale[i]);
            double scaled = scalbn(raw, e);
            if (!isfinite(raw) || !isfinite(scaled) || scalbn(scaled, -e) != raw) {
                bad = 1;
                break;
            }
            w->b[i] = scaled;
        }
        if (w->upper_row[j] != SIZE_MAX) {
            size_t i = position[w->upper_row[j]];
            if (i == SIZE_MAX) {
                bad = 1;
                break;
            }
            double raw = m->col_upper[j] - w->shift[j];
            int e = ilogb(w->equation_scale[i]);
            double scaled = scalbn(raw, e);
            if (!isfinite(raw) || !isfinite(scaled) || scalbn(scaled, -e) != raw) {
                bad = 1;
                break;
            }
            w->b[i] = scaled;
        }
    }
    free(position);
    return bad ? -1 : 0;
}
void lab_lp_session_solve(LabLpSession *s, LabLpResult *result) {
    memset(result, 0, sizeof *result);
    result->status = "numerical_failure";
    double started = lab_now();
    if (!s) {
        snprintf(result->message, sizeof result->message, "Invalid LP session");
        return;
    }
    Work *w = &s->work;
    for (size_t j = 0; j < s->model.cols; j++)
        if (s->model.col_lower[j] > s->model.col_upper[j]) {
            result->status = "infeasible";
            snprintf(result->message, sizeof result->message, "Conflicting column bounds");
            goto done;
        }
    s->running = s->options;
    s->running.time_limit = s->deadline - lab_now();
    if (s->running.time_limit <= 0) {
        result->status = "time_limit";
        snprintf(result->message, sizeof result->message, "Cumulative LP deadline exhausted");
        goto done;
    }
    if (s->ready && !update_rhs(s)) {
        w->options = &s->running;
        w->result = result;
        w->started = lab_now();
        result->warm_start = 1;
        result->standard_rows = w->total_rows;
        result->standard_columns = w->n;
        result->redundant_rows = w->total_rows - w->m;
        int status = dual_simplex(w, 0);
        if (status == 1) {
            s->ready = 0;
            goto done;
        }
        if (status == 0 && !recover(w, &s->model)) {
            message(w, "optimal", "Dual revised simplex repaired the saved basis");
            if (save_basis(w))
                s->ready = 0;
            goto done;
        }
        if (!strcmp(result->status, "time_limit") || !strcmp(result->status, "iteration_limit")) {
            s->ready = 0;
            goto done;
        }
        result->cold_fallbacks++;
        result->warm_start = 0;
        free(result->primal);
        free(result->dual);
        result->primal = result->dual = NULL;
        result->has_primal = 0;
    } else if (s->ready) {
        result->cold_fallbacks++;
    }
    work_free(w);
    s->ready = 0;
    s->running.time_limit = s->deadline - lab_now();
    w->options = &s->running;
    w->result = result;
    w->started = lab_now();
    result->status = "numerical_failure";
    if (s->running.time_limit <= 0) {
        message(w, "time_limit", "Cumulative LP deadline exhausted");
        goto done;
    }
    (void)cold_solve(w, &s->model);
    s->ready = !strcmp(result->status, "optimal");
done:
    result->scaled_rows = w->scaled_rows;
    result->scaled_columns = w->scaled_columns;
    if (w->a)
        for (size_t j = 0; j < w->n; j++)
            result->standard_nnz += w->a[j].count;
    /* build initially records nnz; use one consistent active-matrix count. */
    if (!result->warm_start && w->a) {
        result->standard_nnz = 0;
        for (size_t j = 0; j < w->n; j++)
            result->standard_nnz += w->a[j].count;
    }
    result->elapsed = lab_now() - started;
}
void lab_lp_session_free(LabLpSession *s) {
    if (!s)
        return;
    work_free(&s->work);
    free(s->model.col_lower);
    free(s->model.col_upper);
    free(s);
}
void lab_lp_solve(const LabModel *m, const LabLpOptions *options, LabLpResult *result) {
    LabLpSession *s = lab_lp_session_create(m, options);
    if (!s) {
        memset(result, 0, sizeof *result);
        result->status = "input_error";
        snprintf(result->message, sizeof result->message, "Invalid LP options or basis identity");
        return;
    }
    lab_lp_session_solve(s, result);
    lab_lp_session_free(s);
}

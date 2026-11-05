#include "lp.h"
#include "common.h"
#include "factor.h"
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
    unsigned char *artificial;
    size_t *lower_row, *upper_row;
    unsigned char *transform;
    LabFactor factor;
    double *matrix, *xb, *pi, *cb, *direction, *column;
    const LabLpOptions *options;
    LabLpResult *result;
    double started;
} Work;
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
    long double norm = 0;
    for (size_t k = 0; k < w->a[j].count; k++)
        norm = fmaxl(norm, fabs(w->a[j].value[k]));
    return norm;
}
static int transpose_solve(Work *w, const double *rhs, double *x) {
    if (lab_factor_solve(&w->factor, rhs, x, 1))
        return -1;
    for (size_t j = 0; j < w->m; j++) {
        long double magnitude = fabs(rhs[j]);
        long double sum = dot(w, w->basis[j], x, &magnitude);
        if (!isfinite(sum) || !isfinite(magnitude) ||
            fabsl(sum - rhs[j]) > 1e-10L * (1 + magnitude))
            return -1;
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
static int refactor(Work *w, int primal_required) {
    w->result->refactorizations++;
    int failed;
    if (w->options->dense) {
        memset(w->matrix, 0, w->m * w->m * sizeof(double));
        for (size_t j = 0; j < w->m; j++) {
            const Column *c = &w->a[w->basis[j]];
            for (size_t k = 0; k < c->count; k++)
                w->matrix[c->row[k] * w->m + j] += c->value[k];
        }
        failed = lab_factorize(&w->factor, w->matrix, w->m);
    } else {
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
        failed = lab_factorize_sparse(&w->factor, w->m, start, row, value);
        free(start);
        free(row);
        free(value);
    }
    if (w->factor.peak_nnz > w->result->factor_nnz_peak)
        w->result->factor_nnz_peak = w->factor.peak_nnz;
    if (failed || lab_factor_solve(&w->factor, w->b, w->xb, 0)) {
        message(w, "numerical_failure", "Singular or ill-conditioned basis");
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
        if (primal_required && w->xb[i] < -1e-8) {
            message(w, "numerical_failure", "Basis lost primal feasibility");
            return -1;
        }
        if (primal_required && w->xb[i] < 0)
            w->xb[i] = 0;
    }
    return 0;
}
static int is_basic(const Work *w, size_t j) {
    for (size_t i = 0; i < w->m; i++)
        if (w->basis[i] == j)
            return 1;
    return 0;
}
/* 0 optimum, 1 unbounded, -1 explicit failure/limit. Bland's index rules
 * resolve entering variables and ratio ties deterministically. */
static int simplex(Work *w, int phase_one) {
    for (;;) {
        if (budget(w) || refactor(w, 1))
            return -1;
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
            if (phase_one && worst > 0) {
                double artificial_value = 0;
                for (size_t i = 0; i < w->m; i++)
                    if (w->artificial[w->basis[i]])
                        artificial_value += w->xb[i];
                if (artificial_value > 1e-8) {
                    message(w, "numerical_failure",
                            "Uncertain Phase I dual sign; infeasibility unproved");
                    return -1;
                }
            }
            return 0;
        }
        scatter(w, entering, w->column);
        if (lab_factor_solve(&w->factor, w->column, w->direction, 0)) {
            message(w, "numerical_failure", "Pivot direction solve failed");
            return -1;
        }
        size_t leaving = SIZE_MAX;
        double theta = INFINITY;
        int tiny_positive = 0;
        for (size_t i = 0; i < w->m; i++) {
            double d = w->direction[i];
            if (d > 0 && d <= 1e-12)
                tiny_positive = 1;
            if (d <= 1e-12)
                continue;
            double ratio = w->xb[i] / d;
            if (!isfinite(ratio)) {
                message(w, "numerical_failure", "Ratio test overflow");
                return -1;
            }
            double tie = 1e-12 * (1 + fabs(ratio));
            if (leaving == SIZE_MAX || ratio < theta - tie ||
                (fabs(ratio - theta) <= tie && w->basis[i] < w->basis[leaving])) {
                theta = ratio;
                leaving = i;
            }
        }
        if (leaving == SIZE_MAX) {
            if (phase_one || tiny_positive) {
                message(w, "numerical_failure", "No stable leaving pivot");
                return -1;
            }
            return 1;
        }
        w->basis[leaving] = entering;
        w->result->iterations++;
    }
}
static void remove_row(Work *w, size_t row, size_t basis_position) {
    for (size_t j = 0; j < w->n; j++) {
        Column *c = &w->a[j];
        size_t keep = 0;
        for (size_t k = 0; k < c->count; k++) {
            if (c->row[k] == row)
                continue;
            c->row[keep] = c->row[k] - (size_t)(c->row[k] > row);
            c->value[keep++] = c->value[k];
        }
        c->count = keep;
    }
    for (size_t i = row; i + 1 < w->m; i++) {
        w->b[i] = w->b[i + 1];
        w->row_id[i] = w->row_id[i + 1];
        w->source_row[i] = w->source_row[i + 1];
        w->row_scale[i] = w->row_scale[i + 1];
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
            w->basis[i] = entering;
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
static int build(Work *w, const LabModel *m) {
    size_t dim = w->options->max_dimension, nz = 0, nr = 0;
    if (dim > 32768 || (w->options->dense && dim > 4096)) {
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
    if (m->nnz > (20000000U - 3 * nr) / 4) {
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
    w->row_id = lab_alloc(cap, sizeof(size_t));
    w->artificial = lab_alloc(n, 1);
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
            double flip = rhs < 0 ? -1 : 1;
            w->b[r] = (double)(flip * rhs);
            w->source_row[r] = i;
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
    for (size_t i = 0; i < nr; i++) {
        if (!isfinite(w->b[i])) {
            message(w, "numerical_failure", "Bound transformation overflow");
            return -1;
        }
        if (w->artificial[w->basis[i]] && w->b[i] > 0) {
            int nonzero = 0;
            for (size_t j = 0; j < nz && !nonzero; j++)
                for (size_t k = 0; k < w->a[j].count; k++)
                    if (w->a[j].row[k] == i)
                        nonzero = 1;
            if (!nonzero) {
                message(w, "infeasible", "Contradictory constant row");
                return -1;
            }
        }
    }
    for (size_t j = 0; j < n; j++) {
        w->result->standard_nnz += w->a[j].count;
        for (size_t k = 0; k < w->a[j].count; k++)
            if (!isfinite(w->a[j].value[k])) {
                message(w, "numerical_failure", "Coefficient transformation overflow");
                return -1;
            }
    }
    w->matrix = w->options->dense ? lab_alloc(cap * cap, sizeof(double)) : NULL;
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
static int recover(Work *w, const LabModel *m) {
    LabLpResult *r = w->result;
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
    long double objective = m->offset, variation = 0;
    double residual = 0;
    long double *activity = lab_alloc(m->rows, sizeof(long double));
    for (size_t j = 0; j < m->cols; j++) {
        double x = r->primal[j];
        if (!isfinite(x)) {
            free(activity);
            message(w, "numerical_failure", "Primal reconstruction overflow");
            return -1;
        }
        objective += (long double)m->objective[j] * x;
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
    r->objective = (double)objective;
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
    if (available && (!isfinite(gap) || !isfinite(variation) || gap > 1e-7L * (1 + variation))) {
        message(w, "numerical_failure", "Primal/support-bound gap exceeds stopping tolerance");
        return -1;
    }
    r->has_primal = 1;
    return 0;
}
void lab_lp_result_free(LabLpResult *result) {
    free(result->primal);
    free(result->dual);
    memset(result, 0, sizeof *result);
}
static void work_free(Work *w) {
    lab_factor_free(&w->factor);
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
    free(w->artificial);
    free(w->lower_row);
    free(w->upper_row);
    free(w->transform);
    free(w->matrix);
    free(w->xb);
    free(w->pi);
    free(w->cb);
    free(w->direction);
    free(w->column);
    memset(w, 0, sizeof *w);
}
static int cold_solve(Work *w, const LabModel *m) {
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
    if (status < 0 || recover(w, m))
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
static int dual_simplex(Work *w) {
    double *reduced = lab_alloc(w->n, sizeof(double));
    int outcome = -1;
    for (;;) {
        if (budget(w) || refactor(w, 0))
            break;
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
            if (w->xb[i] < -1e-8 && (leaving == SIZE_MAX || w->basis[i] < w->basis[leaving]))
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
        w->basis[leaving] = entering;
        w->result->iterations++;
        w->result->dual_iterations++;
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
            w->b[i] = w->shift[j] - m->col_lower[j];
            if (!isfinite(w->b[i])) {
                bad = 1;
                break;
            }
        }
        if (w->upper_row[j] != SIZE_MAX) {
            size_t i = position[w->upper_row[j]];
            if (i == SIZE_MAX) {
                bad = 1;
                break;
            }
            w->b[i] = m->col_upper[j] - w->shift[j];
            if (!isfinite(w->b[i])) {
                bad = 1;
                break;
            }
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
        int status = dual_simplex(w);
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

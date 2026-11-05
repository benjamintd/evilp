#include "lp.h"
#include "common.h"
#include "factor.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Revised primal simplex reference engine. The original model stays immutable.
 * The standard form is deliberately dense and resource capped. Every pivot
 * refactors the basis; factor updates and sparse LU belong to the next backend. */
typedef struct {
    size_t m, n, structural, total_rows;
    double *a, *b, *cost, *shift, *sign;
    size_t *basis, *original, *row_id, *art_row;
    size_t *source_row;
    double *row_scale;
    unsigned char *artificial;
    LabFactor factor;
    double *matrix, *xb, *pi, *cb, *direction, *column;
    const LabLpOptions *options;
    LabLpResult *result;
    double started;
} Work;
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
static int refactor(Work *w) {
    for (size_t i = 0; i < w->m; i++)
        for (size_t j = 0; j < w->m; j++)
            w->matrix[i * w->m + j] = w->a[i * w->n + w->basis[j]];
    w->result->refactorizations++;
    if (lab_factorize(&w->factor, w->matrix, w->m) ||
        lab_factor_solve(&w->factor, w->b, w->xb, 0)) {
        message(w, "numerical_failure", "Singular or ill-conditioned basis");
        return -1;
    }
    /* Recompute residuals rather than trusting a successful factorization. */
    double residual = 0, scale = 1;
    for (size_t i = 0; i < w->m; i++) {
        long double sum = 0, mag = fabsl(w->b[i]);
        for (size_t j = 0; j < w->m; j++) {
            long double term = (long double)w->matrix[i * w->m + j] * w->xb[j];
            sum += term;
            mag += fabsl(term);
        }
        residual = fmax(residual, (double)fabsl(sum - w->b[i]));
        scale = fmax(scale, (double)mag);
    }
    if (!isfinite(residual) || residual > 1e-10 * scale) {
        message(w, "numerical_failure", "Basis solve failed residual check");
        return -1;
    }
    for (size_t i = 0; i < w->m; i++) {
        if (w->xb[i] < -1e-8) {
            message(w, "numerical_failure", "Basis lost primal feasibility");
            return -1;
        }
        if (w->xb[i] < 0)
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
    double cost_scale = phase_one ? 1 : 0;
    if (!phase_one)
        for (size_t j = 0; j < w->n; j++)
            cost_scale = fmax(cost_scale, fabs(w->cost[j]));
    for (;;) {
        if (budget(w) || refactor(w))
            return -1;
        for (size_t i = 0; i < w->m; i++)
            w->cb[i] = phase_one ? (double)w->artificial[w->basis[i]] : w->cost[w->basis[i]];
        if (lab_factor_solve(&w->factor, w->cb, w->pi, 1)) {
            message(w, "numerical_failure", "Transpose basis solve failed");
            return -1;
        }
        size_t entering = SIZE_MAX;
        double worst = 0;
        for (size_t j = 0; j < w->n; j++) {
            if ((!phase_one && w->artificial[j]) || is_basic(w, j))
                continue;
            long double reduced = phase_one ? (double)w->artificial[j] : w->cost[j];
            /* Preserve small explicit costs. Zero-cost columns need a floor
             * from the whole objective to avoid pricing transpose-solve noise. */
            long double scale = reduced == 0 ? cost_scale : fabsl(reduced);
            for (size_t i = 0; i < w->m; i++) {
                long double term = (long double)w->a[i * w->n + j] * w->pi[i];
                reduced -= term;
                scale += fabsl(term);
            }
            if (!isfinite((double)reduced) || !isfinite((double)scale)) {
                message(w, "numerical_failure", "Reduced cost overflow");
                return -1;
            }
            worst = fmax(worst, (double)-reduced);
            if (entering == SIZE_MAX && reduced < -1e-10L * fmaxl(scale, DBL_MIN))
                entering = j;
        }
        if (!phase_one)
            w->result->dual_residual = worst;
        if (entering == SIZE_MAX)
            return 0;
        for (size_t i = 0; i < w->m; i++)
            w->column[i] = w->a[i * w->n + entering];
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
    for (size_t i = row; i + 1 < w->m; i++) {
        memmove(w->a + i * w->n, w->a + (i + 1) * w->n, w->n * sizeof(double));
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
        if (budget(w) || refactor(w))
            return -1;
        if (fabs(w->xb[i]) > 1e-9) {
            message(w, "numerical_failure", "Nonzero artificial variable after Phase I");
            return -1;
        }
        memset(w->cb, 0, w->m * sizeof(double));
        w->cb[i] = 1;
        if (lab_factor_solve(&w->factor, w->cb, w->pi, 1))
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
            for (size_t k = 0; k < w->m; k++) {
                long double term = (long double)w->pi[k] * w->a[k * w->n + j];
                v += term;
                magnitude += fabsl(term);
                column_norm = fmaxl(column_norm, fabsl(w->a[k * w->n + j]));
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
    for (size_t j = 0; j < m->cols; j++) {
        if (m->col_lower[j] > m->col_upper[j]) {
            message(w, "infeasible", "Conflicting column bounds");
            return -1;
        }
        if (m->col_lower[j] == m->col_upper[j])
            continue;
        nz += (!isfinite(m->col_lower[j]) && !isfinite(m->col_upper[j])) ? 2U : 1U;
        if (isfinite(m->col_lower[j]) && isfinite(m->col_upper[j]))
            nr++;
    }
    for (size_t i = 0; i < m->rows; i++) {
        if (m->row_lower[i] > m->row_upper[i]) {
            message(w, "infeasible", "Conflicting row bounds");
            return -1;
        }
        if (m->row_lower[i] == m->row_upper[i])
            nr++;
        else
            nr += (size_t)isfinite(m->row_lower[i]) + (size_t)isfinite(m->row_upper[i]);
    }
    /* At most one slack and one artificial per equation. Check dimensions
     * before allocating any dense matrix, including the transformation. */
    if (nr > dim || nz > dim || nr > (dim - nz) / 2) {
        message(w, "unsupported", "Dense reference backend dimension limit exceeded");
        return -1;
    }
    w->structural = nz;
    w->m = nr;
    w->total_rows = nr;
    w->n = nz + 2 * nr;
    size_t n = w->n, cap = nr ? nr : 1;
    w->a = lab_alloc(cap * (n ? n : 1), sizeof(double));
    w->b = lab_alloc(cap, sizeof(double));
    w->cost = lab_alloc(n ? n : 1, sizeof(double));
    w->basis = lab_alloc(cap, sizeof(size_t));
    w->source_row = lab_alloc(cap, sizeof(size_t));
    w->row_scale = lab_alloc(cap, sizeof(double));
    w->row_id = lab_alloc(cap, sizeof(size_t));
    w->artificial = lab_alloc(n ? n : 1, 1);
    w->art_row = lab_alloc(n ? n : 1, sizeof(size_t));
    w->shift = lab_alloc(m->cols, sizeof(double));
    w->original = lab_alloc(nz ? nz : 1, sizeof(size_t));
    w->sign = lab_alloc(nz ? nz : 1, sizeof(double));
    size_t *first = lab_alloc(m->cols, sizeof(size_t));
    size_t *second = lab_alloc(m->cols, sizeof(size_t));
    size_t z = 0;
    double sense = m->maximize ? -1 : 1;
    for (size_t j = 0; j < m->cols; j++) {
        first[j] = second[j] = SIZE_MAX;
        if (m->col_lower[j] == m->col_upper[j]) {
            w->shift[j] = m->col_lower[j];
            continue;
        }
        w->shift[j] = isfinite(m->col_lower[j])   ? m->col_lower[j]
                      : isfinite(m->col_upper[j]) ? m->col_upper[j]
                                                  : 0;
        first[j] = z;
        w->original[z] = j;
        w->sign[z] = isfinite(m->col_lower[j]) || !isfinite(m->col_upper[j]) ? 1 : -1;
        w->cost[z] = sense * m->objective[j] * w->sign[z];
        z++;
        if (!isfinite(m->col_lower[j]) && !isfinite(m->col_upper[j])) {
            second[j] = z;
            w->original[z] = j;
            w->sign[z] = -1;
            w->cost[z++] = -sense * m->objective[j];
        }
    }
    /* Temporary row-major matrix of the shifted structural variables. */
    double *row = lab_alloc(m->rows * (nz ? nz : 1), sizeof(double));
    double *shifted = lab_alloc(m->rows, sizeof(double));
    for (size_t j = 0; j < m->cols; j++)
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            size_t i = m->row_index[k];
            double a = m->value[k];
            shifted[i] += a * w->shift[j];
            if (first[j] != SIZE_MAX)
                row[i * nz + first[j]] += a * w->sign[first[j]];
            if (second[j] != SIZE_MAX)
                row[i * nz + second[j]] -= a;
        }
    size_t r = 0;
    for (size_t i = 0; i < m->rows; i++) {
        int equality = m->row_lower[i] == m->row_upper[i];
        for (int side = 0; side < 2; side++) {
            double bound = side ? m->row_upper[i] : m->row_lower[i];
            if (!isfinite(bound) || (equality && side))
                continue;
            double rhs = bound - shifted[i], flip = rhs < 0 ? -1 : 1;
            w->b[r] = flip * rhs;
            w->source_row[r] = i;
            w->row_scale[r] = flip;
            w->row_id[r] = r;
            for (size_t j = 0; j < nz; j++)
                w->a[r * n + j] = flip * row[i * nz + j];
            double slack = equality ? 0 : flip * (side ? 1 : -1);
            int nonzero = 0;
            for (size_t j = 0; j < nz; j++)
                nonzero |= w->a[r * n + j] != 0;
            if (!nonzero && slack <= 0 && w->b[r] > 0) {
                free(first);
                free(second);
                free(row);
                free(shifted);
                message(w, "infeasible", "Contradictory constant row");
                return -1;
            }
            w->a[r * n + nz + r] = slack;
            w->artificial[nz + nr + r] = 1;
            w->art_row[nz + nr + r] = r;
            w->a[r * n + nz + nr + r] = 1;
            w->basis[r] = slack == 1 ? nz + r : nz + nr + r;
            r++;
        }
    }
    for (size_t j = 0; j < m->cols; j++)
        if (isfinite(m->col_lower[j]) && isfinite(m->col_upper[j]) && first[j] != SIZE_MAX) {
            w->b[r] = m->col_upper[j] - m->col_lower[j];
            w->a[r * n + first[j]] = 1;
            w->a[r * n + nz + r] = 1;
            w->a[r * n + nz + nr + r] = 1;
            w->artificial[nz + nr + r] = 1;
            w->art_row[nz + nr + r] = r;
            w->basis[r] = nz + r;
            w->source_row[r] = SIZE_MAX;
            w->row_id[r] = r;
            r++;
        }
    free(first);
    free(second);
    free(row);
    free(shifted);
    for (size_t i = 0; i < nr; i++)
        if (!isfinite(w->b[i])) {
            message(w, "numerical_failure", "Bound transformation overflow");
            return -1;
        }
    for (size_t i = 0; i < nr * n; i++)
        if (!isfinite(w->a[i])) {
            message(w, "numerical_failure", "Coefficient transformation overflow");
            return -1;
        }
    w->matrix = lab_alloc(cap * cap, sizeof(double));
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
    long double objective = m->offset;
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
    free(activity);
    r->objective = (double)objective;
    r->primal_residual = residual;
    if (!isfinite(r->objective) || !isfinite(residual) || residual > 1e-7) {
        message(w, "numerical_failure", "Original-model primal residual exceeds 1e-7");
        return -1;
    }
    /* A tolerance-small reduced cost can imply a large objective error over
     * wide bounds. This floating support calculation is a rejection safeguard,
     * not a certificate; milpcheck recomputes a rigorous bound independently. */
    long double sense = m->maximize ? -1 : 1;
    long double lower = sense * m->offset;
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
        lower += (long double)y * side;
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
        lower += reduced * side;
    }
    if (available &&
        (!isfinite(lower) || fabsl(sense * objective - lower) > 1e-7L * (1 + fabsl(objective)))) {
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
void lab_lp_solve(const LabModel *m, const LabLpOptions *options, LabLpResult *result) {
    memset(result, 0, sizeof *result);
    result->status = "numerical_failure";
    Work w = {0};
    w.options = options;
    w.result = result;
    w.started = lab_now();
    if (build(&w, m))
        goto cleanup;
    {
        if (simplex(&w, 1))
            goto cleanup;
        double violation = 0;
        for (size_t i = 0; i < w.m; i++)
            if (w.artificial[w.basis[i]])
                violation += w.xb[i];
        if (violation > 1e-8) {
            message(&w, "infeasible", "Positive Phase I artificial objective (numerical claim)");
            goto cleanup;
        }
        if (remove_artificials(&w))
            goto cleanup;
    }
    if (load_basis(&w) < 0)
        goto cleanup;
    int status = simplex(&w, 0);
    if (status < 0)
        goto cleanup;
    if (recover(&w, m))
        goto cleanup;
    if (status) {
        message(&w, "unbounded", "Improving LP ray found (numerical claim)");
    } else {
        message(&w, "optimal", "LP reduced costs satisfy numerical stopping tolerance");
        (void)save_basis(&w);
    }
cleanup:
    result->elapsed = lab_now() - w.started;
    lab_factor_free(&w.factor);
    free(w.a);
    free(w.b);
    free(w.cost);
    free(w.shift);
    free(w.sign);
    free(w.basis);
    free(w.original);
    free(w.row_id);
    free(w.art_row);
    free(w.source_row);
    free(w.row_scale);
    free(w.artificial);
    free(w.matrix);
    free(w.xb);
    free(w.pi);
    free(w.cb);
    free(w.direction);
    free(w.column);
}

#include "factor.h"
#include "common.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#endif

typedef struct {
    size_t column;
    double value;
} SparseEntry;
typedef struct {
    SparseEntry *entry;
    size_t count, capacity;
} SparseRow;
typedef struct {
    SparseRow *row;
    size_t *col_start, *row_index;
    double *value, *diagonal;
    /* Structural singleton elimination keeps the original sparse matrix and
     * factors only its coupled core. No numerical coefficient is dropped. */
    LabFactor core;
    size_t core_n, steps;
    size_t *core_row, *core_col, *step_row, *step_col;
    unsigned char *step_kind, *row_kind, *col_kind;
    double *step_pivot;
    size_t *core_start, *core_index;
    double *core_value;
} SparseFactor;

/* Both limits are independent of LP/model resource policy. There is no
 * silently selected dense fallback when fill exhausts the sparse budget. */
enum { SPARSE_FILL_BUDGET = 64000000, SPARSE_DIMENSION_LIMIT = 8000000 };

static void sparse_free(SparseFactor *s, size_t n) {
    if (!s)
        return;
    if (s->row)
        for (size_t i = 0; i < n; i++)
            free(s->row[i].entry);
    free(s->row);
    free(s->col_start);
    free(s->row_index);
    free(s->value);
    free(s->diagonal);
    lab_factor_free(&s->core);
    free(s->core_row); free(s->core_col);
    free(s->step_row); free(s->step_col); free(s->step_kind); free(s->step_pivot);
    free(s->row_kind); free(s->col_kind);
    free(s->core_start); free(s->core_index); free(s->core_value);
    free(s);
}
void lab_factor_free(LabFactor *f) {
    sparse_free(f->sparse_data, f->n);
    free(f->lu);
    free(f->pivot);
    memset(f, 0, sizeof *f);
}
int lab_factorize(LabFactor *f, const double *matrix, size_t n) {
    lab_factor_free(f);
    if ((n && !matrix) || (n && (n > SIZE_MAX / n || n * n > SIZE_MAX / sizeof(double))))
        return -1;
    f->n = n;
    f->core_dimension = n;
    f->core_refinement = 1;
    f->lu = lab_alloc(n * n, sizeof(double));
    f->pivot = lab_alloc(n, sizeof(size_t));
    if (n)
        memcpy(f->lu, matrix, n * n * sizeof(double));
    f->peak_nnz = n * n;
    double norm = 0;
    for (size_t i = 0; i < n * n; i++) {
        if (!isfinite(f->lu[i]))
            return -1;
        norm = fmax(norm, fabs(f->lu[i]));
    }
    for (size_t k = 0; k < n; k++) {
        size_t p = k;
        for (size_t i = k + 1; i < n; i++)
            if (fabs(f->lu[i * n + k]) > fabs(f->lu[p * n + k]))
                p = i;
        if (fabs(f->lu[p * n + k]) <= fmax(DBL_MIN, 16 * DBL_EPSILON * norm))
            return -1;
        f->pivot[k] = p;
        if (p != k)
            for (size_t j = 0; j < n; j++) {
                double t = f->lu[k * n + j];
                f->lu[k * n + j] = f->lu[p * n + j];
                f->lu[p * n + j] = t;
            }
        for (size_t i = k + 1; i < n; i++) {
            f->lu[i * n + k] /= f->lu[k * n + k];
            for (size_t j = k + 1; j < n; j++)
                f->lu[i * n + j] -= f->lu[i * n + k] * f->lu[k * n + j];
        }
    }
    for (size_t i = 0; i < n * n; i++) {
        if (!isfinite(f->lu[i]))
            return -1;
        f->factor_nnz += f->lu[i] != 0;
    }
    f->valid = 1;
    return 0;
}
static double sparse_get(const SparseRow *row, size_t column) {
    size_t lo = 0, hi = row->count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (row->entry[mid].column < column)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo < row->count && row->entry[lo].column == column ? row->entry[lo].value : 0;
}
static int sparse_append(SparseRow *row, size_t column, double value, size_t *total) {
    if (value == 0)
        return 0;
    if (row->count && row->entry[row->count - 1].column == column) {
        SparseEntry *last = &row->entry[row->count - 1];
        last->value += value;
        if (!isfinite(last->value))
            return -1;
        if (last->value == 0) {
            row->count--;
            (*total)--;
        }
        return 0;
    }
    if (*total >= SPARSE_FILL_BUDGET)
        return -1;
    if (row->count == row->capacity) {
        size_t cap = row->capacity ? 2 * row->capacity : 4;
        if (cap > SPARSE_FILL_BUDGET)
            cap = SPARSE_FILL_BUDGET;
        SparseEntry *next = realloc(row->entry, cap * sizeof *next);
        if (!next)
            return -1;
        row->entry = next;
        row->capacity = cap;
    }
    row->entry[row->count++] = (SparseEntry){column, value};
    (*total)++;
    return 0;
}
static int sparse_eliminate(SparseRow *row, const SparseRow *pivot, size_t column,
                            double multiplier, size_t *total, size_t *peak) {
    size_t p = 0;
    while (p < pivot->count && pivot->entry[p].column <= column)
        p++;
    size_t tail = pivot->count - p;
    if (row->count > SPARSE_FILL_BUDGET - tail)
        return -1;
    size_t capacity = row->count + tail;
    SparseEntry *next = malloc((capacity ? capacity : 1) * sizeof *next);
    if (!next)
        return -1;
    size_t old = 0, count = 0;
    while (old < row->count && row->entry[old].column < column)
        next[count++] = row->entry[old++];
    /* The eliminated coefficient is stored as the L multiplier. */
    if (old == row->count || row->entry[old].column != column) {
        free(next);
        return -1;
    }
    old++;
    next[count++] = (SparseEntry){column, multiplier};
    while (old < row->count || p < pivot->count) {
        size_t j;
        double value;
        if (p == pivot->count ||
            (old < row->count && row->entry[old].column < pivot->entry[p].column)) {
            j = row->entry[old].column;
            value = row->entry[old++].value;
        } else if (old == row->count || pivot->entry[p].column < row->entry[old].column) {
            j = pivot->entry[p].column;
            value = -multiplier * pivot->entry[p++].value;
        } else {
            j = row->entry[old].column;
            value = row->entry[old++].value - multiplier * pivot->entry[p++].value;
        }
        if (!isfinite(value)) {
            free(next);
            return -1;
        }
        if (value != 0)
            next[count++] = (SparseEntry){j, value};
    }
    if (count > SPARSE_FILL_BUDGET - (*total - row->count)) {
        free(next);
        return -1;
    }
    /* Do not retain historical merge capacity after cancellation: storage
     * should follow live fill, rather than accumulate old row high-water marks. */
    if (count < capacity) {
        SparseEntry *compact = realloc(next, count * sizeof *compact);
        if (!compact) {
            free(next);
            return -1;
        }
        next = compact;
        capacity = count;
    }
    *total = *total - row->count + count;
    if (*total > *peak)
        *peak = *total;
    free(row->entry);
    row->entry = next;
    row->count = count;
    row->capacity = capacity;
    return 0;
}
static int sparse_transpose(SparseFactor *s, size_t n, size_t nnz) {
    s->col_start = calloc(n + 1, sizeof *s->col_start);
    s->row_index = malloc((nnz ? nnz : 1) * sizeof *s->row_index);
    s->value = malloc((nnz ? nnz : 1) * sizeof *s->value);
    s->diagonal = malloc((n ? n : 1) * sizeof *s->diagonal);
    size_t *cursor = calloc(n ? n : 1, sizeof *cursor);
    if (!s->col_start || !s->row_index || !s->value || !s->diagonal || !cursor) {
        free(cursor);
        return -1;
    }
    for (size_t i = 0; i < n; i++)
        for (size_t e = 0; e < s->row[i].count; e++)
            s->col_start[s->row[i].entry[e].column + 1]++;
    for (size_t j = 0; j < n; j++) {
        s->col_start[j + 1] += s->col_start[j];
        cursor[j] = s->col_start[j];
    }
    for (size_t i = 0; i < n; i++) {
        s->diagonal[i] = sparse_get(&s->row[i], i);
        if (!isfinite(s->diagonal[i]) || s->diagonal[i] == 0) {
            free(cursor);
            return -1;
        }
        for (size_t e = 0; e < s->row[i].count; e++) {
            SparseEntry a = s->row[i].entry[e];
            size_t k = cursor[a.column]++;
            s->row_index[k] = i;
            s->value[k] = a.value;
        }
    }
    free(cursor);
    return 0;
}
static double smallest_pivot(const LabFactor *f) {
    if (!f->n)
        return INFINITY;
    const SparseFactor *s = f->sparse_data;
    double smallest = INFINITY;
    if (s->steps) {
        for (size_t k = 0; k < s->steps; k++)
            smallest = fmin(smallest, fabs(s->step_pivot[k]));
        return fmin(smallest, smallest_pivot(&s->core));
    }
    for (size_t i = 0; i < f->n; i++)
        smallest = fmin(smallest, fabs(s->diagonal[i]));
    return smallest;
}
/* 1: peeled factor ready; 0: no structural reduction; -1: failure.
 * Singleton rows fix a variable; singleton columns defer one equation. Their
 * removal produces no fill, and reverse substitution recovers deferred values. */
static int peel_factor(LabFactor *f, size_t total, double norm) {
    SparseFactor *s = f->sparse_data;
    size_t n = f->n;
    if (!n)
        return 0;
    size_t *rd = lab_alloc(n, sizeof(size_t)), *cd = lab_alloc(n, sizeof(size_t));
    size_t *start = lab_alloc(n + 1, sizeof(size_t));
    size_t *ri = lab_alloc(total, sizeof(size_t)), *cursor = lab_alloc(n, sizeof(size_t));
    size_t *queue = lab_alloc(2 * n, sizeof(size_t));
    unsigned char *active_r = lab_alloc(n, 1), *active_c = lab_alloc(n, 1);
    for (size_t i = 0; i < n; i++) {
        active_r[i] = active_c[i] = 1;
        rd[i] = s->row[i].count;
        for (size_t k = 0; k < s->row[i].count; k++)
            start[s->row[i].entry[k].column + 1]++;
    }
    for (size_t j = 0; j < n; j++) {
        cd[j] = start[j + 1];
        start[j + 1] += start[j];
        cursor[j] = start[j];
    }
    for (size_t i = 0; i < n; i++)
        for (size_t k = 0; k < s->row[i].count; k++)
            ri[cursor[s->row[i].entry[k].column]++] = i;
    s->step_row = lab_alloc(n, sizeof(size_t));
    s->step_col = lab_alloc(n, sizeof(size_t));
    s->step_kind = lab_alloc(n, 1);
    s->step_pivot = lab_alloc(n, sizeof(double));
    s->row_kind = lab_alloc(n, 1);
    s->col_kind = lab_alloc(n, 1);
    size_t head = 0, tail = 0;
    for (size_t j = 0; j < n; j++)
        if (cd[j] == 1)
            queue[tail++] = n + j;
    for (size_t i = 0; i < n; i++)
        if (rd[i] == 1)
            queue[tail++] = i;
    while (head < tail) {
        size_t token = queue[head++], r = SIZE_MAX, c = SIZE_MAX;
        unsigned char kind = token < n ? 1 : 2;
        if (kind == 1) {
            r = token;
            if (!active_r[r] || rd[r] != 1)
                continue;
            for (size_t k = 0; k < s->row[r].count; k++)
                if (active_c[s->row[r].entry[k].column]) {
                    c = s->row[r].entry[k].column;
                    break;
                }
        } else {
            c = token - n;
            if (!active_c[c] || cd[c] != 1)
                continue;
            for (size_t k = start[c]; k < start[c + 1]; k++)
                if (active_r[ri[k]]) {
                    r = ri[k];
                    break;
                }
        }
        if (r == SIZE_MAX || c == SIZE_MAX)
            continue;
        size_t at = s->steps++;
        s->step_row[at] = r; s->step_col[at] = c; s->step_kind[at] = kind;
        s->step_pivot[at] = sparse_get(&s->row[r], c);
        s->row_kind[r] = s->col_kind[c] = kind;
        active_r[r] = active_c[c] = 0;
        for (size_t k = 0; k < s->row[r].count; k++) {
            size_t j = s->row[r].entry[k].column;
            if (active_c[j] && --cd[j] == 1)
                queue[tail++] = n + j;
        }
        for (size_t k = start[c]; k < start[c + 1]; k++) {
            size_t i = ri[k];
            if (active_r[i] && --rd[i] == 1)
                queue[tail++] = i;
        }
    }
    int outcome = 0;
    if (s->steps) {
        s->core_n = n - s->steps;
        s->core_row = lab_alloc(s->core_n, sizeof(size_t));
        s->core_col = lab_alloc(s->core_n, sizeof(size_t));
        size_t *row_map = lab_alloc(n, sizeof(size_t));
        size_t rcount = 0, ccount = 0, cnz = 0;
        for (size_t i = 0; i < n; i++) {
            row_map[i] = SIZE_MAX;
            if (active_r[i]) {
                row_map[i] = rcount;
                s->core_row[rcount++] = i;
            }
            if (active_c[i])
                s->core_col[ccount++] = i;
        }
        size_t *cs = lab_alloc(s->core_n + 1, sizeof(size_t));
        size_t *cr = lab_alloc(total, sizeof(size_t));
        double *cv = lab_alloc(total, sizeof(double));
        for (size_t j = 0; j < s->core_n; j++) {
            size_t original = s->core_col[j];
            cs[j] = cnz;
            for (size_t k = start[original]; k < start[original + 1]; k++) {
                size_t i = ri[k];
                if (active_r[i]) {
                    cr[cnz] = row_map[i];
                    cv[cnz++] = sparse_get(&s->row[i], original);
                }
            }
        }
        cs[s->core_n] = cnz;
        outcome = -1;
        if (rcount == s->core_n && ccount == s->core_n &&
            !lab_factorize_sparse(&s->core, s->core_n, cs, cr, cv) &&
            fmin(smallest_pivot(&s->core), smallest_pivot(f)) >
                fmax(DBL_MIN, 16 * DBL_EPSILON * norm) &&
            s->core.peak_nnz <= SPARSE_FILL_BUDGET - total) {
            /* Build an original CSC mirror for forward/transpose substitution.
             * Unlike the ordinary LU mirror, no original diagonal is needed. */
            s->col_start = start; start = NULL;
            s->row_index = ri; ri = NULL;
            s->value = lab_alloc(total, sizeof(double));
            for (size_t j = 0; j < n; j++)
                for (size_t k = s->col_start[j]; k < s->col_start[j + 1]; k++)
                    s->value[k] = sparse_get(&s->row[s->row_index[k]], j);
            f->factor_nnz = total + s->core.factor_nnz;
            f->peak_nnz = total + s->core.peak_nnz;
            f->core_dimension = s->core.core_dimension;
            f->eliminated_singletons = s->steps + s->core.eliminated_singletons;
            f->valid = 1;
            s->core_start = cs; cs = NULL;
            s->core_index = cr; cr = NULL;
            s->core_value = cv; cv = NULL;
            outcome = 1;
        }
        free(row_map); free(cs); free(cr); free(cv);
    }
    free(rd); free(cd); free(start); free(ri); free(cursor); free(queue);
    free(active_r); free(active_c);
    return outcome;
}
int lab_factorize_sparse(LabFactor *f, size_t n, const size_t *col_start, const size_t *row_index,
                         const double *value) {
    lab_factor_free(f);
    if (n > SPARSE_DIMENSION_LIMIT || !col_start || col_start[0] != 0 ||
        col_start[n] > SPARSE_FILL_BUDGET || (col_start[n] && (!row_index || !value)))
        return -1;
    for (size_t j = 0; j < n; j++)
        if (col_start[j] > col_start[j + 1])
            return -1;
    f->n = n;
    f->core_dimension = n;
    f->core_refinement = 1;
    f->sparse = 1;
    SparseFactor *s = calloc(1, sizeof *s);
    if (!s)
        return -1;
    f->sparse_data = s;
    s->row = calloc(n ? n : 1, sizeof *s->row);
    f->pivot = malloc((n ? n : 1) * sizeof *f->pivot);
    if (!s->row || !f->pivot)
        return -1;
    size_t total = 0;
    for (size_t j = 0; j < n; j++) {
        for (size_t e = col_start[j]; e < col_start[j + 1]; e++) {
            if (row_index[e] >= n || !isfinite(value[e]) ||
                sparse_append(&s->row[row_index[e]], j, value[e], &total))
                return -1;
            if (total > f->peak_nnz)
                f->peak_nnz = total;
        }
    }
    double norm = 0;
    for (size_t i = 0; i < n; i++) {
        SparseRow *row = &s->row[i];
        if (row->count == 0) {
            free(row->entry);
            row->entry = NULL;
            row->capacity = 0;
        } else if (row->count < row->capacity) {
            SparseEntry *compact = realloc(row->entry, row->count * sizeof *compact);
            if (!compact)
                return -1;
            row->entry = compact;
            row->capacity = row->count;
        }
        for (size_t e = 0; e < s->row[i].count; e++)
            norm = fmax(norm, fabs(s->row[i].entry[e].value));
    }
    int peeled = peel_factor(f, total, norm);
    if (peeled)
        return peeled > 0 ? 0 : -1;
    for (size_t k = 0; k < n; k++) {
        size_t p = k;
        double largest = fabs(sparse_get(&s->row[k], k));
        for (size_t i = k + 1; i < n; i++) {
            double candidate = fabs(sparse_get(&s->row[i], k));
            if (candidate > largest) {
                largest = candidate;
                p = i;
            }
        }
        if (largest <= fmax(DBL_MIN, 16 * DBL_EPSILON * norm))
            return -1;
        f->pivot[k] = p;
        if (p != k) {
            SparseRow tmp = s->row[k];
            s->row[k] = s->row[p];
            s->row[p] = tmp;
        }
        double diagonal = sparse_get(&s->row[k], k);
        for (size_t i = k + 1; i < n; i++) {
            double entry = sparse_get(&s->row[i], k);
            if (entry == 0)
                continue;
            double multiplier = entry / diagonal;
            if (!isfinite(multiplier) || multiplier == 0 ||
                sparse_eliminate(&s->row[i], &s->row[k], k, multiplier, &total, &f->peak_nnz))
                return -1;
        }
    }
    f->factor_nnz = total;
    if (sparse_transpose(s, n, total))
        return -1;
    f->valid = 1;
    return 0;
}
static int compensated_add(double *high, double *low, double term) {
    double next = *high + term;
    double virtual_term = next - *high;
    double error = (*high - (next - virtual_term)) + (term - virtual_term);
    *low += error;
    *high = next;
    return !isfinite(next) || !isfinite(*low) ? -1 : 0;
}
static int core_residual(const SparseFactor *s, const double *b, const double *x, double *r,
                          int transpose, double *error) {
    size_t n = s->core_n;
    double *low = lab_alloc(n, sizeof(double));
    memcpy(r, b, n * sizeof(double));
    int bad = 0;
    for (size_t j = 0; j < n && !bad; j++)
        for (size_t k = s->core_start[j]; k < s->core_start[j + 1]; k++) {
            size_t i = s->core_index[k], destination = transpose ? j : i;
            double a = s->core_value[k], v = x[transpose ? i : j], product = a * v;
            double tail = fma(a, v, -product);
            if (!isfinite(product) || !isfinite(tail) ||
                compensated_add(&r[destination], &low[destination], -product) ||
                compensated_add(&r[destination], &low[destination], -tail)) {
                bad = 1;
                break;
            }
        }
    *error = 0;
    for (size_t i = 0; i < n && !bad; i++) {
        r[i] += low[i];
        if (!isfinite(r[i]))
            bad = 1;
        *error = fmax(*error, fabs(r[i]));
    }
    free(low);
    return bad ? -1 : 0;
}
static int peeled_solve(const LabFactor *f, const double *rhs, double *x, int transpose) {
    const SparseFactor *s = f->sparse_data;
    size_t n = f->n;
    double *work = lab_alloc(n, sizeof(double));
    memcpy(work, rhs, n * sizeof(double));
    memset(x, 0, n * sizeof(double));
    for (size_t k = 0; k < s->steps; k++) {
        if (s->step_kind[k] != (transpose ? 2 : 1))
            continue;
        size_t r = s->step_row[k], c = s->step_col[k];
        double v = work[transpose ? c : r] / s->step_pivot[k];
        if (!isfinite(v)) {
            free(work);
            return -1;
        }
        x[transpose ? r : c] = v;
        if (!transpose) {
            for (size_t e = s->col_start[c]; e < s->col_start[c + 1]; e++) {
                size_t i = s->row_index[e];
                if (i != r)
                    work[i] = fma(-s->value[e], v, work[i]);
            }
        } else {
            for (size_t e = 0; e < s->row[r].count; e++) {
                SparseEntry a = s->row[r].entry[e];
                if (a.column != c)
                    work[a.column] = fma(-a.value, v, work[a.column]);
            }
        }
    }
    double *b = lab_alloc(s->core_n, sizeof(double));
    double *solution = lab_alloc(s->core_n, sizeof(double));
    for (size_t i = 0; i < s->core_n; i++)
        b[i] = work[transpose ? s->core_col[i] : s->core_row[i]];
    int bad = s->core_n && lab_factor_solve(&s->core, b, solution, transpose);
    /* Eliminating exact unit equations can hide forward error in a small
     * ill-conditioned core even when whole-basis backward error is tiny.
     * Correct that core using compensated residuals before recovering slacks. */
    if (!bad && s->core_n && f->core_refinement) {
        double *r = lab_alloc(s->core_n, sizeof(double));
        double *next_r = lab_alloc(s->core_n, sizeof(double));
        double *delta = lab_alloc(s->core_n, sizeof(double));
        double *trial = lab_alloc(s->core_n, sizeof(double));
        double error;
        bad = core_residual(s, b, solution, r, transpose, &error);
        for (size_t pass = 0; pass < 2 && !bad && error > 0; pass++) {
            if (lab_factor_solve(&s->core, r, delta, transpose))
                break;
            for (size_t i = 0; i < s->core_n; i++)
                trial[i] = solution[i] + delta[i];
            double next_error;
            if (core_residual(s, b, trial, next_r, transpose, &next_error) || next_error >= error)
                break;
            memcpy(solution, trial, s->core_n * sizeof(double));
            memcpy(r, next_r, s->core_n * sizeof(double));
            error = next_error;
        }
        free(r); free(next_r); free(delta); free(trial);
    }
    if (!bad)
        for (size_t i = 0; i < s->core_n; i++)
            x[transpose ? s->core_row[i] : s->core_col[i]] = solution[i];
    for (size_t k = s->steps; k-- > 0 && !bad;) {
        if (s->step_kind[k] != (transpose ? 1 : 2))
            continue;
        size_t r = s->step_row[k], c = s->step_col[k];
        double v = work[transpose ? c : r];
        if (!transpose) {
            for (size_t e = 0; e < s->row[r].count; e++) {
                SparseEntry a = s->row[r].entry[e];
                if (a.column != c && s->col_kind[a.column] != 1)
                    v = fma(-a.value, x[a.column], v);
            }
        } else {
            for (size_t e = s->col_start[c]; e < s->col_start[c + 1]; e++) {
                size_t i = s->row_index[e];
                if (i != r && s->row_kind[i] != 2)
                    v = fma(-s->value[e], x[i], v);
            }
        }
        v /= s->step_pivot[k];
        if (!isfinite(v))
            bad = 1;
        else
            x[transpose ? r : c] = v;
    }
    free(work); free(b); free(solution);
    return bad ? -1 : 0;
}
static int sparse_solve(const LabFactor *f, const double *rhs, double *x, int transpose) {
    const SparseFactor *s = f->sparse_data;
    size_t n = f->n;
    if (s->steps)
        return peeled_solve(f, rhs, x, transpose);
    if (!transpose) {
        for (size_t k = 0; k < n; k++) {
            size_t p = f->pivot[k];
            double t = x[k];
            x[k] = x[p];
            x[p] = t;
        }
        for (size_t i = 0; i < n; i++)
            for (size_t e = 0; e < s->row[i].count; e++) {
                SparseEntry a = s->row[i].entry[e];
                if (a.column >= i)
                    break;
                x[i] -= a.value * x[a.column];
            }
        for (size_t i = n; i-- > 0;) {
            for (size_t e = s->row[i].count; e-- > 0;) {
                SparseEntry a = s->row[i].entry[e];
                if (a.column <= i)
                    break;
                x[i] -= a.value * x[a.column];
            }
            x[i] /= s->diagonal[i];
        }
    } else {
        for (size_t i = 0; i < n; i++) {
            for (size_t e = s->col_start[i]; e < s->col_start[i + 1]; e++) {
                size_t row = s->row_index[e];
                if (row >= i)
                    break;
                x[i] -= s->value[e] * x[row];
            }
            x[i] /= s->diagonal[i];
        }
        for (size_t i = n; i-- > 0;)
            for (size_t e = s->col_start[i + 1]; e-- > s->col_start[i];) {
                size_t row = s->row_index[e];
                if (row <= i)
                    break;
                x[i] -= s->value[e] * x[row];
            }
        for (size_t k = n; k-- > 0;) {
            size_t p = f->pivot[k];
            double t = x[k];
            x[k] = x[p];
            x[p] = t;
        }
    }
    (void)rhs;
    return 0;
}
int lab_factor_solve(const LabFactor *f, const double *rhs, double *x, int transpose) {
    size_t n = f->n;
    if (!f->valid || (n && (!rhs || !x)))
        return -1;
    for (size_t i = 0; i < n; i++)
        if (!isfinite(rhs[i]))
            return -1;
    if (n)
        memmove(x, rhs, n * sizeof(double));
    if (f->sparse) {
        if (sparse_solve(f, rhs, x, transpose))
            return -1;
    } else if (!transpose) {
        for (size_t k = 0; k < n; k++) {
            size_t p = f->pivot[k];
            double t = x[k];
            x[k] = x[p];
            x[p] = t;
        }
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < i; j++)
                x[i] -= f->lu[i * n + j] * x[j];
        for (size_t i = n; i-- > 0;) {
            for (size_t j = i + 1; j < n; j++)
                x[i] -= f->lu[i * n + j] * x[j];
            x[i] /= f->lu[i * n + i];
        }
    } else {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < i; j++)
                x[i] -= f->lu[j * n + i] * x[j];
            x[i] /= f->lu[i * n + i];
        }
        for (size_t i = n; i-- > 0;)
            for (size_t j = i + 1; j < n; j++)
                x[i] -= f->lu[j * n + i] * x[j];
        for (size_t k = n; k-- > 0;) {
            size_t p = f->pivot[k];
            double t = x[k];
            x[k] = x[p];
            x[p] = t;
        }
    }
    for (size_t i = 0; i < n; i++)
        if (!isfinite(x[i]))
            return -1;
    return 0;
}

#include "factor.h"
#include "common.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

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
} SparseFactor;

/* Both limits are independent of LP/model resource policy. There is no
 * silently selected dense fallback when fill exhausts the sparse budget. */
enum { SPARSE_FILL_BUDGET = 20000000, SPARSE_DIMENSION_LIMIT = 1000000 };

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
static int sparse_solve(const LabFactor *f, const double *rhs, double *x, int transpose) {
    const SparseFactor *s = f->sparse_data;
    size_t n = f->n;
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
        sparse_solve(f, rhs, x, transpose);
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

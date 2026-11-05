#include "rational.h"
#include "common.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static size_t bits(const mpq_t q) {
    size_t p = mpz_sizeinbase(mpq_numref(q), 2);
    size_t d = mpz_sizeinbase(mpq_denref(q), 2);
    return p > d ? p : d;
}
int lab_rational_solve(size_t n, mpq_t *a, mpq_t *rhs, double deadline) {
    if (n > LAB_RATIONAL_MAX_DIMENSION || !isfinite(deadline) || lab_now() >= deadline ||
        (n && (!a || !rhs)))
        return -1;
    for (size_t i = 0; i < n * n; i++)
        if (bits(a[i]) > LAB_RATIONAL_MAX_BITS)
            return -1;
    for (size_t i = 0; i < n; i++)
        if (bits(rhs[i]) > LAB_RATIONAL_MAX_BITS)
            return -1;
    mpq_t factor, term;
    mpq_inits(factor, term, NULL);
    int bad = 0;
    for (size_t k = 0; k < n && !bad; k++) {
        if (lab_now() >= deadline) {
            bad = 1;
            break;
        }
        size_t pivot = SIZE_MAX, cost = SIZE_MAX;
        for (size_t i = k; i < n; i++)
            if (mpq_sgn(a[i * n + k]) && bits(a[i * n + k]) < cost) {
                pivot = i;
                cost = bits(a[i * n + k]);
            }
        if (pivot == SIZE_MAX) {
            bad = 1;
            break;
        }
        if (pivot != k) {
            for (size_t j = k; j < n; j++)
                mpq_swap(a[k * n + j], a[pivot * n + j]);
            mpq_swap(rhs[k], rhs[pivot]);
        }
        for (size_t i = k + 1; i < n && !bad; i++) {
            if (!mpq_sgn(a[i * n + k]))
                continue;
            mpq_div(factor, a[i * n + k], a[k * n + k]);
            mpq_set_ui(a[i * n + k], 0, 1);
            for (size_t j = k + 1; j < n; j++) {
                mpq_mul(term, factor, a[k * n + j]);
                mpq_sub(a[i * n + j], a[i * n + j], term);
                if (bits(a[i * n + j]) > LAB_RATIONAL_MAX_BITS) {
                    bad = 1;
                    break;
                }
            }
            mpq_mul(term, factor, rhs[k]);
            mpq_sub(rhs[i], rhs[i], term);
            if (bits(rhs[i]) > LAB_RATIONAL_MAX_BITS || lab_now() >= deadline)
                bad = 1;
        }
    }
    for (size_t k = n; k-- > 0 && !bad;) {
        for (size_t j = k + 1; j < n; j++) {
            mpq_mul(term, a[k * n + j], rhs[j]);
            mpq_sub(rhs[k], rhs[k], term);
        }
        mpq_div(rhs[k], rhs[k], a[k * n + k]);
        if (bits(rhs[k]) > LAB_RATIONAL_MAX_BITS || lab_now() >= deadline)
            bad = 1;
    }
    mpq_clears(factor, term, NULL);
    return bad ? -1 : 0;
}

static int acceptable(const mpq_t q) {
    return mpz_sgn(mpq_denref(q)) > 0 && bits(q) <= LAB_RATIONAL_MAX_BITS;
}
static mpq_t *new_vector(size_t n) {
    mpq_t *result = lab_alloc(n, sizeof *result);
    for (size_t i = 0; i < n; i++)
        mpq_init(result[i]);
    return result;
}
static void free_vector(mpq_t *v, size_t n) {
    for (size_t i = 0; i < n; i++)
        mpq_clear(v[i]);
    free(v);
}
int lab_rational_sparse_solve(size_t n, const size_t *start, const size_t *row,
                              mpq_t *value, mpq_t *rhs, mpq_t *x, double deadline) {
    if (n > LAB_RATIONAL_SPARSE_MAX_DIMENSION || !isfinite(deadline) ||
        lab_now() >= deadline || (n && (!start || !rhs || !x)) ||
        (start && start[0]) || (n && (rhs == x || value == x)))
        return -1;
    size_t nnz = n ? start[n] : 0;
    if (nnz > LAB_RATIONAL_SPARSE_MAX_NONZEROS || (nnz && (!row || !value)))
        return -1;
    for (size_t j = 0; j < n; j++)
        if (start[j] > start[j + 1] || start[j + 1] > nnz || !acceptable(rhs[j]))
            return -1;
    size_t *seen = lab_alloc(n, sizeof(size_t));
    size_t *row_start = lab_alloc(n + 1, sizeof(size_t));
    size_t *entry_col = lab_alloc(nnz, sizeof(size_t));
    size_t *row_degree = lab_alloc(n, sizeof(size_t));
    size_t *col_degree = lab_alloc(n, sizeof(size_t));
    for (size_t i = 0; i < n; i++)
        seen[i] = SIZE_MAX;
    int bad = 0;
    for (size_t j = 0; j < n && !bad; j++) {
        for (size_t k = start[j]; k < start[j + 1]; k++) {
            if (row[k] >= n || !acceptable(value[k]) || seen[row[k]] == j) {
                bad = 1;
                break;
            }
            seen[row[k]] = j;
            entry_col[k] = j;
            if (mpq_sgn(value[k])) {
                row_degree[row[k]]++;
                col_degree[j]++;
                row_start[row[k] + 1]++;
            }
            if ((k & 127U) == 0 && lab_now() >= deadline) {
                bad = 1;
                break;
            }
        }
    }
    for (size_t i = 0; i < n; i++)
        row_start[i + 1] += row_start[i];
    size_t *row_entry = lab_alloc(row_start[n], sizeof(size_t));
    for (size_t i = 0; i < n; i++)
        seen[i] = row_start[i];
    if (!bad)
        for (size_t k = 0; k < nnz; k++)
            if (mpq_sgn(value[k]))
                row_entry[seen[row[k]]++] = k;
    unsigned char *active_row = lab_alloc(n, 1), *active_col = lab_alloc(n, 1);
    size_t *removed_at = lab_alloc(n, sizeof(size_t));
    size_t *deferred_col = lab_alloc(n, sizeof(size_t));
    size_t *deferred_row = lab_alloc(n, sizeof(size_t));
    size_t *deferred_pivot = lab_alloc(n, sizeof(size_t));
    mpq_t *b = new_vector(n), *answer = new_vector(n);
    for (size_t i = 0; i < n; i++) {
        active_row[i] = active_col[i] = 1;
        removed_at[i] = SIZE_MAX;
        mpq_set(b[i], rhs[i]);
    }
    mpq_t term, total;
    mpq_inits(term, total, NULL);
    size_t remaining = n, step = 0, deferred = 0;
    while (remaining && !bad) {
        if (lab_now() >= deadline) {
            bad = 1;
            break;
        }
        size_t column = SIZE_MAX, equation = SIZE_MAX;
        for (size_t i = 0; i < n; i++) {
            if ((active_row[i] && !row_degree[i]) || (active_col[i] && !col_degree[i])) {
                bad = 1;
                break;
            }
            if (active_col[i] && col_degree[i] == 1 && column == SIZE_MAX)
                column = i;
            if (active_row[i] && row_degree[i] == 1 && equation == SIZE_MAX)
                equation = i;
        }
        if (bad || (column == SIZE_MAX && equation == SIZE_MAX))
            break;
        if (column != SIZE_MAX) {
            /* This variable appears only in the removed equation. Its value
             * is deferred until all variables still in that equation are known. */
            size_t pivot = SIZE_MAX;
            for (size_t k = start[column]; k < start[column + 1]; k++)
                if (active_row[row[k]] && mpq_sgn(value[k])) {
                    pivot = k;
                    equation = row[k];
                    break;
                }
            if (pivot == SIZE_MAX) {
                bad = 1;
                break;
            }
            deferred_col[deferred] = column;
            deferred_row[deferred] = equation;
            deferred_pivot[deferred++] = pivot;
            active_row[equation] = active_col[column] = 0;
            removed_at[column] = step++;
            for (size_t p = row_start[equation]; p < row_start[equation + 1]; p++) {
                size_t j = entry_col[row_entry[p]];
                if (active_col[j])
                    col_degree[j]--;
            }
        } else {
            /* A singleton equation gives this variable directly. Substitute
             * it only into equations that are still active. */
            size_t pivot = SIZE_MAX;
            for (size_t p = row_start[equation]; p < row_start[equation + 1]; p++) {
                size_t k = row_entry[p];
                if (active_col[entry_col[k]]) {
                    pivot = k;
                    column = entry_col[k];
                    break;
                }
            }
            if (pivot == SIZE_MAX) {
                bad = 1;
                break;
            }
            mpq_div(answer[column], b[equation], value[pivot]);
            if (!acceptable(answer[column])) {
                bad = 1;
                break;
            }
            active_row[equation] = active_col[column] = 0;
            removed_at[column] = step++;
            for (size_t k = start[column]; k < start[column + 1]; k++) {
                size_t i = row[k];
                if (!active_row[i] || !mpq_sgn(value[k]))
                    continue;
                mpq_mul(term, value[k], answer[column]);
                mpq_sub(b[i], b[i], term);
                row_degree[i]--;
                if (!acceptable(term) || !acceptable(b[i])) {
                    bad = 1;
                    break;
                }
            }
        }
        remaining--;
    }
    if (remaining > LAB_RATIONAL_MAX_DIMENSION)
        bad = 1;
    if (!bad && remaining) {
        mpq_t *a = new_vector(remaining * remaining), *core_rhs = new_vector(remaining);
        size_t *columns = lab_alloc(remaining, sizeof(size_t));
        size_t r = 0, c = 0;
        for (size_t i = 0; i < n; i++)
            if (active_row[i]) {
                seen[i] = r;
                mpq_set(core_rhs[r++], b[i]);
            }
        for (size_t j = 0; j < n; j++)
            if (active_col[j]) {
                columns[c] = j;
                for (size_t k = start[j]; k < start[j + 1]; k++)
                    if (active_row[row[k]])
                        mpq_set(a[seen[row[k]] * remaining + c], value[k]);
                c++;
            }
        if (r != remaining || c != remaining ||
            lab_rational_solve(remaining, a, core_rhs, deadline))
            bad = 1;
        if (!bad)
            for (size_t j = 0; j < remaining; j++)
                mpq_set(answer[columns[j]], core_rhs[j]);
        free(columns);
        free_vector(a, remaining * remaining);
        free_vector(core_rhs, remaining);
    }
    for (size_t p = deferred; p-- > 0 && !bad;) {
        size_t column = deferred_col[p], equation = deferred_row[p];
        mpq_set(total, b[equation]);
        for (size_t at = row_start[equation]; at < row_start[equation + 1]; at++) {
            size_t k = row_entry[at], j = entry_col[k];
            if (j == column || removed_at[j] < removed_at[column])
                continue;
            mpq_mul(term, value[k], answer[j]);
            mpq_sub(total, total, term);
            if (!acceptable(term) || !acceptable(total)) {
                bad = 1;
                break;
            }
        }
        if (!bad) {
            mpq_div(answer[column], total, value[deferred_pivot[p]]);
            if (!acceptable(answer[column]) || lab_now() >= deadline)
                bad = 1;
        }
    }
    /* Recheck every original equation exactly; no success is inferred merely
     * from the bookkeeping of the peeled system. */
    if (!bad) {
        for (size_t i = 0; i < n; i++)
            mpq_neg(b[i], rhs[i]);
        for (size_t j = 0; j < n && !bad; j++) {
            for (size_t k = start[j]; k < start[j + 1]; k++) {
                mpq_mul(term, value[k], answer[j]);
                mpq_add(b[row[k]], b[row[k]], term);
                if (!acceptable(term) || !acceptable(b[row[k]])) {
                    bad = 1;
                    break;
                }
            }
            if (lab_now() >= deadline)
                bad = 1;
        }
        for (size_t i = 0; i < n && !bad; i++)
            if (mpq_sgn(b[i]))
                bad = 1;
    }
    if (!bad)
        for (size_t i = 0; i < n; i++)
            mpq_set(x[i], answer[i]);
    mpq_clears(term, total, NULL);
    free_vector(b, n);
    free_vector(answer, n);
    free(seen);
    free(row_start);
    free(entry_col);
    free(row_degree);
    free(col_degree);
    free(row_entry);
    free(active_row);
    free(active_col);
    free(removed_at);
    free(deferred_col);
    free(deferred_row);
    free(deferred_pivot);
    return bad ? -1 : 0;
}

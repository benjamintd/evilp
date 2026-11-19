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
/* Fraction-free elimination. Rational row denominators are cleared once;
 * exact integer divisions replace repeated rational GCD reduction. A failed
 * cap or divisibility check falls back to the original bounded solver. */
static int rational_fraction_free(size_t n, mpq_t *a, mpq_t *rhs, double deadline) {
    size_t width = n + 1;
    mpz_t *z = lab_alloc(n * width, sizeof(mpz_t));
    for (size_t i = 0; i < n * width; i++)
        mpz_init(z[i]);
    mpz_t lcm, factor, previous, pivot_value, term, numerator;
    mpz_inits(lcm, factor, previous, pivot_value, term, numerator, NULL);
    mpz_set_ui(previous, 1);
    int bad = 0;
    for (size_t i = 0; i < n && !bad; i++) {
        mpz_set_ui(lcm, 1);
        for (size_t j = 0; j <= n; j++) {
            mpq_srcptr q = j == n ? rhs[i] : a[i * n + j];
            mpz_lcm(lcm, lcm, mpq_denref(q));
            if (mpz_sizeinbase(lcm, 2) > LAB_RATIONAL_MAX_BITS || lab_now() >= deadline) {
                bad = 1;
                break;
            }
        }
        for (size_t j = 0; j <= n && !bad; j++) {
            mpq_srcptr q = j == n ? rhs[i] : a[i * n + j];
            mpz_divexact(factor, lcm, mpq_denref(q));
            mpz_mul(z[i * width + j], factor, mpq_numref(q));
            if (mpz_sizeinbase(z[i * width + j], 2) > LAB_RATIONAL_MAX_BITS)
                bad = 1;
        }
    }
    for (size_t k = 0; k < n && !bad; k++) {
        size_t pivot = SIZE_MAX, cost = SIZE_MAX;
        for (size_t i = k; i < n; i++)
            if (mpz_sgn(z[i * width + k]) && mpz_sizeinbase(z[i * width + k], 2) < cost) {
                pivot = i;
                cost = mpz_sizeinbase(z[i * width + k], 2);
            }
        if (pivot == SIZE_MAX || lab_now() >= deadline) {
            bad = 1;
            break;
        }
        if (pivot != k)
            for (size_t j = k; j <= n; j++)
                mpz_swap(z[k * width + j], z[pivot * width + j]);
        mpz_set(pivot_value, z[k * width + k]);
        for (size_t i = k + 1; i < n && !bad; i++) {
            for (size_t j = k + 1; j <= n; j++) {
                mpz_mul(numerator, pivot_value, z[i * width + j]);
                mpz_mul(term, z[i * width + k], z[k * width + j]);
                mpz_sub(numerator, numerator, term);
                if (!mpz_divisible_p(numerator, previous)) {
                    bad = 1;
                    break;
                }
                mpz_divexact(z[i * width + j], numerator, previous);
                if (mpz_sizeinbase(z[i * width + j], 2) > LAB_RATIONAL_MAX_BITS ||
                    lab_now() >= deadline) {
                    bad = 1;
                    break;
                }
            }
            mpz_set_ui(z[i * width + k], 0);
        }
        mpz_set(previous, pivot_value);
    }
    mpq_t *answer = lab_alloc(n, sizeof(mpq_t));
    for (size_t i = 0; i < n; i++)
        mpq_init(answer[i]);
    mpq_t coefficient, product, total;
    mpq_inits(coefficient, product, total, NULL);
    for (size_t i = n; i-- > 0 && !bad;) {
        mpq_set_z(answer[i], z[i * width + n]);
        for (size_t j = i + 1; j < n; j++) {
            mpq_set_z(coefficient, z[i * width + j]);
            mpq_mul(product, coefficient, answer[j]);
            mpq_sub(answer[i], answer[i], product);
        }
        mpq_set_z(coefficient, z[i * width + i]);
        mpq_div(answer[i], answer[i], coefficient);
        if (bits(answer[i]) > LAB_RATIONAL_MAX_BITS || lab_now() >= deadline)
            bad = 1;
    }
    /* Verify against the untouched rational input, including RHS. */
    for (size_t i = 0; i < n && !bad; i++) {
        mpq_set_ui(total, 0, 1);
        for (size_t j = 0; j < n; j++) {
            mpq_mul(product, a[i * n + j], answer[j]);
            mpq_add(total, total, product);
            if (bits(total) > LAB_RATIONAL_MAX_BITS || lab_now() >= deadline) {
                bad = 1;
                break;
            }
        }
        if (mpq_cmp(total, rhs[i]))
            bad = 1;
    }
    if (!bad)
        for (size_t i = 0; i < n; i++)
            mpq_set(rhs[i], answer[i]);
    for (size_t i = 0; i < n; i++)
        mpq_clear(answer[i]);
    free(answer);
    mpq_clears(coefficient, product, total, NULL);
    mpz_clears(lcm, factor, previous, pivot_value, term, numerator, NULL);
    for (size_t i = 0; i < n * width; i++)
        mpz_clear(z[i]);
    free(z);
    return bad ? -1 : 0;
}

static int solve_dense(size_t n, mpq_t *a, mpq_t *rhs, double deadline) {
    if (!isfinite(deadline) || lab_now() >= deadline ||
        (n && (!a || !rhs)))
        return -1;
    for (size_t i = 0; i < n * n; i++)
        if (bits(a[i]) > LAB_RATIONAL_MAX_BITS)
            return -1;
    for (size_t i = 0; i < n; i++)
        if (bits(rhs[i]) > LAB_RATIONAL_MAX_BITS)
            return -1;
    if (!rational_fraction_free(n, a, rhs, deadline))
        return 0;
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

int lab_rational_solve(size_t n, mpq_t *a, mpq_t *rhs, double deadline) {
    return n <= LAB_RATIONAL_MAX_DIMENSION ? solve_dense(n, a, rhs, deadline) : -1;
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
    /* Degrees only decrease. A row or column enters its queue at most
     * once, when it becomes a singleton. Rows consumed by column pivots
     * (and conversely) are discarded lazily. */
    size_t *row_queue = lab_alloc(n, sizeof *row_queue);
    size_t *col_queue = lab_alloc(n, sizeof *col_queue);
    size_t row_head = 0, row_tail = 0, col_head = 0, col_tail = 0;
    for (size_t i = 0; i < n; ++i) {
        if (!row_degree[i] || !col_degree[i])
            bad = 1;
        if (row_degree[i] == 1)
            row_queue[row_tail++] = i;
        if (col_degree[i] == 1)
            col_queue[col_tail++] = i;
    }
    size_t remaining = n, step = 0, deferred = 0;
    while (remaining && !bad) {
        if (lab_now() >= deadline) {
            bad = 1;
            break;
        }
        size_t column = SIZE_MAX, equation = SIZE_MAX;
        while (col_head < col_tail && !active_col[col_queue[col_head]])
            ++col_head;
        if (col_head < col_tail) {
            column = col_queue[col_head++];
        } else {
            while (row_head < row_tail && !active_row[row_queue[row_head]])
                ++row_head;
            if (row_head == row_tail)
                break;
            equation = row_queue[row_head++];
        }
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
                if (active_col[j]) {
                    if (--col_degree[j] == 1)
                        col_queue[col_tail++] = j;
                    else if (!col_degree[j])
                        bad = 1;
                }
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
                if (--row_degree[i] == 1)
                    row_queue[row_tail++] = i;
                else if (!row_degree[i])
                    bad = 1;
                if (!acceptable(term) || !acceptable(b[i])) {
                    bad = 1;
                    break;
                }
            }
        }
        remaining--;
    }
    if (remaining > LAB_RATIONAL_SPARSE_MAX_CORE)
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
            solve_dense(remaining, a, core_rhs, deadline))
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
    free(col_queue);
    free(row_queue);
    free(deferred_pivot);
    return bad ? -1 : 0;
}

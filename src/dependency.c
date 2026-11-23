#include "dependency.h"
#include "common.h"
#include "exact.h"
#include "rational.h"
#include <math.h>
#include <stdlib.h>

static int bounded(const mpq_t q) {
    return mpz_sizeinbase(mpq_numref(q), 2) <= LAB_RATIONAL_MAX_BITS &&
           mpz_sizeinbase(mpq_denref(q), 2) <= LAB_RATIONAL_MAX_BITS;
}
int lab_source_row_dependency(const LabModel *m, const double *y, double deadline) {
    if (!m || !y || m->rows > 4096 || m->cols > 16384 || m->nnz > 100000 ||
        lab_now() >= deadline)
        return 0;
    mpq_t sum, a, b, q, term;
    mpq_inits(sum, a, b, q, term, NULL);
    int valid = 1, nonzero = 0;
    for (size_t i = 0; i < m->rows && valid; i++) {
        if (!isfinite(y[i])) {
            valid = 0;
            break;
        }
        if (y[i] == 0)
            continue;
        nonzero = 1;
        if (!isfinite(m->row_lower[i]) || !isfinite(m->row_upper[i]) ||
            lab_source_exact_row_bound(m, i, 0, a) ||
            lab_source_exact_row_bound(m, i, 1, b) || mpq_cmp(a, b)) {
            valid = 0;
            break;
        }
        mpq_set_d(q, y[i]);
        mpq_mul(term, a, q);
        mpq_add(sum, sum, term);
        if (!bounded(term) || !bounded(sum) || lab_now() >= deadline)
            valid = 0;
    }
    valid &= nonzero && mpq_sgn(sum) == 0;
    for (size_t j = 0; j < m->cols && valid; j++) {
        mpq_set_ui(sum, 0, 1);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            size_t i = m->row_index[k];
            if (i >= m->rows) {
                valid = 0;
                break;
            }
            if (y[i] == 0)
                continue;
            if (lab_source_exact_coefficient(m, k, a)) {
                valid = 0;
                break;
            }
            mpq_set_d(q, y[i]);
            mpq_mul(term, a, q);
            mpq_add(sum, sum, term);
            if (!bounded(term) || !bounded(sum) || lab_now() >= deadline) {
                valid = 0;
                break;
            }
        }
        if (mpq_sgn(sum) || lab_now() >= deadline)
            valid = 0;
    }
    mpq_clears(sum, a, b, q, term, NULL);
    return valid;
}

/* A nonsingular equality matrix fixes every participating variable. A
 * disconnected zero-cost column can take zero if its source bounds allow it.
 * Keep the original CSC entries and remove only empty column pointers. */
int lab_recover_square_equalities(const LabModel *m, double *point, double *objective,
                                  double deadline) {
    const LabExactModel *e = m->exact;
    size_t n = m->rows;
    if (!e || !n || n > LAB_RATIONAL_MAX_DIMENSION || lab_now() >= deadline)
        return 0;
    for (size_t i = 0; i < n; ++i)
        if (!isfinite(m->row_lower[i]) || !isfinite(m->row_upper[i]) ||
            mpq_cmp(e->row_lower[i], e->row_upper[i]))
            return 0;

    size_t start[LAB_RATIONAL_MAX_DIMENSION + 1], column[LAB_RATIONAL_MAX_DIMENSION];
    size_t count = 0;
    for (size_t j = 0; j < m->cols; ++j) {
        if (m->col_start[j] == m->col_start[j + 1]) {
            if (mpq_sgn(e->objective[j]) ||
                (isfinite(m->col_lower[j]) && mpq_sgn(e->col_lower[j]) > 0) ||
                (isfinite(m->col_upper[j]) && mpq_sgn(e->col_upper[j]) < 0))
                return 0;
            point[j] = 0;
        } else {
            if (count == n)
                return 0;
            start[count] = m->col_start[j];
            column[count++] = j;
        }
    }
    if (count != n)
        return 0;
    start[n] = m->nnz;
    mpq_t *x = lab_alloc(n, sizeof *x);
    for (size_t i = 0; i < n; ++i)
        mpq_init(x[i]);
    int feasible = !lab_rational_sparse_solve(n, start, m->row_index, e->value,
                                               e->row_lower, x, deadline);
    int representable = 1;
    mpq_t binary;
    mpq_init(binary);
    for (size_t i = 0; i < n && feasible; ++i) {
        size_t j = column[i];
        if ((m->integer[j] && mpz_cmp_ui(mpq_denref(x[i]), 1)) ||
            (isfinite(m->col_lower[j]) && mpq_cmp(x[i], e->col_lower[j]) < 0) ||
            (isfinite(m->col_upper[j]) && mpq_cmp(x[i], e->col_upper[j]) > 0))
            feasible = 0;
        if (lab_exact_nearest(x[i], &point[j])) {
            representable = 0;
        } else {
            mpq_set_d(binary, point[j]);
            if (mpq_cmp(binary, x[i]))
                representable = 0;
        }
    }
    mpq_clear(binary);
    for (size_t i = 0; i < n; ++i)
        mpq_clear(x[i]);
    free(x);
    if (!feasible)
        return 0;
    double residual;
    if (!representable || lab_exact_evaluate(m, point, objective, &residual) || residual != 0)
        return -1;
    return 1;
}

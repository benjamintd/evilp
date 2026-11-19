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

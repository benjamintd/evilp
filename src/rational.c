#include "rational.h"
#include "common.h"
#include <math.h>
#include <stdint.h>

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

#include "factor.h"
#include "common.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void lab_factor_free(LabFactor *f) {
    free(f->lu);
    free(f->pivot);
    memset(f, 0, sizeof *f);
}
int lab_factorize(LabFactor *f, const double *matrix, size_t n) {
    lab_factor_free(f);
    if (n && (n > SIZE_MAX / n || n * n > SIZE_MAX / sizeof(double)))
        return -1;
    f->n = n;
    f->lu = lab_alloc(n * n, sizeof(double));
    f->pivot = lab_alloc(n, sizeof(size_t));
    memcpy(f->lu, matrix, n * n * sizeof(double));
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
    for (size_t i = 0; i < n * n; i++)
        if (!isfinite(f->lu[i]))
            return -1;
    return 0;
}
int lab_factor_solve(const LabFactor *f, const double *rhs, double *x, int transpose) {
    size_t n = f->n;
    memcpy(x, rhs, n * sizeof(double));
    if (!transpose) {
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

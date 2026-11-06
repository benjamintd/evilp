#include "pivot.h"
#include <math.h>
#include <stdint.h>
#include <float.h>

int lab_primal_ratio_scaled(size_t n, const double *xb, const double *direction,
                           const size_t *basis, const double *tolerance,
                           int bland, size_t *leaving) {
    if (!leaving || (n && (!xb || !direction || !basis || !tolerance)))
        return -1;
    *leaving = SIZE_MAX;
    long double permitted = INFINITY, chosen = INFINITY;
    double strongest = 0;
    int positive = 0;
    for (size_t i = 0; i < n; ++i) {
        double d = direction[i], t = tolerance[i];
        if (!isfinite(d) || !isfinite(xb[i]) || xb[i] < 0 ||
            !isfinite(t) || t < 0 || t > LAB_PRIMAL_TOLERANCE)
            return -1;
        if (d > 0) {
            positive = 1;
            long double relaxed = ((long double)xb[i] + .75L * t) / d;
            if (relaxed < permitted)
                permitted = relaxed;
        }
    }
    for (size_t i = 0; i < n; ++i) {
        double d = direction[i];
        if (d <= 1e-12)
            continue;
        long double ratio = (long double)xb[i] / d;
        if (!isfinite(ratio) || ratio > permitted)
            continue;
        if (bland) {
            /* No dimensionless floor: a 1e-12 absolute tie can dwarf a
             * genuine 1e-15 difference between two tiny limiting steps. */
            long double tie = 64 * LDBL_EPSILON * fmaxl(fabsl(ratio), fabsl(chosen));
            if (*leaving == SIZE_MAX || ratio < chosen - tie ||
                (fabsl(ratio - chosen) <= tie && basis[i] < basis[*leaving])) {
                *leaving = i;
                chosen = ratio;
            }
        } else if (*leaving == SIZE_MAX || d > strongest ||
                   (d == strongest && basis[i] < basis[*leaving])) {
            *leaving = i;
            chosen = ratio;
            strongest = d;
        }
    }
    if (*leaving == SIZE_MAX)
        return positive ? -1 : 1;
    /* Use the actual binary step that the floating basis update predicts.
     * Rounded division must also respect the same per-basic headroom. */
    double theta = xb[*leaving] / direction[*leaving];
    if (!isfinite(theta)) {
        *leaving = SIZE_MAX;
        return -1;
    }
    for (size_t i = 0; i < n; ++i)
        if (direction[i] > 0 &&
            fmal(-(long double)theta, direction[i], xb[i]) < -.75L * tolerance[i]) {
            *leaving = SIZE_MAX;
            return -1;
        }
    return 0;
}

int lab_primal_ratio(size_t n, const double *xb, const double *direction,
                     const size_t *basis, int bland, size_t *leaving) {
    if (!leaving || (n && (!xb || !direction || !basis)))
        return -1;
    *leaving = SIZE_MAX;
    double limit = INFINITY, strongest = 0;
    int positive = 0;
    for (size_t i = 0; i < n; i++) {
        double d = direction[i];
        if (!isfinite(d) || !isfinite(xb[i]) || xb[i] < 0)
            return -1;
        if (d <= 0)
            continue;
        positive = 1;
        if (bland) {
            if (d <= 1e-12)
                continue;
            double ratio = xb[i] / d;
            if (!isfinite(ratio))
                return -1;
            double tie = 1e-12 * (1 + fabs(ratio));
            if (*leaving == SIZE_MAX || ratio < limit - tie ||
                (fabs(ratio - limit) <= tie && basis[i] < basis[*leaving])) {
                limit = ratio;
                *leaving = i;
            }
        } else {
            /* Even too-small positive pivots constrain the permitted step.
             * Discarding them here can step through another row's bound. */
            double relaxed = (xb[i] + .75 * LAB_PRIMAL_TOLERANCE) / d;
            limit = fmin(limit, relaxed);
        }
    }
    if (!bland)
        for (size_t i = 0; i < n; i++) {
            double d = direction[i];
            if (d <= 1e-12)
                continue;
            double ratio = xb[i] / d;
            if (!isfinite(ratio) || ratio > limit)
                continue;
            if (*leaving == SIZE_MAX || d > strongest ||
                (d == strongest && basis[i] < basis[*leaving])) {
                strongest = d;
                *leaving = i;
            }
        }
    if (!bland && *leaving != SIZE_MAX) {
        double theta = xb[*leaving] / direction[*leaving];
        for (size_t i = 0; i < n; i++)
            if (direction[i] > 0 &&
                fma(-theta, direction[i], xb[i]) < -.75 * LAB_PRIMAL_TOLERANCE) {
                *leaving = SIZE_MAX;
                return -1;
            }
    }
    if (*leaving != SIZE_MAX)
        return 0;
    return positive ? -1 : 1;
}

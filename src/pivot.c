#include "pivot.h"
#include <math.h>
#include <stdint.h>

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

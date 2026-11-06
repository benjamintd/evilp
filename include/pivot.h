#ifndef LAB_PIVOT_H
#define LAB_PIVOT_H
#include <stddef.h>
/* Existing canonical primal acceptance tolerance. Harris uses 75% of this
 * allowance so refactorization retains headroom for solve roundoff. */
#define LAB_PRIMAL_TOLERANCE 1e-8
/* 0: chosen row; 1: no positive direction (unbounded); -1: unsafe/invalid.
 * Bland mode preserves the earlier ratio rule for controlled comparisons. */
int lab_primal_ratio(size_t n, const double *xb, const double *direction,
                     const size_t *basis, int bland, size_t *leaving);
/* Per-basic feasibility allowances in the same units as xb. All positive
 * directions constrain the step, including directions too weak to pivot on.
 * Bland ties are relative to the step size and remain inside every row's
 * allowance. The legacy entry point above remains a comparison control. */
int lab_primal_ratio_scaled(size_t n, const double *xb, const double *direction,
                           const size_t *basis, const double *tolerance,
                           int bland, size_t *leaving);
#endif

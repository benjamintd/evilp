#ifndef LAB_RATIONAL_H
#define LAB_RATIONAL_H
#include <gmp.h>
#include <stddef.h>

/* A deliberately bounded exact reconstruction, separate from normal sparse
 * basis solves. A and rhs are initialized, row-major rationals, overwritten.
 * On success rhs contains the exact solution. Failure is never a proof. */
#define LAB_RATIONAL_MAX_DIMENSION 128
#define LAB_RATIONAL_MAX_BITS 8192
int lab_rational_solve(size_t n, mpq_t *a, mpq_t *rhs, double deadline);
#endif

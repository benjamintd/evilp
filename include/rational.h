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
/* Read-only canonical-rational CSC input, distinct initialized output x.
 * Singleton peeling precedes a dense remaining core bounded separately below. Duplicate
 * row entries within a column are rejected. x is unchanged on failure. */
#define LAB_RATIONAL_SPARSE_MAX_DIMENSION 8192
#define LAB_RATIONAL_SPARSE_MAX_CORE 512
#define LAB_RATIONAL_SPARSE_MAX_NONZEROS 100000
int lab_rational_sparse_solve(size_t n, const size_t *start, const size_t *row,
                              mpq_t *value, mpq_t *rhs, mpq_t *x, double deadline);
#endif

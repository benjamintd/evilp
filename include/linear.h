#ifndef LAB_LINEAR_H
#define LAB_LINEAR_H
#include "factor.h"

/* Equilibrated basis solve. Original CSC coefficients remain immutable and are
 * retained for compensated residuals; sparse mode never builds a dense matrix.
 * Start from zero and free after success or failure. */
typedef struct {
    LabFactor factor;
    size_t n, nnz;
    size_t *col_start, *row_index;
    double *value;
    int *row_exponent, *col_exponent;
    int valid, scaling_enabled, refinement_enabled, scaling_changed;
    size_t scaled_rows, scaled_columns;
    /* Cumulative solve/correction counts; errors describe the latest call. */
    size_t solve_calls, refinement_steps;
    double initial_backward_error, max_backward_error;
    char error[192];
} LabLinear;

int lab_linear_factorize(LabLinear *s, size_t n, const size_t *col_start, const size_t *row_index,
                         const double *value, int dense);
int lab_linear_factorize_ex(LabLinear *s, size_t n, const size_t *col_start,
                            const size_t *row_index, const double *value, int dense, int scaling,
                            int refinement);
/* Solves B*x=rhs or B^T*x=rhs in original units, including in-place RHS.
 * Success requires original componentwise backward error <=1e-10. Refinement
 * targets 1e-14 with at most eight corrections and returns the best iterate. */
int lab_linear_solve(LabLinear *s, const double *rhs, double *x, int transpose);
void lab_linear_free(LabLinear *s);
#endif

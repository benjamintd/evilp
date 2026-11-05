#ifndef LAB_LINEAR_H
#define LAB_LINEAR_H
#include "factor.h"

#define LAB_LINEAR_UPDATE_LIMIT 32

/* Equilibrated basis solve. Caller CSC inputs remain immutable. Current
 * original-unit coefficients are retained for compensated residuals; sparse
 * mode never builds a dense matrix. Start from zero and free after either
 * success or failure. */
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
    /* B_current = B_base E_0 ... E_(update_count-1). Each eta replaces one
     * column of the identity and is stored in original basis coordinates. */
    size_t update_count, accepted_updates, rejected_updates;
    size_t *eta_column;
    double **eta_value;
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
/* Replace one current basis column. Sparse input requires unique row indices
 * (arbitrary order, explicit zeros allowed); dense input has n entries. An
 * optional direction is B_current^-1 * new_column in ORIGINAL units. Supplied
 * directions are independently checked; NULL computes it with guarded solves.
 * Return 0 when accepted, 1 when a full refactorization is required (including
 * the bounded update limit or a weak/inaccurate pivot), -1 on invalid input.
 * Every nonzero return leaves the old basis and factorization usable. The
 * caller may refactor its NEW basis directly with lab_linear_factorize_ex. */
int lab_linear_update(LabLinear *s, size_t column, size_t nnz, const size_t *row,
                      const double *value, const double *direction);
int lab_linear_update_dense(LabLinear *s, size_t column, const double *new_column,
                            const double *direction);
void lab_linear_free(LabLinear *s);
#endif

#ifndef LAB_LINEAR_H
#define LAB_LINEAR_H
#include "factor.h"

#define LAB_LINEAR_UPDATE_LIMIT 32
#define LAB_LINEAR_DIMENSION_LIMIT 8000000

/* Explicit owned copy of the current original-unit basis. Initialize to zero;
 * free with lab_linear_csc_free. Borrowed column views require no allocation. */
typedef struct {
    size_t n, nnz;
    size_t *col_start, *row_index;
    double *value;
} LabLinearCSC;
typedef struct {
    size_t column, nnz;
    size_t *row_index;
    double *value;
} LabLinearColumn;
typedef struct {
    size_t column, count;
    size_t *row_index;
    double *value;
    double pivot;
    int sparse;
} LabLinearEta;

/* Equilibrated basis solve. Caller CSC inputs remain immutable. Current
 * original-unit coefficients are retained for compensated residuals; sparse
 * mode never builds a dense matrix. Start from zero and free after either
 * success or failure. */
typedef struct {
    LabFactor factor;
    size_t n, nnz, base_nnz;
    /* Immutable original-unit base CSC. These are deliberately named base_*:
     * use the accessors below for the current matrix after column updates. */
    size_t *base_col_start, *base_row_index;
    double *base_value;
    size_t overlay_count;
    unsigned char *overlay_slot; /* zero: base column; otherwise slot+1 */
    LabLinearColumn *overlay;
    int *row_exponent, *col_exponent;
    int valid, scaling_enabled, refinement_enabled, scaling_changed;
    size_t scaled_rows, scaled_columns;
    /* Cumulative solve/correction counts; errors describe the latest call. */
    size_t solve_calls, refinement_steps;
    /* B_current = B_base E_0 ... E_(update_count-1). Each eta replaces one
     * column of the identity and is stored in original basis coordinates. */
    size_t update_count, accepted_updates, rejected_updates;
    /* Sparse etas retain every nonzero off-diagonal coefficient, in row
     * order. Dense etas retain the original full vector. The pivot is never
     * dropped. Packed payloads use at most half a dense vector's bytes. */
    LabLinearEta *eta;
    /* Work counters since factorization. Updates copy only replacement-column
     * entries; current CSC materialization is explicit and never part of a
     * solve/update. Entry counters saturate if accumulated usage overflows. */
    size_t column_entries_copied, csc_materializations, csc_materialized_entries;
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
 * the 32-update/four-million-entry limit or a weak/inaccurate pivot), -1 on invalid input.
 * Every nonzero return leaves the old basis and factorization usable. The
 * caller may refactor its NEW basis directly with lab_linear_factorize_ex. */
int lab_linear_update(LabLinear *s, size_t column, size_t nnz, const size_t *row,
                      const double *value, const double *direction);
int lab_linear_update_dense(LabLinear *s, size_t column, const double *new_column,
                            const double *direction);
/* Current original-unit column, including unsorted rows/explicit zeros from
 * the source input. Views are borrowed, immutable and valid until the next
 * accepted update, refactorization or free. Invalid requests return -1 and
 * clear provided outputs. No allocation or materialization occurs. */
int lab_linear_get_column(const LabLinear *s, size_t column, size_t *nnz, const size_t **row_index,
                          const double **value);
/* Allocate an independent current CSC copy. The destination must be zeroed
 * or freed first; an invalid request leaves it unchanged. No internal pointer
 * is transferred to the caller. This explicit O(n+nnz) operation is optional. */
int lab_linear_copy_csc(LabLinear *s, LabLinearCSC *copy);
void lab_linear_csc_free(LabLinearCSC *copy);
void lab_linear_free(LabLinear *s);
#endif

#ifndef LAB_FACTOR_H
#define LAB_FACTOR_H
#include <stddef.h>
/* Factorization backends share numerical guards, but no LP model semantics.
 * Start from a zero-initialized factor and free after either success or failure. */
typedef struct {
    size_t n;
    double *lu;
    size_t *pivot;
    void *sparse_data;
    /* Logical LU entries (diagonal counted once), and peak live row entries.
     * The sparse transpose index mirrors these entries; it is not counted twice.
     * For dense reference storage, peak_nnz is the allocated n*n capacity. */
    size_t factor_nnz, peak_nnz;
    size_t core_dimension, eliminated_singletons;
    int sparse, valid;
    int core_refinement;
} LabFactor;
int lab_factorize(LabFactor *f, const double *matrix, size_t n);
/* CSC input may contain unsorted rows and duplicate entries. Duplicate values
 * are summed. Exact zeros are discarded; no tolerance-based dropping is used.
 * Storage is proportional to dimension + nonzeros + elimination fill, never n*n.
 * A fixed fill budget causes an explicit failure before unbounded storage growth. */
int lab_factorize_sparse(LabFactor *f, size_t n, const size_t *col_start, const size_t *row_index,
                         const double *value);
int lab_factor_solve(const LabFactor *f, const double *rhs, double *x, int transpose);
void lab_factor_free(LabFactor *f);
#endif

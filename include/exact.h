#ifndef LAB_EXACT_H
#define LAB_EXACT_H
#include <gmp.h>
#include <stddef.h>

/* Immutable original-decimal feasibility data. CSC values align with LabModel.
 * Infinite bounds have a zero placeholder; consult the model's finite sides. */
typedef struct LabExactModel {
    size_t rows, cols, nnz;
    mpq_t *value, *row_lower, *row_upper, *col_lower, *col_upper;
    /* Original parsed binary bounds identify subsequent in-memory overlays. */
    double *col_lower_source, *col_upper_source;
} LabExactModel;

/* Bounded decimal/scientific parser; supports MPS D exponents. Returns -1 for
 * malformed input or resource limits, without changing a model's float reader. */
int lab_exact_decimal(mpq_t value, const char *text);
LabExactModel *lab_exact_model_new(size_t rows, size_t cols, size_t nnz);
void lab_exact_model_free(LabExactModel *model);
#endif

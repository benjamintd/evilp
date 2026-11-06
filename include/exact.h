#ifndef LAB_EXACT_H
#define LAB_EXACT_H
#include <gmp.h>
#include <stddef.h>
#include "model.h"

typedef struct {
    double lower, upper;
} LabInterval;
typedef struct LabSourceEnvelope {
    size_t rows, cols, nnz;
    LabInterval *value, *row_lower, *row_upper, *col_lower, *col_upper;
    double *col_lower_source, *col_upper_source;
    int complete;
} LabSourceEnvelope;
/* Outward binary64 interval operations. A denominator containing zero gives
 * the whole real line; indeterminate extended operations never narrow it. */
LabInterval lab_interval_add(LabInterval a, LabInterval b);
LabInterval lab_interval_subtract(LabInterval a, LabInterval b);
LabInterval lab_interval_multiply(LabInterval a, LabInterval b);
LabInterval lab_interval_divide(LabInterval a, LabInterval b);
/* text has already passed the optimizer's finite decimal grammar/conversion. */
LabInterval lab_source_number_interval(const char *text, double rounded);
LabSourceEnvelope *lab_source_envelope_new(size_t rows, size_t cols, size_t nnz);
void lab_source_envelope_free(LabSourceEnvelope *source);

/* Immutable original-decimal feasibility data. CSC values align with LabModel.
 * Infinite bounds have a zero placeholder; consult the model's finite sides. */
typedef struct LabExactModel {
    size_t rows, cols, nnz;
    mpq_t *value, *row_lower, *row_upper, *col_lower, *col_upper;
    mpq_t *objective;
    mpq_t offset;
    /* Original parsed binary bounds identify subsequent in-memory overlays. */
    double *col_lower_source, *col_upper_source;
} LabExactModel;

/* Bounded decimal/scientific parser; supports MPS D exponents. Returns -1 for
 * malformed input or resource limits, without changing a model's float reader. */
int lab_exact_decimal(mpq_t value, const char *text);
LabExactModel *lab_exact_model_new(size_t rows, size_t cols, size_t nnz);
void lab_exact_model_free(LabExactModel *model);
/* Closest finite double (either nearest value is allowed for an exact tie).
 * Returns -1 when the rational magnitude exceeds the finite double range. */
int lab_exact_nearest(const mpq_t value, double *result);
/* Evaluate a binary-double candidate in original decimal units where source
 * provenance is available; bound overlays use their exact binary values.
 * Otherwise use sparse FMA/compensated binary arithmetic. Violation is a finite
 * nonnegative upper bound; objective is a finite nearest/compensated value.
 * No feasibility tolerance or acceptance policy is applied here. */
int lab_exact_evaluate(const LabModel *model, const double *x, double *objective,
                       double *max_violation);
/* Exact accessors preserve initialized out on unavailability. Parsed source
 * without an MPQ snapshot requires a complete matching point enclosure.
 * Column values distinct from original source doubles are binary overlays. */
int lab_source_exact_coefficient(const LabModel *, size_t entry, mpq_t out);
int lab_source_exact_row_bound(const LabModel *, size_t row, int upper, mpq_t out);
int lab_source_exact_column_value(const LabModel *, size_t column, int upper, double value, mpq_t out);
int lab_source_exact_column_bound(const LabModel *, size_t column, int upper, mpq_t out);
#endif

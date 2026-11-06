#ifndef LAB_NATIVE_TABLEAU_H
#define LAB_NATIVE_TABLEAU_H
#include "lp.h"
#include <gmp.h>

typedef enum {
    LAB_NATIVE_OK, LAB_NATIVE_NO_ROWS, LAB_NATIVE_STALE, LAB_NATIVE_SOURCE_UNAVAILABLE,
    LAB_NATIVE_FREE_SPLIT, LAB_NATIVE_RESOURCE, LAB_NATIVE_DEADLINE, LAB_NATIVE_INVALID,
    LAB_NATIVE_EXACT_FAILURE, LAB_NATIVE_IDENTITY_FAILURE
} LabNativeStatus;
typedef struct {
    size_t max_rows, max_terms, max_basis_dimension, max_basis_nnz, max_canonical_nnz;
    size_t max_bits, max_total_bits;
    double deadline;
} LabNativeLimits;
typedef struct {
    int has_lower, has_upper;
    mpq_t lower, upper;
} LabNativeBox;
typedef enum { LAB_NATIVE_VAR, LAB_NATIVE_ROW, LAB_NATIVE_BOUND } LabNativeTermKind;
typedef struct {
    LabNativeTermKind kind;
    size_t index, canonical_column;
    int upper, is_integer;
    mpq_t coefficient, anchor, canonical_scale, canonical_offset;
} LabNativeTerm;
typedef struct {
    size_t index;
    mpq_t multiplier;
} LabNativeMultiplier;
typedef struct {
    size_t basic_column, basis_position, canonical_basic;
    mpq_t beta, signed_basic_scale, basic_shift;
    size_t terms, equalities, fixed;
    LabNativeTerm *term;
    LabNativeMultiplier *equality, *fix;
    /* Exact canonical inverse-row weights, diagnostic evidence only. */
    mpq_t *weights;
    size_t weight_count;
    int basis_identity_verified, native_identity_verified;
} LabNativeTableau;
typedef struct {
    size_t rows, cols, nnz;
    size_t *start, *row;
    mpq_t *value;
} LabNativeMatrix;
typedef struct {
    size_t count, columns, source_rows;
    LabNativeTableau *tableau;
    LabNativeBox *box, *source_row_box;
    unsigned char *integer;
    LabNativeMatrix source_matrix, canonical_matrix, basis_matrix;
    mpq_t *rhs, *basic_values;
    size_t dimension, total_rows, canonical_columns;
    size_t *basis, *row_id;
    int row_limit_reached;
    /* Local premises only: presolve tightenings and the tree are NOT proved. */
    const char *scope;
} LabNativeBatch;

void lab_native_batch_init(LabNativeBatch *batch);
void lab_native_batch_free(LabNativeBatch *batch);
const char *lab_native_status_name(LabNativeStatus status);
/* Requires a successful solve at the current bound generation. Pure export:
 * no pivots, caches, metrics, result statuses, input data or bounds change.
 * Initialized out is replaced only by a successful, complete owned batch.
 * Unsupported/limited exact reconstruction is abstention, never a proof.
 * Current bounds: 16,384 original rows/columns, 100,000 original nonzeros,
 * 1,024 canonical basis rows, 100,000 basis/canonical nonzeros and four rows.
 * Parsed decimal source requires retained rationals or matching complete point
 * enclosures for every used value. Uncertain data and free splits abstain.
 */
LabNativeStatus lab_lp_session_export_native(LabLpSession *session,
    const LabNativeLimits *limits, LabNativeBatch *out);
#ifdef LAB_NATIVE_TABLEAU_TESTING
void lab_native_test_fail_after(size_t allocations);
#endif
#endif

#ifndef LAB_ROOT_CUT_H
#define LAB_ROOT_CUT_H

#include "native_tableau.h"

typedef enum {
    LAB_ROOT_CUT_OK,
    LAB_ROOT_CUT_NONE,
    LAB_ROOT_CUT_STALE,
    LAB_ROOT_CUT_SOURCE_UNAVAILABLE,
    LAB_ROOT_CUT_RESOURCE,
    LAB_ROOT_CUT_DEADLINE,
    LAB_ROOT_CUT_INVALID,
    LAB_ROOT_CUT_EXACT_FAILURE
} LabRootCutStatus;

typedef struct {
    LabNativeLimits native;
    size_t max_certificate_bytes; /* at most 16 MiB per cut */
} LabRootCutLimits;

typedef struct {
    size_t columns, basic_column;
    double *coefficient;
    double lower;
    char *certificate; /* complete MILPLAB_GMI_1 derivation, owned */
} LabRootCut;

typedef struct {
    size_t count, columns;
    LabRootCut *cut;
    /* Exact premises for every cut: ORIGINAL_MIP intersected with this box.
     * No implication from presolve, global validity or tree proof is claimed.
     * A caller may use cuts only in this box or a proved contained box. */
    LabNativeBox *box;
} LabRootCutBatch;

void lab_root_cut_batch_init(LabRootCutBatch *batch);
void lab_root_cut_batch_free(LabRootCutBatch *batch);
LabRootCutLimits lab_root_cut_limits(double deadline);
const char *lab_root_cut_status_name(LabRootCutStatus status);
/* Export the current optimal, nonstale session and derive up to four exact
 * GMI cuts. Source must have the same matrix, integrality and row premises.
 * Conversion to binary rows is conservatively weakened over the exact box.
 * Success transactionally replaces initialized out. Failure leaves it intact.
 * This constructor does not independently check proofs, add rows or prune. */
/* The native export's bounded source/exact-basis restrictions apply. Exceeding
 * them is abstention and cannot count as benchmark completion. */
LabRootCutStatus lab_root_cut_derive(LabLpSession *session, const LabModel *source,
                                   const LabRootCutLimits *limits,
                                   LabRootCutBatch *out);

#endif

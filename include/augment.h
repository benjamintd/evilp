#ifndef LAB_AUGMENT_H
#define LAB_AUGMENT_H
#include "model.h"

/* Appended rows are already justified by their caller. This module only owns
 * the numerical model and faithfully preserves its source arithmetic. It does
 * not establish that a cut or its recorded box follows from the original MIP. */
typedef struct {
    const char *name;
    const double *coefficient; /* dense, base->cols finite binary64 values */
    double lower, upper;
} LabAddedRow;
typedef struct {
    size_t max_rows, max_nonzeros, max_added_rows;
} LabAugmentLimits;
typedef enum {
    LAB_AUGMENT_OK, LAB_AUGMENT_INVALID, LAB_AUGMENT_RESOURCE
} LabAugmentStatus;

/* base stays borrowed and immutable. out must be a zero-initialized or owned
 * LabModel, distinct from base. Success atomically replaces out with a fully
 * owned copy, including exact/enclosure provenance and names. Failure leaves
 * out untouched. Destroy sessions borrowing out before replacing/freeing it.
 * Added coefficients/bounds have their exact binary64 rational meaning. */
LabAugmentStatus lab_model_append_rows(const LabModel *base, size_t count,
    const LabAddedRow *rows, const LabAugmentLimits *limits, LabModel *out);
#endif

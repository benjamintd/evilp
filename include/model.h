#ifndef LAB_MODEL_H
#define LAB_MODEL_H

#include <stddef.h>
typedef struct LabExactModel LabExactModel;
typedef struct LabSourceEnvelope LabSourceEnvelope;
typedef struct LabSourceObjective LabSourceObjective;

/* Original (unscaled) bounded-row model. The matrix is compressed by column.
 * Bounds use +/-INFINITY; every stored coefficient and objective is finite.
 * Integer flags are metadata: the LP engine must explicitly relax them. */
typedef struct {
    size_t rows, cols, nnz;
    char **row_name, **col_name;
    double *row_lower, *row_upper, *col_lower, *col_upper, *objective;
    unsigned char *integer;
    size_t *col_start, *row_index;
    double *value;
    double offset;
    int maximize;
    char error[512];
    int unsupported;
    /* Optional bounded decimal provenance, owned by this model. In-memory
     * sessions borrow it; manually constructed models leave both fields zero. */
    LabExactModel *exact;
    /* Optional larger, bounded feasibility-only decimal provenance for native
     * cut derivations. Never enables global primal/objective exact paths.
     * Owned independently; sessions/views borrow it. */
    LabExactModel *native_exact;
    /* Bounded objective-only decimal provenance, independent of basis data. */
    LabSourceObjective *objective_exact;
    int parsed_decimal;
    /* Outward source-decimal enclosures, independent of the exact-size cap.
     * Session views borrow them; this model owns the allocation. */
    LabSourceEnvelope *envelope;
} LabModel;

/* Start from a zero-initialized model. On failure the diagnostic stays in
 * model->error; lab_model_free is valid after either success or failure. */
int lab_model_read(LabModel *model, const char *path);
void lab_model_free(LabModel *model);

#endif

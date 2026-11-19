#ifndef LAB_PACKING_H
#define LAB_PACKING_H
#include "model.h"

/* Return a complete point only when its opening cost meets an exact original
 * model lower bound. Attempts include restarts and local search steps. */
int lab_resource_packing(const LabModel *model, double deadline, size_t limit, double *point,
                         size_t *attempts);
#ifdef LAB_GLPK_BASIS
int lab_packing_patterns(const LabModel *model, double deadline, size_t limit, double *point,
                         size_t *work);
#endif
#endif

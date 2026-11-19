#ifndef LAB_DEPENDENCY_H
#define LAB_DEPENDENCY_H
#include "model.h"
/* Verify an equality-row dependence against exact original-source scalars.
 * Binary multipliers are candidates only. Return1 for a nonzero exact null
 * combination (including RHS),0 for rejection/unavailable/budget. */
int lab_source_row_dependency(const LabModel *model, const double *multipliers, double deadline);
#endif

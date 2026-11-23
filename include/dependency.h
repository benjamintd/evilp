#ifndef LAB_DEPENDENCY_H
#define LAB_DEPENDENCY_H
#include "model.h"
/* Verify an equality-row dependence against exact original-source scalars.
 * Binary multipliers are candidates only. Return1 for a nonzero exact null
 * combination (including RHS),0 for rejection/unavailable/budget. */
int lab_source_row_dependency(const LabModel *model, const double *multipliers, double deadline);
/* 1 for a verified optimal point, -1 for an exact feasible point which cannot
 * be exported and verified as binary64, 0 when this recovery does not apply. */
int lab_recover_square_equalities(const LabModel *, double *point, double *objective, double deadline);
#endif

#ifndef LAB_ASSIGNMENT_H
#define LAB_ASSIGNMENT_H
#include "model.h"

/* Recognize binary perfect matchings with optional difference constraints on
 * continuous potentials. Bounds must be valid root-presolve implications. */
int lab_assignment_solve(const LabModel *model, const double *lower, const double *upper,
                         double deadline, size_t limit, double *point, char ***decimal,
                         size_t *trials);
#ifdef LAB_SCIP_INTEGER
#include "scip_mip.h"
int lab_assignment_quotient(const LabModel *model, const double *lower, const double *upper,
                           double deadline, size_t node_limit, size_t iteration_limit,
                           double *point, char ***decimal, LabScipMipStats *stats);
#endif
void lab_assignment_free_decimals(char **decimal, size_t cols);
#endif

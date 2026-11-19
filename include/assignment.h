#ifndef LAB_ASSIGNMENT_H
#define LAB_ASSIGNMENT_H
#include "model.h"

/* Recognize binary perfect matchings with optional difference constraints on
 * continuous potentials. Bounds must be valid root-presolve implications. */
int lab_assignment_solve(const LabModel *model, const double *lower, const double *upper,
                         double deadline, size_t limit, double *point, char ***decimal,
                         size_t *trials);
void lab_assignment_free_decimals(char **decimal, size_t cols);
#endif

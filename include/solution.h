#ifndef LAB_SOLUTION_H
#define LAB_SOLUTION_H
#include "model.h"
#include <stdio.h>

/* Shortest significant-digit decimal that roundtrips to the same finite double. */
void lab_solution_decimal(FILE *file, double value);
/* Serialize a finite point and, when bounded source arithmetic is available,
 * its exact decimal objective. The objective describes the emitted coordinates.
 * A missing objective claim leaves scoring to the independent checker.
 * This performs serialization only, with no feasibility or optimality claim. */
int lab_solution_write_point(FILE *file, const LabModel *model, const double *point);
int lab_solution_write_decimals(FILE *file, const LabModel *model, const double *point,
                                char *const *decimal);
#endif

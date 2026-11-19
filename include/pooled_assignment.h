#ifndef LAB_POOLED_ASSIGNMENT_H
#define LAB_POOLED_ASSIGNMENT_H
#include "model.h"
typedef struct {
    size_t nodes, lp_solves, iterations;
} LabPooledStats;
/* Return 2 for an exhausted search with an original integer witness, 1 for
 * an incumbent, and 0 otherwise. No infeasibility status is returned. */
int lab_pooled_assignment(const LabModel *model, double deadline, size_t nodes, size_t iterations,
                          double *point, LabPooledStats *stats);
#endif

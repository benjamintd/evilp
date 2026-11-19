#ifndef LAB_BINARY_SEARCH_H
#define LAB_BINARY_SEARCH_H

#include "model.h"

/* Exact finite binary subproblem. Returns 1 with a minimizing/maximizing
 * point, 0 to defer to branch-and-bound (including empty or limited searches).
 * Queries is the number of half assignments searched against the other half. */
int lab_binary_search(const LabModel *model, double deadline, size_t query_budget, double *point,
                      size_t *queries);
#endif

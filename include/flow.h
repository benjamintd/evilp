#ifndef LAB_FLOW_H
#define LAB_FLOW_H
#include "model.h"
/* Returns 1 only for a reconstructed optimum. Unsupported, empty, limited
 * or infeasible searches defer to ordinary branch-and-bound. */
int lab_fixed_charge_flow(const LabModel *m, double deadline, size_t limit, double *point,
                          size_t *eliminations, size_t *states);
#endif

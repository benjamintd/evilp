#ifndef LAB_BINARY_SAT_H
#define LAB_BINARY_SAT_H
#include "model.h"

/* Solve an exactly representable binary model. Unknown and unsupported
 * return 0. A checked encoding yields 1 for optimal, 2 for an incumbent
 * only and -1 for infeasible. lo/hi are safe root propagation bounds. */
int lab_binary_sat(const LabModel *model, const double *lo, const double *hi, double deadline,
                   size_t decision_budget, double *point, size_t *calls);
#endif

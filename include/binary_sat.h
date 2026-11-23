#ifndef LAB_BINARY_SAT_H
#define LAB_BINARY_SAT_H
#include "search_result.h"
#include "model.h"

/* The returned point borrows the supplied point buffer. */
LabSearchResult lab_binary_sat(const LabModel *model, const double *lo, const double *hi, double deadline,
                   size_t decision_budget, double *point, size_t *calls);
#endif

#ifndef LAB_INTEGER_MIP_H
#define LAB_INTEGER_MIP_H
#include "search_result.h"
#include "model.h"

typedef struct {
    size_t nodes, lp_solves, lp_iterations;
} LabIntegerMipStats;

/* The returned point borrows the supplied point buffer. */
LabSearchResult lab_integer_mip(const LabModel *, const double *lower, const double *upper, double deadline,
                    size_t node_limit, size_t lp_iteration_limit, double *point,
                    LabIntegerMipStats *);
#endif

#ifndef LAB_INTEGER_MIP_H
#define LAB_INTEGER_MIP_H
#include "model.h"

typedef struct {
    size_t nodes, lp_solves, lp_iterations;
} LabIntegerMipStats;

/* 2 for an original-feasible incumbent at a solver-reported optimum, 1 for
 * a checked incumbent, 0 for an unsupported or unsuccessful search. */
int lab_integer_mip(const LabModel *, const double *lower, const double *upper, double deadline,
                    size_t node_limit, size_t lp_iteration_limit, double *point,
                    LabIntegerMipStats *);
#endif

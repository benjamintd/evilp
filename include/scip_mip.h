#ifndef LAB_SCIP_MIP_H
#define LAB_SCIP_MIP_H
#include "model.h"

typedef struct {
    size_t nodes, lp_solves, lp_iterations;
    int status; /* 0 = abstained/other, 1 = time_limit, 2 = node_limit, 3 = optimal, 4 = infeasible */
    int has_bound;
    double best_bound;
} LabScipMipStats;

/* General-purpose C17 SCIP branch-and-cut search over presolved bounds.
 * Returns:
 *   2: original-model-verified optimal witness
 *  -2: solver-reported infeasibility
 *   1: original-model-verified incumbent (incomplete search)
 *   0: no verified incumbent (see stats->status / stats->has_bound)
 */
int lab_scip_mip(const LabModel *model, const double *lower, const double *upper,
                 const double *row_lower, const double *row_upper, const double *incumbent,
                 double deadline, size_t node_limit, size_t lp_iteration_limit, double *point,
                 LabScipMipStats *stats);
#endif

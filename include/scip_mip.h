#ifndef LAB_SCIP_MIP_H
#define LAB_SCIP_MIP_H
#include "search_result.h"
#include "model.h"

typedef struct {
    size_t nodes, lp_solves, lp_iterations, path_cliques;
    LabMipStatus status;
    const char *fallback_reason;
    int has_bound;
    double best_bound;
} LabScipMipStats;

/* A returned point borrows the supplied point buffer. */
LabSearchResult lab_scip_mip(const LabModel *model, const double *lower, const double *upper,
                 const double *row_lower, const double *row_upper, const double *incumbent,
                 double deadline, size_t node_limit, size_t lp_iteration_limit, double *point,
                 LabScipMipStats *stats);
#endif

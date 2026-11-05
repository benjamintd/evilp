#ifndef LAB_BB_H
#define LAB_BB_H
#include "lp.h"

/* Minimal numerical branch-and-bound. Outcomes and bounds are solver claims;
 * this module does not implement an independently checkable tree proof. */
typedef struct {
    LabLpOptions lp;
    double time_limit;
    size_t node_limit;
    int cold;
    const char *trace_path;
} LabBbOptions;

typedef struct {
    const char *status;
    char message[256];
    double *primal;
    double objective, best_bound, elapsed;
    int has_primal, has_bound, bound_certified;
    size_t nodes, lp_solves, lp_iterations, refactorizations;
    size_t warm_starts, cold_starts, pruned_nodes;
    size_t dual_iterations, cold_fallbacks, factor_nnz_peak, standard_nnz_peak;
    size_t refinement_steps, scaled_rows, scaled_columns;
    double initial_backward_error, max_backward_error;
} LabBbResult;

void lab_bb_solve(const LabModel *model, const LabBbOptions *options, LabBbResult *result);
void lab_bb_result_free(LabBbResult *result);
#endif

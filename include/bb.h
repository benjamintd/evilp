#ifndef LAB_BB_H
#define LAB_BB_H
#include "lp.h"

/* Minimal numerical branch-and-bound. Outcomes and bounds are solver claims;
 * this module does not implement an independently checkable tree proof. */
typedef struct {
    LabLpOptions lp;
    double time_limit;
    size_t node_limit;
    int cold, no_presolve, first_fractional;
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
    size_t exact_attempts, exact_infeasibilities;
    size_t source_infeasibility_attempts, source_infeasibilities, source_projected_rows;
    size_t source_removed_terms;
    size_t rejected_pivots, harris_rescues;
    size_t primal_repairs, primal_repair_iterations;
    size_t exact_primal_attempts, exact_primal_recoveries;
    size_t basis_updates, basis_update_rejections;
    size_t factor_core_peak, singleton_eliminations_peak;
    size_t presolve_passes, presolve_tightenings, presolve_fixed, presolve_rows;
    size_t branch_observations;
    size_t primal_accepted_pivots, primal_zero_steps, primal_short_steps, pricing_bland_fallbacks;
} LabBbResult;

void lab_bb_solve(const LabModel *model, const LabBbOptions *options, LabBbResult *result);
void lab_bb_result_free(LabBbResult *result);
#endif

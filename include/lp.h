#ifndef LAB_LP_H
#define LAB_LP_H
#include "model.h"
#define LAB_LP_SPARSE_DIMENSION_LIMIT 16000000
#define LAB_LP_DENSE_DIMENSION_LIMIT 4096
typedef struct {
    double time_limit;
    size_t iteration_limit, max_dimension;
    const char *basis_in, *basis_out, *model_sha256;
    int dense; /* zero: sparse LU; one: dense reference backend */
    int no_scaling, no_refinement, no_exact, bland_ratio;
    int no_primal_repair;
    int no_updates;
    int bland_pricing;
} LabLpOptions;
typedef struct {
    const char *status;
    char message[256];
    double *primal, *dual;
    double objective, primal_residual, dual_residual, elapsed;
    size_t iterations, refactorizations, standard_rows, standard_columns, redundant_rows;
    int has_primal;
    int warm_start;
    size_t dual_iterations, cold_fallbacks, factor_nnz_peak, standard_nnz;
    size_t refinement_steps, scaled_rows, scaled_columns;
    double initial_backward_error, max_backward_error;
    /* Rational original-row Farkas candidate; external checking is separate.
     * Bound overlays are part of the LP to which this candidate applies.
     * The margin is exact unless farkas_margin_is_lower_bound is set. */
    char **farkas;
    char *farkas_margin;
    size_t farkas_rows, exact_attempts, exact_infeasibilities;
    int farkas_margin_is_lower_bound;
    size_t source_infeasibility_attempts, source_infeasibilities, source_projected_rows;
    size_t source_removed_terms;
    size_t source_reconstructed_terms;
    size_t rejected_pivots, harris_rescues;
    size_t primal_repairs, primal_repair_iterations;
    size_t exact_primal_attempts, exact_primal_recoveries;
    int primal_repair_required;
    size_t basis_updates, basis_update_rejections;
    size_t factor_core_peak, singleton_eliminations_peak;
    size_t primal_accepted_pivots, primal_zero_steps, primal_short_steps, pricing_bland_fallbacks;
} LabLpResult;
void lab_lp_solve(const LabModel *m, const LabLpOptions *options, LabLpResult *result);
void lab_lp_result_free(LabLpResult *result);
#define LAB_LP_SESSION_API 1
typedef struct LabLpSession LabLpSession;
/* The original model must outlive the session and stay immutable. Limits use
 * one cumulative deadline. Bounds can be restored within the original domain. */
LabLpSession *lab_lp_session_create(const LabModel *m, const LabLpOptions *options);
int lab_lp_session_set_bounds(LabLpSession *s, size_t column, double lower, double upper);
void lab_lp_session_solve(LabLpSession *s, LabLpResult *result);
void lab_lp_session_free(LabLpSession *s);
#endif

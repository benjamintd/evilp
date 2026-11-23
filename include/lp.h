#ifndef LAB_LP_H
#define LAB_LP_H
#include "augment.h"
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
    int crash_basis; /* guarded singleton structural starting columns */
    int fast_residual;
    int adaptive_pricing; /* delayed fallback with repeated-basis/stall guard */
    int glpk_basis; /* optional starting basis; requires LAB_GLPK_BASIS */
    int verify_glpk_point; /* standalone solves may verify a proposal before lifting its basis */
    int packing_bounds;
    int keep_bound_rows; /* retain explicit rows for transformed implicit bounds */
    /* Zero selects the established defaults (1e-7 absolute source primal,
     * 1e-10 relative reduced cost). Positive values can only tighten them. */
    double primal_tolerance, dual_relative_tolerance;
} LabLpOptions;
typedef struct {
    const char *status;
    char message[256];
    double *primal, *dual;
    char **dual_rational; /* exact source closure, when available */
    size_t dual_rows;
    double objective, primal_residual, dual_residual, elapsed;
    size_t iterations, refactorizations, standard_rows, standard_columns, redundant_rows;
    size_t implicit_bound_rows, packing_bound_rows;
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
    size_t glpk_basis_iterations, glpk_bases_accepted;
    const char *glpk_basis_status;
    size_t basis_updates, basis_update_rejections;
    size_t factor_core_peak, singleton_eliminations_peak;
    size_t primal_accepted_pivots, primal_zero_steps, primal_short_steps, pricing_bland_fallbacks;
    size_t primal_candidate_rejections;
    size_t source_pricing_attempts, source_pricing_pivots, source_pricing_completions;
    size_t source_dual_closures;
    size_t source_support_scans, source_support_binary_paths;
} LabLpResult;
void lab_lp_solve(const LabModel *m, const LabLpOptions *options, LabLpResult *result);
void lab_lp_result_free(LabLpResult *result);
#define LAB_LP_SESSION_API 1
typedef struct LabLpSession LabLpSession;
/* The original model must outlive the session and stay immutable. Limits use
 * one cumulative deadline. Bounds can be restored within the original domain. */
LabLpSession *lab_lp_session_create(const LabModel *m, const LabLpOptions *options);
/* Update policy on a live session and recheck the basis on its next solve.
 * Zero restores the default for that quantity. An open probe rejects changes. */
int lab_lp_session_set_tolerances(LabLpSession *s, double primal, double dual_relative);
int lab_lp_session_set_bounds(LabLpSession *s, size_t column, double lower, double upper);
/* One reversible bound probe at a time. Rollback restores the parent domain
 * and, when the standard form survives, its saved basis. Deadlines never rewind. */
int lab_lp_session_begin_probe(LabLpSession *s);
int lab_lp_session_rollback_probe(LabLpSession *s);
/* Permanently append a batch to the LP domain. The session owns the new
 * sparse matrix and its exact source metadata. Active bound overlays survive.
 * A structural edit invalidates the old standard form and disk basis identity. */
LabAugmentStatus lab_lp_session_append_rows(LabLpSession *s, size_t count,
    const LabSparseRow *rows, const LabAugmentLimits *limits);
/* Same transaction for cuts already produced as dense coefficient vectors. */
LabAugmentStatus lab_lp_session_append_dense_rows(LabLpSession *s, size_t count,
    const LabAddedRow *rows, const LabAugmentLimits *limits);
/* Borrowed immutable domain, valid until the next structural edit or free. */
const LabModel *lab_lp_session_domain(const LabLpSession *s);
void lab_lp_session_solve(LabLpSession *s, LabLpResult *result);
void lab_lp_session_free(LabLpSession *s);
#endif

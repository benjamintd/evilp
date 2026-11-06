#ifndef LAB_DIVE_H
#define LAB_DIVE_H
#include "lp.h"

/* A primal heuristic only. A failed or infeasible dive establishes no node
 * bound or original-model infeasibility. The caller must validate any offered
 * point against its original model and node domain before accepting it. */
typedef enum {
    LAB_DIVE_CANDIDATE, LAB_DIVE_NOT_FOUND, LAB_DIVE_TIME_LIMIT,
    LAB_DIVE_STEP_LIMIT, LAB_DIVE_LP_ABSTAIN, LAB_DIVE_INVALID,
    LAB_DIVE_RESOURCE
} LabDiveStatus;
typedef struct {
    double time_limit;
    size_t max_steps, max_backtracks, max_rows, max_columns, max_nonzeros;
    int cold;
} LabDiveOptions;
typedef struct {
    LabDiveStatus status;
    char message[160];
    double *primal, objective, elapsed;
    size_t steps, backtracks, lp_solves, lp_iterations, warm_starts, rounded_checks;
    size_t final_resolves;
    const char *last_lp_status;
} LabDiveResult;

LabDiveOptions lab_dive_options(double time_limit);
const char *lab_dive_status_name(LabDiveStatus status);
/* out starts zero-initialized. All input models, bounds, points and sessions
 * stay immutable; the heuristic owns its separate domain, state and session.
 * lp_model may contain valid root cuts. The caller establishes their scope;
 * root bounds must be contained in lp_model's declared column domain. */
void lab_primal_dive(const LabModel *original, const LabModel *lp_model,
                     const double *root_lower, const double *root_upper,
                     const double *initial, const LabLpOptions *lp_options,
                     const LabDiveOptions *options, LabDiveResult *out);
void lab_dive_result_free(LabDiveResult *result);
#endif

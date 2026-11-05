#ifndef LAB_PRESOLVE_H
#define LAB_PRESOLVE_H
#include "model.h"
typedef struct {
    double time_limit;
    size_t max_passes;
    int relax; /* ignore integrality when presolving an LP relaxation */
} LabPresolveOptions;
typedef struct {
    const char *status; /* unresolved, infeasible, or input_error; solver claim */
    char message[256];
    double *lower, *upper;
    unsigned char *redundant;
    size_t passes, tightenings, fixed, redundant_rows;
    int source_available, limited;
    double elapsed;
} LabPresolveResult;
/* Original model is immutable. Optional input bounds are numerical tightenings
 * of its current box; NULL selects that side of the model. Output bounds stay
 * inside this numerical box. Proof calculations instead use an OUTER source
 * box, so rounded input values alone cannot justify an infeasibility claim.
 * Redundant-row masks are valid only with the returned box and integrality.
 * No independent certificate is emitted by this optimizer module. */
void lab_presolve(const LabModel *model, const double *lower, const double *upper,
                  const LabPresolveOptions *options, LabPresolveResult *result);
void lab_presolve_result_free(LabPresolveResult *result);
#endif

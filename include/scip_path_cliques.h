#ifndef LAB_SCIP_PATH_CLIQUES_H
#define LAB_SCIP_PATH_CLIQUES_H
#include "model.h"
#include "scip/scip.h"

/* Optional source-checked clique inequalities for conditional differences.
 * Partial separation under the deadline leaves previously added cuts valid. */
SCIP_RETCODE lab_scip_path_cliques(SCIP *scip, const LabModel *model, SCIP_VAR **vars,
                                   const size_t *start, const int *entry, const double *lower,
                                   const double *upper, const double *row_lower,
                                   const double *row_upper, double deadline, size_t *added);
#endif

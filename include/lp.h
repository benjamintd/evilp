#ifndef LAB_LP_H
#define LAB_LP_H
#include "model.h"
typedef struct {
    double time_limit;
    size_t iteration_limit, max_dimension;
    const char *basis_in, *basis_out, *model_sha256;
} LabLpOptions;
typedef struct {
    const char *status;
    char message[256];
    double *primal, *dual;
    double objective, primal_residual, dual_residual, elapsed;
    size_t iterations, refactorizations, standard_rows, standard_columns, redundant_rows;
    int has_primal;
} LabLpResult;
void lab_lp_solve(const LabModel *m, const LabLpOptions *options, LabLpResult *result);
void lab_lp_result_free(LabLpResult *result);
#endif

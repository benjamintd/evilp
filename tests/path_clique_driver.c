#include "common.h"
#include "model.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include "scip_path_cliques.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    LabModel model = {0};
    if (argc != 2 || lab_model_read_mip(&model, argv[1]))
        return 2;
    SCIP *scip = NULL;
    if (SCIPcreate(&scip) != SCIP_OKAY || SCIPincludeDefaultPlugins(scip) != SCIP_OKAY ||
        SCIPcreateProbBasic(scip, "path_cliques") != SCIP_OKAY)
        return 2;
    SCIP_VAR **vars = lab_alloc(model.cols, sizeof *vars);
    size_t *start = lab_alloc(model.rows + 1, sizeof *start);
    size_t *cursor = lab_alloc(model.rows, sizeof *cursor);
    int *entry = lab_alloc(model.nnz, sizeof *entry);
    for (size_t j = 0; j < model.cols; ++j) {
        SCIP_VARTYPE type = model.integer[j] ? SCIP_VARTYPE_INTEGER : SCIP_VARTYPE_CONTINUOUS;
        if (SCIPcreateVarBasic(scip, &vars[j], model.col_name[j], model.col_lower[j],
                               model.col_upper[j], model.objective[j], type) != SCIP_OKAY ||
            SCIPaddVar(scip, vars[j]) != SCIP_OKAY)
            return 2;
    }
    for (size_t k = 0; k < model.nnz; ++k)
        ++start[model.row_index[k] + 1];
    for (size_t i = 0; i < model.rows; ++i)
        start[i + 1] += start[i];
    memcpy(cursor, start, model.rows * sizeof *start);
    for (size_t k = 0; k < model.nnz; ++k)
        entry[cursor[model.row_index[k]]++] = (int)k;
    size_t count = 0;
    if (lab_scip_path_cliques(scip, &model, vars, start, entry, model.col_lower, model.col_upper,
                              model.row_lower, model.row_upper, lab_now() + 10,
                              &count) != SCIP_OKAY ||
        count != (size_t)SCIPgetNOrigConss(scip))
        return 2;
    SCIP_CONS **constraints = SCIPgetOrigConss(scip);
    for (size_t i = 0; i < count; ++i) {
        SCIP_CONS *c = constraints[i];
        if (SCIPgetLhsLinear(scip, c) != -SCIPinfinity(scip) || SCIPgetRhsLinear(scip, c) != 1)
            return 3;
        for (int k = 0; k < SCIPgetNVarsLinear(scip, c); ++k) {
            if (SCIPgetValsLinear(scip, c)[k] != 1)
                return 3;
            printf("%s%s", k ? " " : "", SCIPvarGetName(SCIPgetVarsLinear(scip, c)[k]));
        }
        putchar('\n');
    }
    free(entry);
    free(cursor);
    free(start);
    for (size_t j = 0; j < model.cols; ++j)
        if (SCIPreleaseVar(scip, &vars[j]) != SCIP_OKAY)
            return 3;
    free(vars);
    lab_model_free(&model);
    return SCIPfree(&scip) != SCIP_OKAY;
}

#include "scip_mip.h"
#include "scip_feasjump.h"
#include "scip_presol_parallel.h"
#include "common.h"
#include "exact.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int extract_verified(const LabModel *m, const double *lower, const double *upper,
                            SCIP *s, SCIP_VAR *const *vars, double *point, double *out_obj) {
    int nsols = SCIPgetNSols(s);
    if (nsols <= 0)
        return 0;
    SCIP_SOL **sols = SCIPgetSols(s);
    int limit = nsols < 8 ? nsols : 8;
    double *cand = lab_alloc(m->cols ? m->cols : 1, sizeof *cand);
    LabModel view = *m;
    view.col_lower = (double *)lower;
    view.col_upper = (double *)upper;
    int found = 0;
    for (int idx = 0; idx < limit; ++idx) {
        SCIP_SOL *sol = sols[idx];
        if (!sol)
            continue;
        int ok = 1;
        for (size_t j = 0; j < m->cols; ++j) {
            double v = SCIPgetSolVal(s, sol, vars[j]);
            if (!isfinite(v)) {
                ok = 0;
                break;
            }
            if (m->integer && m->integer[j]) {
                if (fabs(v) > 9007199254740991.0 || fabs(v - round(v)) > 1e-5) {
                    ok = 0;
                    break;
                }
                v = round(v);
            }
            if (v < lower[j] - 1e-7 || v > upper[j] + 1e-7) {
                ok = 0;
                break;
            }
            if (v < lower[j])
                v = lower[j];
            if (v > upper[j])
                v = upper[j];
            cand[j] = v;
        }
        if (!ok)
            continue;
        double obj = 0, viol = 0;
        if (!lab_exact_evaluate(&view, cand, &obj, &viol) && viol <= 1e-7) {
            memcpy(point, cand, m->cols * sizeof *point);
            *out_obj = obj;
            found = 1;
            break;
        }
    }
    free(cand);
    return found;
}

int lab_scip_mip(const LabModel *m, const double *lower, const double *upper,
                 const double *row_lower, const double *row_upper, const double *incumbent,
                 double deadline, size_t node_limit, size_t lp_iteration_limit, double *point,
                 LabScipMipStats *stats) {
    *stats = (LabScipMipStats){0};
    if (!m || !m->integer || !m->rows || !m->cols || !node_limit || !lp_iteration_limit ||
        node_limit > (size_t)LLONG_MAX || lp_iteration_limit > (size_t)LLONG_MAX ||
        m->rows > (size_t)INT_MAX || m->cols > (size_t)INT_MAX || m->nnz > (size_t)INT_MAX ||
        lab_now() >= deadline)
        return 0;

    int has_int = 0;
    if (!isfinite(m->offset) || fabs(m->offset) > 1e12)
        return 0;
    for (size_t j = 0; j < m->cols; ++j) {
        if (lower[j] > upper[j])
            return 0;
        if ((isfinite(lower[j]) && fabs(lower[j]) > 1e12) ||
            (isfinite(upper[j]) && fabs(upper[j]) > 1e12) ||
            fabs(m->objective[j]) > 1e12)
            return 0;
        if (m->integer[j])
            has_int = 1;
    }
    if (!has_int)
        return 0;

    SCIP *s = NULL;
    SCIP_VAR **vars = NULL;
    SCIP_CONS **cons = NULL;
    size_t *rstart = NULL;
    size_t *rcur = NULL;
    SCIP_VAR **rvars = NULL;
    SCIP_Real *rvals = NULL;
    unsigned char *active = NULL;
    int result = 0;

    if (SCIPcreate(&s) != SCIP_OKAY)
        goto done;
    SCIPsetMessagehdlrQuiet(s, TRUE);
    if (SCIPincludeDefaultPlugins(s) != SCIP_OKAY ||
        lab_scip_include_feasjump(s) != SCIP_OKAY ||
        lab_scip_include_presol_parallel(s) != SCIP_OKAY ||
        SCIPcreateProbBasic(s, "milp_lab") != SCIP_OKAY ||
        SCIPsetIntParam(s, "timing/clocktype", 2) != SCIP_OKAY ||
        SCIPsetIntParam(s, "parallel/maxnthreads", 1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "randomization/randomseedshift", 1) != SCIP_OKAY ||
        SCIPsetLongintParam(s, "limits/totalnodes", (SCIP_Longint)node_limit) != SCIP_OKAY ||
        SCIPsetLongintParam(s, "lp/iterlim", (SCIP_Longint)lp_iteration_limit) != SCIP_OKAY ||
        SCIPsetLongintParam(s, "lp/rootiterlim", (SCIP_Longint)lp_iteration_limit) != SCIP_OKAY ||
        SCIPsetRealParam(s, "limits/gap", 0.0) != SCIP_OKAY ||
        SCIPsetRealParam(s, "limits/absgap", 0.0) != SCIP_OKAY ||
        /* Enable built-in C presolvers that default to 0 when PaPILO is expected,
         * and cap root aggregation-cut rounds to avoid dense LP blowup. */
        SCIPsetIntParam(s, "presolving/dualagg/maxrounds", -1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "presolving/dualinfer/maxrounds", -1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "presolving/stuffing/maxrounds", -1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "presolving/tworowbnd/maxrounds", -1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "presolving/redvub/maxrounds", -1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "presolving/sparsify/maxnonzeros", 64) != SCIP_OKAY ||
        SCIPsetRealParam(s, "presolving/restartfac", 0.02) != SCIP_OKAY ||
        SCIPsetRealParam(s, "presolving/immrestartfac", 0.08) != SCIP_OKAY ||
        SCIPsetRealParam(s, "presolving/restartminred", 0.05) != SCIP_OKAY ||
        SCIPsetIntParam(s, "propagating/probing/maxprerounds", 2) != SCIP_OKAY ||
        SCIPsetIntParam(s, "propagating/probing/maxuseless", 200) != SCIP_OKAY ||
        SCIPsetIntParam(s, "propagating/probing/maxtotaluseless", 400) != SCIP_OKAY ||
        SCIPsetIntParam(s, "separating/rapidlearning/freq", 0) != SCIP_OKAY ||
        SCIPsetIntParam(s, "separating/aggregation/maxroundsroot", 10) != SCIP_OKAY ||
        SCIPsetRealParam(s, "heuristics/rens/minfixingrate", 0.3) != SCIP_OKAY ||
        SCIPsetIntParam(s, "heuristics/farkasdiving/freq", 20) != SCIP_OKAY ||
        SCIPsetIntParam(s, "heuristics/conflictdiving/freq", 20) != SCIP_OKAY)
        goto done;
    if ((m->rows >= 2000 || m->nnz >= 100000) &&
        SCIPsetBoolParam(s, "constraints/knapsack/disaggregation", FALSE) != SCIP_OKAY)
        goto done;
    if (m->maximize && SCIPsetObjsense(s, SCIP_OBJSENSE_MAXIMIZE) != SCIP_OKAY)
        goto done;
    if (m->offset != 0.0 &&
        SCIPaddOrigObjoffset(s, m->offset) != SCIP_OKAY)
        goto done;

    vars = lab_alloc(m->cols, sizeof *vars);
    double inf = SCIPinfinity(s);
    for (size_t j = 0; j < m->cols; ++j) {
        if ((j & 1023U) == 0 && lab_now() >= deadline)
            goto done;
        double lo = isfinite(lower[j]) ? lower[j] : -inf;
        double hi = isfinite(upper[j]) ? upper[j] : inf;
        if (lo < -inf)
            lo = -inf;
        if (hi > inf)
            hi = inf;
        SCIP_VARTYPE vt = SCIP_VARTYPE_CONTINUOUS;
        if (m->integer[j])
            vt = (lo >= 0.0 && hi <= 1.0) ? SCIP_VARTYPE_BINARY : SCIP_VARTYPE_INTEGER;
        char name[32];
        snprintf(name, sizeof name, "C%zu", j);
        if (SCIPcreateVarBasic(s, &vars[j], name, lo, hi, m->objective[j], vt) != SCIP_OKAY ||
            SCIPaddVar(s, vars[j]) != SCIP_OKAY)
            goto done;
    }

    active = lab_alloc(m->rows, sizeof *active);
    rstart = lab_alloc(m->rows + 1, sizeof *rstart);
    const double *rl = row_lower ? row_lower : m->row_lower;
    const double *ru = row_upper ? row_upper : m->row_upper;
    for (size_t i = 0; i < m->rows; ++i)
        if (isfinite(rl[i]) || isfinite(ru[i]))
            active[i] = 1;
    for (size_t j = 0; j < m->cols; ++j) {
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            if (active[i] && m->value[k] != 0.0)
                rstart[i + 1]++;
        }
    }
    for (size_t i = 0; i < m->rows; ++i)
        rstart[i + 1] += rstart[i];
    size_t total_nz = rstart[m->rows];
    rcur = lab_alloc(m->rows, sizeof *rcur);
    memcpy(rcur, rstart, m->rows * sizeof *rcur);
    rvars = lab_alloc(total_nz ? total_nz : 1, sizeof *rvars);
    rvals = lab_alloc(total_nz ? total_nz : 1, sizeof *rvals);
    for (size_t j = 0; j < m->cols; ++j) {
        if ((j & 1023U) == 0 && lab_now() >= deadline)
            goto done;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            double v = m->value[k];
            if (active[i] && v != 0.0) {
                size_t pos = rcur[i]++;
                rvars[pos] = vars[j];
                rvals[pos] = v;
            }
        }
    }
    cons = lab_alloc(m->rows, sizeof *cons);
    for (size_t i = 0; i < m->rows; ++i) {
        if ((i & 1023U) == 0 && lab_now() >= deadline)
            goto done;
        if (!active[i])
            continue;
        double l = isfinite(rl[i]) ? rl[i] : -inf;
        double u = isfinite(ru[i]) ? ru[i] : inf;
        if (l < -inf)
            l = -inf;
        if (u > inf)
            u = inf;
        int len = (int)(rstart[i + 1] - rstart[i]);
        char name[32];
        snprintf(name, sizeof name, "R%zu", i);
        if (SCIPcreateConsBasicLinear(s, &cons[i], name, len, len ? &rvars[rstart[i]] : NULL,
                                      len ? &rvals[rstart[i]] : NULL, l, u) != SCIP_OKAY ||
            SCIPaddCons(s, cons[i]) != SCIP_OKAY)
            goto done;
    }

    if (incumbent) {
        SCIP_SOL *init_sol = NULL;
        if (SCIPcreateOrigSol(s, &init_sol, NULL) == SCIP_OKAY && init_sol) {
            int ok = 1;
            for (size_t j = 0; j < m->cols && ok; ++j)
                if (SCIPsetSolVal(s, init_sol, vars[j], incumbent[j]) != SCIP_OKAY)
                    ok = 0;
            SCIP_Bool stored = FALSE;
            if (ok)
                (void)SCIPaddSolFree(s, &init_sol, &stored);
            else
                (void)SCIPfreeSol(s, &init_sol);
        }
    }

    double remaining = deadline - lab_now();
    if (remaining <= 0.01) {
        stats->status = 1;
        goto done;
    }
    if (SCIPsetRealParam(s, "limits/time", remaining * 0.96) != SCIP_OKAY)
        goto done;

    SCIP_RETCODE rc = SCIPsolve(s);
    if (rc != SCIP_OKAY)
        goto done;

    stats->nodes = (size_t)SCIPgetNTotalNodes(s);
    stats->lp_iterations = (size_t)SCIPgetNLPIterations(s);
    stats->lp_solves = (size_t)SCIPgetNLPs(s);

    SCIP_STATUS st = SCIPgetStatus(s);
    double dual = SCIPgetDualbound(s);
    if (isfinite(dual) && fabs(dual) < 0.5 * inf) {
        stats->has_bound = 1;
        stats->best_bound = dual;
    }

    if (st == SCIP_STATUS_INFEASIBLE) {
        stats->status = 4;
        result = -2;
        goto done;
    }

    double orig_obj = 0.0;
    int verified = extract_verified(m, lower, upper, s, vars, point, &orig_obj);
    if (st == SCIP_STATUS_OPTIMAL && verified) {
        double primal = SCIPgetPrimalbound(s);
        double scale = 1.0 + fabs(orig_obj);
        if (isfinite(primal) && isfinite(dual) &&
            fabs(primal - orig_obj) <= 1e-4 * scale &&
            fabs(dual - orig_obj) <= 1e-4 * scale) {
            stats->status = 3;
            stats->has_bound = 1;
            stats->best_bound = orig_obj;
            result = 2;
            goto done;
        }
    }
    if (verified)
        result = 1;
    if (st == SCIP_STATUS_TIMELIMIT)
        stats->status = 1;
    else if (st == SCIP_STATUS_NODELIMIT || st == SCIP_STATUS_TOTALNODELIMIT)
        stats->status = 2;

done:
    if (s) {
        if (cons)
            for (size_t i = 0; i < m->rows; ++i)
                if (cons[i])
                    (void)SCIPreleaseCons(s, &cons[i]);
        if (vars)
            for (size_t j = 0; j < m->cols; ++j)
                if (vars[j])
                    (void)SCIPreleaseVar(s, &vars[j]);
        if (SCIPfree(&s) != SCIP_OKAY && result == 2)
            result = 0;
    }
    free(cons);
    free(rvars);
    free(rvals);
    free(rcur);
    free(rstart);
    free(active);
    free(vars);
    return result;
}

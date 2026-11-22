#include "scip_presol_coefftight.h"
#include "scip/cons_linear.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define COEFF_MAX_CONSS 300000
#define COEFF_MAX_ROWLEN 4096
#define COEFF_MAX_MAG 1e10

static SCIP_DECL_PRESOLEXEC(presolExecCoefftight) {
    (void)presol;
    (void)nrounds;
    (void)presoltiming;
    (void)nnewfixedvars;
    (void)nnewaggrvars;
    (void)nnewchgvartypes;
    (void)nnewchgbds;
    (void)nnewholes;
    (void)nnewdelconss;
    (void)nnewaddconss;
    (void)nnewupgdconss;
    (void)nnewchgcoefs;
    (void)nnewchgsides;
    (void)nfixedvars;
    (void)naggrvars;
    (void)nchgvartypes;
    (void)nchgbds;
    (void)naddholes;
    (void)ndelconss;
    (void)naddconss;
    (void)nupgdconss;

    *result = SCIP_DIDNOTRUN;
    if (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING || SCIPisStopped(scip))
        return SCIP_OKAY;

    SCIP_CONSHDLR *linhdlr = SCIPfindConshdlr(scip, "linear");
    if (!linhdlr)
        return SCIP_OKAY;

    int nconss = SCIPconshdlrGetNConss(linhdlr);
    if (nconss <= 0 || nconss > COEFF_MAX_CONSS)
        return SCIP_OKAY;

    SCIP_CONS **conss = SCIPconshdlrGetConss(linhdlr);
    if (!conss)
        return SCIP_OKAY;

    /* Snapshot the current constraint list since adding replacement constraints
     * appends to the handler's internal array. */
    SCIP_CONS **snap = malloc((size_t)nconss * sizeof *snap);
    SCIP_VAR **buf_vars = malloc((size_t)COEFF_MAX_ROWLEN * sizeof *buf_vars);
    SCIP_Real *buf_vals = malloc((size_t)COEFF_MAX_ROWLEN * sizeof *buf_vals);
    if (!snap || !buf_vars || !buf_vals) {
        free(snap);
        free(buf_vars);
        free(buf_vals);
        return SCIP_OKAY;
    }
    memcpy(snap, conss, (size_t)nconss * sizeof *snap);

    *result = SCIP_DIDNOTFIND;
    double inf = SCIPinfinity(scip);
    double feastol = SCIPfeastol(scip);
    int local_chgcoefs = 0;
    int local_chgsides = 0;
    SCIP_RETCODE rc = SCIP_OKAY;

    for (int i = 0; i < nconss; ++i) {
        if ((i & 1023) == 0 && SCIPisStopped(scip))
            break;
        SCIP_CONS *cons = snap[i];
        if (!cons || !SCIPconsIsActive(cons) || SCIPconsIsLocal(cons) ||
            SCIPconsIsModifiable(cons))
            continue;

        double lhs = SCIPgetLhsLinear(scip, cons);
        double rhs = SCIPgetRhsLinear(scip, cons);
        int has_lhs = (lhs > -inf) ? 1 : 0;
        int has_rhs = (rhs < inf) ? 1 : 0;
        if (has_lhs == has_rhs)
            continue; /* skip free, ranged, and equality rows */

        int len = SCIPgetNVarsLinear(scip, cons);
        if (len < 2 || len > COEFF_MAX_ROWLEN)
            continue;
        SCIP_VAR **vars = SCIPgetVarsLinear(scip, cons);
        SCIP_Real *vals = SCIPgetValsLinear(scip, cons);
        if (!vars || !vals)
            continue;

        int row_chgcoefs = 0;
        int row_chgsides = 0;

        if (has_rhs) {
            if (fabs(rhs) > COEFF_MAX_MAG)
                continue;
            double maxact = 0.0;
            int has_bin = 0;
            int ok = 1;
            for (int k = 0; k < len; ++k) {
                SCIP_VAR *v = vars[k];
                double a = vals[k];
                if (!v || !SCIPvarIsActive(v) || !isfinite(a) || fabs(a) > COEFF_MAX_MAG) {
                    ok = 0;
                    break;
                }
                double lb = SCIPvarGetLbGlobal(v);
                double ub = SCIPvarGetUbGlobal(v);
                if (a > 0.0) {
                    if (ub >= inf || fabs(ub) > COEFF_MAX_MAG) {
                        ok = 0;
                        break;
                    }
                    maxact += a * ub;
                } else if (a < 0.0) {
                    if (lb <= -inf || fabs(lb) > COEFF_MAX_MAG) {
                        ok = 0;
                        break;
                    }
                    maxact += a * lb;
                }
                if (SCIPvarIsBinary(v) && lb == 0.0 && ub == 1.0)
                    has_bin = 1;
            }
            if (!ok || !has_bin || !isfinite(maxact) || fabs(maxact) > COEFF_MAX_MAG)
                continue;
            if (maxact <= rhs + feastol)
                continue;

            memcpy(buf_vars, vars, (size_t)len * sizeof *buf_vars);
            memcpy(buf_vals, vals, (size_t)len * sizeof *buf_vals);

            /* Conservative rounding safety margin so valid integer points are never cut off. */
            double eps = 1e-9 * (1.0 + fabs(maxact) + fabs(rhs));
            for (int k = 0; k < len; ++k) {
                SCIP_VAR *v = buf_vars[k];
                double a = buf_vals[k];
                if (!SCIPvarIsBinary(v) || SCIPvarGetLbGlobal(v) != 0.0 ||
                    SCIPvarGetUbGlobal(v) != 1.0)
                    continue;
                double min_d = fmax(1e-6, 1e-6 * fabs(a));
                if (a > 1e-9) {
                    double rem = maxact - a;
                    double d = rhs - rem - eps;
                    if (d > min_d && a - d > 1e-7) {
                        buf_vals[k] = a - d;
                        rhs -= d;
                        maxact -= d;
                        row_chgcoefs++;
                        row_chgsides++;
                    }
                } else if (a < -1e-9) {
                    double rem = maxact + a;
                    double d = rhs - rem - eps;
                    if (d > min_d && a + d < -1e-7) {
                        buf_vals[k] = a + d;
                        row_chgcoefs++;
                    }
                }
            }
        } else {
            if (fabs(lhs) > COEFF_MAX_MAG)
                continue;
            double minact = 0.0;
            int has_bin = 0;
            int ok = 1;
            for (int k = 0; k < len; ++k) {
                SCIP_VAR *v = vars[k];
                double a = vals[k];
                if (!v || !SCIPvarIsActive(v) || !isfinite(a) || fabs(a) > COEFF_MAX_MAG) {
                    ok = 0;
                    break;
                }
                double lb = SCIPvarGetLbGlobal(v);
                double ub = SCIPvarGetUbGlobal(v);
                if (a > 0.0) {
                    if (lb <= -inf || fabs(lb) > COEFF_MAX_MAG) {
                        ok = 0;
                        break;
                    }
                    minact += a * lb;
                } else if (a < 0.0) {
                    if (ub >= inf || fabs(ub) > COEFF_MAX_MAG) {
                        ok = 0;
                        break;
                    }
                    minact += a * ub;
                }
                if (SCIPvarIsBinary(v) && lb == 0.0 && ub == 1.0)
                    has_bin = 1;
            }
            if (!ok || !has_bin || !isfinite(minact) || fabs(minact) > COEFF_MAX_MAG)
                continue;
            if (minact >= lhs - feastol)
                continue;

            memcpy(buf_vars, vars, (size_t)len * sizeof *buf_vars);
            memcpy(buf_vals, vals, (size_t)len * sizeof *buf_vals);

            double eps = 1e-9 * (1.0 + fabs(minact) + fabs(lhs));
            for (int k = 0; k < len; ++k) {
                SCIP_VAR *v = buf_vars[k];
                double a = buf_vals[k];
                if (!SCIPvarIsBinary(v) || SCIPvarGetLbGlobal(v) != 0.0 ||
                    SCIPvarGetUbGlobal(v) != 1.0)
                    continue;
                double min_d = fmax(1e-6, 1e-6 * fabs(a));
                if (a > 1e-9) {
                    double rem = minact + a;
                    double d = rem - lhs - eps;
                    if (d > min_d && a - d > 1e-7) {
                        buf_vals[k] = a - d;
                        row_chgcoefs++;
                    }
                } else if (a < -1e-9) {
                    double rem = minact - a;
                    double d = rem - lhs - eps;
                    if (d > min_d && a + d < -1e-7) {
                        buf_vals[k] = a + d;
                        lhs += d;
                        minact += d;
                        row_chgcoefs++;
                        row_chgsides++;
                    }
                }
            }
        }

        if (row_chgcoefs > 0) {
            SCIP_CONS *newcons = NULL;
            rc = SCIPcreateConsLinear(
                scip, &newcons, SCIPconsGetName(cons), len, buf_vars, buf_vals, lhs, rhs,
                SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                SCIPconsIsStickingAtNode(cons));
            if (rc != SCIP_OKAY)
                break;
            rc = SCIPdelCons(scip, cons);
            if (rc != SCIP_OKAY) {
                (void)SCIPreleaseCons(scip, &newcons);
                break;
            }
            rc = SCIPaddCons(scip, newcons);
            (void)SCIPreleaseCons(scip, &newcons);
            if (rc != SCIP_OKAY)
                break;
            local_chgcoefs += row_chgcoefs;
            local_chgsides += row_chgsides;
        }
    }

    free(snap);
    free(buf_vars);
    free(buf_vals);
    if (rc != SCIP_OKAY)
        return rc;

    if (local_chgcoefs > 0 || local_chgsides > 0) {
        *nchgcoefs += local_chgcoefs;
        *nchgsides += local_chgsides;
        *result = SCIP_SUCCESS;
    }
    return SCIP_OKAY;
}

SCIP_RETCODE lab_scip_include_presol_coefftight(SCIP *scip) {
    SCIP_PRESOL *presol = NULL;
    SCIP_CALL(SCIPincludePresolBasic(
        scip, &presol, "coefftight", "pure-C17 PaPILO-style binary coefficient tightening",
        -9000, 3, SCIP_PRESOLTIMING_FAST | SCIP_PRESOLTIMING_MEDIUM,
        presolExecCoefftight, NULL));
    return SCIP_OKAY;
}

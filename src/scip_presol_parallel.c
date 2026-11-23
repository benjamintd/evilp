#include "scip_presol_parallel.h"
#include "scip/cons_linear.h"
#include "scip/pub_matrix.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define PARROW_MAX_ROWS 300000

static uint64_t mix64(uint64_t x) {
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return x;
}

static uint64_t hash_row_support(const int *cols, const SCIP_Real *vals, int len) {
    uint64_t h = mix64((uint64_t)(unsigned int)len);
    double v0 = vals[0];
    if (fabs(v0) <= 1e-15)
        return 0;
    for (int k = 0; k < len; ++k) {
        h ^= mix64((uint64_t)(unsigned int)cols[k] + (uint64_t)(unsigned int)k * UINT64_C(0x9e3779b97f4a7c15));
        double r = vals[k] / v0;
        /* Coarse quantized bucket so proportional rows with 1e-13 noise still collide */
        int64_t q = (fabs(r) < 1e9) ? (int64_t)llround(r * 1024.0) : (r > 0 ? 1 : -1);
        h ^= mix64((uint64_t)q + UINT64_C(0x517cc1b727220a95));
    }
    return h ? h : 1;
}

/* fma rounds the product difference only once. Away from underflow,
 * a nonzero difference between binary64 products cannot round to zero. */
static int exact_product(double a, double b, double value) {
    if (!isfinite(a) || !isfinite(b) || !isfinite(value))
        return 0;
    if (value == 0.0)
        return a == 0.0 || b == 0.0;
    return fabs(value) >= DBL_MIN / (DBL_EPSILON * DBL_EPSILON) &&
           fma(a, b, -value) == 0.0;
}

static int rows_proportional(const int *c1, const SCIP_Real *v1, const int *c2,
                             const SCIP_Real *v2, int len, double *out_alpha) {
    if (fabs(v1[0]) <= 1e-15 || fabs(v2[0]) <= 1e-15)
        return 0;
    for (int k = 0; k < len; ++k)
        if (c1[k] != c2[k])
            return 0;
    double alpha = v2[0] / v1[0];
    if (!isfinite(alpha) || fabs(alpha) < 1e-9 || fabs(alpha) > 1e9)
        return 0;
    for (int k = 0; k < len; ++k) {
        if (!exact_product(alpha, v1[k], v2[k]))
            return 0;
    }
    *out_alpha = alpha;
    return 1;
}

static int flags_dominate(SCIP_CONS *keep, SCIP_CONS *drop) {
    if (SCIPconsIsLocal(keep) || SCIPconsIsLocal(drop))
        return 0;
    if (SCIPconsIsModifiable(keep) || SCIPconsIsModifiable(drop))
        return 0;
    if (!SCIPconsIsChecked(keep) && SCIPconsIsChecked(drop))
        return 0;
    if (!SCIPconsIsEnforced(keep) && SCIPconsIsEnforced(drop))
        return 0;
    if (!SCIPconsIsInitial(keep) && SCIPconsIsInitial(drop))
        return 0;
    if (SCIPconsIsDynamic(keep) && !SCIPconsIsDynamic(drop))
        return 0;
    return 1;
}

static int is_linear_exact_sides(SCIP *scip, SCIP_CONS *cons, int len, double mlhs, double mrhs) {
    SCIP_CONSHDLR *hdlr = SCIPconsGetHdlr(cons);
    if (!hdlr || strcmp(SCIPconshdlrGetName(hdlr), "linear") != 0)
        return 0;
    if (SCIPgetNVarsLinear(scip, cons) != len)
        return 0;
    double inf = SCIPinfinity(scip);
    double clhs = SCIPgetLhsLinear(scip, cons);
    double crhs = SCIPgetRhsLinear(scip, cons);
    int lhs_ok = (clhs <= -inf && mlhs <= -inf) ||
                 (clhs > -inf && mlhs > -inf && clhs == mlhs);
    int rhs_ok = (crhs >= inf && mrhs >= inf) ||
                 (crhs < inf && mrhs < inf && crhs == mrhs);
    return lhs_ok && rhs_ok;
}

static SCIP_DECL_PRESOLEXEC(presolExecParallel) {
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
    (void)naggrvars;
    (void)nchgvartypes;
    (void)naddholes;
    (void)nupgdconss;

    *result = SCIP_DIDNOTRUN;
    if (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING || SCIPisStopped(scip))
        return SCIP_OKAY;
    if (SCIPgetNConss(scip) < 2 || SCIPgetNConss(scip) > PARROW_MAX_ROWS)
        return SCIP_OKAY;

    SCIP_MATRIX *matrix = NULL;
    SCIP_Bool initialized = FALSE;
    SCIP_Bool complete = FALSE;
    SCIP_Bool infeasible = FALSE;

    SCIP_CALL(SCIPmatrixCreate(scip, &matrix, TRUE, &initialized, &complete, &infeasible,
                               naddconss, ndelconss, nchgcoefs, nchgbds, nfixedvars));
    if (infeasible) {
        if (initialized)
            SCIPmatrixFree(scip, &matrix);
        *result = SCIP_CUTOFF;
        return SCIP_OKAY;
    }
    if (!initialized)
        return SCIP_OKAY;

    int nrows = SCIPmatrixGetNRows(matrix);
    if (nrows < 2 || nrows > PARROW_MAX_ROWS) {
        SCIPmatrixFree(scip, &matrix);
        return SCIP_OKAY;
    }

    *result = SCIP_DIDNOTFIND;
    size_t hcap = 16;
    while (hcap < (size_t)nrows * 2)
        hcap <<= 1;
    int *head = malloc(hcap * sizeof *head);
    int *next = malloc((size_t)nrows * sizeof *next);
    uint64_t *rhash = malloc((size_t)nrows * sizeof *rhash);
    unsigned char *deleted = calloc((size_t)nrows, 1);
    double *cur_lhs = malloc((size_t)nrows * sizeof *cur_lhs);
    double *cur_rhs = malloc((size_t)nrows * sizeof *cur_rhs);
    if (!head || !next || !rhash || !deleted || !cur_lhs || !cur_rhs) {
        free(head);
        free(next);
        free(rhash);
        free(deleted);
        free(cur_lhs);
        free(cur_rhs);
        SCIPmatrixFree(scip, &matrix);
        return SCIP_OKAY;
    }
    for (size_t i = 0; i < hcap; ++i)
        head[i] = -1;

    double inf = SCIPinfinity(scip);
    double feastol = SCIPfeastol(scip);

    for (int r = 0; r < nrows; ++r) {
        int len = SCIPmatrixGetRowNNonzs(matrix, r);
        SCIP_CONS *cons = SCIPmatrixGetCons(matrix, r);
        cur_lhs[r] = SCIPmatrixGetRowLhs(matrix, r);
        cur_rhs[r] = SCIPmatrixIsRowRhsInfinity(matrix, r) ? inf : SCIPmatrixGetRowRhs(matrix, r);
        if (len < 2 || !cons || !SCIPconsIsActive(cons) || SCIPconsIsModifiable(cons)) {
            rhash[r] = 0;
            continue;
        }
        const int *cols = SCIPmatrixGetRowIdxPtr(matrix, r);
        const SCIP_Real *vals = SCIPmatrixGetRowValPtr(matrix, r);
        /* Verify ascending column indices for fast exact comparison */
        int sorted = 1;
        for (int k = 1; k < len; ++k) {
            if (cols[k] <= cols[k - 1]) {
                sorted = 0;
                break;
            }
        }
        if (!sorted) {
            rhash[r] = 0;
            continue;
        }
        uint64_t h = hash_row_support(cols, vals, len);
        rhash[r] = h;
        if (!h)
            continue;
        size_t b = (size_t)(h & (uint64_t)(hcap - 1));
        next[r] = head[b];
        head[b] = r;
    }

    int local_del = 0;
    int local_chg = 0;
    SCIP_Bool cutoff = FALSE;

    for (size_t b = 0; b < hcap && !cutoff; ++b) {
        for (int r1 = head[b]; r1 >= 0 && !cutoff; r1 = next[r1]) {
            if (deleted[r1])
                continue;
            SCIP_CONS *c1 = SCIPmatrixGetCons(matrix, r1);
            if (!c1 || !SCIPconsIsActive(c1))
                continue;
            int len = SCIPmatrixGetRowNNonzs(matrix, r1);
            const int *cols1 = SCIPmatrixGetRowIdxPtr(matrix, r1);
            const SCIP_Real *vals1 = SCIPmatrixGetRowValPtr(matrix, r1);
            int c1_lin = is_linear_exact_sides(scip, c1, len, cur_lhs[r1], cur_rhs[r1]);

            for (int r2 = next[r1]; r2 >= 0 && !cutoff; r2 = next[r2]) {
                if (deleted[r2] || rhash[r2] != rhash[r1])
                    continue;
                if (SCIPmatrixGetRowNNonzs(matrix, r2) != len)
                    continue;
                SCIP_CONS *c2 = SCIPmatrixGetCons(matrix, r2);
                if (!c2 || !SCIPconsIsActive(c2))
                    continue;
                const int *cols2 = SCIPmatrixGetRowIdxPtr(matrix, r2);
                const SCIP_Real *vals2 = SCIPmatrixGetRowValPtr(matrix, r2);
                double alpha = 0.0;
                if (!rows_proportional(cols1, vals1, cols2, vals2, len, &alpha))
                    continue;

                /* Scale row r2 bounds into row r1's coordinate system */
                double l2_raw = cur_lhs[r2];
                double u2_raw = cur_rhs[r2];
                double l2, u2;
                if (alpha > 0.0) {
                    l2 = (l2_raw > -inf) ? (l2_raw / alpha) : -inf;
                    u2 = (u2_raw < inf) ? (u2_raw / alpha) : inf;
                } else {
                    l2 = (u2_raw < inf) ? (u2_raw / alpha) : -inf;
                    u2 = (l2_raw > -inf) ? (l2_raw / alpha) : inf;
                }

                if ((l2_raw > -inf && (fabs(l2_raw / alpha) >= inf ||
                        !exact_product(l2_raw / alpha, alpha, l2_raw))) ||
                    (u2_raw < inf && (fabs(u2_raw / alpha) >= inf ||
                        !exact_product(u2_raw / alpha, alpha, u2_raw))))
                    continue;
                double l1 = cur_lhs[r1];
                double u1 = cur_rhs[r1];
                double new_l = (l1 > l2) ? l1 : l2;
                double new_u = (u1 < u2) ? u1 : u2;

                if (new_l > -inf && new_u < inf && new_l > new_u + feastol) {
                    cutoff = TRUE;
                    break;
                }

                if (new_l > new_u)
                    continue;
                /* Case 1: r1 already dominates r2 (no side change needed on r1) */
                int l1_dom = (l2 <= -inf) || (l1 > -inf && l1 >= l2);
                int u1_dom = (u2 >= inf) || (u1 < inf && u1 <= u2);
                if (l1_dom && u1_dom && flags_dominate(c1, c2)) {
                    SCIP_CALL(SCIPdelCons(scip, c2));
                    deleted[r2] = 1;
                    local_del++;
                    continue;
                }

                /* Case 2: r2 dominates r1 */
                int l2_dom = (l1 <= -inf) || (l2 > -inf && l2 >= l1);
                int u2_dom = (u1 >= inf) || (u2 < inf && u2 <= u1);
                if (l2_dom && u2_dom && flags_dominate(c2, c1)) {
                    SCIP_CALL(SCIPdelCons(scip, c1));
                    deleted[r1] = 1;
                    local_del++;
                    break;
                }

                /* Case 3: merge strict sides into linear c1 and delete c2 */
                if (c1_lin && flags_dominate(c1, c2)) {
                    int changed = 0;
                    if (new_l > -inf && (l1 <= -inf || new_l > l1)) {
                        SCIP_CALL(SCIPchgLhsLinear(scip, c1, new_l));
                        cur_lhs[r1] = SCIPgetLhsLinear(scip, c1);
                        changed = 1;
                    }
                    if (new_u < inf && (u1 >= inf || new_u < u1)) {
                        SCIP_CALL(SCIPchgRhsLinear(scip, c1, new_u));
                        cur_rhs[r1] = SCIPgetRhsLinear(scip, c1);
                        changed = 1;
                    }
                    if (changed)
                        local_chg++;
                    if (cur_lhs[r1] == new_l && cur_rhs[r1] == new_u) {
                        SCIP_CALL(SCIPdelCons(scip, c2));
                        deleted[r2] = 1;
                        local_del++;
                    }
                }
            }
        }
    }

    free(head);
    free(next);
    free(rhash);
    free(deleted);
    free(cur_lhs);
    free(cur_rhs);
    SCIPmatrixFree(scip, &matrix);

    if (cutoff) {
        *result = SCIP_CUTOFF;
    } else if (local_del > 0 || local_chg > 0) {
        *ndelconss += local_del;
        *nchgsides += local_chg;
        *result = SCIP_SUCCESS;
    }
    return SCIP_OKAY;
}

SCIP_RETCODE lab_scip_include_presol_parallel(SCIP *scip) {
    SCIP_PRESOL *presol = NULL;
    SCIP_CALL(SCIPincludePresolBasic(
        scip, &presol, "parrow", "pure-C17 PaPILO-style parallel-row presolver",
        -10000, 2, SCIP_PRESOLTIMING_MEDIUM, presolExecParallel, NULL));
    return SCIP_OKAY;
}

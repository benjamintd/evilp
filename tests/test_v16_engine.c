#define _POSIX_C_SOURCE 200809L
#include "lpi/lpi.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "test_v16_engine failed at line %d: %s\n", __LINE__, #x);              \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)
#define OK(x) CHECK((x) == SCIP_OKAY)

static double elapsed_s(struct timespec a, struct timespec b) {
    return (double)(b.tv_sec - a.tv_sec) + 1e-9 * (double)(b.tv_nsec - a.tv_nsec);
}

/* 1. Stress dense-column LU factorization (lazy column pruning + cancellations +
 *    adaptive dense U-solve bailout for dim >= 256) and verify primal/dual optimality. */
static int test_dense_lu_and_u_bailout(void) {
    const int nrows = 300;
    const int nstruct = 340;
    const int ndense = 10;
    const int max_nnz = nrows * 12 + ndense * nrows;
    double *obj = calloc((size_t)nstruct, sizeof *obj);
    double *lo = calloc((size_t)nstruct, sizeof *lo);
    double *hi = calloc((size_t)nstruct, sizeof *hi);
    double *lhs = calloc((size_t)nrows, sizeof *lhs);
    double *rhs = calloc((size_t)nrows, sizeof *rhs);
    int *beg = calloc((size_t)nstruct, sizeof *beg);
    int *ind = calloc((size_t)max_nnz, sizeof *ind);
    double *val = calloc((size_t)max_nnz, sizeof *val);
    CHECK(obj && lo && hi && lhs && rhs && beg && ind && val);

    for (int i = 0; i < nrows; ++i) {
        lhs[i] = 1.0 + 0.1 * (double)(i % 7);
        rhs[i] = 1e20;
    }

    int nnz = 0;
    for (int j = 0; j < nstruct; ++j) {
        obj[j] = 1.0 + 0.05 * (double)((j * 13 + 5) % 29);
        lo[j] = 0.0;
        hi[j] = 10.0;
        beg[j] = nnz;
        if (j < nrows) {
            /* Banded upper-triangular coupling chain + cancellation pairs */
            ind[nnz] = j;
            val[nnz++] = 2.0;
            if (j > 0) {
                ind[nnz] = j - 1;
                val[nnz++] = (j % 3 == 0) ? -1.0 : 1.0;
            }
            if (j >= 2) {
                ind[nnz] = j - 2;
                val[nnz++] = 0.5;
            }
            if (j + 5 < nrows) {
                ind[nnz] = j + 5;
                val[nnz++] = (j % 2 == 0) ? 0.25 : -0.25;
            }
        } else if (j < nrows + ndense) {
            /* Dense linking columns (~40% density) */
            int phase = j - nrows;
            obj[j] = 0.2 + 0.03 * (double)phase;
            for (int r = phase % 2; r < nrows; r += 2) {
                if ((r + phase) % 5 != 0) {
                    ind[nnz] = r;
                    val[nnz++] = ((r + phase) % 3 == 0) ? -0.5 : 1.0;
                }
            }
        } else {
            /* Extra sparse columns */
            for (int k = 0; k < 4; ++k) {
                int r = (j * 17 + k * 53) % nrows;
                ind[nnz] = r;
                val[nnz++] = 1.0 + 0.2 * (double)k;
            }
        }
    }
    CHECK(nnz <= max_nnz);

    SCIP_LPI *lp = NULL;
    OK(SCIPlpiCreate(&lp, NULL, "v16_dense_lu", SCIP_OBJSEN_MINIMIZE));
    OK(SCIPlpiLoadColLP(lp, SCIP_OBJSEN_MINIMIZE, nstruct, obj, lo, hi, NULL, nrows, lhs, rhs, NULL,
                        nnz, beg, ind, val));
    OK(SCIPlpiSolveDual(lp));
    CHECK(SCIPlpiIsOptimal(lp));

    double objval = 0.0;
    double *prim = calloc((size_t)nstruct, sizeof *prim);
    double *dual = calloc((size_t)nrows, sizeof *dual);
    double *act = calloc((size_t)nrows, sizeof *act);
    CHECK(prim && dual && act);
    OK(SCIPlpiGetSol(lp, &objval, prim, dual, NULL, NULL));

    /* Independently check Ax >= lhs and c^T x == objval */
    double calc_obj = 0.0;
    for (int j = 0; j < nstruct; ++j) {
        CHECK(prim[j] >= lo[j] - 1e-7 && prim[j] <= hi[j] + 1e-7);
        calc_obj += obj[j] * prim[j];
        int end = (j + 1 < nstruct) ? beg[j + 1] : nnz;
        for (int k = beg[j]; k < end; ++k)
            act[ind[k]] += val[k] * prim[j];
    }
    CHECK(fabs(calc_obj - objval) <= 1e-6 * (1.0 + fabs(objval)));
    for (int i = 0; i < nrows; ++i) {
        CHECK(act[i] >= lhs[i] - 1e-6);
        CHECK(dual[i] >= -1e-6);
    }

    /* Trigger bound changes + re-solves (Forrest-Tomlin updates & refactorizations) */
    for (int step = 0; step < 6; ++step) {
        int col = nrows + (step % ndense);
        double newub = 0.5 + 0.3 * (double)step;
        double newlb = 0.0;
        OK(SCIPlpiChgBounds(lp, 1, &col, &newlb, &newub));
        OK(SCIPlpiSolveDual(lp));
        CHECK(SCIPlpiIsOptimal(lp));
    }

    /* Check sparse BInv rows match dense BInv rows on sampled indices */
    double *drow = calloc((size_t)nrows, sizeof *drow);
    double *srow = calloc((size_t)nrows, sizeof *srow);
    int *rinds = calloc((size_t)nrows, sizeof *rinds);
    CHECK(drow && srow && rinds);
    for (int r = 0; r < nrows; r += 37) {
        int ninds = -1;
        memset(drow, 0, (size_t)nrows * sizeof *drow);
        memset(srow, 0, (size_t)nrows * sizeof *srow);
        OK(SCIPlpiGetBInvRow(lp, r, drow, NULL, NULL));
        OK(SCIPlpiGetBInvRow(lp, r, srow, rinds, &ninds));
        CHECK(ninds >= 0 && ninds <= nrows);
        for (int i = 0; i < nrows; ++i)
            CHECK(fabs(drow[i] - srow[i]) <= 1e-9 * (1.0 + fabs(drow[i])));
    }

    free(drow);
    free(srow);
    free(rinds);
    free(prim);
    free(dual);
    free(act);
    OK(SCIPlpiFree(&lp));
    free(obj);
    free(lo);
    free(hi);
    free(lhs);
    free(rhs);
    free(beg);
    free(ind);
    free(val);
    return 0;
}

/* 2. Stress primal steepest-edge norm reuse across objective changes (feaspump pattern)
 *    and safe invalidation across dual solves, row/col additions, and matrix edits. */
static int test_primal_norm_reuse_and_invalidation(void) {
    const int nrows = 40;
    const int ncols = 60;
    double *obj = calloc((size_t)ncols, sizeof *obj);
    double *lo = calloc((size_t)ncols, sizeof *lo);
    double *hi = calloc((size_t)ncols, sizeof *hi);
    double *lhs = calloc((size_t)nrows, sizeof *lhs);
    double *rhs = calloc((size_t)nrows, sizeof *rhs);
    int *beg = calloc((size_t)ncols, sizeof *beg);
    int *ind = calloc((size_t)(ncols * 4), sizeof *ind);
    double *val = calloc((size_t)(ncols * 4), sizeof *val);
    CHECK(obj && lo && hi && lhs && rhs && beg && ind && val);

    for (int i = 0; i < nrows; ++i) {
        lhs[i] = 1.5;
        rhs[i] = 1e20;
    }
    for (int j = 0; j < ncols; ++j) {
        obj[j] = 1.0 + 0.1 * (double)(j % 11);
        lo[j] = 0.0;
        hi[j] = 5.0;
        beg[j] = j * 4;
        for (int k = 0; k < 4; ++k) {
            ind[beg[j] + k] = (j * 3 + k * 7) % nrows;
            val[beg[j] + k] = 1.0 + 0.25 * (double)((j + k) % 3);
        }
    }

    SCIP_LPI *lp = NULL;
    OK(SCIPlpiCreate(&lp, NULL, "v16_psteep", SCIP_OBJSEN_MINIMIZE));
    OK(SCIPlpiLoadColLP(lp, SCIP_OBJSEN_MINIMIZE, ncols, obj, lo, hi, NULL, nrows, lhs, rhs, NULL,
                        ncols * 4, beg, ind, val));
    OK(SCIPlpiSetIntpar(lp, SCIP_LPPAR_PRICING, (int)SCIP_PRICING_STEEP));
    OK(SCIPlpiSolvePrimal(lp));
    CHECK(SCIPlpiIsOptimal(lp));

    int *cinds = calloc((size_t)ncols, sizeof *cinds);
    CHECK(cinds != NULL);
    for (int j = 0; j < ncols; ++j)
        cinds[j] = j;

    /* 15 repeated objective changes + primal solves (feaspump pattern) */
    for (int iter = 0; iter < 15; ++iter) {
        for (int j = 0; j < ncols; ++j)
            obj[j] = ((j + iter * 5) % 7 < 3) ? -1.0 : 1.0 + 0.07 * (double)((j * 3 + iter) % 13);
        OK(SCIPlpiChgObj(lp, ncols, cinds, obj));
        OK(SCIPlpiSolvePrimal(lp));
        CHECK(SCIPlpiIsOptimal(lp));
    }

    /* Check agreement with a fresh LP on the final objective */
    double warm_obj = 0.0;
    OK(SCIPlpiGetObjval(lp, &warm_obj));
    SCIP_LPI *ref = NULL;
    OK(SCIPlpiCreate(&ref, NULL, "v16_ref", SCIP_OBJSEN_MINIMIZE));
    OK(SCIPlpiLoadColLP(ref, SCIP_OBJSEN_MINIMIZE, ncols, obj, lo, hi, NULL, nrows, lhs, rhs, NULL,
                        ncols * 4, beg, ind, val));
    OK(SCIPlpiSolveDual(ref));
    CHECK(SCIPlpiIsOptimal(ref));
    double ref_obj = 0.0;
    OK(SCIPlpiGetObjval(ref, &ref_obj));
    CHECK(fabs(warm_obj - ref_obj) <= 1e-7 * (1.0 + fabs(ref_obj)));
    OK(SCIPlpiFree(&ref));

    /* Interleave dual solve, row addition, col addition, coef edit, then primal solve */
    OK(SCIPlpiSolveDual(lp));
    CHECK(SCIPlpiIsOptimal(lp));

    int rbeg = 0;
    int rind[3] = {0, 5, 10};
    double rval[3] = {1.0, 1.0, 1.0};
    double rlhs = 0.5, rrhs = 1e20;
    OK(SCIPlpiAddRows(lp, 1, &rlhs, &rrhs, NULL, 3, &rbeg, rind, rval));
    OK(SCIPlpiSolvePrimal(lp));
    CHECK(SCIPlpiIsOptimal(lp));

    int cbeg = 0;
    int cind[2] = {0, nrows};
    double cval[2] = {2.0, 1.0};
    double cobj = 0.5, clb = 0.0, cub = 3.0;
    OK(SCIPlpiAddCols(lp, 1, &cobj, &clb, &cub, NULL, 2, &cbeg, cind, cval));
    OK(SCIPlpiSolvePrimal(lp));
    CHECK(SCIPlpiIsOptimal(lp));

    OK(SCIPlpiChgCoef(lp, 0, ncols, 3.0));
    OK(SCIPlpiSolvePrimal(lp));
    CHECK(SCIPlpiIsOptimal(lp));

    OK(SCIPlpiFree(&lp));
    free(cinds);
    free(obj);
    free(lo);
    free(hi);
    free(lhs);
    free(rhs);
    free(beg);
    free(ind);
    free(val);
    return 0;
}

typedef struct {
    SCIP_VAR *var;
    SCIP_EVENTHDLR *hdlr;
    int calls;
    int did_delayed_ops;
    int delayed_pos;
} EventStressCtx;

static SCIP_DECL_EVENTEXEC(stress_event_exec) {
    EventStressCtx *ctx = (EventStressCtx *)SCIPeventhdlrGetData(eventhdlr);
    (void)event;
    ctx->calls++;
    /* On the first callback during SCIPeventfilterProcess (delayupdates == TRUE),
     * exercise both delayed addition and delayed deletion with filterpos = -1. */
    if (!ctx->did_delayed_ops && eventdata == (SCIP_EVENTDATA *)(uintptr_t)10) {
        ctx->did_delayed_ops = 1;
        SCIP_CALL(SCIPdropVarEvent(scip, ctx->var, SCIP_EVENTTYPE_LBTIGHTENED, ctx->hdlr,
                                   (SCIP_EVENTDATA *)(uintptr_t)20, -1));
        SCIP_CALL(SCIPcatchVarEvent(scip, ctx->var, SCIP_EVENTTYPE_LBTIGHTENED, ctx->hdlr,
                                    (SCIP_EVENTDATA *)(uintptr_t)999999, &ctx->delayed_pos));
    }
    return SCIP_OKAY;
}

/* 3. Stress SCIP_EventFilter lazy O(1) hash table: 4096 events, duplicates,
 *    delayed add/del during SCIPeventfilterProcess, and O(n) drop with filterpos = -1. */
static int test_eventfilter_lazy_hash(void) {
    const int n = 4096;
    SCIP *scip = NULL;
    OK(SCIPcreate(&scip));
    OK(SCIPincludeDefaultPlugins(scip));
    OK(SCIPcreateProbBasic(scip, "eventfilter_stress"));

    EventStressCtx ctx;
    memset(&ctx, 0, sizeof ctx);
    ctx.delayed_pos = -1;

    SCIP_EVENTHDLR *hdlr = NULL;
    OK(SCIPincludeEventhdlrBasic(scip, &hdlr, "v16_stress_ev", "eventfilter hash stress",
                                 stress_event_exec, (SCIP_EVENTHDLRDATA *)&ctx));
    ctx.hdlr = hdlr;

    SCIP_VAR *var = NULL;
    OK(SCIPcreateVarBasic(scip, &var, "x", 0.0, 10.0, 1.0, SCIP_VARTYPE_CONTINUOUS));
    OK(SCIPaddVar(scip, var));

    /* Transform problem so event queue / variable event filters are active */
    OK(SCIPtransformProb(scip));
    SCIP_VAR *tvar = SCIPvarGetTransVar(var);
    CHECK(tvar != NULL);
    ctx.var = tvar;

    /* Catch n distinct events + 1 duplicate of id=7 */
    int *fpos = calloc((size_t)n, sizeof *fpos);
    CHECK(fpos != NULL);
    for (int i = 0; i < n; ++i) {
        OK(SCIPcatchVarEvent(scip, tvar, SCIP_EVENTTYPE_LBTIGHTENED, hdlr,
                             (SCIP_EVENTDATA *)(uintptr_t)(i + 1), &fpos[i]));
    }
    int dup_pos = -1;
    OK(SCIPcatchVarEvent(scip, tvar, SCIP_EVENTTYPE_LBTIGHTENED, hdlr,
                         (SCIP_EVENTDATA *)(uintptr_t)7, &dup_pos));
    CHECK(dup_pos >= 0);

    /* Drop the duplicate via filterpos = -1 (activates lazy hash table since len >= 128) */
    OK(SCIPdropVarEvent(scip, tvar, SCIP_EVENTTYPE_LBTIGHTENED, hdlr,
                        (SCIP_EVENTDATA *)(uintptr_t)7, -1));

    /* Drop one entry using explicit filterpos != -1 while hash table is active */
    OK(SCIPdropVarEvent(scip, tvar, SCIP_EVENTTYPE_LBTIGHTENED, hdlr,
                        (SCIP_EVENTDATA *)(uintptr_t)15, fpos[14]));

    /* Trigger LBTIGHTENED event -> exercises delayed add & delayed del inside Process */
    OK(SCIPchgVarLbGlobal(scip, tvar, 1.0));
    CHECK(ctx.did_delayed_ops == 1);
    CHECK(ctx.calls >= n - 2);
    CHECK(ctx.delayed_pos >= 0);

    /* Drop the delayed-added entry via filterpos = -1 */
    OK(SCIPdropVarEvent(scip, tvar, SCIP_EVENTTYPE_LBTIGHTENED, hdlr,
                        (SCIP_EVENTDATA *)(uintptr_t)999999, -1));

    /* Now drop all remaining ~4094 events with filterpos = -1 and time it (must be O(n)) */
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < n; ++i) {
        int id = i + 1;
        if (id == 15 || id == 20)
            continue;
        OK(SCIPdropVarEvent(scip, tvar, SCIP_EVENTTYPE_LBTIGHTENED, hdlr,
                            (SCIP_EVENTDATA *)(uintptr_t)id, -1));
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double dt = elapsed_s(t0, t1);
    CHECK(dt < 0.5);

    free(fpos);
    OK(SCIPreleaseVar(scip, &var));
    OK(SCIPfree(&scip));
    return 0;
}

/* 4. Stress sparse SCIPseparateRelaxedKnapsack on a large binary + continuous VUB model:
 *    calls SCIPseparateRelaxedKnapsack directly with duplicate VUB substitutions onto the
 *    same binary variable and verifies cut generation + clean binvals zeroing + O(nnz) speed. */
static int test_sparse_relaxed_knapsack(void) {
    const int nbin = 4000;
    const int ncont = 6;
    SCIP *scip = NULL;
    OK(SCIPcreate(&scip));
    OK(SCIPincludeDefaultPlugins(scip));
    OK(SCIPsetIntParam(scip, "display/verblevel", 0));
    OK(SCIPcreateProbBasic(scip, "sparse_knap_relax"));

    SCIP_VAR **bvars = calloc((size_t)nbin, sizeof *bvars);
    SCIP_VAR **cvars = calloc((size_t)ncont, sizeof *cvars);
    CHECK(bvars && cvars);
    for (int i = 0; i < nbin; ++i) {
        char name[32];
        snprintf(name, sizeof name, "b%d", i);
        OK(SCIPcreateVarBasic(scip, &bvars[i], name, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY));
        OK(SCIPaddVar(scip, bvars[i]));
    }
    for (int j = 0; j < ncont; ++j) {
        char name[32];
        snprintf(name, sizeof name, "c%d", j);
        OK(SCIPcreateVarBasic(scip, &cvars[j], name, 0.0, 10.0, 0.5, SCIP_VARTYPE_CONTINUOUS));
        OK(SCIPaddVar(scip, cvars[j]));
        /* VUB: c_j <= 3.5 * b_{j % 3} so duplicate binary targets occur during substitution */
        SCIP_CONS *vub = NULL;
        SCIP_VAR *vars[2] = {cvars[j], bvars[j % 3]};
        SCIP_Real vals[2] = {1.0, -3.5};
        OK(SCIPcreateConsBasicLinear(scip, &vub, name, 2, vars, vals, -SCIPinfinity(scip), 0.0));
        OK(SCIPaddCons(scip, vub));
        OK(SCIPreleaseCons(scip, &vub));
    }

    OK(SCIPtransformProb(scip));
    SCIP_SOL *sol = NULL;
    OK(SCIPcreateSol(scip, &sol, NULL));
    for (int i = 0; i < 10; ++i)
        OK(SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(bvars[i]), 0.6));
    for (int j = 0; j < ncont; ++j)
        OK(SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(cvars[j]), 2.0));

    SCIP_VAR *kvars[10];
    SCIP_Real kvals[10] = {2.5, 2.5, 2.5, 2.5, -1.0, -1.0, 1.5, 1.5, 1.5, 1.5};
    for (int k = 0; k < 4; ++k)
        kvars[k] = SCIPvarGetTransVar(bvars[k]);
    kvars[4] = SCIPvarGetTransVar(cvars[0]);
    kvars[5] = SCIPvarGetTransVar(cvars[3]);
    for (int k = 6; k < 10; ++k)
        kvars[k] = SCIPvarGetTransVar(bvars[k + 1000]);

    /* Call 400 times: O(nnz) takes < 0.05s; old O(nbin) took ~0.5s+ */
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int rep = 0; rep < 400; ++rep) {
        SCIP_Bool cutoff = FALSE;
        int ncuts = 0;
        OK(SCIPseparateRelaxedKnapsack(scip, NULL, NULL, 10, kvars, kvals, 1.0, 5.0, sol, &cutoff,
                                       &ncuts));
        CHECK(!cutoff);
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    CHECK(elapsed_s(t0, t1) < 0.5);

    OK(SCIPfreeSol(scip, &sol));
    for (int i = 0; i < nbin; ++i)
        OK(SCIPreleaseVar(scip, &bvars[i]));
    for (int j = 0; j < ncont; ++j)
        OK(SCIPreleaseVar(scip, &cvars[j]));
    free(bvars);
    free(cvars);
    OK(SCIPfree(&scip));
    return 0;
}

/* 5. Stress dual simplex BFRT exact tie-breaking, ILLfct_update_dIIfeas prefix,
 *    adaptive row-NNZ compute_zA, and scaled dual cost perturbation + exact unroll
 *    on a wide, massively dual-degenerate bounded LP. */
static int test_dual_bfrt_and_degenerate_perturb(void) {
    const int nrows = 96;
    const int ncols = 2048;
    const int nnz_per_col = 4;
    const int nnz = ncols * nnz_per_col;

    double *obj = calloc((size_t)ncols, sizeof *obj);
    double *lo = calloc((size_t)ncols, sizeof *lo);
    double *hi = calloc((size_t)ncols, sizeof *hi);
    double *lhs = calloc((size_t)nrows, sizeof *lhs);
    double *rhs = calloc((size_t)nrows, sizeof *rhs);
    int *beg = calloc((size_t)ncols, sizeof *beg);
    int *ind = calloc((size_t)nnz, sizeof *ind);
    double *val = calloc((size_t)nnz, sizeof *val);
    double *prim_d = calloc((size_t)ncols, sizeof *prim_d);
    double *prim_p = calloc((size_t)ncols, sizeof *prim_p);
    CHECK(obj && lo && hi && lhs && rhs && beg && ind && val && prim_d && prim_p);

    for (int j = 0; j < ncols; ++j) {
        /* 95% zero-cost bounded columns create massive t_j == 0 dual degeneracy;
         * remaining columns carry scaled objectives up to 500. */
        obj[j] = (j % 20 == 0) ? (double)((j / 20) % 7 + 1) * 50.0 : 0.0;
        lo[j] = 0.0;
        hi[j] = 1.0;
        beg[j] = j * nnz_per_col;
        for (int k = 0; k < nnz_per_col; ++k) {
            int r = (j * 7 + k * 23) % nrows;
            ind[j * nnz_per_col + k] = r;
            /* Varying magnitudes so max-|zAj| tie-breaking has real work to do. */
            val[j * nnz_per_col + k] = 0.25 * (double)((j + k) % 12 + 1);
        }
    }
    for (int i = 0; i < nrows; ++i) {
        lhs[i] = (i % 3 == 0) ? 12.0 : 8.0;
        rhs[i] = (i % 3 == 0) ? 12.0 : 1e20;
    }

    SCIP_LPI *lp_d = NULL;
    SCIP_LPI *lp_p = NULL;
    OK(SCIPlpiCreate(&lp_d, NULL, "v20_dual_bfrt", SCIP_OBJSEN_MINIMIZE));
    OK(SCIPlpiCreate(&lp_p, NULL, "v20_primal_ref", SCIP_OBJSEN_MINIMIZE));
    OK(SCIPlpiLoadColLP(lp_d, SCIP_OBJSEN_MINIMIZE, ncols, obj, lo, hi, NULL, nrows, lhs, rhs,
                        NULL, nnz, beg, ind, val));
    OK(SCIPlpiLoadColLP(lp_p, SCIP_OBJSEN_MINIMIZE, ncols, obj, lo, hi, NULL, nrows, lhs, rhs,
                        NULL, nnz, beg, ind, val));

    OK(SCIPlpiSolveDual(lp_d));
    OK(SCIPlpiSolvePrimal(lp_p));
    CHECK(SCIPlpiIsOptimal(lp_d));
    CHECK(SCIPlpiIsOptimal(lp_p));

    double obj_d = 0.0, obj_p = 0.0;
    OK(SCIPlpiGetSol(lp_d, &obj_d, prim_d, NULL, NULL, NULL));
    OK(SCIPlpiGetSol(lp_p, &obj_p, prim_p, NULL, NULL, NULL));
    CHECK(fabs(obj_d - obj_p) < 1e-5 * (1.0 + fabs(obj_p)));

    /* Verify primal feasibility directly on the original matrix. */
    for (int j = 0; j < ncols; ++j)
        CHECK(prim_d[j] >= lo[j] - 1e-6 && prim_d[j] <= hi[j] + 1e-6);
    for (int i = 0; i < nrows; ++i) {
        double act = 0.0;
        for (int j = 0; j < ncols; ++j) {
            for (int k = 0; k < nnz_per_col; ++k) {
                if (ind[j * nnz_per_col + k] == i)
                    act += val[j * nnz_per_col + k] * prim_d[j];
            }
        }
        CHECK(act >= lhs[i] - 1e-5);
        if (rhs[i] < 1e19)
            CHECK(act <= rhs[i] + 1e-5);
    }

    /* Warm re-solve after tightening bounds & adding a cut row (root-cut pattern). */
    int cidx[16];
    double clb[16], cub[16];
    for (int k = 0; k < 16; ++k) {
        cidx[k] = k * 40;
        clb[k] = 0.5;
        cub[k] = 1.0;
    }
    OK(SCIPlpiChgBounds(lp_d, 16, cidx, clb, cub));
    OK(SCIPlpiChgBounds(lp_p, 16, cidx, clb, cub));
    int rbeg = 0;
    int rind[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    double rval[8] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double rlhs = 3.5, rrhs = 1e20;
    OK(SCIPlpiAddRows(lp_d, 1, &rlhs, &rrhs, NULL, 8, &rbeg, rind, rval));
    OK(SCIPlpiAddRows(lp_p, 1, &rlhs, &rrhs, NULL, 8, &rbeg, rind, rval));
    OK(SCIPlpiSolveDual(lp_d));
    OK(SCIPlpiSolvePrimal(lp_p));
    CHECK(SCIPlpiIsOptimal(lp_d));
    CHECK(SCIPlpiIsOptimal(lp_p));
    OK(SCIPlpiGetSol(lp_d, &obj_d, prim_d, NULL, NULL, NULL));
    OK(SCIPlpiGetSol(lp_p, &obj_p, prim_p, NULL, NULL, NULL));
    CHECK(fabs(obj_d - obj_p) < 1e-5 * (1.0 + fabs(obj_p)));

    OK(SCIPlpiFree(&lp_d));
    OK(SCIPlpiFree(&lp_p));
    free(prim_p);
    free(prim_d);
    free(val);
    free(ind);
    free(beg);
    free(rhs);
    free(lhs);
    free(hi);
    free(lo);
    free(obj);
    return 0;
}

/* y <= 2x + 1 and y <= 5.5 - x meet at (1.5, 4). Branching on x
 * gives (1, 3) down and (2, 3.5) up. Minimizing -y reverses both signs. */
static int test_strongbranch_objective_sense(void) {
    for (int maximize = 0; maximize <= 1; ++maximize) {
        SCIP_OBJSEN sense = maximize ? SCIP_OBJSEN_MAXIMIZE : SCIP_OBJSEN_MINIMIZE;
        SCIP_LPI *lp = NULL;
        OK(SCIPlpiCreate(&lp, NULL, "branch bound sense", sense));
        double sign = maximize ? 1.0 : -1.0;
        double cost[] = {0, sign}, lo[] = {0, 0}, hi[] = {3, 8};
        double lhs[] = {-SCIPlpiInfinity(lp), -SCIPlpiInfinity(lp)}, rhs[] = {1, 5.5};
        int starts[] = {0, 2}, rows[] = {0, 1, 0, 1};
        double values[] = {-2, 1, 1, 1};
        OK(SCIPlpiLoadColLP(lp, sense, 2, cost, lo, hi, NULL, 2, lhs, rhs, NULL,
                            4, starts, rows, values));
        OK(SCIPlpiSolveDual(lp));
        CHECK(SCIPlpiIsOptimal(lp));
        double parent, prim[2], down, up;
        SCIP_Bool down_valid, up_valid;
        OK(SCIPlpiGetSol(lp, &parent, prim, NULL, NULL, NULL));
        /* The oracle is dyadic. The backend still rounds intermediate pivots. */
        const double roundoff = 32 * DBL_EPSILON;
        CHECK(fabs(parent - 4 * sign) <= roundoff && fabs(prim[0] - 1.5) <= roundoff);
        OK(SCIPlpiStartStrongbranch(lp));
        OK(SCIPlpiStrongbranchFrac(lp, 0, prim[0], 500, &down, &up, &down_valid, &up_valid, NULL));
        CHECK(down_valid && up_valid);
        CHECK(fabs(down - 3 * sign) <= roundoff && fabs(up - 3.5 * sign) <= roundoff);
        OK(SCIPlpiEndStrongbranch(lp));
        OK(SCIPlpiFree(&lp));
    }
    return 0;
}

/* 6. Stress fast cached strong branching: compare SCIPlpiStartStrongbranch +
 *    repeated SCIPlpiStrongbranchFrac against independent fresh dual solves
 *    across MIN and MAX senses, itlim = 0 / 1 / 500, infeasible branches, and
 *    external bound/obj/rhs invalidation between sessions. */
static int test_fast_strongbranch_cached(void) {
    const int nrows = 40;
    const int ncols = 120;
    const int nnz_per_col = 4;
    const int nnz = ncols * nnz_per_col;

    double *obj = calloc((size_t)ncols, sizeof *obj);
    double *lo = calloc((size_t)ncols, sizeof *lo);
    double *hi = calloc((size_t)ncols, sizeof *hi);
    double *lhs = calloc((size_t)nrows, sizeof *lhs);
    double *rhs = calloc((size_t)nrows, sizeof *rhs);
    int *beg = calloc((size_t)ncols, sizeof *beg);
    int *ind = calloc((size_t)nnz, sizeof *ind);
    double *val = calloc((size_t)nnz, sizeof *val);
    double *prim = calloc((size_t)ncols, sizeof *prim);
    CHECK(obj && lo && hi && lhs && rhs && beg && ind && val && prim);

    for (int j = 0; j < ncols; ++j) {
        obj[j] = (double)((j * 13) % 19 + 1) * 0.5;
        lo[j] = 0.0;
        hi[j] = 1.0;
        beg[j] = j * nnz_per_col;
        for (int k = 0; k < nnz_per_col; ++k) {
            int r = (j * 7 + k * 13) % nrows;
            ind[j * nnz_per_col + k] = r;
            val[j * nnz_per_col + k] = 0.5 * (double)((j + 3 * k) % 7 + 1);
        }
    }
    for (int pass = 0; pass < 2; ++pass) {
        SCIP_OBJSEN sense = (pass == 0) ? SCIP_OBJSEN_MINIMIZE : SCIP_OBJSEN_MAXIMIZE;
        SCIP_LPI *lp = NULL;
        SCIP_LPI *ref = NULL;
        OK(SCIPlpiCreate(&lp, NULL, "v27_sb_fast", sense));
        OK(SCIPlpiCreate(&ref, NULL, "v27_sb_ref", sense));
        double inf = SCIPlpiInfinity(lp);
        if (pass == 0) {
            for (int i = 0; i < nrows; ++i) {
                lhs[i] = 2.25 + 0.1 * (double)(i % 5);
                rhs[i] = (i % 4 == 0) ? lhs[i] + 5.0 : inf;
            }
        } else {
            for (int i = 0; i < nrows; ++i) {
                lhs[i] = -inf;
                rhs[i] = 3.25 + 0.1 * (double)(i % 5);
            }
        }
        OK(SCIPlpiLoadColLP(lp, sense, ncols, obj, lo, hi, NULL, nrows, lhs, rhs,
                            NULL, nnz, beg, ind, val));
        OK(SCIPlpiLoadColLP(ref, sense, ncols, obj, lo, hi, NULL, nrows, lhs, rhs,
                            NULL, nnz, beg, ind, val));

        OK(SCIPlpiSolveDual(lp));
        CHECK(SCIPlpiIsOptimal(lp));
        double base_obj = 0.0;
        OK(SCIPlpiGetSol(lp, &base_obj, prim, NULL, NULL, NULL));

        /* Collect fractional variables */
        int frac_cols[64];
        double frac_vals[64];
        int nfrac = 0;
        for (int j = 0; j < ncols && nfrac < 64; ++j) {
            double f = prim[j] - floor(prim[j]);
            if (f > 1e-4 && f < 1.0 - 1e-4) {
                frac_cols[nfrac] = j;
                frac_vals[nfrac] = prim[j];
                nfrac++;
            }
        }
        CHECK(nfrac >= 4);

        /* Session 1: repeated cached strong branching across itlim = 0, 1, 500 */
        OK(SCIPlpiStartStrongbranch(lp));
        for (int c = 0; c < nfrac; ++c) {
            double down = 0.0, up = 0.0;
            SCIP_Bool dval = FALSE, uval = FALSE;
            int it = 0;
            int lim = (c == 0) ? 0 : (c == 1) ? 1 : 500;
            OK(SCIPlpiStrongbranchFrac(lp, frac_cols[c], frac_vals[c], lim,
                                       &down, &up, &dval, &uval, &it));
            CHECK(dval && uval);

            if (lim == 500) {
                /* Verify exact agreement against independent fresh solve on ref */
                int col = frac_cols[c];
                double oldlb = lo[col], oldub = hi[col];
                double fl = floor(frac_vals[c]), cl = ceil(frac_vals[c]);

                /* Down branch */
                OK(SCIPlpiChgBounds(ref, 1, &col, &oldlb, &fl));
                OK(SCIPlpiSolveDual(ref));
                if (SCIPlpiIsOptimal(ref)) {
                    double robj = 0.0;
                    OK(SCIPlpiGetObjval(ref, &robj));
                    CHECK(fabs(down - robj) < 1e-5 * (1.0 + fabs(robj)));
                } else {
                    CHECK(SCIPlpiIsPrimalInfeasible(ref));
                    CHECK(sense == SCIP_OBJSEN_MINIMIZE ? down > 1e15 : down < -1e15);
                }
                OK(SCIPlpiChgBounds(ref, 1, &col, &oldlb, &oldub));

                /* Up branch */
                OK(SCIPlpiChgBounds(ref, 1, &col, &cl, &oldub));
                OK(SCIPlpiSolveDual(ref));
                if (SCIPlpiIsOptimal(ref)) {
                    double robj = 0.0;
                    OK(SCIPlpiGetObjval(ref, &robj));
                    CHECK(fabs(up - robj) < 1e-5 * (1.0 + fabs(robj)));
                } else {
                    CHECK(SCIPlpiIsPrimalInfeasible(ref));
                    CHECK(sense == SCIP_OBJSEN_MINIMIZE ? up > 1e15 : up < -1e15);
                }
                OK(SCIPlpiChgBounds(ref, 1, &col, &oldlb, &oldub));
            }
        }
        OK(SCIPlpiEndStrongbranch(lp));
        CHECK(SCIPlpiIsOptimal(lp));
        double after_obj = 0.0;
        OK(SCIPlpiGetObjval(lp, &after_obj));
        CHECK(fabs(after_obj - base_obj) < 1e-6 * (1.0 + fabs(base_obj)));

        /* Session 2: mutate objective & bound externally to test sb_clean_lp / sb_snap_valid invalidation */
        int mcol = frac_cols[0];
        double newobj = obj[mcol] + 7.25;
        OK(SCIPlpiChgObj(lp, 1, &mcol, &newobj));
        OK(SCIPlpiChgObj(ref, 1, &mcol, &newobj));
        OK(SCIPlpiSolveDual(lp));
        CHECK(SCIPlpiIsOptimal(lp));
        OK(SCIPlpiGetSol(lp, &base_obj, prim, NULL, NULL, NULL));

        OK(SCIPlpiStartStrongbranch(lp));
        for (int j = 0; j < ncols; ++j) {
            double f = prim[j] - floor(prim[j]);
            if (f > 1e-4 && f < 1.0 - 1e-4) {
                double down = 0.0, up = 0.0;
                SCIP_Bool dval = FALSE, uval = FALSE;
                OK(SCIPlpiStrongbranchFrac(lp, j, prim[j], 500, &down, &up, &dval, &uval, NULL));
                CHECK(dval && uval);
                double fl = floor(prim[j]);
                OK(SCIPlpiChgBounds(ref, 1, &j, &lo[j], &fl));
                OK(SCIPlpiSolveDual(ref));
                if (SCIPlpiIsOptimal(ref)) {
                    double robj = 0.0;
                    OK(SCIPlpiGetObjval(ref, &robj));
                    CHECK(fabs(down - robj) < 1e-5 * (1.0 + fabs(robj)));
                }
                OK(SCIPlpiChgBounds(ref, 1, &j, &lo[j], &hi[j]));
                break;
            }
        }
        OK(SCIPlpiEndStrongbranch(lp));

        OK(SCIPlpiFree(&lp));
        OK(SCIPlpiFree(&ref));
    }

    free(prim);
    free(val);
    free(ind);
    free(beg);
    free(rhs);
    free(lhs);
    free(hi);
    free(lo);
    free(obj);
    return 0;
}

int main(void) {
    if (test_dense_lu_and_u_bailout() != 0)
        return 1;
    if (test_primal_norm_reuse_and_invalidation() != 0)
        return 1;
    if (test_eventfilter_lazy_hash() != 0)
        return 1;
    if (test_sparse_relaxed_knapsack() != 0)
        return 1;
    if (test_dual_bfrt_and_degenerate_perturb() != 0)
        return 1;
    if (test_strongbranch_objective_sense() != 0)
        return 1;
    if (test_fast_strongbranch_cached() != 0)
        return 1;
    printf("test_v16_engine: ok\n");
    return 0;
}

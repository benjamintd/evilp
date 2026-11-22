#include "lpi/lpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "failed at line %d: %s\n", __LINE__, #x);                              \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)
#define OK(x) CHECK((x) == SCIP_OKAY)

static int check_sparse_vs_dense(SCIP_LPI *lp, int nrows, int ncols) {
    double *drow = calloc((size_t)nrows, sizeof *drow);
    double *srow = calloc((size_t)nrows, sizeof *srow);
    int *rinds = calloc((size_t)nrows, sizeof *rinds);
    unsigned char *rmark = calloc((size_t)nrows, 1);

    double *darow = calloc((size_t)ncols, sizeof *darow);
    double *sarow = calloc((size_t)ncols, sizeof *sarow);
    int *cinds = calloc((size_t)ncols, sizeof *cinds);
    unsigned char *cmark = calloc((size_t)ncols, 1);
    CHECK(drow && srow && rinds && rmark && darow && sarow && cinds && cmark);

    for (int r = 0; r < nrows; ++r) {
        int rninds = -1;
        /* Poison buffers before each call to confirm complete zeroing. */
        for (int i = 0; i < nrows; ++i) {
            drow[i] = 123.456;
            srow[i] = -987.654;
            rmark[i] = 0;
        }
        OK(SCIPlpiGetBInvRow(lp, r, drow, NULL, NULL));
        OK(SCIPlpiGetBInvRow(lp, r, srow, rinds, &rninds));
        CHECK(rninds >= 0 && rninds <= nrows);
        for (int k = 0; k < rninds; ++k) {
            int idx = rinds[k];
            CHECK(idx >= 0 && idx < nrows);
            CHECK(!rmark[idx]);
            rmark[idx] = 1;
        }
        for (int i = 0; i < nrows; ++i) {
            CHECK(fabs(srow[i] - drow[i]) <= 1e-12 * (1.0 + fabs(drow[i])));
            if (!rmark[i])
                CHECK(srow[i] == 0.0);
        }

        int cninds = -1;
        for (int j = 0; j < ncols; ++j) {
            darow[j] = 111.222;
            sarow[j] = -333.444;
            cmark[j] = 0;
        }
        OK(SCIPlpiGetBInvARow(lp, r, drow, darow, NULL, NULL));
        OK(SCIPlpiGetBInvARow(lp, r, srow, sarow, cinds, &cninds));
        CHECK(cninds >= 0 && cninds <= ncols);
        for (int k = 0; k < cninds; ++k) {
            int idx = cinds[k];
            CHECK(idx >= 0 && idx < ncols);
            CHECK(!cmark[idx]);
            cmark[idx] = 1;
        }
        for (int j = 0; j < ncols; ++j) {
            CHECK(fabs(sarow[j] - darow[j]) <= 1e-11 * (1.0 + fabs(darow[j])));
            if (!cmark[j])
                CHECK(sarow[j] == 0.0);
        }
    }

    free(drow);
    free(srow);
    free(rinds);
    free(rmark);
    free(darow);
    free(sarow);
    free(cinds);
    free(cmark);
    return 0;
}

static int build_knapsack_lp(SCIP_LPI **out, int ncols, int nrows) {
    SCIP_LPI *lp = NULL;
    OK(SCIPlpiCreate(&lp, NULL, "v15_test", SCIP_OBJSEN_MINIMIZE));
    double *obj = calloc((size_t)ncols, sizeof *obj);
    double *lo = calloc((size_t)ncols, sizeof *lo);
    double *hi = calloc((size_t)ncols, sizeof *hi);
    double *lhs = calloc((size_t)nrows, sizeof *lhs);
    double *rhs = calloc((size_t)nrows, sizeof *rhs);
    int *beg = calloc((size_t)ncols, sizeof *beg);
    int nnz_per_col = 4;
    int *ind = calloc((size_t)(ncols * nnz_per_col), sizeof *ind);
    double *val = calloc((size_t)(ncols * nnz_per_col), sizeof *val);
    CHECK(obj && lo && hi && lhs && rhs && beg && ind && val);

    for (int i = 0; i < nrows; ++i) {
        lhs[i] = 2.5 + 0.3 * (double)(i % 5);
        rhs[i] = 1e20;
    }
    for (int j = 0; j < ncols; ++j) {
        obj[j] = 1.0 + 0.17 * (double)((j * 7 + 3) % 19);
        lo[j] = 0.0;
        hi[j] = 1.0;
        beg[j] = j * nnz_per_col;
        for (int k = 0; k < nnz_per_col; ++k) {
            ind[beg[j] + k] = (j * 3 + k * 7) % nrows;
            val[beg[j] + k] = 1.0 + 0.25 * (double)((j + k) % 4);
        }
    }
    OK(SCIPlpiLoadColLP(lp, SCIP_OBJSEN_MINIMIZE, ncols, obj, lo, hi, NULL, nrows, lhs, rhs, NULL,
                        ncols * nnz_per_col, beg, ind, val));
    free(obj);
    free(lo);
    free(hi);
    free(lhs);
    free(rhs);
    free(beg);
    free(ind);
    free(val);
    *out = lp;
    return 0;
}

int main(void) {
    const int ncols = 60;
    const int nrows = 24;
    SCIP_LPI *lp = NULL;
    SCIP_LPI *ref = NULL;
    CHECK(build_knapsack_lp(&lp, ncols, nrows) == 0);
    CHECK(build_knapsack_lp(&ref, ncols, nrows) == 0);

    OK(SCIPlpiSolveDual(lp));
    CHECK(SCIPlpiIsOptimal(lp));
    OK(SCIPlpiSolveDual(ref));
    CHECK(SCIPlpiIsOptimal(ref));

    double obj0 = 0.0, ref0 = 0.0;
    double x0[60], xref[60];
    OK(SCIPlpiGetSol(lp, &obj0, x0, NULL, NULL, NULL));
    OK(SCIPlpiGetSol(ref, &ref0, xref, NULL, NULL, NULL));
    CHECK(fabs(obj0 - ref0) < 1e-9);

    /* 1. Verify sparse BInvRow & BInvARow match dense bit-for-bit. */
    CHECK(check_sparse_vs_dense(lp, nrows, ncols) == 0);

    /* 2. Collect fractional basic candidates and compare cached PF-rollback strong
     *    branching against independent fresh bound-changed LP solves. */
    int cands[60];
    double psols[60];
    int ncands = 0;
    for (int j = 0; j < ncols; ++j) {
        if (x0[j] > 1e-4 && x0[j] < 1.0 - 1e-4) {
            cands[ncands] = j;
            psols[ncands] = x0[j];
            ncands++;
        }
    }
    CHECK(ncands >= 4);

    double sb_down[60], sb_up[60];
    SCIP_Bool dv[60], uv[60];
    OK(SCIPlpiStartStrongbranch(lp));
    for (int k = 0; k < ncands; ++k) {
        int iters = 0;
        OK(SCIPlpiStrongbranchFrac(lp, cands[k], psols[k], 200, &sb_down[k], &sb_up[k], &dv[k],
                                   &uv[k], &iters));
        CHECK(dv[k] && uv[k]);
    }
    OK(SCIPlpiEndStrongbranch(lp));

    /* Verify parent LP is still intact after EndStrongbranch. */
    double obj_after = 0.0;
    OK(SCIPlpiGetSol(lp, &obj_after, NULL, NULL, NULL, NULL));
    CHECK(fabs(obj_after - obj0) < 1e-9);

    /* Compare every single strongbranch dual bound against exact dual solves on ref. */
    for (int k = 0; k < ncands; ++k) {
        int col = cands[k];
        double zero = 0.0, one = 1.0;

        OK(SCIPlpiChgBounds(ref, 1, &col, &zero, &zero));
        OK(SCIPlpiSolveDual(ref));
        if (SCIPlpiIsOptimal(ref)) {
            double dref = 0.0;
            OK(SCIPlpiGetObjval(ref, &dref));
            CHECK(fabs(sb_down[k] - dref) <= 1e-7 * (1.0 + fabs(dref)));
        } else {
            CHECK(SCIPlpiIsPrimalInfeasible(ref));
            CHECK(sb_down[k] > 1e20);
        }
        OK(SCIPlpiChgBounds(ref, 1, &col, &zero, &one));

        OK(SCIPlpiChgBounds(ref, 1, &col, &one, &one));
        OK(SCIPlpiSolveDual(ref));
        if (SCIPlpiIsOptimal(ref)) {
            double uref = 0.0;
            OK(SCIPlpiGetObjval(ref, &uref));
            CHECK(fabs(sb_up[k] - uref) <= 1e-7 * (1.0 + fabs(uref)));
        } else {
            CHECK(SCIPlpiIsPrimalInfeasible(ref));
            CHECK(sb_up[k] > 1e20);
        }
        OK(SCIPlpiChgBounds(ref, 1, &col, &zero, &one));
    }

    /* 3. Verify reuse_internal_lpinfo coherence across bound/obj/rhs changes AND
     *    structural mutations (add row, chg coef, del row, add col, del col). */
    int cidx[2] = {0, 1};
    double nobj[2] = {3.5, 0.4};
    OK(SCIPlpiChgObj(lp, 2, cidx, nobj));
    OK(SCIPlpiChgObj(ref, 2, cidx, nobj));
    OK(SCIPlpiSolveDual(lp));
    OK(SCIPlpiSolveDual(ref));
    CHECK(SCIPlpiIsOptimal(lp) && SCIPlpiIsOptimal(ref));
    double o1 = 0.0, o2 = 0.0;
    OK(SCIPlpiGetObjval(lp, &o1));
    OK(SCIPlpiGetObjval(ref, &o2));
    CHECK(fabs(o1 - o2) <= 1e-8 * (1.0 + fabs(o2)));

    int rbeg = 0;
    int rind[3] = {0, 5, 10};
    double rval[3] = {1.0, 2.0, 1.5};
    double rlhs = 1.2, rrhs = 1e20;
    OK(SCIPlpiAddRows(lp, 1, &rlhs, &rrhs, NULL, 3, &rbeg, rind, rval));
    OK(SCIPlpiAddRows(ref, 1, &rlhs, &rrhs, NULL, 3, &rbeg, rind, rval));
    OK(SCIPlpiChgCoef(lp, 0, 2, 2.75));
    OK(SCIPlpiChgCoef(ref, 0, 2, 2.75));
    OK(SCIPlpiSolveDual(lp));
    OK(SCIPlpiSolveDual(ref));
    CHECK(SCIPlpiIsOptimal(lp) && SCIPlpiIsOptimal(ref));
    OK(SCIPlpiGetObjval(lp, &o1));
    OK(SCIPlpiGetObjval(ref, &o2));
    CHECK(fabs(o1 - o2) <= 1e-8 * (1.0 + fabs(o2)));
    CHECK(check_sparse_vs_dense(lp, nrows + 1, ncols) == 0);

    OK(SCIPlpiDelRows(lp, nrows, nrows));
    OK(SCIPlpiDelRows(ref, nrows, nrows));
    OK(SCIPlpiSolveDual(lp));
    OK(SCIPlpiSolveDual(ref));
    CHECK(SCIPlpiIsOptimal(lp) && SCIPlpiIsOptimal(ref));
    OK(SCIPlpiGetObjval(lp, &o1));
    OK(SCIPlpiGetObjval(ref, &o2));
    CHECK(fabs(o1 - o2) <= 1e-8 * (1.0 + fabs(o2)));

    /* 4. Verify SCIP_LPPAR_FROMSCRATCH get/set and solve. */
    int fs = -1;
    OK(SCIPlpiGetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, &fs));
    CHECK(fs == FALSE);
    OK(SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, TRUE));
    OK(SCIPlpiGetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, &fs));
    CHECK(fs == TRUE);
    OK(SCIPlpiSolveDual(lp));
    CHECK(SCIPlpiIsOptimal(lp));
    OK(SCIPlpiGetObjval(lp, &o1));
    CHECK(fabs(o1 - o2) <= 1e-8 * (1.0 + fabs(o2)));
    OK(SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, FALSE));

    OK(SCIPlpiFree(&lp));
    OK(SCIPlpiFree(&ref));
    return 0;
}

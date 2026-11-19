#include "scip_feasjump.h"
#include "common.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define FJ_MAX_VARS 300000
#define FJ_MAX_ROWS 300000
#define FJ_MAX_NNZ 3500000
#define FJ_TOL 1e-6
#define FJ_INF 1e19

typedef struct {
    int nvars;
    int nrows; /* includes final objective cutoff row at index nrows - 1 */
    int obj_row;
    size_t nnz;
    SCIP_VAR **vars;
    double *lb;
    double *ub;
    double *obj;
    unsigned char *is_int;
    unsigned char *is_bin;

    size_t *rstart;
    int *rcol;
    double *rval;
    double *rlhs;
    double *rrhs;

    size_t *cstart;
    int *crow;
    double *cval;

    double *x;
    double *best_x;
    double *act;
    double *w;
    int *viol_list;
    int *viol_pos;
    int nviol;

    int *cand_list;
    int *cand_pos;
    int ncand;
    int *var_viol_cnt;
    int *last_move;
    int moves_since_recompute;
    uint64_t rng;
} FjState;

static uint64_t fj_rand(uint64_t *s) {
    uint64_t x = *s;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *s = x;
    return x;
}

static double row_viol(double a, double l, double u) {
    if (l > -FJ_INF && a < l - FJ_TOL)
        return l - a;
    if (u < FJ_INF && a > u + FJ_TOL)
        return a - u;
    return 0.0;
}

static double clamp_val(double v, double lo, double hi) {
    if (v < lo)
        v = lo;
    if (v > hi)
        v = hi;
    return v;
}

static void fj_free(FjState *st) {
    free(st->lb);
    free(st->ub);
    free(st->obj);
    free(st->is_int);
    free(st->is_bin);
    free(st->rstart);
    free(st->rcol);
    free(st->rval);
    free(st->rlhs);
    free(st->rrhs);
    free(st->cstart);
    free(st->crow);
    free(st->cval);
    free(st->x);
    free(st->best_x);
    free(st->act);
    free(st->w);
    free(st->viol_list);
    free(st->viol_pos);
    free(st->cand_list);
    free(st->cand_pos);
    free(st->var_viol_cnt);
    free(st->last_move);
    memset(st, 0, sizeof *st);
}

static void cand_add(FjState *st, int j) {
    if (st->cand_pos[j] < 0) {
        st->cand_pos[j] = st->ncand;
        st->cand_list[st->ncand++] = j;
    }
}

static void cand_remove(FjState *st, int j) {
    int pos = st->cand_pos[j];
    if (pos >= 0) {
        int last = st->cand_list[--st->ncand];
        st->cand_list[pos] = last;
        st->cand_pos[last] = pos;
        st->cand_pos[j] = -1;
    }
}

static void viol_add(FjState *st, int r) {
    if (st->viol_pos[r] >= 0)
        return;
    st->viol_pos[r] = st->nviol;
    st->viol_list[st->nviol++] = r;
    for (size_t k = st->rstart[r]; k < st->rstart[r + 1]; ++k) {
        int j = st->rcol[k];
        if (st->lb[j] < st->ub[j] - 1e-12) {
            if (++st->var_viol_cnt[j] == 1)
                cand_add(st, j);
        }
    }
}

static void viol_remove(FjState *st, int r) {
    int pos = st->viol_pos[r];
    if (pos < 0)
        return;
    int last = st->viol_list[--st->nviol];
    st->viol_list[pos] = last;
    st->viol_pos[last] = pos;
    st->viol_pos[r] = -1;
    for (size_t k = st->rstart[r]; k < st->rstart[r + 1]; ++k) {
        int j = st->rcol[k];
        if (st->var_viol_cnt[j] > 0 && --st->var_viol_cnt[j] == 0)
            cand_remove(st, j);
    }
}

static void fj_recompute_all(FjState *st) {
    st->nviol = 0;
    st->ncand = 0;
    st->moves_since_recompute = 0;
    for (int r = 0; r < st->nrows; ++r)
        st->viol_pos[r] = -1;
    for (int j = 0; j < st->nvars; ++j) {
        st->cand_pos[j] = -1;
        st->var_viol_cnt[j] = 0;
    }
    for (int r = 0; r < st->nrows; ++r) {
        double a = 0.0;
        for (size_t k = st->rstart[r]; k < st->rstart[r + 1]; ++k)
            a += st->rval[k] * st->x[st->rcol[k]];
        st->act[r] = a;
        if (row_viol(a, st->rlhs[r], st->rrhs[r]) > 0.0)
            viol_add(st, r);
    }
}

static void fj_apply_move(FjState *st, int j, double new_val, int step) {
    double delta = new_val - st->x[j];
    if (fabs(delta) <= 1e-15)
        return;
    st->x[j] = new_val;
    st->last_move[j] = step;
    st->moves_since_recompute++;
    for (size_t k = st->cstart[j]; k < st->cstart[j + 1]; ++k) {
        int r = st->crow[k];
        st->act[r] += st->cval[k] * delta;
        double v = row_viol(st->act[r], st->rlhs[r], st->rrhs[r]);
        if (v > 0.0)
            viol_add(st, r);
        else
            viol_remove(st, r);
    }
}

/* Greedy coordinate-descent objective polish on a feasible point (nviol == 0). */
static void fj_polish_feasible(FjState *st, double deadline) {
    for (int pass = 0; pass < 4; ++pass) {
        int improved = 0;
        for (int j = 0; j < st->nvars; ++j) {
            if ((j & 1023) == 0 && lab_now() >= deadline)
                return;
            double cj = st->obj[j];
            if (fabs(cj) <= 1e-15)
                continue;
            double lo = st->lb[j];
            double hi = st->ub[j];
            if (lo >= hi - 1e-12)
                continue;
            double cur = st->x[j];
            double dlo = (lo > -FJ_INF) ? (lo - cur) : -1e20;
            double dhi = (hi < FJ_INF) ? (hi - cur) : 1e20;

            for (size_t k = st->cstart[j]; k < st->cstart[j + 1]; ++k) {
                int r = st->crow[k];
                if (r == st->obj_row)
                    continue;
                double aij = st->cval[k];
                if (fabs(aij) <= 1e-12)
                    continue;
                double a = st->act[r];
                if (aij > 0.0) {
                    if (st->rlhs[r] > -FJ_INF) {
                        double need = st->rlhs[r] - a;
                        if (need > 0.0)
                            need = 0.0;
                        double lim = need / aij;
                        if (lim > dlo)
                            dlo = lim;
                    }
                    if (st->rrhs[r] < FJ_INF) {
                        double slack = st->rrhs[r] - a;
                        if (slack < 0.0)
                            slack = 0.0;
                        double lim = slack / aij;
                        if (lim < dhi)
                            dhi = lim;
                    }
                } else {
                    if (st->rlhs[r] > -FJ_INF) {
                        double need = st->rlhs[r] - a;
                        if (need > 0.0)
                            need = 0.0;
                        double lim = need / aij;
                        if (lim < dhi)
                            dhi = lim;
                    }
                    if (st->rrhs[r] < FJ_INF) {
                        double slack = st->rrhs[r] - a;
                        if (slack < 0.0)
                            slack = 0.0;
                        double lim = slack / aij;
                        if (lim > dlo)
                            dlo = lim;
                    }
                }
                if (dlo >= -1e-12 && cj > 0.0)
                    break;
                if (dhi <= 1e-12 && cj < 0.0)
                    break;
            }

            double target = cur;
            if (cj > 0.0) {
                if (dlo <= -1e18 || dlo >= -1e-12)
                    continue;
                target = cur + dlo;
                if (st->is_int[j]) {
                    target = ceil(target - 1e-9);
                    if (target < cur + dlo - 1e-9)
                        target += 1.0;
                }
                if (target > cur)
                    target = cur;
            } else {
                if (dhi >= 1e18 || dhi <= 1e-12)
                    continue;
                target = cur + dhi;
                if (st->is_int[j]) {
                    target = floor(target + 1e-9);
                    if (target > cur + dhi + 1e-9)
                        target -= 1.0;
                }
                if (target < cur)
                    target = cur;
            }
            target = clamp_val(target, lo, hi);
            double delta = target - cur;
            if (fabs(delta) <= 1e-12 || cj * delta >= -1e-12)
                continue;

            st->x[j] = target;
            st->moves_since_recompute++;
            for (size_t k = st->cstart[j]; k < st->cstart[j + 1]; ++k) {
                int r = st->crow[k];
                st->act[r] += st->cval[k] * delta;
            }
            improved = 1;
        }
        if (!improved)
            break;
    }
}

static double eval_candidate_delta(const FjState *st, int j, double cand, double *out_obj_delta) {
    double delta = cand - st->x[j];
    double dv = 0.0;
    for (size_t k = st->cstart[j]; k < st->cstart[j + 1]; ++k) {
        int r = st->crow[k];
        double old_a = st->act[r];
        double new_a = old_a + st->cval[k] * delta;
        double old_v = row_viol(old_a, st->rlhs[r], st->rrhs[r]);
        double new_v = row_viol(new_a, st->rlhs[r], st->rrhs[r]);
        dv += st->w[r] * (new_v - old_v);
    }
    *out_obj_delta = st->obj[j] * delta;
    return -dv;
}

static void push_cand(double *cands, int *nc, int maxc, double v, double cur, double lo, double hi) {
    if (!isfinite(v))
        return;
    v = clamp_val(v, lo, hi);
    if (fabs(v - cur) <= 1e-12)
        return;
    for (int i = 0; i < *nc; ++i)
        if (fabs(cands[i] - v) <= 1e-12)
            return;
    if (*nc < maxc)
        cands[(*nc)++] = v;
}

static double fj_best_jump(const FjState *st, int j, double *best_val) {
    double cur = st->x[j];
    double lo = st->lb[j];
    double hi = st->ub[j];
    *best_val = cur;
    if (lo >= hi - 1e-12)
        return -1e30;

    if (st->is_bin[j]) {
        double cand = (cur < 0.5) ? 1.0 : 0.0;
        double od = 0.0;
        double gain = eval_candidate_delta(st, j, cand, &od);
        *best_val = cand;
        return gain - 1e-9 * od;
    }

    double cands[24];
    int nc = 0;
    size_t c0 = st->cstart[j];
    size_t c1 = st->cstart[j + 1];
    int max_cands = (c1 - c0 > 256) ? 10 : 20;

    for (size_t k = c0; k < c1 && nc < max_cands; ++k) {
        int r = st->crow[k];
        if (st->viol_pos[r] < 0)
            continue;
        double aij = st->cval[k];
        if (fabs(aij) <= 1e-12)
            continue;
        double a = st->act[r];
        double target = cur;
        if (st->rlhs[r] > -FJ_INF && a < st->rlhs[r] - FJ_TOL)
            target = cur + (st->rlhs[r] - a) / aij;
        else if (st->rrhs[r] < FJ_INF && a > st->rrhs[r] + FJ_TOL)
            target = cur + (st->rrhs[r] - a) / aij;
        else
            continue;
        if (st->is_int[j]) {
            push_cand(cands, &nc, max_cands, floor(target), cur, lo, hi);
            push_cand(cands, &nc, max_cands, ceil(target), cur, lo, hi);
        } else {
            push_cand(cands, &nc, max_cands, target, cur, lo, hi);
        }
    }
    if (lo > -FJ_INF)
        push_cand(cands, &nc, 24, st->is_int[j] ? ceil(lo - 1e-9) : lo, cur, lo, hi);
    if (hi < FJ_INF)
        push_cand(cands, &nc, 24, st->is_int[j] ? floor(hi + 1e-9) : hi, cur, lo, hi);
    if (lo <= 0.0 && hi >= 0.0)
        push_cand(cands, &nc, 24, 0.0, cur, lo, hi);

    if (nc == 0)
        return -1e30;

    double best_score = -1e30;
    for (int i = 0; i < nc; ++i) {
        double od = 0.0;
        double g = eval_candidate_delta(st, j, cands[i], &od);
        double score = g - 1e-9 * od;
        if (score > best_score) {
            best_score = score;
            *best_val = cands[i];
        }
    }
    return best_score;
}

typedef struct {
    int *col;
    double *val;
    size_t cap;
    size_t len;
    size_t *rstart;
    double *rlhs;
    double *rrhs;
    int rcap;
    int nrows;
} RowBuf;

static int rb_ensure_row(RowBuf *rb) {
    if (rb->nrows + 2 <= rb->rcap)
        return 0;
    int ncap = rb->rcap ? rb->rcap * 2 : 256;
    if (ncap > FJ_MAX_ROWS + 4)
        ncap = FJ_MAX_ROWS + 4;
    if (rb->nrows + 2 > ncap)
        return 1;
    size_t *ns = realloc(rb->rstart, (size_t)ncap * sizeof *ns);
    double *nl = realloc(rb->rlhs, (size_t)ncap * sizeof *nl);
    double *nu = realloc(rb->rrhs, (size_t)ncap * sizeof *nu);
    if (!ns || !nl || !nu) {
        if (ns) rb->rstart = ns;
        if (nl) rb->rlhs = nl;
        if (nu) rb->rrhs = nu;
        return 1;
    }
    rb->rstart = ns;
    rb->rlhs = nl;
    rb->rrhs = nu;
    rb->rcap = ncap;
    return 0;
}

static int rb_ensure_nz(RowBuf *rb, size_t add) {
    if (rb->len + add <= rb->cap)
        return 0;
    size_t ncap = rb->cap ? rb->cap * 2 : 1024;
    while (ncap < rb->len + add)
        ncap *= 2;
    if (ncap > FJ_MAX_NNZ)
        ncap = FJ_MAX_NNZ;
    if (rb->len + add > ncap)
        return 1;
    int *nc = realloc(rb->col, ncap * sizeof *nc);
    double *nv = realloc(rb->val, ncap * sizeof *nv);
    if (!nc || !nv) {
        if (nc) rb->col = nc;
        if (nv) rb->val = nv;
        return 1;
    }
    rb->col = nc;
    rb->val = nv;
    rb->cap = ncap;
    return 0;
}

static void rb_free(RowBuf *rb) {
    free(rb->col);
    free(rb->val);
    free(rb->rstart);
    free(rb->rlhs);
    free(rb->rrhs);
    memset(rb, 0, sizeof *rb);
}

static int rb_push_term(RowBuf *rb, SCIP *scip, SCIP_VAR *var, double coef,
                        double *lhs, double *rhs, double inf) {
    if (fabs(coef) <= 1e-15 || !var)
        return 0;
    SCIP_Real scalar = coef;
    SCIP_Real constant = 0.0;
    SCIP_VAR *v = var;
    if (SCIPgetProbvarSum(scip, &v, &scalar, &constant) != SCIP_OKAY)
        return 1;
    if (constant != 0.0) {
        if (*lhs > -inf)
            *lhs -= constant;
        if (*rhs < inf)
            *rhs -= constant;
    }
    if (!v || !SCIPvarIsActive(v) || fabs(scalar) <= 1e-15)
        return 0;
    int idx = SCIPvarGetProbindex(v);
    if (idx < 0)
        return 0;
    if (rb_ensure_nz(rb, 1))
        return 1;
    rb->col[rb->len] = idx;
    rb->val[rb->len] = scalar;
    rb->len++;
    return 0;
}

static int build_rows_from_lp(SCIP *scip, RowBuf *rb, int nvars, double inf) {
    SCIP_ROW **rows = NULL;
    int nrows = 0;
    if (SCIPgetLPRowsData(scip, &rows, &nrows) != SCIP_OKAY || nrows <= 0 || nrows > FJ_MAX_ROWS)
        return 1;
    for (int i = 0; i < nrows; ++i) {
        SCIP_ROW *row = rows[i];
        if (!row)
            continue;
        double c = SCIProwGetConstant(row);
        double lhs = SCIProwGetLhs(row);
        double rhs = SCIProwGetRhs(row);
        if (lhs <= -inf && rhs >= inf)
            continue;
        if (lhs > -inf)
            lhs -= c;
        else
            lhs = -1e20;
        if (rhs < inf)
            rhs -= c;
        else
            rhs = 1e20;
        int len = SCIProwGetNNonz(row);
        if (len <= 0)
            continue;
        if (rb_ensure_row(rb) || rb_ensure_nz(rb, (size_t)len))
            return 1;
        SCIP_COL **cols = SCIProwGetCols(row);
        SCIP_Real *vals = SCIProwGetVals(row);
        rb->rstart[rb->nrows] = rb->len;
        for (int k = 0; k < len; ++k) {
            SCIP_VAR *v = SCIPcolGetVar(cols[k]);
            int idx = v ? SCIPvarGetProbindex(v) : -1;
            if (idx >= 0 && idx < nvars && fabs(vals[k]) > 1e-15) {
                rb->col[rb->len] = idx;
                rb->val[rb->len] = vals[k];
                rb->len++;
            }
        }
        if (rb->len == rb->rstart[rb->nrows])
            continue;
        rb->rlhs[rb->nrows] = lhs;
        rb->rrhs[rb->nrows] = rhs;
        rb->nrows++;
    }
    return rb->nrows == 0;
}

static int build_rows_from_conss(SCIP *scip, RowBuf *rb, double inf) {
    int nconss = SCIPgetNConss(scip);
    SCIP_CONS **conss = SCIPgetConss(scip);
    if (nconss <= 0 || !conss || nconss > FJ_MAX_ROWS)
        return 1;
    for (int i = 0; i < nconss; ++i) {
        SCIP_CONS *cons = conss[i];
        if (!cons || !SCIPconsIsActive(cons) || !SCIPconsIsChecked(cons))
            continue;
        SCIP_CONSHDLR *hdlr = SCIPconsGetHdlr(cons);
        const char *hname = hdlr ? SCIPconshdlrGetName(hdlr) : "";
        if (rb_ensure_row(rb))
            return 1;
        size_t start = rb->len;
        double lhs = -inf;
        double rhs = inf;
        if (!strcmp(hname, "linear")) {
            int n = SCIPgetNVarsLinear(scip, cons);
            SCIP_VAR **vs = SCIPgetVarsLinear(scip, cons);
            SCIP_Real *vals = SCIPgetValsLinear(scip, cons);
            lhs = SCIPgetLhsLinear(scip, cons);
            rhs = SCIPgetRhsLinear(scip, cons);
            for (int k = 0; k < n; ++k)
                if (rb_push_term(rb, scip, vs[k], vals[k], &lhs, &rhs, inf))
                    return 1;
        } else if (!strcmp(hname, "setppc")) {
            int n = SCIPgetNVarsSetppc(scip, cons);
            SCIP_VAR **vs = SCIPgetVarsSetppc(scip, cons);
            SCIP_SETPPCTYPE tp = SCIPgetTypeSetppc(scip, cons);
            lhs = (tp == SCIP_SETPPCTYPE_PACKING) ? -inf : 1.0;
            rhs = (tp == SCIP_SETPPCTYPE_COVERING) ? inf : 1.0;
            for (int k = 0; k < n; ++k)
                if (rb_push_term(rb, scip, vs[k], 1.0, &lhs, &rhs, inf))
                    return 1;
        } else if (!strcmp(hname, "knapsack")) {
            int n = SCIPgetNVarsKnapsack(scip, cons);
            SCIP_VAR **vs = SCIPgetVarsKnapsack(scip, cons);
            SCIP_Longint *wts = SCIPgetWeightsKnapsack(scip, cons);
            lhs = -inf;
            rhs = (double)SCIPgetCapacityKnapsack(scip, cons);
            for (int k = 0; k < n; ++k)
                if (rb_push_term(rb, scip, vs[k], (double)wts[k], &lhs, &rhs, inf))
                    return 1;
        } else if (!strcmp(hname, "logicor")) {
            int n = SCIPgetNVarsLogicor(scip, cons);
            SCIP_VAR **vs = SCIPgetVarsLogicor(scip, cons);
            lhs = 1.0;
            rhs = inf;
            for (int k = 0; k < n; ++k)
                if (rb_push_term(rb, scip, vs[k], 1.0, &lhs, &rhs, inf))
                    return 1;
        } else if (!strcmp(hname, "varbound")) {
            lhs = SCIPgetLhsVarbound(scip, cons);
            rhs = SCIPgetRhsVarbound(scip, cons);
            if (rb_push_term(rb, scip, SCIPgetVarVarbound(scip, cons), 1.0, &lhs, &rhs, inf) ||
                rb_push_term(rb, scip, SCIPgetVbdvarVarbound(scip, cons),
                             SCIPgetVbdcoefVarbound(scip, cons), &lhs, &rhs, inf))
                return 1;
        } else {
            continue;
        }
        if (rb->len == start || (lhs <= -inf && rhs >= inf)) {
            rb->len = start;
            continue;
        }
        rb->rstart[rb->nrows] = start;
        rb->rlhs[rb->nrows] = (lhs > -inf) ? lhs : -1e20;
        rb->rrhs[rb->nrows] = (rhs < inf) ? rhs : 1e20;
        rb->nrows++;
    }
    return rb->nrows == 0;
}

static int fj_build_state(SCIP *scip, FjState *st, int from_lp) {
    memset(st, 0, sizeof *st);
    int nvars = SCIPgetNVars(scip);
    SCIP_VAR **vars = SCIPgetVars(scip);
    if (nvars <= 0 || !vars || nvars > FJ_MAX_VARS)
        return 1;
    double inf = SCIPinfinity(scip);

    RowBuf rb = {0};
    int rc = from_lp ? build_rows_from_lp(scip, &rb, nvars, inf)
                     : build_rows_from_conss(scip, &rb, inf);
    if (rc) {
        rb_free(&rb);
        return 1;
    }

    /* Append objective cutoff row */
    if (rb_ensure_row(&rb)) {
        rb_free(&rb);
        return 1;
    }
    int obj_nnz = 0;
    for (int j = 0; j < nvars; ++j)
        if (fabs(SCIPvarGetObj(vars[j])) > 1e-15)
            obj_nnz++;
    if (rb_ensure_nz(&rb, (size_t)(obj_nnz > 0 ? obj_nnz : 1))) {
        rb_free(&rb);
        return 1;
    }
    int obj_row = rb.nrows;
    rb.rstart[obj_row] = rb.len;
    for (int j = 0; j < nvars; ++j) {
        double c = SCIPvarGetObj(vars[j]);
        if (fabs(c) > 1e-15) {
            rb.col[rb.len] = j;
            rb.val[rb.len] = c;
            rb.len++;
        }
    }
    rb.rlhs[obj_row] = -1e20;
    rb.rrhs[obj_row] = 1e20; /* activated once a feasible point or incumbent exists */
    rb.nrows++;
    rb.rstart[rb.nrows] = rb.len;

    st->nvars = nvars;
    st->nrows = rb.nrows;
    st->obj_row = obj_row;
    st->nnz = rb.len;
    st->vars = vars;
    st->rstart = rb.rstart;
    st->rcol = rb.col;
    st->rval = rb.val;
    st->rlhs = rb.rlhs;
    st->rrhs = rb.rrhs;

    st->lb = calloc((size_t)nvars, sizeof *st->lb);
    st->ub = calloc((size_t)nvars, sizeof *st->ub);
    st->obj = calloc((size_t)nvars, sizeof *st->obj);
    st->is_int = calloc((size_t)nvars, sizeof *st->is_int);
    st->is_bin = calloc((size_t)nvars, sizeof *st->is_bin);
    st->cstart = calloc((size_t)nvars + 1, sizeof *st->cstart);
    st->crow = calloc(st->nnz ? st->nnz : 1, sizeof *st->crow);
    st->cval = calloc(st->nnz ? st->nnz : 1, sizeof *st->cval);
    st->x = calloc((size_t)nvars, sizeof *st->x);
    st->best_x = calloc((size_t)nvars, sizeof *st->best_x);
    st->act = calloc((size_t)st->nrows, sizeof *st->act);
    st->w = calloc((size_t)st->nrows, sizeof *st->w);
    st->viol_list = calloc((size_t)st->nrows, sizeof *st->viol_list);
    st->viol_pos = calloc((size_t)st->nrows, sizeof *st->viol_pos);
    st->cand_list = calloc((size_t)nvars, sizeof *st->cand_list);
    st->cand_pos = calloc((size_t)nvars, sizeof *st->cand_pos);
    st->var_viol_cnt = calloc((size_t)nvars, sizeof *st->var_viol_cnt);
    st->last_move = calloc((size_t)nvars, sizeof *st->last_move);

    if (!st->lb || !st->ub || !st->obj || !st->is_int || !st->is_bin || !st->cstart ||
        !st->crow || !st->cval || !st->x || !st->best_x || !st->act || !st->w || !st->viol_list ||
        !st->viol_pos || !st->cand_list || !st->cand_pos || !st->var_viol_cnt || !st->last_move) {
        fj_free(st);
        return 1;
    }

    for (int j = 0; j < nvars; ++j) {
        double l = SCIPvarGetLbLocal(vars[j]);
        double u = SCIPvarGetUbLocal(vars[j]);
        st->lb[j] = (l > -inf) ? l : -1e20;
        st->ub[j] = (u < inf) ? u : 1e20;
        st->obj[j] = SCIPvarGetObj(vars[j]);
        SCIP_VARTYPE vt = SCIPvarGetType(vars[j]);
        st->is_int[j] = (vt == SCIP_VARTYPE_BINARY || vt == SCIP_VARTYPE_INTEGER) ? 1U : 0U;
        st->is_bin[j] = (st->is_int[j] && fabs(st->lb[j]) <= 1e-12 && fabs(st->ub[j] - 1.0) <= 1e-12) ? 1U : 0U;
    }

    /* Build CSC from CSR */
    for (int r = 0; r < st->nrows; ++r) {
        for (size_t k = st->rstart[r]; k < st->rstart[r + 1]; ++k)
            st->cstart[st->rcol[k] + 1]++;
    }
    for (int j = 0; j < nvars; ++j)
        st->cstart[j + 1] += st->cstart[j];
    size_t *cur = calloc((size_t)nvars, sizeof *cur);
    if (!cur) {
        fj_free(st);
        return 1;
    }
    memcpy(cur, st->cstart, (size_t)nvars * sizeof *cur);
    for (int r = 0; r < st->nrows; ++r) {
        for (size_t k = st->rstart[r]; k < st->rstart[r + 1]; ++k) {
            int j = st->rcol[k];
            size_t pos = cur[j]++;
            st->crow[pos] = r;
            st->cval[pos] = st->rval[k];
        }
    }
    free(cur);

    for (int r = 0; r < st->nrows; ++r)
        st->w[r] = 1.0;
    st->rng = 0x9e3779b97f4a7c15ULL ^ (uint64_t)nvars ^ ((uint64_t)st->nrows << 20);
    return 0;
}

static void fj_init_point(SCIP *scip, FjState *st, int mode) {
    for (int j = 0; j < st->nvars; ++j) {
        double lo = st->lb[j];
        double hi = st->ub[j];
        double v = 0.0;
        if (mode == 0 && SCIPhasCurrentNodeLP(scip) && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL) {
            v = SCIPvarGetLPSol(st->vars[j]);
            if (st->is_int[j])
                v = round(v);
        } else if (mode == 1) {
            if (st->obj[j] > 0.0 && lo > -FJ_INF)
                v = lo;
            else if (st->obj[j] < 0.0 && hi < FJ_INF)
                v = hi;
            else if (lo <= 0.0 && hi >= 0.0)
                v = 0.0;
            else if (lo > -FJ_INF)
                v = lo;
            else if (hi < FJ_INF)
                v = hi;
        } else {
            if (lo <= 0.0 && hi >= 0.0)
                v = 0.0;
            else if (lo > -FJ_INF && hi < FJ_INF)
                v = 0.5 * (lo + hi);
            else if (lo > -FJ_INF)
                v = lo;
            else if (hi < FJ_INF)
                v = hi;
            if (st->is_int[j])
                v = round(v);
        }
        v = clamp_val(v, lo, hi);
        if (st->is_int[j])
            v = clamp_val(round(v), lo, hi);
        st->x[j] = v;
        st->last_move[j] = -1000;
    }
    for (int r = 0; r < st->nrows; ++r)
        st->w[r] = 1.0;
    fj_recompute_all(st);
}

static SCIP_RETCODE fj_submit_vec(SCIP *scip, SCIP_HEUR *heur, const FjState *st,
                                  const double *vec, SCIP_Bool *stored) {
    *stored = FALSE;
    SCIP_SOL *sol = NULL;
    SCIP_CALL(SCIPcreateSol(scip, &sol, heur));
    for (int j = 0; j < st->nvars; ++j) {
        SCIP_CALL(SCIPsetSolVal(scip, sol, st->vars[j], vec[j]));
    }
    SCIP_CALL(SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, TRUE, stored));
    return SCIP_OKAY;
}

static SCIP_DECL_HEUREXEC(heurExecFeasjump) {
    (void)nodeinfeasible;
    *result = SCIP_DIDNOTRUN;

    if (SCIPgetDepth(scip) > 0 || SCIPgetNVars(scip) <= 0)
        return SCIP_OKAY;
    if (SCIPisStopped(scip))
        return SCIP_OKAY;

    int from_lp = (heurtiming == SCIP_HEURTIMING_AFTERLPNODE) ? 1 : 0;
    if (from_lp) {
        if (!SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL)
            return SCIP_OKAY;
        if (SCIPgetNSols(scip) > 0 && SCIPisEQ(scip, SCIPgetPrimalbound(scip), SCIPgetLPObjval(scip)))
            return SCIP_OKAY;
    }

    SCIP_Real tlim = 0.0;
    SCIP_CALL(SCIPgetRealParam(scip, "limits/time", &tlim));
    double rem = tlim - SCIPgetSolvingTime(scip);
    if (rem <= 0.1)
        return SCIP_OKAY;
    double budget = fmin(0.15, 0.03 * rem);
    if (budget < 0.02)
        budget = 0.02;
    double deadline = lab_now() + budget;

    FjState st;
    if (fj_build_state(scip, &st, from_lp))
        return SCIP_OKAY;
    *result = SCIP_DIDNOTFIND;

    int has_obj = (st.rstart[st.obj_row + 1] > st.rstart[st.obj_row]) ? 1 : 0;

    /* If SCIP already has a primal incumbent, immediately try to beat it */
    int has_best = 0;
    int have_local_best = 0;
    double best_obj = 1e20;
    double last_sub_obj = 1e20;
    int n_submitted = 0;
    int best_submitted = 0;
    int rejects = 0;

    if (SCIPgetNSols(scip) > 0) {
        SCIP_SOL *bs = SCIPgetBestSol(scip);
        if (bs) {
            double o = 0.0;
            for (int j = 0; j < st.nvars; ++j)
                o += st.obj[j] * SCIPgetSolVal(scip, bs, st.vars[j]);
            if (isfinite(o)) {
                has_best = 1;
                best_obj = o;
                last_sub_obj = o;
                double gap = fmax(1e-4 * (1.0 + fabs(best_obj)), 1e-4);
                st.rrhs[st.obj_row] = best_obj - gap;
            }
        }
    }

    int modes[2];
    int nmodes = 0;
    if (from_lp) {
        modes[nmodes++] = 0;
        modes[nmodes++] = 1;
    } else {
        modes[nmodes++] = 1;
        modes[nmodes++] = 2;
    }

    double max_w = 1.0;
    int max_steps = 80000;
    for (int m = 0; m < nmodes && lab_now() < deadline && rejects < 3; ++m) {
        double mode_dl = (m + 1 < nmodes) ? (lab_now() + 0.65 * (deadline - lab_now())) : deadline;
        fj_init_point(scip, &st, modes[m]);
        max_w = 1.0;

        for (int step = 0; step < max_steps; ++step) {
            if ((step & 255) == 0 && lab_now() >= mode_dl)
                break;

            if (st.nviol == 0) {
                /* Eliminate any floating-point drift before polishing/accepting */
                if (st.moves_since_recompute > 64) {
                    fj_recompute_all(&st);
                    if (st.nviol > 0)
                        continue;
                }
                if (has_obj)
                    fj_polish_feasible(&st, mode_dl);

                double cur_obj = st.act[st.obj_row];
                int is_better = (!has_best || cur_obj < best_obj - 1e-6);

                if (n_submitted == 0) {
                    SCIP_Bool stored = FALSE;
                    SCIP_CALL(fj_submit_vec(scip, heur, &st, st.x, &stored));
                    if (stored) {
                        *result = SCIP_FOUNDSOL;
                        n_submitted++;
                        last_sub_obj = cur_obj;
                        has_best = 1;
                        have_local_best = 1;
                        best_obj = cur_obj;
                        memcpy(st.best_x, st.x, (size_t)st.nvars * sizeof *st.best_x);
                        best_submitted = 1;
                    } else if (is_better) {
                        /* Rejected despite beating best_obj (e.g. unmodeled nonlinear/SOS constraint) */
                        rejects++;
                        if (rejects >= 3)
                            break;
                    }
                } else if (is_better) {
                    has_best = 1;
                    have_local_best = 1;
                    best_obj = cur_obj;
                    memcpy(st.best_x, st.x, (size_t)st.nvars * sizeof *st.best_x);
                    best_submitted = 0;
                    double sub_thresh = fmax(0.005 * (1.0 + fabs(last_sub_obj)), 1e-4);
                    if (n_submitted < 12 && best_obj < last_sub_obj - sub_thresh) {
                        SCIP_Bool stored = FALSE;
                        SCIP_CALL(fj_submit_vec(scip, heur, &st, st.best_x, &stored));
                        if (stored) {
                            *result = SCIP_FOUNDSOL;
                            n_submitted++;
                            last_sub_obj = best_obj;
                            best_submitted = 1;
                        }
                    }
                }

                if (!has_obj)
                    break; /* pure feasibility model: done! */
                if (lab_now() >= mode_dl)
                    break;

                double gap = fmax(1e-4 * (1.0 + fabs(best_obj)), 1e-4);
                st.rrhs[st.obj_row] = best_obj - gap;
                st.w[st.obj_row] = 1.0;
                if (row_viol(st.act[st.obj_row], st.rlhs[st.obj_row], st.rrhs[st.obj_row]) > 0.0)
                    viol_add(&st, st.obj_row);
                else
                    break;
                continue;
            }

            if (st.ncand == 0)
                break;

            /* Sample up to 40 active candidate variables for best improving jump */
            int best_j = -1;
            double best_score = 1e-9;
            double best_val = 0.0;
            int sample_k = st.ncand <= 40 ? st.ncand : 40;
            for (int s = 0; s < sample_k; ++s) {
                int j = (st.ncand <= 40) ? st.cand_list[s]
                                         : st.cand_list[fj_rand(&st.rng) % (uint64_t)st.ncand];
                if (step - st.last_move[j] < 2)
                    continue;
                double cv = 0.0;
                double sc = fj_best_jump(&st, j, &cv);
                if (sc > best_score) {
                    best_score = sc;
                    best_j = j;
                    best_val = cv;
                }
            }

            if (best_j >= 0) {
                fj_apply_move(&st, best_j, best_val, step);
                continue;
            }

            /* Local minimum: bump weights on violated rows and escape via a random violated row */
            for (int vi = 0; vi < st.nviol; ++vi) {
                int r = st.viol_list[vi];
                st.w[r] += 1.0;
                if (st.w[r] > max_w)
                    max_w = st.w[r];
            }
            if (max_w > 1e5) {
                for (int r = 0; r < st.nrows; ++r)
                    st.w[r] = 0.5 * st.w[r] + 0.5;
                max_w = 0.5 * max_w + 0.5;
            }

            int r = st.viol_list[fj_rand(&st.rng) % (uint64_t)st.nviol];
            size_t r0 = st.rstart[r];
            size_t r1 = st.rstart[r + 1];
            size_t rlen = r1 - r0;
            if (rlen == 0)
                break;
            int esc_j = -1;
            double esc_score = -1e30;
            double esc_val = 0.0;
            size_t scan = rlen <= 48 ? rlen : 48;
            for (size_t t = 0; t < scan; ++t) {
                size_t idx = (rlen <= 48) ? (r0 + t) : (r0 + (size_t)(fj_rand(&st.rng) % (uint64_t)rlen));
                int j = st.rcol[idx];
                if (st.lb[j] >= st.ub[j] - 1e-12)
                    continue;
                double cv = 0.0;
                double sc = fj_best_jump(&st, j, &cv);
                if (step - st.last_move[j] < 3)
                    sc -= 1.0;
                if (sc > esc_score) {
                    esc_score = sc;
                    esc_j = j;
                    esc_val = cv;
                }
            }
            if (esc_j >= 0)
                fj_apply_move(&st, esc_j, esc_val, step);
            else
                break;
        }
        /* If pure feasibility succeeded, no need for second restart */
        if (has_best && !has_obj)
            break;
    }

    if (have_local_best && !best_submitted) {
        SCIP_Bool stored = FALSE;
        SCIP_CALL(fj_submit_vec(scip, heur, &st, st.best_x, &stored));
        if (stored)
            *result = SCIP_FOUNDSOL;
    }

    fj_free(&st);
    return SCIP_OKAY;
}

SCIP_RETCODE lab_scip_include_feasjump(SCIP *scip) {
    SCIP_HEUR *heur = NULL;
    SCIP_CALL(SCIPincludeHeurBasic(
        scip, &heur, "feasjump", "pure-C17 Feasibility Jump primal heuristic", 'J',
        50000, 0, 0, 0,
        SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_AFTERLPNODE,
        FALSE, heurExecFeasjump, NULL));
    return SCIP_OKAY;
}

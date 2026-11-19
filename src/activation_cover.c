#include "activation_cover.h"
#include "common.h"
#include "exact.h"
#include <glpk.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

enum {
    MAX_ROWS = 10000,
    MAX_COLS = 400000,
    MAX_CELLS = 5000000,
    MAX_SCALAR = 10000,
    PRICE_SCALE = 100000
};

typedef struct {
    size_t gate, demand;
    int capacity, members;
} SourceRow;

typedef struct {
    const LabModel *m;
    SourceRow *row;
    size_t *owner;
    int *cost, *rhs;
    unsigned char *coeff;
    size_t columns, demands, stride;
    double deadline;
} Cover;

static int integer(mpq_srcptr q, int *out) {
    if (mpz_cmp_ui(mpq_denref(q), 1) || mpz_cmpabs_ui(mpq_numref(q), MAX_SCALAR) > 0)
        return 0;
    *out = (int)mpz_get_si(mpq_numref(q));
    return 1;
}

/* A priced column opens a row of the form sum x <= M*y. All auxiliary
 * columns are binary, costless, and occur positively in covering rows.
 * When M covers every member, setting all members equal to y preserves
 * feasibility in both directions and preserves the objective. */
static int project(Cover *c, mpq_t q) {
    const LabModel *m = c->m;
    for (size_t i = 0; i < m->rows; ++i) {
        if (m->row_lower[i] == -INFINITY && m->row_upper[i] == 0 &&
            !lab_source_exact_row_bound(m, i, 1, q) && !mpq_sgn(q))
            continue;
        int rhs;
        if (m->row_upper[i] != INFINITY || !isfinite(m->row_lower[i]) ||
            lab_source_exact_row_bound(m, i, 0, q) || !integer(q, &rhs) || rhs < 0)
            return 0;
        c->row[i].demand = ++c->demands;
        c->rhs[c->demands] = rhs;
    }
    for (size_t j = 0; j < m->cols; ++j) {
        if ((j % 256 == 0 && lab_now() >= c->deadline) || !m->integer[j] || m->col_lower[j] != 0 ||
            m->col_upper[j] != 1 || lab_source_exact_column_bound(m, j, 0, q) || mpq_sgn(q) ||
            lab_source_exact_column_bound(m, j, 1, q) || mpq_cmp_ui(q, 1, 1) ||
            lab_source_exact_objective(m, j, q))
            return 0;
        if (!mpq_sgn(q))
            continue;
        int cost, coefficient;
        size_t k = m->col_start[j];
        if (!integer(q, &cost) || (m->maximize ? -cost : cost) <= 0 ||
            m->col_start[j + 1] != k + 1 || lab_source_exact_coefficient(m, k, q) ||
            !integer(q, &coefficient) || coefficient >= 0)
            return 0;
        SourceRow *r = &c->row[m->row_index[k]];
        if (r->gate || r->demand)
            return 0;
        r->gate = ++c->columns;
        r->capacity = -coefficient;
        c->owner[j] = r->gate;
        c->cost[r->gate] = m->maximize ? -cost : cost;
    }
    if (!c->columns || !c->demands || c->columns > MAX_CELLS / c->demands)
        return 0;
    c->stride = c->demands + 1;
    c->coeff = lab_alloc((c->columns + 1) * c->stride, sizeof *c->coeff);
    for (size_t j = 0; j < m->cols; ++j) {
        if (j % 256 == 0 && lab_now() >= c->deadline)
            return 0;
        if (c->owner[j])
            continue;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            SourceRow *r = &c->row[m->row_index[k]];
            if (m->value[k] != 1 || lab_source_exact_coefficient(m, k, q) || mpq_cmp_ui(q, 1, 1))
                return 0;
            if (!r->demand) {
                if (c->owner[j] || !r->gate)
                    return 0;
                c->owner[j] = r->gate;
                ++r->members;
            }
        }
        if (!c->owner[j])
            return 0;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = c->row[m->row_index[k]].demand;
            if (i) {
                unsigned char *entry = &c->coeff[c->owner[j] * c->stride + i];
                if (*entry == UCHAR_MAX)
                    return 0;
                ++*entry;
            }
        }
    }
    for (size_t i = 0; i < m->rows; ++i)
        if (!c->row[i].demand && (!c->row[i].gate || c->row[i].members > c->row[i].capacity))
            return 0;
    return 1;
}

/* For any y >= 0, b*y + sum_j min(0, c_j - A_j*y) is a lower bound.
 * All terms use integer arithmetic, including the nonbasic upper bounds.
 * MAX_CELLS*UCHAR_MAX*MAX_SCALAR*PRICE_SCALE fits in signed 64 bits on
 * the supported 8-bit-byte platforms. An integral objective closes the
 * gap whenever this bound is strictly greater than the incumbent minus one. */
static int certify(const Cover *c, glp_prob *lp, int *price, double *point) {
    int64_t bound = 0, objective = 0;
    for (size_t i = 1; i <= c->demands; ++i) {
        double y = glp_get_row_dual(lp, (int)i);
        if (!isfinite(y) || y > MAX_SCALAR)
            return 0;
        price[i] = (int)llround(fmax(0, y) * PRICE_SCALE);
        bound += (int64_t)price[i] * c->rhs[i];
    }
    int *activity = lab_alloc(c->stride, sizeof *activity);
    unsigned char *chosen = lab_alloc(c->columns + 1, sizeof *chosen);
    int ok = 0;
    for (size_t j = 1; j <= c->columns; ++j) {
        double x = glp_get_col_prim(lp, (int)j);
        if (!isfinite(x) || (j % 64 == 0 && lab_now() >= c->deadline))
            goto done;
        chosen[j] = x >= .5;
        int64_t reduced = (int64_t)c->cost[j] * PRICE_SCALE;
        objective += chosen[j] * c->cost[j];
        for (size_t i = 1; i <= c->demands; ++i) {
            int a = c->coeff[j * c->stride + i];
            reduced -= (int64_t)price[i] * a;
            activity[i] += chosen[j] * a;
        }
        if (reduced < 0)
            bound += reduced;
    }
    if (bound <= (objective - 1) * PRICE_SCALE)
        goto done;
    for (size_t i = 1; i <= c->demands; ++i)
        if (activity[i] < c->rhs[i])
            goto done;
    for (size_t j = 0; j < c->m->cols; ++j)
        point[j] = chosen[c->owner[j]];
    ok = lab_now() < c->deadline;
done:
    free(chosen);
    free(activity);
    return ok;
}

int lab_activation_cover(const LabModel *m, double deadline, size_t limit, double *point,
                         size_t *iterations) {
    *iterations = 0;
    if (!limit || !m->integer || !m->rows || !m->cols || m->rows > MAX_ROWS || m->cols > MAX_COLS ||
        m->nnz > 1000000 || lab_now() >= deadline)
        return 0;
    /* A priced column has one negative gate coefficient. Check the cheap
     * conditions first to avoid rational scans on unrelated formulations. */
    for (size_t j = 0; j < m->cols; ++j) {
        size_t k = m->col_start[j];
        if (!m->integer[j] || m->col_lower[j] != 0 || m->col_upper[j] != 1 ||
            (m->objective[j] != 0 && (m->col_start[j + 1] != k + 1 || m->value[k] >= 0)))
            return 0;
    }
    Cover c = {.m = m, .deadline = deadline};
    c.row = lab_alloc(m->rows, sizeof *c.row);
    c.owner = lab_alloc(m->cols, sizeof *c.owner);
    c.cost = lab_alloc(m->cols + 1, sizeof *c.cost);
    c.rhs = lab_alloc(m->rows + 1, sizeof *c.rhs);
    mpq_t q;
    mpq_init(q);
    int ok = 0;
    if (project(&c, q)) {
        glp_prob *lp = glp_create_prob();
        glp_add_rows(lp, (int)c.demands);
        glp_add_cols(lp, (int)c.columns);
        int *index = lab_alloc(c.stride, sizeof *index);
        double *value = lab_alloc(c.stride, sizeof *value);
        for (size_t i = 1; i <= c.demands; ++i)
            glp_set_row_bnds(lp, (int)i, GLP_LO, c.rhs[i], 0);
        for (size_t j = 1; j <= c.columns; ++j) {
            int nz = 0;
            for (size_t i = 1; i <= c.demands; ++i)
                if (c.coeff[j * c.stride + i]) {
                    index[++nz] = (int)i;
                    value[nz] = c.coeff[j * c.stride + i];
                }
            glp_set_mat_col(lp, (int)j, nz, index, value);
            glp_set_col_bnds(lp, (int)j, GLP_DB, 0, 1);
            glp_set_obj_coef(lp, (int)j, c.cost[j]);
        }
        glp_smcp params;
        glp_init_smcp(&params);
        params.msg_lev = GLP_MSG_OFF;
        params.meth = GLP_DUALP;
        params.it_lim = (int)(limit > INT_MAX ? INT_MAX : limit);
        double ms = (deadline - lab_now()) * 900;
        if (ms >= 1) {
            params.tm_lim = (int)fmin(INT_MAX, ms);
            int status = glp_simplex(lp, &params);
            *iterations = (size_t)glp_get_it_cnt(lp);
            if (!status && glp_get_status(lp) == GLP_OPT && lab_now() < deadline)
                ok = certify(&c, lp, index, point);
        }
        free(value);
        free(index);
        glp_delete_prob(lp);
    }
    mpq_clear(q);
    free(c.coeff);
    free(c.rhs);
    free(c.cost);
    free(c.owner);
    free(c.row);
    return ok;
}

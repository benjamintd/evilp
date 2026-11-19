#include "common.h"
#include "exact.h"
#include "packing.h"
#include <glpk.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Uniform bins and bounded integer multiplicities. A pattern describes one
 * bin, independent of its label. Enumerate every allowed nonempty pattern. */
enum { TYPES = 8, BINS = 4096, PATTERNS = 4096, MAX_COUNT = 1000000 };
typedef struct {
    int types, bins, capacity, demand[TYPES], weight[TYPES], upper[TYPES];
    int *column, *gate;
    int (*pattern)[TYPES], count;
    size_t spent, limit;
    double deadline;
} Patterns;

static int exact_small(const mpq_t q, int lo, int hi, int *out) {
    if (mpz_cmp_ui(mpq_denref(q), 1) || mpz_cmp_si(mpq_numref(q), lo) < 0 ||
        mpz_cmp_si(mpq_numref(q), hi) > 0)
        return 0;
    *out = (int)mpz_get_si(mpq_numref(q));
    return 1;
}

static int read_types(Patterns *p, const LabModel *m) {
    int types = 0, bins = 0;
    for (size_t r = 0; r < m->rows; ++r) {
        if (m->row_lower[r] == m->row_upper[r] && m->row_lower[r] > 0)
            ++types;
        else if (m->row_lower[r] == -INFINITY && m->row_upper[r] == 0)
            ++bins;
        else
            return 0;
    }
    if (!types || types > TYPES || !bins || bins > BINS ||
        (size_t)(types + 1) * (size_t)bins != m->cols ||
        (size_t)(2 * types + 1) * (size_t)bins != m->nnz)
        return 0;
    int *id = lab_alloc(m->rows, sizeof *id);
    p->column = lab_alloc((size_t)bins * (size_t)types, sizeof *p->column);
    p->gate = lab_alloc((size_t)bins, sizeof *p->gate);
    mpq_t q, charge;
    mpq_inits(q, charge, NULL);
    int ok = 0;
    for (size_t r = 0; r < m->rows; ++r) {
        if (lab_source_exact_row_bound(m, r, 1, q))
            goto done;
        if (m->row_lower[r] == -INFINITY) {
            if (mpq_sgn(q))
                goto done;
            id[r] = -1 - p->bins++;
        } else {
            int d;
            if (!exact_small(q, 1, MAX_COUNT, &d) || lab_source_exact_row_bound(m, r, 0, q) ||
                mpq_cmp_si(q, d, 1))
                goto done;
            id[r] = p->types;
            p->demand[p->types++] = d;
        }
    }
    for (size_t j = 0; j < m->cols; ++j) {
        if (((j & 127U) == 0 && lab_now() >= p->deadline) || !m->integer[j] ||
            lab_source_exact_column_bound(m, j, 0, q) || mpq_sgn(q) ||
            lab_source_exact_column_bound(m, j, 1, q))
            goto done;
        int hi;
        if (!exact_small(q, 1, MAX_COUNT, &hi) || lab_source_exact_objective(m, j, q))
            goto done;
        size_t first = m->col_start[j], end = m->col_start[j + 1];
        if (mpq_sgn(q)) {
            if (m->maximize)
                mpq_neg(q, q);
            if (mpq_sgn(q) <= 0 || (mpq_sgn(charge) && mpq_cmp(charge, q)) || hi != 1 ||
                end - first != 1)
                goto done;
            mpq_set(charge, q);
            int b = -1 - id[m->row_index[first]], neg;
            if (b < 0 || p->gate[b] || lab_source_exact_coefficient(m, first, q) ||
                !exact_small(q, -MAX_COUNT, -1, &neg) || (p->capacity && p->capacity != -neg))
                goto done;
            p->capacity = -neg;
            p->gate[b] = (int)j + 1;
        } else {
            if (end - first != 2)
                goto done;
            size_t k = first, l = first + 1;
            if (id[m->row_index[k]] < 0) {
                k = first + 1;
                l = first;
            }
            int t = id[m->row_index[k]], b = -1 - id[m->row_index[l]], weight;
            if (t < 0 || b < 0 || p->column[b * types + t] ||
                lab_source_exact_coefficient(m, k, q) || mpq_cmp_ui(q, 1, 1) ||
                lab_source_exact_coefficient(m, l, q) || !exact_small(q, 1, MAX_COUNT, &weight) ||
                (p->weight[t] && (p->weight[t] != weight || p->upper[t] != hi)))
                goto done;
            p->weight[t] = weight;
            p->upper[t] = hi;
            p->column[b * types + t] = (int)j + 1;
        }
    }
    /* Uniqueness and the column count imply every bin has every type and gate. */
    ok = p->capacity > 0;
done:
    mpq_clears(q, charge, NULL);
    free(id);
    return ok;
}

static int enumerate(Patterns *p, int t, int capacity, int *counts) {
    if (t == p->types) {
        if (capacity == p->capacity)
            return 1;
        if (p->count == PATTERNS || p->spent >= p->limit || lab_now() >= p->deadline)
            return 0;
        memcpy(p->pattern[p->count++], counts, (size_t)p->types * sizeof *counts);
        ++p->spent;
        return 1;
    }
    int hi = capacity / p->weight[t];
    if (hi > p->upper[t])
        hi = p->upper[t];
    if (hi > p->demand[t])
        hi = p->demand[t];
    for (int v = 0; v <= hi; ++v) {
        counts[t] = v;
        if (!enumerate(p, t + 1, capacity - v * p->weight[t], counts))
            return 0;
    }
    return 1;
}

static int remaining_ms(const Patterns *p) {
    return (int)fmax(1, fmin(1000000, (p->deadline - lab_now()) * 1000));
}

typedef struct {
    int nodes;
    size_t limit;
} PatternBudget;
static void check_tree(glp_tree *tree, void *data) {
    PatternBudget *budget = data;
    glp_ios_tree_size(tree, NULL, NULL, &budget->nodes);
    if ((size_t)budget->nodes >= budget->limit)
        glp_ios_terminate(tree);
}

static int fit_patterns(Patterns *p, const LabModel *m, double *point) {
    glp_prob *lp = glp_create_prob();
    glp_add_rows(lp, p->types);
    glp_add_cols(lp, p->count);
    for (int t = 0; t < p->types; ++t)
        glp_set_row_bnds(lp, t + 1, GLP_FX, p->demand[t], p->demand[t]);
    for (int k = 0; k < p->count; ++k) {
        int indices[TYPES + 1], n = 0;
        double values[TYPES + 1];
        for (int t = 0; t < p->types; ++t)
            if (p->pattern[k][t]) {
                indices[++n] = t + 1;
                values[n] = p->pattern[k][t];
            }
        glp_set_mat_col(lp, k + 1, n, indices, values);
        glp_set_col_bnds(lp, k + 1, GLP_LO, 0, 0);
        glp_set_col_kind(lp, k + 1, GLP_IV);
        glp_set_obj_coef(lp, k + 1, 1);
    }
    glp_smcp sm;
    glp_init_smcp(&sm);
    sm.msg_lev = GLP_MSG_OFF;
    sm.tm_lim = remaining_ms(p);
    int ok = 0;
    if (glp_simplex(lp, &sm) || glp_get_status(lp) != GLP_OPT || lab_now() >= p->deadline)
        goto done;
    /* Turn floating prices into an integer inequality valid for every
     * pattern. Even inaccurate prices can only weaken this lower bound. */
    int64_t price[TYPES], demand_price = 0, denominator = 1000000;
    for (int t = 0; t < p->types; ++t) {
        double y = glp_get_row_dual(lp, t + 1);
        if (!isfinite(y) || fabs(y) > 1000)
            goto done;
        price[t] = (int64_t)floor(y * 1000000);
        demand_price += price[t] * p->demand[t];
    }
    for (int k = 0; k < p->count; ++k) {
        int64_t value = 0;
        for (int t = 0; t < p->types; ++t)
            value += price[t] * p->pattern[k][t];
        if (value > denominator)
            denominator = value;
    }
    if (demand_price <= 0)
        goto done;
    int64_t lower = (demand_price + denominator - 1) / denominator;
    if (lower > p->bins || p->spent >= p->limit)
        goto done;
    /* Ask for a witness at the proved lower bound. Only verified integer
     * counts with exactly this objective can leave this routine. */
    int r = glp_add_rows(lp, 1);
    glp_set_row_bnds(lp, r, GLP_FX, (double)lower, (double)lower);
    int *ind = lab_alloc((size_t)p->count + 1, sizeof *ind);
    double *val = lab_alloc((size_t)p->count + 1, sizeof *val);
    for (int k = 1; k <= p->count; ++k) {
        ind[k] = k;
        val[k] = 1;
        glp_set_obj_coef(lp, k, 0);
    }
    glp_set_mat_row(lp, r, p->count, ind, val);
    free(ind);
    free(val);
    glp_iocp ip;
    glp_init_iocp(&ip);
    ip.msg_lev = GLP_MSG_OFF;
    ip.presolve = GLP_ON;
    ip.tm_lim = remaining_ms(p);
    PatternBudget budget = {.limit = p->limit - p->spent};
    ip.cb_func = check_tree;
    ip.cb_info = &budget;
    (void)glp_intopt(lp, &ip);
    p->spent += (size_t)budget.nodes < budget.limit ? (size_t)budget.nodes : budget.limit;
    int status = glp_mip_status(lp);
    if ((status != GLP_OPT && status != GLP_FEAS) || lab_now() >= p->deadline)
        goto done;
    int sums[TYPES] = {0}, used = 0;
    memset(point, 0, m->cols * sizeof *point);
    for (int k = 0; k < p->count; ++k) {
        double copies = round(glp_mip_col_val(lp, k + 1));
        if (!isfinite(copies) || copies < 0 || copies > lower - used)
            goto done;
        for (int i = 0; i < (int)copies; ++i) {
            point[p->gate[used] - 1] = 1;
            for (int t = 0; t < p->types; ++t) {
                int a = p->pattern[k][t];
                if (a > p->demand[t] - sums[t])
                    goto done;
                sums[t] += a;
                point[p->column[used * p->types + t] - 1] = a;
            }
            ++used;
        }
    }
    ok = used == lower && !memcmp(sums, p->demand, (size_t)p->types * sizeof *sums);
done:
    glp_delete_prob(lp);
    return ok;
}

int lab_packing_patterns(const LabModel *m, double deadline, size_t limit, double *point,
                         size_t *work) {
    *work = 0;
    if (!limit || m->rows > TYPES + BINS || !m->cols || lab_now() >= deadline)
        return 0;
    Patterns p = {.deadline = deadline, .limit = limit};
    int ok = 0, counts[TYPES] = {0};
    if (read_types(&p, m)) {
        p.pattern = lab_alloc(PATTERNS, sizeof *p.pattern);
        if (enumerate(&p, 0, p.capacity, counts) && p.count && lab_now() < deadline)
            ok = fit_patterns(&p, m, point);
    }
    free(p.pattern);
    free(p.column);
    free(p.gate);
    *work = p.spent;
    return ok;
}

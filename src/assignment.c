#include "assignment.h"
#include "common.h"
#include "exact.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { MAX_SIDE = 256, MAX_POTENTIALS = 256, MAX_TRIALS = 1000 };
typedef struct {
    int from, to;
    double lag;
    size_t row;
    int direction;
} BoundArc;
typedef struct {
    const LabModel *m;
    size_t *offset, *entry, *column;
    int *cover, *left, *right, *potential, *edge;
    int n, nc;
    double deadline;
} Assignment;

static int is_one(const LabModel *m, size_t k, mpq_t q) {
    return m->value[k] == 1 && !lab_source_exact_coefficient(m, k, q) && !mpq_cmp_ui(q, 1, 1);
}

/* Cover equations define vertices and each binary column joins exactly two.
 * Costs of every other column must be source-exact zero. Duplicate choices
 * and unsupported side constraints defer to the general solver. */
static int recognize(Assignment *g, const double *upper) {
    const LabModel *m = g->m;
    size_t covers = 0;
    mpq_t q;
    mpq_init(q);
    for (size_t i = 0; i < m->rows; i++) {
        if (m->row_lower[i] == 1 && m->row_upper[i] == 1 &&
            !lab_source_exact_row_bound(m, i, 0, q) && !mpq_cmp_ui(q, 1, 1) &&
            !lab_source_exact_row_bound(m, i, 1, q) && !mpq_cmp_ui(q, 1, 1))
            g->cover[i] = 1;
    }
    int ok = 0;
    for (size_t j = 0; j < m->cols; j++) {
        if ((j & 127U) == 0 && lab_now() >= g->deadline)
            goto done;
        if (!m->integer[j]) {
            if (++g->nc > MAX_POTENTIALS || lab_source_exact_objective(m, j, q) || mpq_sgn(q))
                goto done;
            g->potential[j] = g->nc;
        } else if (m->col_lower[j] != 0 || m->col_upper[j] != 1 ||
                   lab_source_exact_column_bound(m, j, 0, q) || mpq_sgn(q) ||
                   lab_source_exact_column_bound(m, j, 1, q) || mpq_cmp_ui(q, 1, 1))
            goto done;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            size_t i = m->row_index[k];
            g->offset[i + 1]++;
            if (g->cover[i] && (!m->integer[j] || !is_one(m, k, q)))
                g->cover[i] = 0;
            g->column[k] = j;
        }
    }
    size_t *cursor = lab_alloc(m->rows, sizeof *cursor);
    for (size_t i = 0; i < m->rows; i++) {
        g->offset[i + 1] += g->offset[i];
        cursor[i] = g->offset[i];
        covers += (size_t)g->cover[i];
    }
    for (size_t k = 0; k < m->nnz; k++)
        g->entry[cursor[m->row_index[k]]++] = k;
    free(cursor);
    if (!covers || covers % 2 || covers > 2 * MAX_SIDE)
        goto done;
    for (size_t j = 0; j < m->cols; j++)
        if (m->integer[j]) {
            for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
                size_t i = m->row_index[k];
                if (!g->cover[i])
                    continue;
                if (!g->left[j])
                    g->left[j] = (int)i + 1;
                else if (!g->right[j])
                    g->right[j] = (int)i + 1;
                else
                    goto done;
            }
            if (!g->right[j] || g->right[j] == g->left[j])
                goto done;
        }
    int *color = lab_alloc(m->rows, sizeof *color), *slot = lab_alloc(m->rows, sizeof *slot);
    size_t *queue = lab_alloc(covers, sizeof *queue);
    int sides[2] = {0};
    for (size_t i = 0; i < m->rows; i++)
        if (g->cover[i] && !color[i]) {
            size_t head = 0, tail = 1;
            queue[0] = i;
            color[i] = 1;
            while (head < tail) {
                size_t v = queue[head++];
                slot[v] = sides[color[v] < 0]++;
                for (size_t t = g->offset[v]; t < g->offset[v + 1]; t++) {
                    size_t j = g->column[g->entry[t]];
                    size_t other =
                        (size_t)(g->left[j] == (int)v + 1 ? g->right[j] : g->left[j]) - 1;
                    if (color[other] == color[v])
                        goto colored;
                    if (!color[other]) {
                        color[other] = -color[v];
                        queue[tail++] = other;
                    }
                }
            }
        }
    if (sides[0] != sides[1])
        goto colored;
    g->n = sides[0];
    g->edge = lab_alloc((size_t)g->n * (size_t)g->n, sizeof *g->edge);
    for (size_t j = 0; j < m->cols; j++)
        if (m->integer[j] && upper[j] >= 1) {
            int a = g->left[j] - 1, b = g->right[j] - 1;
            if (color[a] < 0) {
                int swap = a;
                a = b;
                b = swap;
            }
            int e = slot[a] * g->n + slot[b];
            if (g->edge[e])
                goto colored;
            g->edge[e] = (int)j + 1;
        }
    ok = 1;
colored:
    free(color);
    free(slot);
    free(queue);
done:
    mpq_clear(q);
    return ok;
}

static unsigned random_word(unsigned *state) {
    *state = *state * 1664525U + 1013904223U;
    return *state;
}
/* Shortest augmenting paths for a minimum-cost bipartite matching. Prices
 * used here order proposals only. The final cost test uses source rationals. */
static int propose(const Assignment *g, const double *penalty, int trial, unsigned *seed,
                   int *match) {
    int n = g->n;
    double *cost = lab_alloc((size_t)n * (size_t)n, sizeof *cost);
    for (int e = 0; e < n * n; e++)
        if (g->edge[e]) {
            int j = g->edge[e] - 1;
            double jitter = trial ? 1e-6 * (random_word(seed) % 100000) / 100000 : 0;
            cost[e] =
                (g->m->maximize ? -g->m->objective[j] : g->m->objective[j]) + penalty[j] + jitter;
        }
    double u[MAX_SIDE + 1] = {0}, v[MAX_SIDE + 1] = {0}, dist[MAX_SIDE + 1];
    int basic[MAX_SIDE + 1] = {0}, pred[MAX_SIDE + 1];
    int ok = 0;
    for (int i = 1; i <= n; i++) {
        if (lab_now() >= g->deadline)
            goto done;
        basic[0] = i;
        int j0 = 0;
        unsigned char used[MAX_SIDE + 1] = {0};
        for (int j = 1; j <= n; j++)
            dist[j] = INFINITY;
        do {
            used[j0] = 1;
            int r = basic[j0], next = 0;
            double delta = INFINITY;
            for (int j = 1; j <= n; j++)
                if (!used[j]) {
                    int e = (r - 1) * n + j - 1;
                    if (g->edge[e]) {
                        double reduced = cost[e] - u[r] - v[j];
                        if (reduced < dist[j]) {
                            dist[j] = reduced;
                            pred[j] = j0;
                        }
                    }
                    if (dist[j] < delta) {
                        delta = dist[j];
                        next = j;
                    }
                }
            if (!isfinite(delta))
                goto done;
            for (int j = 0; j <= n; j++) {
                if (used[j]) {
                    u[basic[j]] += delta;
                    v[j] -= delta;
                } else
                    dist[j] -= delta;
            }
            j0 = next;
        } while (basic[j0]);
        do {
            int next = pred[j0];
            basic[j0] = basic[next];
            j0 = next;
        } while (j0);
    }
    for (int j = 1; j <= n; j++)
        match[basic[j] - 1] = j - 1;
    ok = 1;
done:
    free(cost);
    return ok;
}

/* Given a perfect matching, eliminate its tight column prices from the
 * assignment dual. The remaining inequalities are difference constraints
 * on row prices. A feasible dual with equality on the matching proves its
 * cost is minimum, including arbitrarily small source-decimal differences. */
static int minimum_cost(const Assignment *g, const int *match) {
    const LabModel *m = g->m;
    int n = g->n, owner[MAX_SIDE];
    for (int i = 0; i < n; i++)
        owner[match[i]] = i;
    mpq_t *cost = lab_alloc((size_t)n * (size_t)n, sizeof *cost);
    mpq_t price[MAX_SIDE], bound;
    mpq_init(bound);
    for (int i = 0; i < n; i++)
        mpq_init(price[i]);
    for (int e = 0; e < n * n; e++)
        mpq_init(cost[e]);
    int ok = 0;
    for (int e = 0; e < n * n; e++)
        if (g->edge[e]) {
            if (((e & 255) == 0 && lab_now() >= g->deadline) ||
                lab_source_exact_objective(m, (size_t)g->edge[e] - 1, cost[e]))
                goto done;
            if (m->maximize)
                mpq_neg(cost[e], cost[e]);
        }
    for (int j = 0; j < n; j++) {
        mpq_set(bound, cost[owner[j] * n + j]);
        for (int i = 0; i < n; i++)
            if (g->edge[i * n + j])
                mpq_sub(cost[i * n + j], cost[i * n + j], bound);
    }
    for (int pass = 0; pass < n; pass++) {
        int changed = 0;
        for (int i = 0; i < n; i++) {
            if (lab_now() >= g->deadline)
                goto done;
            for (int j = 0; j < n; j++)
                if (g->edge[i * n + j]) {
                    mpq_add(bound, price[owner[j]], cost[i * n + j]);
                    if (mpq_cmp(price[i], bound) > 0) {
                        mpq_set(price[i], bound);
                        changed = 1;
                    }
                }
        }
        if (!changed) {
            ok = 1;
            break;
        }
    }
done:
    for (int e = 0; e < n * n; e++)
        mpq_clear(cost[e]);
    free(cost);
    for (int i = 0; i < n; i++)
        mpq_clear(price[i]);
    mpq_clear(bound);
    return ok;
}

/* After fixing the selected edges every remaining row has at most one
 * positive and one negative potential. Node zero anchors absolute bounds. */
static int build_arcs(const Assignment *g, const double *x, BoundArc *arcs, int *count) {
    const LabModel *m = g->m;
    int n = 0;
    for (size_t i = 0; i < m->rows; i++)
        if (!g->cover[i]) {
            int plus = 0, minus = 0;
            double fixed = 0;
            for (size_t t = g->offset[i]; t < g->offset[i + 1]; t++) {
                size_t k = g->entry[t], j = g->column[k];
                if (!g->potential[j])
                    fixed += m->value[k] * x[j];
                else if (m->value[k] == 1 && !plus)
                    plus = g->potential[j];
                else if (m->value[k] == -1 && !minus)
                    minus = g->potential[j];
                else
                    return 0;
            }
            if (!isfinite(fixed))
                return 0;
            if (isfinite(m->row_lower[i]))
                arcs[n++] = (BoundArc){minus, plus, m->row_lower[i] - fixed, i, 1};
            if (isfinite(m->row_upper[i]))
                arcs[n++] = (BoundArc){plus, minus, fixed - m->row_upper[i], i, -1};
        }
    for (size_t j = 0; j < m->cols; j++)
        if (g->potential[j]) {
            if (isfinite(m->col_lower[j]))
                arcs[n++] = (BoundArc){0, g->potential[j], m->col_lower[j], j, 0};
            if (isfinite(m->col_upper[j]))
                arcs[n++] = (BoundArc){g->potential[j], 0, -m->col_upper[j], j, 0};
        }
    *count = n;
    return 1;
}
static int check_cycle(const Assignment *g, const BoundArc *arc, int ne, const double *x,
                       double *penalty) {
    double y[MAX_POTENTIALS + 1] = {0};
    int pred[MAX_POTENTIALS + 1];
    for (int i = 0; i <= g->nc; i++)
        pred[i] = -1;
    int last = -1;
    for (int pass = 0; pass <= g->nc; pass++) {
        if (lab_now() >= g->deadline)
            return -1;
        last = -1;
        for (int e = 0; e < ne; e++) {
            double cand = y[arc[e].from] + arc[e].lag;
            if (cand > y[arc[e].to] + 1e-8) {
                y[arc[e].to] = cand;
                pred[arc[e].to] = e;
                last = arc[e].to;
            }
        }
        if (last < 0)
            return 1;
    }
    for (int i = 0; i <= g->nc; i++) {
        if (pred[last] < 0)
            return -1;
        last = arc[pred[last]].from;
    }
    int begin = last;
    for (size_t j = 0; j < g->m->cols; j++)
        penalty[j] *= .995;
    do {
        BoundArc z = arc[pred[last]];
        if (z.direction)
            for (size_t p = g->offset[z.row]; p < g->offset[z.row + 1]; p++) {
                size_t k = g->entry[p], j = g->column[k];
                double c = z.direction * g->m->value[k];
                if (g->m->integer[j] && ((c < 0 && x[j] == 1) || (c > 0 && x[j] == 0)))
                    penalty[j] -= copysign(.001, c);
            }
        last = z.from;
    } while (last != begin);
    return 0;
}
/* Recover a source-decimal coordinate without passing through binary64.
 * Decimal inputs and unit potential coefficients guarantee a denominator
 * with factors two and five. Keep the SOL reader's scalar size limits. */
static char *decimal_coordinate(const mpq_t value) {
    mpz_t d, num;
    mpz_inits(d, num, NULL);
    mpz_set(d, mpq_denref(value));
    mp_bitcnt_t twos = mpz_scan1(d, 0);
    mpz_fdiv_q_2exp(d, d, twos);
    unsigned long fives = 0;
    while (mpz_divisible_ui_p(d, 5)) {
        mpz_divexact_ui(d, d, 5);
        ++fives;
    }
    unsigned long places = twos > fives ? twos : fives;
    char *out = NULL;
    if (mpz_cmp_ui(d, 1) == 0 && places <= 1000) {
        mpz_ui_pow_ui(d, 5, places - fives);
        mpz_mul(num, mpq_numref(value), d);
        mpz_mul_2exp(num, num, places - twos);
        size_t len = mpz_sizeinbase(num, 10) + 2;
        if (len < 1000) {
            out = lab_alloc(len + 32, 1);
            mpz_get_str(out, 10, num);
            if (places)
                snprintf(out + strlen(out), 32, "e-%lu", places);
        }
    }
    mpz_clears(d, num, NULL);
    return out;
}
void lab_assignment_free_decimals(char **decimal, size_t cols) {
    if (decimal)
        for (size_t j = 0; j < cols; j++)
            free(decimal[j]);
    free(decimal);
}
static int recover(const Assignment *g, const BoundArc *arcs, int ne, double *x, char ***out) {
    const LabModel *m = g->m;
    mpq_t *lag = lab_alloc((size_t)ne, sizeof *lag);
    for (int e = 0; e < ne; e++)
        mpq_init(lag[e]);
    mpq_t y[MAX_POTENTIALS + 1], sum, q;
    for (int i = 0; i <= g->nc; i++)
        mpq_init(y[i]);
    mpq_inits(sum, q, NULL);
    char **decimal = lab_alloc(m->cols, sizeof *decimal);
    int ok = 0;
    for (int e = 0; e < ne; e++) {
        if ((e & 63) == 0 && lab_now() >= g->deadline)
            goto done;
        const BoundArc a = arcs[e];
        if (!a.direction) {
            if (lab_source_exact_column_bound(m, a.row, a.to == 0, lag[e]))
                goto done;
            if (a.to == 0)
                mpq_neg(lag[e], lag[e]);
            continue;
        }
        if (lab_source_exact_row_bound(m, a.row, a.direction < 0, lag[e]))
            goto done;
        for (size_t t = g->offset[a.row]; t < g->offset[a.row + 1]; t++) {
            size_t k = g->entry[t], j = g->column[k];
            if (g->potential[j]) {
                if (lab_source_exact_coefficient(m, k, q) ||
                    mpq_cmp_si(q, m->value[k] > 0 ? 1 : -1, 1))
                    goto done;
            } else if (x[j] != 0) {
                if (lab_source_exact_coefficient(m, k, q))
                    goto done;
                mpq_sub(lag[e], lag[e], q);
            }
        }
        if (a.direction < 0)
            mpq_neg(lag[e], lag[e]);
    }
    for (int pass = 0; pass <= g->nc; pass++) {
        int changed = 0;
        for (int e = 0; e < ne; e++) {
            if ((e & 255) == 0 && lab_now() >= g->deadline)
                goto done;
            mpq_add(sum, y[arcs[e].from], lag[e]);
            if (mpq_cmp(sum, y[arcs[e].to]) > 0) {
                mpq_set(y[arcs[e].to], sum);
                changed = 1;
            }
        }
        if (!changed)
            break;
        if (pass == g->nc)
            goto done;
    }
    for (size_t j = 0; j < m->cols; j++)
        if (g->potential[j]) {
            mpq_sub(q, y[g->potential[j]], y[0]);
            decimal[j] = decimal_coordinate(q);
            if (!decimal[j])
                goto done;
            x[j] = strtod(decimal[j], NULL);
            if (!isfinite(x[j]))
                goto done;
        }
    *out = decimal;
    decimal = NULL;
    ok = 1;
done:
    lab_assignment_free_decimals(decimal, m->cols);
    for (int e = 0; e < ne; e++)
        mpq_clear(lag[e]);
    free(lag);
    for (int i = 0; i <= g->nc; i++)
        mpq_clear(y[i]);
    mpq_clears(sum, q, NULL);
    return ok;
}
static Assignment assignment_graph(const LabModel *m, double deadline) {
    Assignment g = {.m = m, .deadline = deadline};
    g.offset = lab_alloc(m->rows + 1, sizeof *g.offset);
    g.cover = lab_alloc(m->rows, sizeof *g.cover);
    g.left = lab_alloc(m->cols, sizeof *g.left);
    g.right = lab_alloc(m->cols, sizeof *g.right);
    g.potential = lab_alloc(m->cols, sizeof *g.potential);
    g.entry = lab_alloc(m->nnz, sizeof *g.entry);
    g.column = lab_alloc(m->nnz, sizeof *g.column);
    return g;
}

static void free_assignment_graph(Assignment *g) {
    free(g->offset);
    free(g->cover);
    free(g->left);
    free(g->right);
    free(g->potential);
    free(g->entry);
    free(g->column);
    free(g->edge);
}

int lab_assignment_solve(const LabModel *m, const double *lower, const double *upper,
                         double deadline, size_t limit, double *point, char ***decimal,
                         size_t *trials) {
    *trials = 0;
    *decimal = NULL;
    if (!limit || !m->integer || !m->rows || !m->cols || m->rows > 8192 || m->cols > 16384 ||
        m->nnz > 100000 || lab_now() >= deadline)
        return 0;
    /* Each edge must touch two cover equations. This cheap filter avoids
     * building secondary indices for unrelated models. */
    size_t nr = 0, ni = 0;
    for (size_t i = 0; i < m->rows; i++)
        nr += m->row_lower[i] == 1 && m->row_upper[i] == 1;
    if (nr < 2)
        return 0;
    for (size_t j = 0; j < m->cols; j++) {
        if (m->integer[j]) {
            if (m->col_lower[j] != 0 || m->col_upper[j] != 1)
                return 0;
            ++ni;
        } else if (m->objective[j] != 0)
            return 0;
    }
    if (!ni || m->cols - ni > MAX_POTENTIALS || m->nnz < 2 * ni)
        return 0;
    Assignment g = assignment_graph(m, deadline);
    int ok = 0;
    if (recognize(&g, upper)) {
        double *penalty = lab_alloc(m->cols, sizeof *penalty);
        BoundArc *arcs = lab_alloc(2 * (m->rows + m->cols), sizeof *arcs);
        int match[MAX_SIDE];
        unsigned seed = 7;
        while (*trials < limit && *trials < MAX_TRIALS && lab_now() < deadline) {
            int trial = (int)(*trials)++;
            if (!propose(&g, penalty, trial, &seed, match))
                break;
            memset(point, 0, m->cols * sizeof *point);
            for (int i = 0; i < g.n; i++)
                point[g.edge[i * g.n + match[i]] - 1] = 1;
            int ne = 0;
            if (!build_arcs(&g, point, arcs, &ne))
                break;
            int feasible = check_cycle(&g, arcs, ne, point, penalty);
            if (feasible < 0)
                break;
            if (feasible) {
                /* Fixed-one presolve bounds can only restrict feasibility.
                 * The matching dual on the larger edge set is still valid. */
                int inside = 1;
                for (size_t j = 0; j < m->cols; j++)
                    if (m->integer[j] && (point[j] < lower[j] || point[j] > upper[j]))
                        inside = 0;
                if (inside && minimum_cost(&g, match) && recover(&g, arcs, ne, point, decimal))
                    ok = 1;
                break;
            }
        }
        free(penalty);
        free(arcs);
    }
    free_assignment_graph(&g);
    return ok;
}

#ifdef LAB_SCIP_INTEGER
#include "assignment_quotient.inc"
#endif

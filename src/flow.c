#include "flow.h"
#include "common.h"
#include "exact.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* A curve is the minimum opening cost as a function of signed integer
 * transfer between two terminals. For any fixed gate assignment the bounded
 * incidence system has an integral flow if it has a real one (augment with
 * a super source/sink and apply integral max-flow). */
enum { MAX_VERTICES = 512, MAX_ARCS = 4096, MAX_PAIRS = 150000, MAX_SEGMENTS = 600000 };
typedef struct {
    int64_t lo, hi, cost;
} Band;
typedef struct Curve Curve;
struct Curve {
    int a, b, kind, flow, gate, depth;
    int64_t demand;
    Curve *left, *right;
    Band *band;
    size_t n;
};
typedef struct {
    Curve *pool, **edge;
    size_t used, segments;
    size_t nv;
    double deadline;
} Network;

static int exact_int(mpq_srcptr q, int64_t *x) {
    if (mpz_cmp_ui(mpq_denref(q), 1) || !mpz_fits_slong_p(mpq_numref(q)))
        return 0;
    long v = mpz_get_si(mpq_numref(q));
    if (v < -1000000000 || v > 1000000000)
        return 0;
    *x = v;
    return 1;
}
static void connect(Network *g, int a, int b, Curve *c) {
    g->edge[(size_t)a * g->nv + (size_t)b] = g->edge[(size_t)b * g->nv + (size_t)a] = c;
}
static Band oriented(const Curve *c, int from, size_t k) {
    Band v = c->band[k];
    if (c->a != from) {
        int64_t lo = v.lo;
        v.lo = -v.hi;
        v.hi = -lo;
    }
    return v;
}
static int cost_order(const void *a, const void *b) {
    int64_t x = ((const Band *)a)->cost, y = ((const Band *)b)->cost;
    return (x > y) - (x < y);
}
static int coord_order(const void *a, const void *b) {
    int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    return (x > y) - (x < y);
}
static size_t first(const int64_t *a, size_t n, int64_t x) {
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t k = lo + (hi - lo) / 2;
        if (a[k] < x)
            lo = k + 1;
        else
            hi = k;
    }
    return lo;
}
static size_t unfilled(size_t *skip, size_t i) {
    while (skip[i] != i) {
        skip[i] = skip[skip[i]];
        i = skip[i];
    }
    return i;
}

/* Resolve overlapping intervals in increasing cost order. Endpoints alone
 * determine the result, independently of the magnitudes of the capacities. */
static int envelope(Network *g, Curve *c, Band *all, size_t n) {
    if (!n)
        return 0;
    int64_t *coord = lab_alloc(2 * n, sizeof *coord);
    for (size_t i = 0; i < n; ++i) {
        coord[2 * i] = all[i].lo;
        coord[2 * i + 1] = all[i].hi + 1;
    }
    qsort(coord, 2 * n, sizeof *coord, coord_order);
    size_t unique = 0;
    for (size_t i = 0; i < 2 * n; ++i)
        if (!unique || coord[i] != coord[unique - 1])
            coord[unique++] = coord[i];
    qsort(all, n, sizeof *all, cost_order);
    size_t *skip = lab_alloc(unique, sizeof *skip);
    int64_t *cost = lab_alloc(unique, sizeof *cost);
    for (size_t k = 0; k < unique; ++k) {
        skip[k] = k;
        cost[k] = INT64_MAX;
    }
    for (size_t k = 0; k < n; ++k) {
        size_t j = unfilled(skip, first(coord, unique, all[k].lo));
        while (coord[j] <= all[k].hi) {
            cost[j] = all[k].cost;
            skip[j] = unfilled(skip, j + 1);
            j = skip[j];
        }
    }
    Band *out = lab_alloc(unique, sizeof *out);
    size_t count = 0;
    for (size_t j = 0; j + 1 < unique; ++j)
        if (cost[j] != INT64_MAX) {
            if (count && out[count - 1].hi + 1 == coord[j] && out[count - 1].cost == cost[j])
                out[count - 1].hi = coord[j + 1] - 1;
            else
                out[count++] = (Band){coord[j], coord[j + 1] - 1, cost[j]};
        }
    c->band = lab_alloc(count, sizeof *c->band);
    memcpy(c->band, out, count * sizeof *out);
    c->n = count;
    g->segments += count;
    free(out);
    free(cost);
    free(skip);
    free(coord);
    return g->segments <= MAX_SEGMENTS && lab_now() < g->deadline;
}

static Curve *combine(Network *g, int a, int b, Curve *left, Curve *right, int middle,
                      int64_t demand, int series) {
    int depth = 1 + (left->depth > right->depth ? left->depth : right->depth);
    if (depth > 512 || left->n > MAX_PAIRS / right->n || lab_now() >= g->deadline)
        return NULL;
    Band *all = lab_alloc(left->n * right->n, sizeof *all);
    size_t n = 0;
    for (size_t i = 0; i < left->n; ++i) {
        Band x = oriented(left, a, i);
        for (size_t j = 0; j < right->n; ++j) {
            Band y = oriented(right, series ? middle : a, j);
            int64_t lo, hi;
            if (series) {
                lo = x.lo > y.lo + demand ? x.lo : y.lo + demand;
                hi = x.hi < y.hi + demand ? x.hi : y.hi + demand;
            } else {
                lo = x.lo + y.lo;
                hi = x.hi + y.hi;
            }
            if (lo <= hi)
                all[n++] = (Band){lo, hi, x.cost + y.cost};
        }
    }
    Curve *c = &g->pool[g->used++];
    *c = (Curve){.a = a,
                 .b = b,
                 .kind = series ? 1 : 2,
                 .depth = depth,
                 .left = left,
                 .right = right,
                 .demand = demand};
    int ok = envelope(g, c, all, n);
    free(all);
    return ok ? c : NULL;
}

static int lookup(const Curve *c, int from, int64_t flow, int64_t *cost) {
    int64_t v = c->a == from ? flow : -flow;
    size_t lo = 0, hi = c->n;
    while (lo < hi) {
        size_t k = lo + (hi - lo) / 2;
        Band b = c->band[k];
        if (v < b.lo)
            hi = k;
        else if (v > b.hi)
            lo = k + 1;
        else {
            *cost = b.cost;
            return 1;
        }
    }
    return 0;
}
static int restore(const Curve *c, int from, int64_t flow, double *point, double deadline) {
    if (lab_now() >= deadline)
        return 0;
    int64_t v = from == c->a ? flow : -flow, cost;
    if (!lookup(c, c->a, v, &cost))
        return 0;
    if (!c->kind) {
        point[c->flow] = (double)v;
        point[c->gate] = v != 0;
        return 1;
    }
    Curve *l = c->left, *r = c->right;
    if (c->kind == 1) {
        int mid = l->a == c->a ? l->b : l->a;
        return restore(l, c->a, v, point, deadline) &&
               restore(r, mid, v - c->demand, point, deadline);
    }
    for (size_t i = 0; i < l->n; ++i) {
        Band x = oriented(l, c->a, i);
        for (size_t j = 0; j < r->n; ++j) {
            Band y = oriented(r, c->a, j);
            if (x.cost + y.cost != cost || v < x.lo + y.lo || v > x.hi + y.hi)
                continue;
            int64_t transfer = x.lo > v - y.hi ? x.lo : v - y.hi;
            return restore(l, c->a, transfer, point, deadline) &&
                   restore(r, c->a, v - transfer, point, deadline);
        }
    }
    return 0;
}

/* With one supply and nonbinding arc capacities, feasibility means reaching
 * every demand. The cheapest directed connecting tree either advances along
 * one arc or splits its remaining terminals into two disjoint subsets.
 * Reversing all arcs gives the corresponding single-demand problem. */
typedef struct {
    int from, to, flow, gate, next_in, next_out;
    int64_t cost;
} TerminalArc;
typedef struct {
    unsigned mask;
    int vertex;
} TerminalState;

static int terminal_tree(const Network *g, const int64_t *d, size_t budget, double *point,
                         size_t cols, size_t *states) {
    enum { MAX_TERMINALS = 10 };
    int supplies = 0, demands = 0, root = -1, term[MAX_TERMINALS], nt = 0;
    int64_t supply = 0, total = 0;
    for (size_t v = 0; v < g->nv; ++v) {
        supplies += d[v] < 0;
        demands += d[v] > 0;
        total += d[v];
        if (d[v] > 0)
            supply += d[v];
    }
    if (total || !supply || (supplies != 1 && demands != 1))
        return 0;
    int reverse = supplies != 1;
    if ((reverse ? supplies : demands) > MAX_TERMINALS)
        return 0;
    for (size_t v = 0; v < g->nv; ++v) {
        int64_t b = reverse ? -d[v] : d[v];
        if (b < 0)
            root = (int)v;
        if (b > 0)
            term[nt++] = (int)v;
    }
    unsigned masks = 1u << nt;
    size_t nstates = (size_t)masks * g->nv;
    if (nstates - g->nv > budget)
        return 0;
    int head[MAX_VERTICES], out[MAX_VERTICES];
    for (size_t v = 0; v < g->nv; ++v)
        head[v] = out[v] = -1;
    TerminalArc *arc = lab_alloc(g->used, sizeof *arc);
    int ne = 0, ok = 0;
    for (size_t i = 0; i < g->used; ++i) {
        const Curve *c = &g->pool[i];
        if (c->kind || c->n == 1)
            continue; /* composed curves and zero capacities */
        if (c->band[1].hi < supply) {
            free(arc);
            return 0;
        }
        int a = reverse ? c->b : c->a, b = reverse ? c->a : c->b;
        arc[ne] = (TerminalArc){a, b, c->flow, c->gate, head[b], out[a], c->band[1].cost};
        head[b] = out[a] = ne++;
    }
    const int64_t inf = INT64_MAX / 4;
    int64_t *dist = lab_alloc(nstates, sizeof *dist);
    int *choice = lab_alloc(nstates, sizeof *choice);
    unsigned char *selected = lab_alloc((size_t)ne, 1);
    TerminalState *todo = NULL;
    for (unsigned mask = 1; mask < masks; ++mask) {
        if (lab_now() >= g->deadline)
            goto done;
        size_t base = (size_t)mask * g->nv;
        for (size_t v = 0; v < g->nv; ++v)
            dist[base + v] = inf;
        if (!(mask & (mask - 1))) {
            int bit = 0;
            while ((1u << bit) != mask)
                ++bit;
            dist[base + (size_t)term[bit]] = 0;
        }
        for (unsigned sub = (mask - 1) & mask; sub; sub = (sub - 1) & mask) {
            unsigned rest = mask ^ sub;
            if (sub > rest)
                continue;
            size_t l = (size_t)sub * g->nv, r = (size_t)rest * g->nv;
            for (size_t v = 0; v < g->nv; ++v) {
                int64_t cost = dist[l + v] + dist[r + v];
                if (cost < dist[base + v]) {
                    dist[base + v] = cost;
                    choice[base + v] = (int)sub;
                }
            }
        }
        unsigned char used[MAX_VERTICES] = {0};
        for (size_t i = 0; i < g->nv; ++i) {
            int v = -1;
            int64_t best = inf;
            for (size_t j = 0; j < g->nv; ++j)
                if (!used[j] && dist[base + j] < best) {
                    v = (int)j;
                    best = dist[base + j];
                }
            if (v < 0)
                break;
            used[v] = 1;
            for (int e = head[v]; e >= 0; e = arc[e].next_in) {
                size_t u = (size_t)arc[e].from;
                if (best + arc[e].cost < dist[base + u]) {
                    dist[base + u] = best + arc[e].cost;
                    choice[base + u] = -e - 1;
                }
            }
        }
        *states += g->nv;
    }
    unsigned full = masks - 1;
    int64_t optimum = dist[(size_t)full * g->nv + (size_t)root];
    if (optimum == inf)
        goto done;
    /* Each terminal traces at most one simple path per merge level. */
    todo = lab_alloc(2 * nstates + 1, sizeof *todo);
    size_t pending = 0;
    todo[pending++] = (TerminalState){full, root};
    while (pending) {
        TerminalState t = todo[--pending];
        size_t pos = (size_t)t.mask * g->nv + (size_t)t.vertex;
        int step = choice[pos];
        choice[pos] = 0; /* already expanded states share all their arcs */
        if (step < 0) {
            int e = -step - 1;
            selected[e] = 1;
            todo[pending++] = (TerminalState){t.mask, arc[e].to};
        } else if (step) {
            todo[pending++] = (TerminalState){(unsigned)step, t.vertex};
            todo[pending++] = (TerminalState){t.mask ^ (unsigned)step, t.vertex};
        }
    }
    int parent[MAX_VERTICES], queue[MAX_VERTICES], tail = 0;
    int64_t load[MAX_VERTICES];
    for (size_t v = 0; v < g->nv; ++v) {
        parent[v] = -1;
        load[v] = reverse ? -d[v] : d[v];
    }
    queue[tail++] = root;
    parent[root] = -2;
    for (int i = 0; i < tail; ++i) {
        int v = queue[i];
        for (int e = out[v]; e >= 0; e = arc[e].next_out) {
            int u = arc[e].to;
            if (selected[e] && parent[u] == -1) {
                parent[u] = e;
                queue[tail++] = u;
            }
        }
    }
    for (int j = 0; j < nt; ++j)
        if (parent[term[j]] < 0)
            goto done;
    memset(point, 0, cols * sizeof *point);
    int64_t cost = 0;
    while (tail > 1) {
        int v = queue[--tail], e = parent[v];
        if (!load[v])
            continue;
        point[arc[e].flow] = (double)load[v];
        point[arc[e].gate] = 1;
        load[arc[e].from] += load[v];
        cost += arc[e].cost;
    }
    ok = !load[root] && cost == optimum && lab_now() < g->deadline;
done:
    free(todo);
    free(selected);
    free(choice);
    free(dist);
    free(arc);
    return ok;
}

int lab_fixed_charge_flow(const LabModel *m, double deadline, size_t limit, double *point,
                          size_t *eliminations, size_t *states) {
    *eliminations = *states = 0;
    if (!m->rows || m->cols < 2 || m->cols > 2 * MAX_ARCS || m->nnz != 2 * m->cols ||
        m->rows > MAX_ARCS + MAX_VERTICES || !limit || lab_now() >= deadline)
        return 0;
    Network g = {.deadline = deadline};
    int *vertex = lab_alloc(m->rows, sizeof *vertex), *gate = lab_alloc(m->rows, sizeof *gate);
    int *linked = lab_alloc(m->rows, sizeof *linked);
    int64_t *cap = lab_alloc(m->rows, sizeof *cap), *price = lab_alloc(m->rows, sizeof *price);
    int64_t d[MAX_VERTICES] = {0};
    mpq_t q;
    mpq_init(q);
    int ok = 0;
    size_t arcs = 0;
    for (size_t i = 0; i < m->rows; ++i) {
        vertex[i] = -1;
        gate[i] = -1;
        if (m->row_lower[i] == m->row_upper[i]) {
            if (g.nv >= MAX_VERTICES || lab_source_exact_row_bound(m, i, 0, q) ||
                !exact_int(q, &d[g.nv]) || lab_source_exact_row_bound(m, i, 1, q) ||
                mpq_cmp_si(q, (long)d[g.nv], 1))
                goto done;
            vertex[i] = (int)g.nv++;
        } else if (m->row_lower[i] != -INFINITY || m->row_upper[i] != 0 ||
                   lab_source_exact_row_bound(m, i, 1, q) || mpq_sgn(q))
            goto done;
    }
    if (!g.nv)
        goto done;
    /* No column is discarded. Every gate controls exactly one flow arc. */
    for (size_t j = 0; j < m->cols; ++j) {
        if (m->col_lower[j] != 0 || lab_source_exact_column_bound(m, j, 0, q) || mpq_sgn(q))
            goto done;
        if (!m->integer[j])
            continue;
        size_t k = m->col_start[j];
        if (m->col_start[j + 1] - k != 1 || m->col_upper[j] != 1 ||
            lab_source_exact_column_bound(m, j, 1, q) || mpq_cmp_ui(q, 1, 1))
            goto done;
        size_t row = m->row_index[k];
        if (vertex[row] >= 0 || gate[row] >= 0 || lab_source_exact_coefficient(m, k, q) ||
            !exact_int(q, &cap[row]) || cap[row] >= 0 || lab_source_exact_objective(m, j, q) ||
            !exact_int(q, &price[row]))
            goto done;
        cap[row] = -cap[row];
        if (m->maximize)
            price[row] = -price[row];
        if (price[row] < 0)
            goto done;
        gate[row] = (int)j;
    }
    g.pool = lab_alloc(3 * m->cols + 2 * (size_t)g.nv, sizeof *g.pool);
    g.edge = lab_alloc((size_t)g.nv * (size_t)g.nv, sizeof *g.edge);
    for (size_t j = 0; j < m->cols; ++j)
        if (!m->integer[j]) {
            if (lab_now() >= deadline || m->col_start[j + 1] - m->col_start[j] != 3 ||
                lab_source_exact_objective(m, j, q) || mpq_sgn(q))
                goto done;
            int u = -1, v = -1, link = -1;
            for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
                size_t row = m->row_index[k];
                if (lab_source_exact_coefficient(m, k, q))
                    goto done;
                if (vertex[row] < 0 && !mpq_cmp_ui(q, 1, 1) && link < 0)
                    link = (int)row;
                else if (vertex[row] >= 0 && !mpq_cmp_ui(q, 1, 1) && v < 0)
                    v = vertex[row];
                else if (vertex[row] >= 0 && !mpq_cmp_si(q, -1, 1) && u < 0)
                    u = vertex[row];
                else
                    goto done;
            }
            if (link < 0 || u < 0 || v < 0 || u == v || gate[link] < 0 || linked[link]++)
                goto done;
            int64_t upper = cap[link];
            if (isfinite(m->col_upper[j])) {
                int64_t bound;
                if (lab_source_exact_column_bound(m, j, 1, q) || !exact_int(q, &bound) || bound < 0)
                    goto done;
                if (bound < upper)
                    upper = bound;
            }
            Curve *c = &g.pool[g.used++];
            *c = (Curve){.a = u, .b = v, .flow = (int)j, .gate = gate[link], .n = upper ? 2 : 1};
            c->band = lab_alloc(c->n, sizeof *c->band);
            c->band[0] = (Band){0, 0, 0};
            if (upper)
                c->band[1] = (Band){1, upper, price[link]};
            Curve *other = g.edge[(size_t)u * g.nv + (size_t)v];
            if (other && !(c = combine(&g, u, v, other, c, 0, 0, 0)))
                goto done;
            connect(&g, u, v, c);
            ++arcs;
        }
    if (2 * arcs != m->cols)
        goto done;
    for (size_t i = 0; i < m->rows; ++i)
        if (vertex[i] < 0 && !linked[i])
            goto done;
    int64_t original_d[MAX_VERTICES];
    memcpy(original_d, d, g.nv * sizeof *d);
    unsigned char removed[MAX_VERTICES] = {0};
    memset(point, 0, m->cols * sizeof *point);
    for (int t = 0; t < (int)g.nv; ++t) {
        if (*eliminations >= limit || lab_now() >= deadline)
            goto done;
        int w = -1, nb[3] = {0}, degree = 3;
        for (int i = 0; i < (int)g.nv && degree > 1; ++i)
            if (!removed[i]) {
                int adj[3] = {0}, n = 0;
                for (int j = 0; j < (int)g.nv && n < 3; ++j)
                    if (g.edge[(size_t)i * g.nv + (size_t)j])
                        adj[n++] = j;
                if (n < degree) {
                    w = i;
                    degree = n;
                    memcpy(nb, adj, sizeof nb);
                }
            }
        if (w < 0) {
            ok = terminal_tree(&g, original_d, limit - *eliminations, point, m->cols, states);
            goto done;
        }
        ++*eliminations;
        removed[w] = 1;
        if (!degree) {
            if (d[w])
                goto done;
            continue;
        }
        int u = nb[0];
        Curve *left = g.edge[(size_t)u * g.nv + (size_t)w];
        connect(&g, u, w, NULL);
        if (degree == 1) {
            if (!restore(left, u, d[w], point, deadline))
                goto done;
            d[u] += d[w];
            continue;
        }
        int v = nb[1];
        Curve *right = g.edge[(size_t)w * g.nv + (size_t)v];
        connect(&g, w, v, NULL);
        Curve *c = combine(&g, u, v, left, right, w, d[w], 1);
        if (!c)
            goto done;
        Curve *other = g.edge[(size_t)u * g.nv + (size_t)v];
        if (other && !(c = combine(&g, u, v, other, c, 0, 0, 0)))
            goto done;
        connect(&g, u, v, c);
        d[v] += d[w];
    }
    ok = lab_now() < deadline;
done:
    for (size_t i = 0; i < g.used; ++i)
        free(g.pool[i].band);
    free(g.edge);
    free(g.pool);
    free(vertex);
    free(linked);
    free(gate);
    free(cap);
    free(price);
    mpq_clear(q);
    return ok;
}

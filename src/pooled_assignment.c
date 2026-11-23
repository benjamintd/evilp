#include "pooled_assignment.h"
#include "common.h"
#include "exact.h"
#include <glpk.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { MAX_SITES = 256, MAX_CLIENTS = 4096, MAX_EDGES = 160000 };
typedef struct {
    int gate, group, capacity, type;
} Site;
typedef struct {
    const LabModel *m;
    Site *site;
    int *client_row, *group_row, *quota;
    int *edge, *type_site;
    int sites, types, groups, clients;
    int64_t *cost;
    double deadline;
} Problem;

static int scaled(mpq_t q, int failed, long scale, int64_t *out) {
    if (failed || !mpz_fits_slong_p(mpq_numref(q)) ||
        mpz_cmp_ui(mpq_denref(q), (unsigned long)scale) > 0)
        return 0;
    long den = mpz_get_si(mpq_denref(q)), num = mpz_get_si(mpq_numref(q));
    if (scale % den)
        return 0;
    long factor = scale / den;
    int64_t cap = INT64_C(1000000000000) / factor;
    if (num < -cap || num > cap)
        return 0;
    *out = num * factor;
    return 1;
}
/* Each assignment uses one unit of a site and serves one client. An opening
 * supplies the site's capacity and belongs to one quota group. Equal costs
 * must match in original decimal units for every client before pooling. */
static int extract(Problem *p, mpq_t q) {
    const LabModel *m = p->m;
    if (m->rows < 3 || m->rows > 8192 || m->cols > MAX_EDGES + MAX_SITES || m->nnz != 2 * m->cols)
        return 0;
    int *site = lab_alloc(m->rows, sizeof *site);
    p->site = lab_alloc(MAX_SITES, sizeof *p->site);
    p->client_row = lab_alloc(m->rows, sizeof *p->client_row);
    p->group_row = lab_alloc(m->rows, sizeof *p->group_row);
    p->quota = lab_alloc(m->rows, sizeof *p->quota);
    int ok = 0;
    for (size_t i = 0; i < m->rows; i++) {
        site[i] = p->client_row[i] = p->group_row[i] = -1;
        if (m->row_lower[i] == -INFINITY && m->row_upper[i] == 0) {
            if (p->sites == MAX_SITES || lab_source_exact_row_bound(m, i, 1, q) || mpq_sgn(q))
                goto done;
            site[i] = p->sites++;
        } else if (m->row_lower[i] != m->row_upper[i] || !isfinite(m->row_lower[i]))
            goto done;
    }
    if (!p->sites)
        goto done;
    for (size_t j = 0; j < m->cols; j++) {
        if ((j & 127U) == 0 && lab_now() >= p->deadline)
            goto done;
        int64_t v;
        if (m->col_start[j + 1] - m->col_start[j] != 2 ||
            !scaled(q, lab_source_exact_column_bound(m, j, 0, q), 1, &v) || v ||
            !scaled(q, lab_source_exact_column_bound(m, j, 1, q), 1, &v) || v != 1)
            goto done;
        size_t a = m->col_start[j], b = a + 1;
        if (site[m->row_index[a]] < 0) {
            size_t z = a;
            a = b;
            b = z;
        }
        size_t r = m->row_index[b];
        int f = site[m->row_index[a]];
        if (f < 0 || site[r] >= 0 || !scaled(q, lab_source_exact_coefficient(m, b, q), 1, &v) ||
            v != 1 || !scaled(q, lab_source_exact_coefficient(m, a, q), 1, &v))
            goto done;
        if (v == 1) {
            if (p->group_row[r] >= 0)
                goto done;
            if (p->client_row[r] < 0)
                p->client_row[r] = p->clients++;
        } else if (v < 0 && v >= -MAX_CLIENTS) {
            if (!m->integer[j] || p->client_row[r] >= 0 || p->site[f].gate)
                goto done;
            if (p->group_row[r] < 0)
                p->group_row[r] = p->groups++;
            p->site[f] = (Site){(int)j + 1, p->group_row[r], (int)-v, 0};
        } else
            goto done;
    }
    if (!p->clients || p->clients > MAX_CLIENTS || !p->groups ||
        (size_t)(p->clients + 1) * (size_t)p->sites != m->cols)
        goto done;
    for (size_t i = 0; i < m->rows; i++) {
        if (site[i] >= 0)
            continue;
        int64_t lo, hi;
        if (!scaled(q, lab_source_exact_row_bound(m, i, 0, q), 1, &lo) ||
            !scaled(q, lab_source_exact_row_bound(m, i, 1, q), 1, &hi) || lo != hi)
            goto done;
        if (p->client_row[i] >= 0) {
            if (lo != 1)
                goto done;
        } else if (p->group_row[i] >= 0) {
            if (lo < 0 || lo > p->sites)
                goto done;
            p->quota[p->group_row[i]] = (int)lo;
        } else
            goto done;
    }
    p->edge = lab_alloc(m->cols, sizeof *p->edge);
    p->cost = lab_alloc(m->cols, sizeof *p->cost);
    for (size_t j = 0; j < m->cols; j++) {
        if ((j & 127U) == 0 && lab_now() >= p->deadline)
            goto done;
        size_t a = m->col_start[j], b = a + 1;
        if (site[m->row_index[a]] < 0) {
            size_t z = a;
            a = b;
            b = z;
        }
        int f = site[m->row_index[a]], d = p->client_row[m->row_index[b]];
        if (!p->site[f].gate ||
            !scaled(q, lab_source_exact_objective(m, j, q), 1000000000L, &p->cost[j]))
            goto done;
        if (d < 0) {
            if (p->cost[j])
                goto done;
            continue;
        }
        int *e = &p->edge[f * p->clients + d];
        if (*e)
            goto done;
        *e = (int)j + 1;
    }
    p->type_site = lab_alloc((size_t)p->sites, sizeof *p->type_site);
    for (int f = 0; f < p->sites; f++) {
        if (lab_now() >= p->deadline)
            goto done;
        int t = 0;
        for (; t < p->types; t++) {
            int g = p->type_site[t];
            if (p->site[f].capacity != p->site[g].capacity)
                continue;
            int d = 0;
            for (; d < p->clients; d++)
                if (p->cost[p->edge[f * p->clients + d] - 1] !=
                    p->cost[p->edge[g * p->clients + d] - 1])
                    break;
            if (d == p->clients)
                break;
        }
        if (t == p->types)
            p->type_site[p->types++] = f;
        p->site[f].type = t;
    }
    ok = p->types * 2 <= p->sites && lab_now() < p->deadline;
done:
    free(site);
    return ok;
}
typedef struct {
    int a, b, lo, hi, flow;
    int64_t cost;
} Arc;
typedef struct {
    const Problem *p;
    glp_prob *lp;
    Arc *arc;
    int arcs, vertices, *supply, *net, *load_arc, cap, *lower, *upper;
    int *path, *best;
    int64_t *price, incumbent;
    size_t nodes, iterations, iteration_limit, node_limit, solves;
    double deadline;
    int stopped;
} Network;
static int arc(Network *s, int a, int b, int lo, int hi, int64_t cost) {
    int j = ++s->arcs;
    s->arc[j] = (Arc){a, b, lo, hi, 0, cost};
    int idx[] = {0, a, b};
    double val[] = {0, 1, -1};
    glp_set_mat_col(s->lp, j, 2, idx, val);
    glp_set_col_bnds(s->lp, j, hi == lo ? GLP_FX : GLP_DB, lo, hi);
    glp_set_obj_coef(s->lp, j, (double)cost);
    return j;
}
/* Ship one unit per client and route unused capacity from one dummy
 * supplier. A class-throughput arc counts units of open capacity. Every
 * original solution has a feasible flow here, even for fractional clients.
 * Branching imposes integrality on the number of sites in each class. */
static int build_network(Network *s) {
    const Problem *p = s->p;
    s->cap = p->site[0].capacity;
    int capacity = 0;
    for (int i = 0; i < p->sites; i++)
        if (p->site[i].capacity != s->cap)
            return 0;
    for (int i = 0; i < p->groups; i++)
        capacity += p->quota[i] * s->cap;
    /* An integer opening consumes a distinct original site. This also
     * bounds all dummy and group supplies by MAX_SITES * MAX_CLIENTS. */
    if (capacity > p->sites * s->cap)
        return 0;
    int excess = capacity - p->clients;
    if (excess < 0)
        return 0;
    int dummy = p->clients + 1, types = dummy + 1, outs = types + p->types,
        groups = outs + p->types;
    s->vertices = groups + p->groups - 1;
    int arcs = (p->clients + 2) * p->types + p->sites;
    s->arc = lab_alloc((size_t)arcs + 1, sizeof *s->arc);
    s->load_arc = lab_alloc((size_t)p->types, sizeof *s->load_arc);
    s->lower = lab_alloc((size_t)p->types, sizeof *s->lower);
    s->upper = lab_alloc((size_t)p->types, sizeof *s->upper);
    s->supply = lab_alloc((size_t)s->vertices + 1, sizeof *s->supply);
    s->net = lab_alloc((size_t)s->vertices + 1, sizeof *s->net);
    s->price = lab_alloc((size_t)s->vertices + 1, sizeof *s->price);
    s->best = lab_alloc(p->m->cols, sizeof *s->best);
    s->path = lab_alloc(p->m->cols, sizeof *s->path);
    s->lp = glp_create_prob();
    glp_add_rows(s->lp, s->vertices);
    glp_add_cols(s->lp, arcs);
    for (int d = 0; d < p->clients; d++) {
        s->supply[d + 1] = 1;
        for (int t = 0; t < p->types; t++) {
            int j = p->edge[p->type_site[t] * p->clients + d] - 1;
            arc(s, d + 1, types + t, 0, 1, (p->m->maximize ? -1 : 1) * p->cost[j]);
        }
    }
    s->supply[dummy] = excess;
    for (int f = 0; f < p->sites; f++)
        s->upper[p->site[f].type]++;
    for (int t = 0; t < p->types; t++) {
        arc(s, dummy, types + t, 0, excess, 0);
        s->load_arc[t] = arc(s, types + t, outs + t, 0, s->cap * s->upper[t], 0);
    }
    for (int f = 0; f < p->sites; f++)
        arc(s, outs + p->site[f].type, groups + p->site[f].group, 0, s->cap, 0);
    for (int g = 0; g < p->groups; g++)
        s->supply[groups + g] = -p->quota[g] * s->cap;
    for (int i = 1; i <= s->vertices; i++)
        glp_set_row_bnds(s->lp, i, GLP_FX, s->supply[i], s->supply[i]);
    glp_std_basis(s->lp);
    return 1;
}
/* Unit augmentations realize the lower and upper opening counts in the
 * group-by-type bipartite graph. Every arc denotes one original site.
 * The first pass fills mandatory type counts. The second fills the group
 * quotas while preserving those minima. An alternating path exchanges open
 * sites without changing the degree of an intermediate vertex. */
static int open_sites(Network *s, const int *lo, const int *hi) {
    const Problem *p = s->p;
    int *on = s->path;
    int degree[2 * MAX_SITES] = {0};
    memset(on, 0, p->m->cols * sizeof *on);
    for (int pass = 0; pass < 2; pass++) {
        int roots = pass ? p->groups : p->types;
        for (int i = 0; i < roots; i++) {
            int root = pass ? i : p->groups + i, target = pass ? p->quota[i] : lo[i];
            while (degree[root] < target) {
                int parent[2 * MAX_SITES], queue[2 * MAX_SITES];
                for (int k = 0; k < p->types + p->groups; k++)
                    parent[k] = -1;
                parent[root] = -2;
                queue[0] = root;
                int end = -1, tail = 1;
                for (int k = 0; k < tail && end < 0; k++) {
                    int v = queue[k];
                    for (int f = 0; f < p->sites; f++) {
                        const Site *a = &p->site[f];
                        int group = a->group, type = p->groups + a->type;
                        int from = pass ? group : type, to = pass ? type : group;
                        int w = -1;
                        if (!on[a->gate - 1] && v == from)
                            w = to;
                        if (on[a->gate - 1] && v == to)
                            w = from;
                        if (w < 0 || parent[w] != -1)
                            continue;
                        parent[w] = f;
                        queue[tail++] = w;
                        int cap = pass ? hi[a->type] : p->quota[group];
                        if (w == to && degree[w] < cap) {
                            end = w;
                            break;
                        }
                    }
                }
                if (end < 0)
                    return 0;
                ++degree[root];
                ++degree[end];
                while (end != root) {
                    const Site *a = &p->site[parent[end]];
                    on[a->gate - 1] = !on[a->gate - 1];
                    end = end == a->group ? p->groups + a->type : a->group;
                }
            }
        }
    }
    return 1;
}
static int assign_clients(Network *s) {
    const Problem *p = s->p;
    int need[MAX_SITES], load[MAX_SITES] = {0};
    for (int t = 0; t < p->types; t++)
        need[t] = s->arc[s->load_arc[t]].flow / s->cap;
    if (!open_sites(s, need, need))
        return 0;
    for (int d = 0; d < p->clients; d++) {
        int t = 0;
        while (t < p->types && !s->arc[1 + d * p->types + t].flow)
            t++;
        if (t == p->types)
            return 0;
        int f = 0;
        while (f < p->sites &&
               (p->site[f].type != t || !s->path[p->site[f].gate - 1] || load[f] == s->cap))
            f++;
        if (f == p->sites)
            return 0;
        load[f]++;
        s->path[p->edge[f * p->clients + d] - 1] = 1;
    }
    return 1;
}
/* Exact node lower bound from integer flow conservation, original
 * integer costs, and complementary slackness on every bounded arc. */
static int certify(Network *s, int64_t *cost) {
    *cost = 0;
    memset(s->net, 0, ((size_t)s->vertices + 1) * sizeof *s->net);
    for (int i = 1; i <= s->vertices; i++) {
        double y = glp_get_row_dual(s->lp, i);
        if (!isfinite(y) || fabs(y) > 0x1p51)
            return 0;
        s->price[i] = (int64_t)round(y);
    }
    for (int j = 1; j <= s->arcs; j++) {
        Arc *a = &s->arc[j];
        double z = glp_get_col_prim(s->lp, j);
        if (!isfinite(z) || z < a->lo - .1 || z > a->hi + .1)
            return 0;
        a->flow = (int)round(z);
        int64_t reduced = a->cost - s->price[a->a] + s->price[a->b];
        if ((a->flow < a->hi && reduced < 0) || (a->flow > a->lo && reduced > 0))
            return 0;
        s->net[a->a] += a->flow;
        s->net[a->b] -= a->flow;
        *cost += a->cost * a->flow;
    }
    for (int i = 1; i <= s->vertices; i++)
        if (s->net[i] != s->supply[i])
            return 0;
    return 1;
}
static void explore_network(Network *s) {
    if (s->stopped)
        return;
    if (s->nodes == s->node_limit || lab_now() >= s->deadline) {
        s->stopped = 1;
        return;
    }
    s->nodes++;
    const Problem *p = s->p;
    for (int t = 0; t < p->types; t++) {
        int j = s->load_arc[t];
        Arc *a = &s->arc[j];
        a->lo = s->lower[t] * s->cap;
        a->hi = s->upper[t] * s->cap;
        glp_set_col_bnds(s->lp, j, a->lo == a->hi ? GLP_FX : GLP_DB, a->lo, a->hi);
    }
    glp_smcp o;
    glp_init_smcp(&o);
    o.msg_lev = GLP_MSG_OFF;
    o.meth = GLP_DUALP;
    o.it_lim = (int)(s->iteration_limit > INT_MAX ? INT_MAX : s->iteration_limit);
    o.tm_lim = (int)fmin(INT_MAX, fmax(1, 900 * (s->deadline - lab_now())));
    ++s->solves;
    int before = glp_get_it_cnt(s->lp);
    int rc = glp_simplex(s->lp, &o);
    s->iterations += (size_t)(glp_get_it_cnt(s->lp) - before);
    if (rc || glp_get_status(s->lp) != GLP_OPT) {
        /* Every client can use every site, and total capacity is sufficient.
         * Only the opening quotas can make the current count box infeasible.
         * Confirm that obstruction with exact integral augmentations. */
        if (!rc && glp_get_status(s->lp) == GLP_NOFEAS && !open_sites(s, s->lower, s->upper))
            return;
        s->stopped = 1;
        return;
    }
    int64_t bound;
    if (!certify(s, &bound)) {
        s->stopped = 1;
        return;
    }
    if (bound >= s->incumbent)
        return;
    int branch = -1, frac = 0, value = 0;
    for (int t = 0; t < p->types; t++) {
        int v = s->arc[s->load_arc[t]].flow, rem = v % s->cap;
        int f = rem < s->cap - rem ? rem : s->cap - rem;
        if (f > frac) {
            branch = t;
            frac = f;
            value = v;
        }
    }
    if (branch < 0) {
        if (!assign_clients(s)) {
            s->stopped = 1;
            return;
        }
        memcpy(s->best, s->path, p->m->cols * sizeof *s->best);
        s->incumbent = bound;
        return;
    }
    int lo = s->lower[branch], hi = s->upper[branch], mid = value / s->cap;
    s->upper[branch] = mid;
    explore_network(s);
    s->upper[branch] = hi;
    s->lower[branch] = mid + 1;
    explore_network(s);
    s->lower[branch] = lo;
}
LabSearchResult lab_pooled_assignment(const LabModel *m, double deadline, size_t node_limit,
                          size_t iteration_limit, double *point, LabPooledStats *stats) {
    LabSearchResult result = {.route = "pooled_assignment"};
    if (!stats)
        return result;
    *stats = (LabPooledStats){0};
    if (!m || !m->integer || !m->rows || !m->cols || m->rows > 8192 ||
        m->cols > MAX_EDGES + MAX_SITES || m->nnz != 2 * m->cols || !point || !node_limit ||
        !iteration_limit || !isfinite(deadline) || deadline <= lab_now())
        return result;
    Problem p = {.m = m, .deadline = deadline};
    Network s = {.p = &p,
                 .deadline = deadline,
                 .incumbent = INT64_MAX,
                 .node_limit = node_limit,
                 .iteration_limit = iteration_limit};
    mpq_t q;
    mpq_init(q);
    if (extract(&p, q) && build_network(&s)) {
        explore_network(&s);
        if (s.incumbent != INT64_MAX) {
            for (size_t j = 0; j < m->cols; j++)
                point[j] = s.best[j];
            result.point = point;
            if (!s.stopped && lab_now() < deadline) {
                result.status = LAB_MIP_OPTIMAL;
                result.closure = LAB_CLOSE_EXHAUSTED;
            }
        }
    }
    stats->nodes = s.nodes;
    stats->lp_solves = s.solves;
    stats->iterations = s.iterations;
    if (s.lp)
        glp_delete_prob(s.lp);
    free(s.best);
    free(s.path);
    free(s.price);
    free(s.net);
    free(s.supply);
    free(s.lower);
    free(s.upper);
    free(s.load_arc);
    free(s.arc);
    free(p.site);
    free(p.client_row);
    free(p.group_row);
    free(p.quota);
    free(p.edge);
    free(p.cost);
    free(p.type_site);
    mpq_clear(q);
    return result;
}

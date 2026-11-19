#include "lot_sizing.h"
#include "common.h"
#include "exact.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { MAX_LINES = 16, MAX_SITES = 64, MAX_TOTAL = MAX_LINES * MAX_SITES, MAX_SCALAR = 1000000 };
#define NO_PLAN (INT64_MAX / 4)

typedef struct {
    int demand, group, supply, gate, unit_cost, setup_cost, supply_capacity;
    int degree, neighbor[2], incoming[2], outgoing[2];
} Site;

typedef struct {
    int count, site[MAX_SITES];
    /* Cost of supplying the first i demand sites entirely from site j. */
    int64_t prefix[MAX_SITES][MAX_SITES + 1];
} Line;

typedef struct {
    uint64_t open, forbidden;
    int64_t cost;
} Plan;

typedef struct {
    const LabModel *m;
    Site *sites;
    Line *lines;
    int line_count, site_count, group_count;
    int *capacity, *cost;
    double deadline;
    size_t nodes, limit;
    int interrupted;
    int64_t incumbent;
    Plan best[MAX_LINES];
    int owner[MAX_TOTAL], position[MAX_TOTAL];
} Routing;

static int small_integer(mpq_srcptr q, int *out) {
    if (mpz_cmp_ui(mpq_denref(q), 1) || mpz_cmpabs_ui(mpq_numref(q), MAX_SCALAR) > 0)
        return 0;
    *out = (int)mpz_get_si(mpq_numref(q));
    return 1;
}

/* Every continuous column is a supply at one node or a directed movement
 * between neighboring nodes. Each binary opens one supply and belongs to
 * one at-most-one row. Every matrix entry and every finite side is checked. */
static int read_graph(Routing *s) {
    const LabModel *m = s->m;
    int *balance = lab_alloc(m->rows, sizeof *balance);
    int *group = lab_alloc(m->rows, sizeof *group);
    int *supply = lab_alloc(m->rows, sizeof *supply);
    int *gate = lab_alloc(m->rows, sizeof *gate);
    int *gate_size = lab_alloc(m->rows, sizeof *gate_size);
    int *supply_site = lab_alloc(m->rows, sizeof *supply_site);
    int *gate_group = lab_alloc(m->rows, sizeof *gate_group);
    mpq_t q;
    mpq_init(q);
    int ok = 0;
    for (size_t i = 0; i < m->rows; ++i) {
        balance[i] = group[i] = supply[i] = gate[i] = -1;
        if (m->row_lower[i] == m->row_upper[i]) {
            if (s->site_count == MAX_TOTAL || !isfinite(m->row_lower[i]))
                goto done;
            int demand;
            if (lab_source_exact_row_bound(m, i, 0, q) || !small_integer(q, &demand) ||
                demand <= 0 || lab_source_exact_row_bound(m, i, 1, q) || mpq_cmp_si(q, demand, 1))
                goto done;
            balance[i] = s->site_count++;
            Site *v = &s->sites[balance[i]];
            v->demand = demand;
            v->supply = -1;
        } else {
            if (m->row_lower[i] != -INFINITY || !isfinite(m->row_upper[i]) ||
                lab_source_exact_row_bound(m, i, 1, q))
                goto done;
            if (mpq_cmp_ui(q, 1, 1) == 0) {
                if (s->group_count == MAX_TOTAL)
                    goto done;
                group[i] = s->group_count++;
            } else if (mpq_sgn(q))
                goto done;
        }
    }
    if (!s->site_count || !s->group_count)
        goto done;
    for (size_t j = 0; j < m->cols; ++j) {
        size_t k = m->col_start[j];
        int a, b;
        if ((j % 128 == 0 && lab_now() >= s->deadline) || m->col_start[j + 1] - k != 2 ||
            m->col_lower[j] != 0 || lab_source_exact_column_bound(m, j, 0, q) || mpq_sgn(q) ||
            !isfinite(m->col_upper[j]) || lab_source_exact_column_bound(m, j, 1, q) ||
            !small_integer(q, &s->capacity[j]) || s->capacity[j] <= 0 ||
            lab_source_exact_objective(m, j, q) || !small_integer(q, &s->cost[j]) ||
            lab_source_exact_coefficient(m, k, q) || !small_integer(q, &a) ||
            lab_source_exact_coefficient(m, k + 1, q) || !small_integer(q, &b))
            goto done;
        if (m->maximize)
            s->cost[j] = -s->cost[j];
        if (s->cost[j] < 0)
            goto done;
        size_t r = m->row_index[k], t = m->row_index[k + 1];
        if (r == t)
            goto done;
        /* Put the balance row first for a supply, the exclusive row first
         * for a gate, and the negative coefficient first for a movement. */
        if ((balance[t] >= 0 && balance[r] < 0) || group[t] >= 0 ||
            (balance[t] >= 0 && balance[r] >= 0 && a > 0)) {
            size_t swap = r;
            r = t;
            t = swap;
            int c = a;
            a = b;
            b = c;
        }
        if (m->integer[j]) {
            if (group[r] < 0 || balance[t] >= 0 || group[t] >= 0 || gate[t] >= 0 || a != 1 ||
                b >= 0 || s->capacity[j] != 1)
                goto done;
            gate[t] = (int)j;
            gate_size[t] = -b;
            gate_group[t] = group[r];
        } else if (balance[r] >= 0 && balance[t] >= 0) {
            if (a != -1 || b != 1)
                goto done;
            for (int end = 0; end < 2; ++end) {
                int other = end ? balance[r] : balance[t];
                Site *v = &s->sites[end ? balance[t] : balance[r]];
                int slot = 0;
                while (slot < v->degree && v->neighbor[slot] != other)
                    ++slot;
                if (slot == 2)
                    goto done;
                if (slot == v->degree) {
                    ++v->degree;
                    v->neighbor[slot] = other;
                    v->incoming[slot] = v->outgoing[slot] = -1;
                }
                int *arc = end ? &v->incoming[slot] : &v->outgoing[slot];
                if (*arc >= 0)
                    goto done;
                *arc = (int)j;
            }
        } else {
            if (balance[r] < 0 || group[t] >= 0 || a != 1 || b != 1 || supply[t] >= 0)
                goto done;
            supply[t] = (int)j;
            supply_site[t] = balance[r];
        }
    }
    for (size_t i = 0; i < m->rows; ++i) {
        if (balance[i] >= 0 || group[i] >= 0)
            continue;
        if (gate[i] < 0 || supply[i] < 0)
            goto done;
        Site *v = &s->sites[supply_site[i]];
        if (v->supply >= 0)
            goto done;
        v->supply = supply[i];
        v->gate = gate[i];
        v->group = gate_group[i];
        v->setup_cost = s->cost[v->gate];
        v->unit_cost = s->cost[v->supply];
        v->supply_capacity =
            s->capacity[v->supply] < gate_size[i] ? s->capacity[v->supply] : gate_size[i];
    }
    ok = 1;
done:
    mpq_clear(q);
    free(balance);
    free(group);
    free(supply);
    free(gate);
    free(gate_size);
    free(supply_site);
    free(gate_group);
    return ok;
}

/* Each component must be a bidirectional path. A capacity at least the
 * component's total demand never restricts a direct shipment plan. */
static int read_lines(Routing *s) {
    unsigned char seen[MAX_TOTAL] = {0};
    int visited = 0;
    for (int root = 0; root < s->site_count; ++root) {
        if (seen[root] || s->sites[root].degree > 1)
            continue;
        if (s->line_count == MAX_LINES || lab_now() >= s->deadline)
            return 0;
        Line *l = &s->lines[s->line_count++];
        int prev = -1, total = 0;
        for (int node = root; node >= 0;) {
            Site *v = &s->sites[node];
            if (v->supply < 0 || l->count == MAX_SITES || seen[node])
                return 0;
            seen[node] = 1;
            ++visited;
            l->site[l->count++] = node;
            total += v->demand;
            int next = -1;
            for (int h = 0; h < v->degree; ++h)
                if (v->neighbor[h] != prev)
                    next = v->neighbor[h];
            prev = node;
            node = next;
        }
        if (total > MAX_SCALAR)
            return 0;
        int right[MAX_SITES] = {0}, left[MAX_SITES] = {0};
        for (int j = 0; j < l->count; ++j) {
            Site *v = &s->sites[l->site[j]];
            if (v->supply_capacity < total)
                return 0;
            for (int h = 0; h < v->degree; ++h) {
                int out = v->outgoing[h], in = v->incoming[h];
                if (out < 0 || in < 0 || s->capacity[out] < total || s->capacity[in] < total)
                    return 0;
                if (j + 1 < l->count && v->neighbor[h] == l->site[j + 1]) {
                    right[j + 1] = right[j] + s->cost[out];
                    left[j + 1] = left[j] + s->cost[in];
                }
            }
        }
        for (int j = 0; j < l->count; ++j)
            for (int i = 0; i < l->count; ++i) {
                int transport = j < i ? right[i] - right[j] : left[j] - left[i];
                int delivered = s->sites[l->site[j]].unit_cost + transport;
                l->prefix[j][i + 1] =
                    l->prefix[j][i] + (int64_t)delivered * s->sites[l->site[i]].demand;
            }
    }
    return visited == s->site_count;
}

/* For any open sites, a cheapest supplier can serve a contiguous interval:
 * the difference in delivered unit cost for two ordered sites is monotone
 * along the path, including with asymmetric nonnegative transport costs.
 * Partition demand into intervals and minimize setup plus delivered costs.
 * Repeatedly paying a nonnegative setup can never improve the true minimum.
 * Prefix sums and a running best interval start give quadratic work. */
static Plan plan_line(const Routing *s, int k, uint64_t forbidden) {
    const Line *l = &s->lines[k];
    int64_t dp[MAX_SITES + 1] = {0}, best[MAX_SITES];
    uint64_t choices[MAX_SITES + 1] = {0};
    int start[MAX_SITES] = {0};
    for (int j = 0; j < l->count; ++j)
        best[j] = NO_PLAN;
    for (int end = 1; end <= l->count; ++end) {
        dp[end] = NO_PLAN;
        for (int j = 0; j < l->count; ++j) {
            uint64_t bit = UINT64_C(1) << j;
            if (forbidden & bit)
                continue;
            int begin = end - 1;
            int64_t offer = dp[begin] + s->sites[l->site[j]].setup_cost - l->prefix[j][begin];
            if (offer < best[j]) {
                best[j] = offer;
                start[j] = begin;
            }
            int64_t cost = best[j] + l->prefix[j][end];
            if (cost < dp[end]) {
                dp[end] = cost;
                choices[end] = choices[start[j]] | bit;
            }
        }
    }
    return (Plan){choices[l->count], forbidden, dp[l->count]};
}

static void search(Routing *s, Plan *plan, int64_t bound) {
    if (bound >= s->incumbent || s->interrupted)
        return;
    if (s->nodes >= s->limit || lab_now() >= s->deadline) {
        s->interrupted = 1;
        return;
    }
    ++s->nodes;
    for (int g = 0; g < s->group_count; ++g)
        s->owner[g] = -1;
    for (int k = 0; k < s->line_count; ++k) {
        const Line *l = &s->lines[k];
        for (int j = 0; j < l->count; ++j) {
            if (!(plan[k].open & (UINT64_C(1) << j)))
                continue;
            int g = s->sites[l->site[j]].group;
            if (s->owner[g] < 0) {
                s->owner[g] = k;
                s->position[g] = j;
                continue;
            }
            /* Every feasible solution forbids at least one of this pair. */
            int ks[] = {s->owner[g], k}, js[] = {s->position[g], j};
            Plan child[2];
            int64_t lower[2];
            for (int i = 0; i < 2; ++i) {
                Plan old = plan[ks[i]];
                child[i] = plan_line(s, ks[i], old.forbidden | (UINT64_C(1) << js[i]));
                lower[i] = bound - old.cost + child[i].cost;
            }
            int first = lower[1] < lower[0];
            for (int t = 0; t < 2; ++t) {
                int i = first ^ t;
                Plan old = plan[ks[i]];
                plan[ks[i]] = child[i];
                search(s, plan, lower[i]);
                plan[ks[i]] = old;
            }
            return;
        }
    }
    s->incumbent = bound;
    memcpy(s->best, plan, (size_t)s->line_count * sizeof *plan);
}

static int reconstruct(const Routing *s, double *x) {
    memset(x, 0, s->m->cols * sizeof *x);
    for (int k = 0; k < s->line_count; ++k) {
        const Line *l = &s->lines[k];
        for (int i = 0; i < l->count; ++i) {
            int site = -1;
            int64_t least = NO_PLAN;
            for (int j = 0; j < l->count; ++j) {
                int64_t cost = l->prefix[j][i + 1] - l->prefix[j][i];
                if ((s->best[k].open & (UINT64_C(1) << j)) && cost < least) {
                    least = cost;
                    site = j;
                }
            }
            if (site < 0)
                return 0;
            const Site *supply = &s->sites[l->site[site]];
            int amount = s->sites[l->site[i]].demand;
            x[supply->gate] = 1;
            x[supply->supply] += amount;
            for (int cur = site; cur != i;) {
                int next = cur + (cur < i ? 1 : -1);
                const Site *v = &s->sites[l->site[cur]];
                int arc = v->neighbor[0] == l->site[next] ? v->outgoing[0] : v->outgoing[1];
                x[arc] += amount;
                cur = next;
            }
        }
    }
    int64_t cost = 0;
    for (size_t j = 0; j < s->m->cols; ++j)
        cost += (int64_t)x[j] * s->cost[j];
    return cost == s->incumbent;
}

int lab_lot_sizing(const LabModel *m, double deadline, size_t limit, double *point, size_t *nodes) {
    *nodes = 0;
    if (!limit || !m->integer || !m->rows || m->rows > 3 * MAX_TOTAL || !m->cols ||
        m->cols > 4 * MAX_TOTAL || m->nnz != 2 * m->cols || lab_now() >= deadline)
        return 0;
    Routing s = {.m = m, .limit = limit, .deadline = deadline, .incumbent = NO_PLAN};
    s.sites = lab_alloc(m->rows, sizeof *s.sites);
    s.cost = lab_alloc(m->cols, sizeof *s.cost);
    s.capacity = lab_alloc(m->cols, sizeof *s.capacity);
    int ok = 0;
    if (read_graph(&s)) {
        s.lines = lab_alloc(MAX_LINES, sizeof *s.lines);
        if (read_lines(&s)) {
            Plan plan[MAX_LINES];
            int64_t root = 0;
            for (int k = 0; k < s.line_count; ++k) {
                plan[k] = plan_line(&s, k, 0);
                root += plan[k].cost;
            }
            search(&s, plan, root);
            ok = !s.interrupted && s.incumbent != NO_PLAN && reconstruct(&s, point) &&
                 lab_now() < deadline;
        }
    }
    *nodes = s.nodes;
    free(s.lines);
    free(s.capacity);
    free(s.cost);
    free(s.sites);
    return ok;
}

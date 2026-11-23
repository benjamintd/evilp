#include "integer_mip.h"
#include "common.h"
#include "exact.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { ROW_CAP = 65536, COL_CAP = 65536, ENTRY_CAP = 400000 };
#define MAX_DEN INT64_C(1000000000000)
#define MAX_EXACT (INT64_C(1) << 46)
#define MAX_SCALAR INT64_C(100000000)
#define MAX_ACTIVITY (INT64_C(1) << 43)

typedef struct {
    int64_t scale, divisor, lo, hi, first_value, largest;
    size_t entries, first, second;
    int has_lo, has_hi, equal_columns, eq, shift;
} Row;
typedef struct {
    size_t parent, group, next;
    int64_t lo, hi, cost;
    int has_lo, has_hi, capped;
} Column;
typedef struct {
    size_t head;
    int64_t lo, hi, cost, chosen;
    int has_lo, has_hi;
    SCIP_VAR *var;
} Group;
typedef struct {
    const LabModel *model;
    const double *presolve_lower, *presolve_upper;
    Row *row;
    Column *col;
    Group *group;
    int64_t *coefficient;
    int64_t obj_scale, obj_gcd;
    size_t groups, merged, path_fixed, presolve_fixed;
    double deadline;
} Projection;

static int64_t gcd64(int64_t a, int64_t b) {
    a = llabs(a);
    b = llabs(b);
    while (b) {
        int64_t r = a % b;
        a = b;
        b = r;
    }
    return a;
}
static int add64_ok(int64_t a, int64_t b, int64_t *out) {
    if ((b > 0 && a > INT64_MAX - b) || (b < 0 && a < INT64_MIN - b))
        return 0;
    *out = a + b;
    return 1;
}
static int sub64_ok(int64_t a, int64_t b, int64_t *out) {
    if ((b < 0 && a > INT64_MAX + b) || (b > 0 && a < INT64_MIN + b))
        return 0;
    *out = a - b;
    return 1;
}
static int mul64_ok(int64_t a, int64_t b, int64_t *out) {
    if (!a || !b) {
        *out = 0;
        return 1;
    }
    if ((a == INT64_MIN && b == -1) || (b == INT64_MIN && a == -1))
        return 0;
    if (a > 0) {
        if (b > 0 ? a > INT64_MAX / b : b < INT64_MIN / a)
            return 0;
    } else {
        if (b > 0 ? a < INT64_MIN / b : a < INT64_MAX / b)
            return 0;
    }
    *out = a * b;
    return 1;
}
static int include_scale(int64_t *scale, mpq_srcptr q, int error) {
    if (error || !mpz_fits_slong_p(mpq_denref(q)) || !mpz_fits_slong_p(mpq_numref(q)))
        return 0;
    int64_t d = (int64_t)mpz_get_si(mpq_denref(q));
    if (d <= 0)
        return 0;
    int64_t g = gcd64(d, *scale), mult = d / g;
    if (*scale > MAX_DEN / mult)
        return 0;
    *scale *= mult;
    return 1;
}
static int exact_units(mpq_srcptr q, int error, int64_t scale, int64_t *out) {
    if (error || !mpz_fits_slong_p(mpq_numref(q)) || !mpz_fits_slong_p(mpq_denref(q)))
        return 0;
    int64_t num = (int64_t)mpz_get_si(mpq_numref(q));
    int64_t den = (int64_t)mpz_get_si(mpq_denref(q));
    if (den <= 0 || scale % den)
        return 0;
    int64_t factor = scale / den;
    int64_t value = 0;
    if (!mul64_ok(num, factor, &value) || llabs(value) > MAX_EXACT)
        return 0;
    *out = value;
    return 1;
}
static size_t representative(Column *columns, size_t j) {
    while (columns[j].parent != j) {
        columns[j].parent = columns[columns[j].parent].parent;
        j = columns[j].parent;
    }
    return j;
}

/* Exhaust the simple paths of each small, isolated unit-flow component.
 * Nonnegative normalized cost and one-sided nonnegative capacity usage make
 * cycle deletion valid for the source model. All searches fail open on size. */
typedef struct {
    size_t u, v, j, next;
    int ok;
} FlowEdge;
typedef struct {
    size_t *parent, *head;
    FlowEdge *e;
    unsigned char *visited, *chosen, *some, *all;
    size_t *member, n, steps, paths, target;
    double deadline;
    int stopped;
} PathState;

static size_t flow_root(size_t *a, size_t v) {
    while (a[v] != v) {
        a[v] = a[a[v]];
        v = a[v];
    }
    return v;
}
static void enumerate_paths(PathState *p, size_t v) {
    if (++p->steps > 150000 || p->paths == 16384 ||
        ((p->steps & 1023U) == 0 && lab_now() >= p->deadline)) {
        p->stopped = 1;
        return;
    }
    if (v == p->target) {
        for (size_t t = 0; t < p->n; ++t) {
            size_t j = p->member[t];
            p->some[j] |= p->chosen[j];
            p->all[j] &= p->chosen[j];
        }
        p->paths++;
        return;
    }
    p->visited[v] = 1;
    for (size_t j = p->head[v]; j && !p->stopped; j = p->e[j - 1].next) {
        FlowEdge *a = &p->e[j - 1];
        if (p->visited[a->v])
            continue;
        p->chosen[j - 1] = 1;
        enumerate_paths(p, a->v);
        p->chosen[j - 1] = 0;
    }
    p->visited[v] = 0;
}
static size_t prune_flow_paths(Projection *proj) {
    const LabModel *m = proj->model;
    PathState p = {.deadline = proj->deadline};
    p.parent = lab_alloc(m->rows, sizeof *p.parent);
    p.head = lab_alloc(m->rows, sizeof *p.head);
    p.e = lab_alloc(m->cols, sizeof *p.e);
    p.visited = lab_alloc(m->rows, 1);
    p.chosen = lab_alloc(m->cols, 1);
    p.some = lab_alloc(m->cols, 1);
    p.all = lab_alloc(m->cols, 1);
    p.member = lab_alloc(m->cols, sizeof *p.member);
    unsigned char *bad = lab_alloc(m->rows, 1);
    size_t *count = lab_alloc(m->rows, sizeof *count);
    size_t *ec = lab_alloc(m->rows, sizeof *ec);
    size_t *member_head = lab_alloc(m->rows, sizeof *member_head);
    size_t *next = lab_alloc(m->cols, sizeof *next);
    size_t *source = lab_alloc(m->rows, sizeof *source);
    size_t *sink = lab_alloc(m->rows, sizeof *sink);
    for (size_t i = 0; i < m->rows; ++i)
        p.parent[i] = i;
    for (size_t j = 0; j < m->cols; ++j) {
        size_t first = SIZE_MAX, positive = SIZE_MAX, negative = SIZE_MAX;
        Column *c = &proj->col[j];
        int ok = c->has_lo && c->has_hi && c->lo == 0 && c->hi == 1 && c->cost >= 0;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            Row *r = &proj->row[i];
            int64_t z = proj->coefficient[k] / r->divisor;
            if (!z)
                continue;
            if (r->eq) {
                if (first == SIZE_MAX)
                    first = i;
                else
                    p.parent[flow_root(p.parent, first)] = flow_root(p.parent, i);
                if (z == 1 && positive == SIZE_MAX)
                    positive = i;
                else if (z == -1 && negative == SIZE_MAX)
                    negative = i;
                else
                    ok = 0;
            } else if ((z > 0 && r->has_lo) || (z < 0 && r->has_hi)) {
                ok = 0;
            }
        }
        if (first != SIZE_MAX) {
            if (!ok || positive == SIZE_MAX || negative == SIZE_MAX) {
                bad[first] = 1;
                continue;
            }
            p.e[j] = (FlowEdge){positive, negative, j, p.head[positive], 1};
            p.head[positive] = j + 1;
        }
    }
    for (size_t i = 0; i < m->rows; ++i)
        if (proj->row[i].eq) {
            size_t r = flow_root(p.parent, i);
            count[r]++;
            bad[r] |= bad[i];
            int64_t rhs = proj->row[i].lo / proj->row[i].divisor;
            if (rhs == 1 && !source[r])
                source[r] = i + 1;
            else if (rhs == -1 && !sink[r])
                sink[r] = i + 1;
            else if (rhs)
                bad[r] = 1;
        }
    for (size_t j = 0; j < m->cols; ++j)
        if (p.e[j].ok) {
            size_t r = flow_root(p.parent, p.e[j].u);
            ec[r]++;
            next[j] = member_head[r];
            member_head[r] = j + 1;
        }
    size_t fixed = 0;
    for (size_t r = 0; r < m->rows; ++r)
        if (!bad[r] && count[r] <= 64 && ec[r] <= 512 && source[r] && sink[r]) {
            if (lab_now() >= proj->deadline)
                break;
            p.n = p.steps = p.paths = 0;
            p.stopped = 0;
            p.target = sink[r] - 1;
            for (size_t j = member_head[r]; j; j = next[j - 1]) {
                p.member[p.n++] = j - 1;
                p.all[j - 1] = 1;
            }
            enumerate_paths(&p, source[r] - 1);
            if (p.stopped || !p.paths)
                continue;
            for (size_t t = 0; t < p.n; ++t) {
                size_t j = p.member[t];
                if (!p.some[j]) {
                    proj->col[j].hi = 0;
                    fixed++;
                } else if (p.all[j]) {
                    proj->col[j].lo = 1;
                    fixed++;
                }
            }
        }
    free(sink);
    free(source);
    free(next);
    free(member_head);
    free(ec);
    free(count);
    free(bad);
    free(p.parent);
    free(p.head);
    free(p.e);
    free(p.visited);
    free(p.chosen);
    free(p.some);
    free(p.all);
    free(p.member);
    return fixed;
}

static void tighten_dominance_caps(Projection *p) {
    const LabModel *m = p->model;
    unsigned char *min_bad = lab_alloc(m->rows, 1);
    unsigned char *max_bad = lab_alloc(m->rows, 1);
    int64_t *minimum = lab_alloc(m->rows, sizeof *minimum);
    int64_t *maximum = lab_alloc(m->rows, sizeof *maximum);
    for (size_t j = 0; j < m->cols; ++j) {
        Column *c = &p->col[j];
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            int64_t a = p->coefficient[k];
            size_t i = m->row_index[k];
            if (!a)
                continue;
            int64_t term = 0;
            if (!(a > 0 ? c->has_lo : c->has_hi) ||
                !mul64_ok(a, a > 0 ? c->lo : c->hi, &term) ||
                !add64_ok(minimum[i], term, &minimum[i]))
                min_bad[i] = 1;
            if (!(a > 0 ? c->has_hi : c->has_lo) ||
                !mul64_ok(a, a > 0 ? c->hi : c->lo, &term) ||
                !add64_ok(maximum[i], term, &maximum[i]))
                max_bad[i] = 1;
        }
    }
    for (size_t j = 0; j < m->cols; ++j) {
        Column *c = &p->col[j];
        if (!c->has_lo || c->has_hi || c->cost < 0)
            continue;
        int valid = 1;
        int64_t cap = c->lo;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            int64_t a = p->coefficient[k], delta = 0;
            if (!a || (a > 0 ? !p->row[i].has_lo : !p->row[i].has_hi))
                continue;
            if ((a > 0 ? min_bad[i] : max_bad[i]) ||
                !sub64_ok(a > 0 ? p->row[i].lo : p->row[i].hi,
                          a > 0 ? minimum[i] : maximum[i], &delta) ||
                (delta == INT64_MIN && a == -1)) {
                valid = 0;
                break;
            }
            int64_t inc = delta / a + (delta % a && ((delta > 0) == (a > 0)));
            int64_t cand = 0;
            if (!add64_ok(c->lo, inc, &cand) || cand > MAX_SCALAR) {
                valid = 0;
                break;
            }
            if (cand > cap)
                cap = cand;
        }
        if (valid) {
            c->hi = cap;
            c->has_hi = 1;
            c->capped = 1;
        }
    }
    free(minimum);
    free(maximum);
    free(min_bad);
    free(max_bad);
}

static int assemble(Projection *p, mpq_t q) {
    const LabModel *m = p->model;
    p->obj_scale = 1;
    p->obj_gcd = 0;
    for (size_t i = 0; i < m->rows; ++i) {
        Row *r = &p->row[i];
        r->scale = 1;
        r->has_lo = isfinite(m->row_lower[i]);
        r->has_hi = isfinite(m->row_upper[i]);
        if ((r->has_lo && !include_scale(&r->scale, q, lab_source_exact_row_bound(m, i, 0, q))) ||
            (r->has_hi && !include_scale(&r->scale, q, lab_source_exact_row_bound(m, i, 1, q))))
            return 0;
    }
    for (size_t j = 0; j < m->cols; ++j) {
        if ((j % 256 == 0 && lab_now() >= p->deadline) ||
            !include_scale(&p->obj_scale, q, lab_source_exact_objective(m, j, q)))
            return 0;
        Column *c = &p->col[j];
        c->parent = j;
        c->has_lo = isfinite(m->col_lower[j]);
        c->has_hi = isfinite(m->col_upper[j]);
        if ((c->has_lo &&
             !exact_units(q, lab_source_exact_column_bound(m, j, 0, q), 1, &c->lo)) ||
            (c->has_hi &&
             !exact_units(q, lab_source_exact_column_bound(m, j, 1, q), 1, &c->hi)))
            return 0;
        if ((c->has_lo && llabs(c->lo) > MAX_SCALAR) || (c->has_hi && llabs(c->hi) > MAX_SCALAR) ||
            (c->has_lo && c->has_hi && c->lo > c->hi))
            return 0;
        int orig_fixed = c->has_lo && c->has_hi && c->lo == c->hi;
        if (p->presolve_lower && isfinite(p->presolve_lower[j]) &&
            trunc(p->presolve_lower[j]) == p->presolve_lower[j] &&
            fabs(p->presolve_lower[j]) <= MAX_SCALAR) {
            int64_t pl = (int64_t)p->presolve_lower[j];
            if (!c->has_lo || pl > c->lo) {
                c->lo = pl;
                c->has_lo = 1;
            }
        }
        if (p->presolve_upper && isfinite(p->presolve_upper[j]) &&
            trunc(p->presolve_upper[j]) == p->presolve_upper[j] &&
            fabs(p->presolve_upper[j]) <= MAX_SCALAR) {
            int64_t pu = (int64_t)p->presolve_upper[j];
            if (!c->has_hi || pu < c->hi) {
                c->hi = pu;
                c->has_hi = 1;
            }
        }
        if (c->has_lo && c->has_hi && c->lo > c->hi)
            return 0;
        if (!orig_fixed && c->has_lo && c->has_hi && c->lo == c->hi)
            ++p->presolve_fixed;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k)
            if (!include_scale(&p->row[m->row_index[k]].scale, q,
                               lab_source_exact_coefficient(m, k, q)))
                return 0;
    }
    for (size_t j = 0; j < m->cols; ++j) {
        if (j % 256 == 0 && lab_now() >= p->deadline)
            return 0;
        Column *c = &p->col[j];
        if (!exact_units(q, lab_source_exact_objective(m, j, q), p->obj_scale, &c->cost))
            return 0;
        if (m->maximize)
            c->cost = -c->cost;
        p->obj_gcd = gcd64(p->obj_gcd, c->cost);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            Row *r = &p->row[i];
            int64_t a = 0;
            if (!exact_units(q, lab_source_exact_coefficient(m, k, q), r->scale, &a))
                return 0;
            p->coefficient[k] = a;
            if (!a)
                continue;
            if (!r->entries) {
                r->first = j;
                r->first_value = a;
            } else if (r->entries == 1) {
                r->second = j;
                r->equal_columns = a == -r->first_value;
            }
            ++r->entries;
            r->divisor = gcd64(r->divisor, a);
            if (llabs(a) > r->largest)
                r->largest = llabs(a);
        }
    }
    if (!p->obj_gcd)
        p->obj_gcd = 1;
    for (size_t i = 0; i < m->rows; ++i) {
        Row *r = &p->row[i];
        if ((r->has_lo && !exact_units(q, lab_source_exact_row_bound(m, i, 0, q), r->scale, &r->lo)) ||
            (r->has_hi && !exact_units(q, lab_source_exact_row_bound(m, i, 1, q), r->scale, &r->hi)))
            return 0;
        r->divisor = gcd64(r->divisor, r->has_lo ? r->lo : 0);
        r->divisor = gcd64(r->divisor, r->has_hi ? r->hi : 0);
        if (!r->divisor)
            r->divisor = 1;
        r->eq = r->has_lo && r->has_hi && r->lo == r->hi;
        r->equal_columns = r->equal_columns && r->entries == 2 && r->eq && !r->lo;
        int64_t peak = r->largest / r->divisor;
        while (peak > 2048) {
            peak = (peak + 1) / 2;
            ++r->shift;
        }
    }
    p->path_fixed = prune_flow_paths(p);
    if (p->path_fixed)
        tighten_dominance_caps(p);
    for (size_t i = 0; i < m->rows; ++i) {
        Row *r = &p->row[i];
        if (!r->equal_columns)
            continue;
        size_t a = representative(p->col, r->first), b = representative(p->col, r->second);
        if (a != b) {
            p->col[b].parent = a;
            ++p->merged;
        }
    }
    if (p->merged < m->cols / 16 && p->path_fixed < m->cols / 16 &&
        p->presolve_fixed < m->cols / 2)
        return 0;
    for (size_t j = 0; j < m->cols; ++j) {
        Column *c = &p->col[j], *rep = &p->col[representative(p->col, j)];
        if (!rep->group) {
            rep->group = ++p->groups;
            Group *g = &p->group[rep->group - 1];
            g->has_lo = c->has_lo;
            g->has_hi = c->has_hi;
            g->lo = c->lo;
            g->hi = c->hi;
        }
        c->group = rep->group;
        Group *g = &p->group[c->group - 1];
        if (c->has_lo) {
            if (!g->has_lo || c->lo > g->lo)
                g->lo = c->lo;
            g->has_lo = 1;
        }
        if (c->has_hi) {
            if (!g->has_hi || c->hi < g->hi)
                g->hi = c->hi;
            g->has_hi = 1;
        }
        if (g->has_lo && g->has_hi && g->lo > g->hi)
            return 0;
        int64_t unit_cost = c->cost / p->obj_gcd;
        if (!add64_ok(g->cost, unit_cost, &g->cost) || llabs(g->cost) > MAX_EXACT)
            return 0;
        c->next = g->head;
        g->head = j + 1;
    }
    return 1;
}

static int build_columns(Projection *p, SCIP *s, SCIP_CONS **rows, int64_t *sum) {
    const LabModel *m = p->model;
    size_t *stamp = lab_alloc(m->rows, sizeof *stamp);
    size_t *touched = lab_alloc(m->rows, sizeof *touched);
    int ok = 0;
    for (size_t t = 0; t < p->groups; ++t) {
        if ((t & 127U) == 0 && lab_now() >= p->deadline)
            goto done;
        Group *g = &p->group[t];
        size_t count = 0;
        for (size_t at = g->head; at; at = p->col[at - 1].next) {
            size_t j = at - 1;
            for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
                size_t i = m->row_index[k];
                if (!rows[i] || !p->coefficient[k])
                    continue;
                if (stamp[i] != t + 1) {
                    stamp[i] = t + 1;
                    touched[count++] = i;
                    sum[i] = 0;
                }
                int64_t term = p->coefficient[k] / p->row[i].divisor;
                if (!add64_ok(sum[i], term, &sum[i]) || llabs(sum[i]) > MAX_EXACT)
                    goto done;
            }
        }
        for (size_t at = 0; at < count; ++at) {
            size_t i = touched[at];
            if (!sum[i])
                continue;
            double coef = ldexp((double)sum[i], -p->row[i].shift);
            if (!isfinite(coef) || SCIPaddCoefLinear(s, rows[i], g->var, coef) != SCIP_OKAY)
                goto done;
        }
    }
    ok = 1;
done:
    free(touched);
    free(stamp);
    return ok;
}

/* Validate the expanded integer witness against every source side, including
 * the deleted equations and original bounds. Only copy after all checks succeed. */
static int recover(const Projection *p, SCIP *s, SCIP_SOL *sol, int64_t *sum, double *point,
                   int64_t *out_cost) {
    const LabModel *m = p->model;
    for (size_t t = 0; t < p->groups; ++t) {
        Group *g = &p->group[t];
        double x = SCIPgetSolVal(s, sol, g->var);
        if (!isfinite(x) || fabs(x) > MAX_SCALAR || fabs(x - round(x)) > 1e-6)
            return 0;
        g->chosen = (int64_t)round(x);
        if ((g->has_lo && g->chosen < g->lo) || (g->has_hi && g->chosen > g->hi))
            return 0;
    }
    memset(sum, 0, m->rows * sizeof *sum);
    for (size_t j = 0; j < m->cols; ++j) {
        int64_t x = p->group[p->col[j].group - 1].chosen;
        if ((isfinite(m->col_lower[j]) && x < (int64_t)ceil(m->col_lower[j])) ||
            (isfinite(m->col_upper[j]) && x > (int64_t)floor(m->col_upper[j])))
            return 0;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            int64_t term = 0;
            if (!mul64_ok(x, p->coefficient[k], &term) || !add64_ok(sum[i], term, &sum[i]))
                return 0;
        }
    }
    for (size_t i = 0; i < m->rows; ++i)
        if ((p->row[i].has_lo && sum[i] < p->row[i].lo) ||
            (p->row[i].has_hi && sum[i] > p->row[i].hi))
            return 0;
    int64_t total_cost = 0;
    for (size_t t = 0; t < p->groups; ++t) {
        int64_t term = 0;
        if (!mul64_ok(p->group[t].cost, p->group[t].chosen, &term) ||
            !add64_ok(total_cost, term, &total_cost))
            return 0;
    }
    if (lab_now() >= p->deadline)
        return 0;
    for (size_t j = 0; j < m->cols; ++j)
        point[j] = (double)p->group[p->col[j].group - 1].chosen;
    *out_cost = total_cost;
    return 1;
}

static int search(Projection *p, size_t nodes, size_t iterations, double *point,
                  LabIntegerMipStats *stats) {
    const LabModel *m = p->model;
    SCIP *s = NULL;
    SCIP_CONS **rows = lab_alloc(m->rows, sizeof *rows);
    int64_t *sum = lab_alloc(m->rows, sizeof *sum);
    int result = 0;
    if (SCIPcreate(&s) != SCIP_OKAY)
        goto done;
    SCIPsetMessagehdlrQuiet(s, TRUE);
    if (SCIPincludeDefaultPlugins(s) != SCIP_OKAY ||
        SCIPcreateProbBasic(s, "integer model") != SCIP_OKAY ||
        SCIPsetIntParam(s, "timing/clocktype", 2) != SCIP_OKAY ||
        SCIPsetIntParam(s, "parallel/maxnthreads", 1) != SCIP_OKAY ||
        /* QSopt_ex cannot temporarily free a row during LP graph conflict analysis. */
        SCIPsetCharParam(s, "conflict/useinflp", 'd') != SCIP_OKAY ||
        SCIPsetCharParam(s, "conflict/useboundlp", 'd') != SCIP_OKAY ||
        SCIPsetBoolParam(s, "conflict/usesb", FALSE) != SCIP_OKAY ||
        SCIPsetIntParam(s, "randomization/randomseedshift", 1) != SCIP_OKAY ||
        SCIPsetLongintParam(s, "limits/totalnodes", (SCIP_Longint)nodes) != SCIP_OKAY ||
        SCIPsetLongintParam(s, "lp/iterlim", (SCIP_Longint)iterations) != SCIP_OKAY ||
        SCIPsetLongintParam(s, "lp/rootiterlim", (SCIP_Longint)iterations) != SCIP_OKAY ||
        SCIPsetRealParam(s, "limits/gap", 0) != SCIP_OKAY ||
        SCIPsetRealParam(s, "limits/absgap", 0) != SCIP_OKAY)
        goto done;
    for (size_t t = 0; t < p->groups; ++t) {
        Group *g = &p->group[t];
        char name[32];
        snprintf(name, sizeof name, "X%zu", t);
        double lo = g->has_lo ? (double)g->lo : -SCIPinfinity(s);
        double hi = g->has_hi ? (double)g->hi : SCIPinfinity(s);
        SCIP_VARTYPE type = g->has_lo && g->has_hi && g->lo >= 0 && g->hi <= 1
                                ? SCIP_VARTYPE_BINARY
                                : SCIP_VARTYPE_INTEGER;
        if (SCIPcreateVarBasic(s, &g->var, name, lo, hi, (double)g->cost, type) != SCIP_OKAY ||
            SCIPaddVar(s, g->var) != SCIP_OKAY)
            goto done;
    }
    for (size_t i = 0; i < m->rows; ++i) {
        Row *r = &p->row[i];
        if (r->equal_columns || (!r->has_hi && !r->has_lo))
            continue;
        double lower = r->has_lo ? ldexp((double)(r->lo / r->divisor), -r->shift) : -SCIPinfinity(s);
        double upper = r->has_hi ? ldexp((double)(r->hi / r->divisor), -r->shift) : SCIPinfinity(s);
        if (SCIPcreateConsBasicLinear(s, &rows[i], m->row_name ? m->row_name[i] : "", 0, NULL, NULL,
                                      lower, upper) != SCIP_OKAY)
            goto done;
    }
    if (!build_columns(p, s, rows, sum))
        goto done;
    for (size_t i = 0; i < m->rows; ++i)
        if (rows[i] && SCIPaddCons(s, rows[i]) != SCIP_OKAY)
            goto done;
    double remaining = p->deadline - lab_now();
    if (remaining <= .01 || SCIPsetRealParam(s, "limits/time", remaining * .95) != SCIP_OKAY)
        goto done;
    SCIP_RETCODE rc = SCIPsolve(s);
    if (rc != SCIP_OKAY)
        goto done;
    stats->nodes = (size_t)SCIPgetNTotalNodes(s);
    stats->lp_iterations = (size_t)SCIPgetNLPIterations(s);
    stats->lp_solves = (size_t)SCIPgetNLPs(s);
    int64_t integer_cost = 0;
    if (SCIPgetBestSol(s) && recover(p, s, SCIPgetBestSol(s), sum, point, &integer_cost)) {
        result = 1;
        double witness = (double)integer_cost;
        double dual = SCIPgetDualbound(s);
        if (SCIPgetStatus(s) == SCIP_STATUS_OPTIMAL && isfinite(dual) &&
            fabs(SCIPgetPrimalbound(s) - witness) <= 1e-5 && dual > witness - 1 + 1e-5 &&
            dual <= witness + 1e-5)
            result = 2;
    }
done:
    if (s) {
        for (size_t i = 0; i < m->rows; ++i)
            if (rows[i])
                (void)SCIPreleaseCons(s, &rows[i]);
        for (size_t t = 0; t < p->groups; ++t)
            if (p->group[t].var)
                (void)SCIPreleaseVar(s, &p->group[t].var);
        if (SCIPfree(&s) != SCIP_OKAY)
            result = 0;
    }
    free(rows);
    free(sum);
    return result;
}

LabSearchResult lab_integer_mip(const LabModel *m, const double *lower, const double *upper, double deadline,
                    size_t nodes, size_t iterations, double *point, LabIntegerMipStats *stats) {
    LabSearchResult r = {.route = "contracted_scip"};
    *stats = (LabIntegerMipStats){0};
    if (!m->integer || !nodes || !iterations || nodes > LLONG_MAX || iterations > LLONG_MAX ||
        !m->rows || m->cols < 32 || m->rows > ROW_CAP || m->cols > COL_CAP || m->nnz > ENTRY_CAP ||
        lab_now() >= deadline)
        return r;
    int general = 0;
    for (size_t j = 0; j < m->cols; ++j) {
        if (!m->integer[j])
            return r;
        if (!isfinite(m->col_lower[j]) || !isfinite(m->col_upper[j]) || m->col_lower[j] < 0 ||
            m->col_upper[j] > 1)
            general = 1;
    }
    if (!general)
        return r;
    Projection p = {.model = m, .presolve_lower = lower, .presolve_upper = upper,
                    .deadline = deadline};
    p.row = lab_alloc(m->rows, sizeof *p.row);
    p.col = lab_alloc(m->cols, sizeof *p.col);
    p.group = lab_alloc(m->cols, sizeof *p.group);
    p.coefficient = lab_alloc(m->nnz, sizeof *p.coefficient);
    mpq_t q;
    mpq_init(q);
    int result = 0;
    if (assemble(&p, q))
        result = search(&p, nodes, iterations, point, stats);
    mpq_clear(q);
    free(p.coefficient);
    free(p.group);
    free(p.col);
    free(p.row);
    if (result > 0)
        r.point = point;
    if (result == 2) {
        r.status = LAB_MIP_OPTIMAL;
        r.closure = LAB_CLOSE_INTEGER_LATTICE;
    }
    return r;
}

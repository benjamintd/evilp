#include "binary_search.h"
#include "common.h"
#include "exact.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Meet-in-the-middle for small multidimensional 0-1 programs. A nonbinary
 * singleton in an equality with coefficient +/-1 may be substituted exactly.
 * Ceil/floor, approximate integer recognition and row relaxation are avoided.
 * These caps bound working memory and keep all integer arithmetic in range. */
enum { MAX_ROWS = 8, MIN_BINARY = 16, MAX_BINARY = 32, MAX_ABS_SCALAR = 1000000, LEAF = 12 };
typedef struct {
    size_t rows, count, variable[MAX_BINARY], slack[MAX_ROWS];
    int coefficient[MAX_BINARY][MAX_ROWS], lo[MAX_ROWS], hi[MAX_ROWS];
    int rhs[MAX_ROWS], sign[MAX_ROWS];
    int64_t cost[MAX_BINARY];
} Problem;

static int small_integer(double x) {
    return fabs(x) <= MAX_ABS_SCALAR && trunc(x) == x;
}
static int source_integer(int status, const mpq_t source, double x) {
    return !status && small_integer(x) && !mpq_cmp_si(source, (long)x, 1);
}
static int prepare(const LabModel *m, Problem *p, mpq_t q) {
    int64_t slack_cost[MAX_ROWS] = {0};
    p->rows = m->rows;
    for (size_t i = 0; i < p->rows; i++) {
        double lo = m->row_lower[i], hi = m->row_upper[i];
        if ((isfinite(lo) && !source_integer(lab_source_exact_row_bound(m, i, 0, q), q, lo)) ||
            (isfinite(hi) && !source_integer(lab_source_exact_row_bound(m, i, 1, q), q, hi)))
            return 0;
        p->lo[i] = isfinite(lo) ? (int)lo : INT_MIN;
        p->hi[i] = isfinite(hi) ? (int)hi : INT_MAX;
        p->slack[i] = SIZE_MAX;
    }
    for (size_t j = 0; j < m->cols; j++) {
        double lo = m->col_lower[j], hi = m->col_upper[j], cost = m->objective[j];
        if (!source_integer(lab_source_exact_objective(m, j, q), q, cost) ||
            (isfinite(lo) && !source_integer(lab_source_exact_column_bound(m, j, 0, q), q, lo)) ||
            (isfinite(hi) && !source_integer(lab_source_exact_column_bound(m, j, 1, q), q, hi)))
            return 0;
        size_t first = m->col_start[j], last = m->col_start[j + 1];
        if (m->integer[j] && lo == 0 && hi == 1) {
            size_t v = p->count;
            if (v == MAX_BINARY)
                return 0;
            p->variable[v] = j;
            p->cost[v] = (int64_t)cost;
            for (size_t k = first; k < last; k++) {
                size_t i = m->row_index[k];
                if (p->coefficient[v][i] ||
                    !source_integer(lab_source_exact_coefficient(m, k, q), q, m->value[k]))
                    return 0;
                p->coefficient[v][i] = (int)m->value[k];
            }
            p->count++;
        } else {
            if (last - first != 1)
                return 0;
            size_t i = m->row_index[first];
            double a = m->value[first];
            if (p->slack[i] != SIZE_MAX || p->lo[i] != p->hi[i] || fabs(a) != 1 ||
                !source_integer(lab_source_exact_coefficient(m, first, q), q, a))
                return 0;
            int b = p->lo[i];
            p->slack[i] = j;
            p->rhs[i] = b;
            p->sign[i] = (int)a;
            slack_cost[i] = (int64_t)a * (int64_t)cost;
            /* lo <= a*(b - activity) <= hi */
            p->lo[i] = a > 0 ? (isfinite(hi) ? b - (int)hi : INT_MIN)
                             : (isfinite(lo) ? b + (int)lo : INT_MIN);
            p->hi[i] = a > 0 ? (isfinite(lo) ? b - (int)lo : INT_MAX)
                             : (isfinite(hi) ? b + (int)hi : INT_MAX);
        }
    }
    for (size_t v = 0; v < p->count; v++) {
        for (size_t i = 0; i < p->rows; i++)
            p->cost[v] -= slack_cost[i] * p->coefficient[v][i];
        if (m->maximize)
            p->cost[v] = -p->cost[v];
    }
    return p->count >= MIN_BINARY;
}

typedef struct {
    int activity[MAX_ROWS];
    uint32_t mask;
    int64_t cost;
} Assignment;
typedef struct {
    int lo[MAX_ROWS], hi[MAX_ROWS];
    size_t begin, end, left, right;
    Assignment best;
} Box;
typedef struct {
    const Problem *p;
    Assignment *right;
    Box *box;
    size_t used;
    int64_t best;
    uint32_t mask;
    size_t split;
    double deadline;
    int expired;
    size_t visits;
} Search;

static void add_variable(const Problem *p, Assignment *a, size_t v, int sign) {
    for (size_t i = 0; i < p->rows; i++)
        a->activity[i] += sign * p->coefficient[v][i];
    a->cost += sign * p->cost[v];
}
/* Bisect the widest coordinate at its midpoint. Each non-leaf has two
 * nonempty children. Identical activities need only the cheapest assignment. */
static size_t build_tree(Search *s, size_t begin, size_t end) {
    size_t id = ++s->used;
    Box *box = &s->box[id];
    box->begin = begin;
    box->end = end;
    box->best = s->right[begin];
    memcpy(box->lo, box->best.activity, sizeof box->lo);
    memcpy(box->hi, box->lo, sizeof box->hi);
    for (size_t k = begin + 1; k < end; k++) {
        const Assignment *a = &s->right[k];
        if (a->cost < box->best.cost)
            box->best = *a;
        for (size_t i = 0; i < s->p->rows; i++) {
            if (a->activity[i] < box->lo[i])
                box->lo[i] = a->activity[i];
            if (a->activity[i] > box->hi[i])
                box->hi[i] = a->activity[i];
        }
    }
    size_t axis = 0;
    for (size_t i = 1; i < s->p->rows; i++)
        if (box->hi[i] - box->lo[i] > box->hi[axis] - box->lo[axis])
            axis = i;
    if (end - begin > LEAF && box->lo[axis] < box->hi[axis]) {
        if (lab_now() >= s->deadline) {
            s->expired = 1;
            return id;
        }
        int mid = box->lo[axis] + (box->hi[axis] - box->lo[axis]) / 2;
        size_t split = begin;
        for (size_t k = begin; k < end; k++)
            if (s->right[k].activity[axis] <= mid) {
                Assignment t = s->right[k];
                s->right[k] = s->right[split];
                s->right[split++] = t;
            }
        box->left = build_tree(s, begin, split);
        if (!s->expired)
            box->right = build_tree(s, split, end);
    }
    return id;
}
static int fits(const Problem *p, const Assignment *left, const Assignment *right) {
    for (size_t i = 0; i < p->rows; i++) {
        int sum = left->activity[i] + right->activity[i];
        if (sum < p->lo[i] || sum > p->hi[i])
            return 0;
    }
    return 1;
}
static void keep(Search *s, const Assignment *left, const Assignment *right) {
    s->best = left->cost + right->cost;
    s->mask = left->mask | (right->mask << s->split);
}
static void query(Search *s, size_t id, const Assignment *left) {
    const Box *box = &s->box[id];
    if (s->expired || left->cost + box->best.cost >= s->best)
        return;
    if ((++s->visits & 1023U) == 0 && lab_now() >= s->deadline) {
        s->expired = 1;
        return;
    }
    for (size_t i = 0; i < s->p->rows; i++)
        if (left->activity[i] + box->lo[i] > s->p->hi[i] ||
            left->activity[i] + box->hi[i] < s->p->lo[i])
            return;
    if (fits(s->p, left, &box->best)) {
        keep(s, left, &box->best);
        return;
    }
    if (box->left) {
        size_t a = box->left, b = box->right;
        if (s->box[b].best.cost < s->box[a].best.cost) {
            a = box->right;
            b = box->left;
        }
        query(s, a, left);
        query(s, b, left);
    } else {
        for (size_t k = box->begin; k < box->end; k++) {
            const Assignment *r = &s->right[k];
            if (left->cost + r->cost < s->best && fits(s->p, left, r))
                keep(s, left, r);
        }
    }
}
int lab_binary_search(const LabModel *m, double deadline, size_t query_budget, double *point,
                      size_t *queries) {
    *queries = 0;
    if (!m->rows || m->rows > MAX_ROWS || m->cols < MIN_BINARY || m->cols > MAX_BINARY + MAX_ROWS ||
        !m->integer || lab_now() >= deadline)
        return 0;
    Problem p = {0};
    mpq_t q;
    mpq_init(q);
    int valid = prepare(m, &p, q);
    mpq_clear(q);
    if (!valid)
        return 0;
    Search s = {.p = &p, .split = p.count / 2, .best = INT64_MAX, .deadline = deadline};
    uint32_t left_count = UINT32_C(1) << s.split;
    uint32_t right_count = UINT32_C(1) << (p.count - s.split);
    if (left_count > query_budget)
        return 0;
    s.right = lab_alloc(right_count, sizeof *s.right);
    s.box = lab_alloc(2U * right_count, sizeof *s.box);
    for (uint32_t base = 1, v = 0; base < right_count; base *= 2, v++) {
        for (uint32_t k = 0; k < base; k++) {
            Assignment *a = &s.right[base + k];
            *a = s.right[k];
            a->mask |= base;
            add_variable(&p, a, s.split + v, 1);
        }
    }
    size_t root = build_tree(&s, 0, right_count);
    Assignment left = {0};
    for (uint32_t k = 0; k < left_count && !s.expired; k++) {
        if ((k & 255U) == 0 && lab_now() >= deadline) {
            s.expired = 1;
            break;
        }
        if (k) {
            uint32_t bit = k & (0U - k);
            size_t v = 0;
            for (uint32_t x = bit; x > 1; x >>= 1)
                v++;
            left.mask ^= bit;
            add_variable(&p, &left, v, (left.mask & bit) ? 1 : -1);
        }
        query(&s, root, &left);
        (*queries)++;
    }
    int found = !s.expired && s.best != INT64_MAX && lab_now() < deadline;
    if (found) {
        Assignment sum = {0};
        for (size_t v = 0; v < p.count; v++) {
            point[p.variable[v]] = (s.mask >> v) & 1U;
            if (point[p.variable[v]] != 0)
                add_variable(&p, &sum, v, 1);
        }
        for (size_t i = 0; i < p.rows; i++)
            if (p.slack[i] != SIZE_MAX)
                point[p.slack[i]] = p.sign[i] * (p.rhs[i] - sum.activity[i]);
    }
    free(s.box);
    free(s.right);
    return found;
}

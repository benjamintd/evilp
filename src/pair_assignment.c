#include "pair_assignment.h"
#include "common.h"
#include "exact.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { MAX_COVERS = 32, MAX_CHOICES = 256, MAX_SELECTED = MAX_COVERS / 2 };

typedef struct {
    uint32_t covers;
    int64_t cost;
} Choice;

typedef struct {
    int owner, last_aux;
    uint32_t common_cover;
    uint32_t cover;
} Row;

typedef struct {
    int first, second;
} Endpoints;

typedef struct {
    Choice choices[MAX_CHOICES];
    int by_cover[MAX_COVERS][MAX_CHOICES], count[MAX_COVERS];
    int selected[MAX_SELECTED], best[MAX_SELECTED];
    int64_t pair[MAX_CHOICES][MAX_CHOICES], minimum;
    int nchoices, ncovers;
    size_t nodes, limit;
    double deadline;
} Search;

/* Every non-cover row must have the form sum(y) = x. A y appearing
 * under two different owners can be positive only if both are selected.
 * Conversely, if all other endpoints in a row share one cover equation,
 * at most one can be selected. When both owners are selected that y is
 * forced to one. Distinct terms per row and unique owner pairs matter. */
static int recognize(const LabModel *m, Row *rows, Endpoints *ends, Search *s, mpq_t q) {
    for (size_t r = 0; r < m->rows; ++r) {
        double rhs = m->row_upper[r];
        if ((rhs != 0 && rhs != 1) || m->row_lower[r] != rhs ||
            lab_source_exact_row_bound(m, r, 0, q) || mpq_cmp_ui(q, (unsigned long)rhs, 1) ||
            lab_source_exact_row_bound(m, r, 1, q) || mpq_cmp_ui(q, (unsigned long)rhs, 1))
            return 0;
        if (rhs) {
            if (s->ncovers == MAX_COVERS)
                return 0;
            rows[r].cover = UINT32_C(1) << s->ncovers++;
        }
        rows[r].owner = rows[r].last_aux = -1;
        rows[r].common_cover = UINT32_MAX;
    }
    if (s->ncovers < 4 || s->ncovers % 2)
        return 0;

    for (size_t j = 0; j < m->cols; ++j) {
        if ((j % 128 == 0 && lab_now() >= s->deadline) || m->col_lower[j] != 0 ||
            lab_source_exact_column_bound(m, j, 0, q) || mpq_sgn(q) ||
            (m->col_upper[j] != INFINITY &&
             (m->col_upper[j] < 1 || lab_source_exact_column_bound(m, j, 1, q) ||
              mpq_cmp_ui(q, 1, 1) < 0)))
            return 0;
        uint32_t mask = 0;
        int nc = 0;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            Row *row = &rows[m->row_index[k]];
            double a = m->value[k];
            if ((a != 1 && a != -1) || lab_source_exact_coefficient(m, k, q) ||
                mpq_cmp_si(q, (long)a, 1))
                return 0;
            if (row->cover) {
                if (a != 1 || (mask & row->cover) || ++nc > 2)
                    return 0;
                mask |= row->cover;
            }
        }
        ends[j] = (Endpoints){-1, -1};
        if (!mask)
            continue;
        if (nc != 2 || !m->integer[j] || s->nchoices == MAX_CHOICES)
            return 0;
        int i = s->nchoices++;
        s->choices[i] = (Choice){.covers = mask};
        ends[j].first = i;
        for (int r = 0; r < s->ncovers; ++r)
            if (mask & (UINT32_C(1) << r))
                s->by_cover[r][s->count[r]++] = i;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            Row *row = &rows[m->row_index[k]];
            if (!row->cover) {
                if (m->value[k] != -1 || row->owner >= 0)
                    return 0;
                row->owner = i;
            }
        }
    }
    if ((size_t)s->nchoices == m->cols)
        return 0;
    for (size_t r = 0; r < m->rows; ++r)
        if (!rows[r].cover && rows[r].owner < 0)
            return 0;

    /* INT64_MAX temporarily marks pairs with no product column. */
    for (int i = 0; i < s->nchoices; ++i)
        for (int j = 0; j < s->nchoices; ++j)
            s->pair[i][j] = INT64_MAX;
    int64_t min_pair = 0, min_linear = 0;
    for (size_t j = 0; j < m->cols; ++j) {
        if ((j % 128 == 0 && lab_now() >= s->deadline) || lab_source_exact_objective(m, j, q) ||
            mpz_cmp_ui(mpq_denref(q), 1) || mpq_cmp_si(q, -1000000000, 1) < 0 ||
            mpq_cmp_si(q, 1000000000, 1) > 0)
            return 0;
        int64_t cost = mpz_get_si(mpq_numref(q));
        if (m->maximize)
            cost = -cost;
        if (ends[j].first >= 0) {
            s->choices[ends[j].first].cost = cost;
            if (cost < min_linear)
                min_linear = cost;
            continue;
        }
        int a = -1, b = -1;
        size_t begin = m->col_start[j], end = m->col_start[j + 1];
        /* Detect repeated entries without assuming CSC row order. */
        for (size_t k = begin; k < end; ++k) {
            Row *row = &rows[m->row_index[k]];
            if (row->last_aux == (int)j)
                return 0;
            row->last_aux = (int)j;
            int owner = row->owner;
            if (m->value[k] != 1 || owner < 0)
                return 0;
            if (a < 0)
                a = owner;
            else if (owner != a && b < 0)
                b = owner;
            else if (owner != a && owner != b)
                return 0;
        }
        if (b < 0 || s->pair[a][b] != INT64_MAX)
            return 0;
        s->pair[a][b] = s->pair[b][a] = cost;
        ends[j] = (Endpoints){a, b};
        if (cost < min_pair)
            min_pair = cost;
        for (size_t k = begin; k < end; ++k) {
            Row *row = &rows[m->row_index[k]];
            int other = row->owner == a ? b : a;
            row->common_cover &= s->choices[other].covers;
            if (!row->common_cover)
                return 0;
        }
    }
    /* All full covers select ncovers/2 columns, so shifting all linear
     * and all pair costs adds the same constant to every leaf. With
     * nonnegative increments the partial cost is a valid lower bound. */
    for (int i = 0; i < s->nchoices; ++i) {
        s->choices[i].cost -= min_linear;
        for (int j = 0; j < s->nchoices; ++j)
            s->pair[i][j] = (s->pair[i][j] == INT64_MAX ? 0 : s->pair[i][j]) - min_pair;
    }
    return 1;
}

static int enumerate(Search *s, uint32_t open, int depth, int64_t cost) {
    if (s->nodes >= s->limit || (!(s->nodes % 256) && lab_now() >= s->deadline))
        return 0;
    ++s->nodes;
    if (cost >= s->minimum)
        return 1;
    if (!open) {
        s->minimum = cost;
        memcpy(s->best, s->selected, (size_t)depth * sizeof *s->best);
        return 1;
    }
    int r = 0;
    while (!(open & (UINT32_C(1) << r)))
        ++r;
    for (int k = 0; k < s->count[r]; ++k) {
        int i = s->by_cover[r][k];
        Choice c = s->choices[i];
        if ((open & c.covers) != c.covers)
            continue;
        int64_t next = cost + c.cost;
        for (int t = 0; t < depth; ++t)
            next += s->pair[i][s->selected[t]];
        s->selected[depth] = i;
        if (!enumerate(s, open ^ c.covers, depth + 1, next))
            return 0;
    }
    return 1;
}

int lab_pair_assignment(const LabModel *m, double deadline, size_t limit, double *point,
                        size_t *nodes) {
    *nodes = 0;
    if (!limit || !m->integer || m->rows < 4 || m->cols < 3 || m->cols > 32768 || m->rows > 32768 ||
        m->nnz > 262144 || lab_now() >= deadline)
        return 0;
    for (size_t r = 0; r < m->rows; ++r)
        if (m->row_lower[r] != m->row_upper[r] || (m->row_upper[r] != 1 && m->row_upper[r] != 0))
            return 0;
    Search *s = lab_alloc(1, sizeof *s);
    s->limit = limit;
    s->deadline = deadline;
    s->minimum = INT64_MAX;
    Row *rows = lab_alloc(m->rows, sizeof *rows);
    Endpoints *ends = lab_alloc(m->cols, sizeof *ends);
    mpq_t q;
    mpq_init(q);
    int ok = 0;
    if (!recognize(m, rows, ends, s, q))
        goto done;
    uint32_t open = UINT32_MAX >> (MAX_COVERS - s->ncovers);
    if (!enumerate(s, open, 0, 0) || s->minimum == INT64_MAX)
        goto done;

    /* A missing product may make a matching infeasible. Only close when
     * the minimum over the relaxation satisfies every original equality. */
    unsigned char selected[MAX_CHOICES] = {0};
    for (int k = 0; k < s->ncovers / 2; ++k)
        selected[s->best[k]] = 1;
    int *activity = lab_alloc(m->rows, sizeof *activity);
    for (size_t j = 0; j < m->cols; ++j) {
        Endpoints e = ends[j];
        point[j] = selected[e.first] && (e.second < 0 || selected[e.second]);
        if (point[j])
            for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k)
                activity[m->row_index[k]] += (int)m->value[k];
    }
    ok = lab_now() < deadline;
    for (size_t r = 0; r < m->rows && ok; ++r)
        ok = activity[r] == (rows[r].cover != 0);
    free(activity);
done:
    *nodes = s->nodes;
    mpq_clear(q);
    free(ends);
    free(rows);
    free(s);
    return ok;
}

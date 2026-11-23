/* Source-exact quotient of a bipartite perfect matching. Each nonempty
 * class pair is complete bipartite with identical costs and side columns.
 * Its integer count has a bijective lift up to vertex permutations. */
#include "scip_mip.h"
enum { QUOTIENT_ROW_LIMIT = 16 };
typedef struct {
    int count, representative[MAX_SIDE], size[MAX_SIDE], of[MAX_SIDE];
} VertexClasses;
typedef struct {
    int left, right;
    size_t original;
} QuotientColumn;
typedef struct {
    size_t row, col;
    mpq_t value;
} QuotientEntry;
typedef struct {
    int equality;
    mpq_t side;
} QuotientRow;
typedef struct {
    const Assignment *graph;
    VertexClasses left, right;
    QuotientColumn *column;
    size_t *index, cols, rows, nz;
    QuotientRow *row;
    QuotientEntry *entries;
    LabModel reduced;
} Quotient;

static int equivalent_edges(const Assignment *g, int a, int b, mpq_t x, mpq_t y) {
    if (a == b)
        return 1;
    if (!a || !b)
        return 0;
    const LabModel *m = g->m;
    size_t j = (size_t)a - 1, h = (size_t)b - 1;
    if (m->objective[j] != m->objective[h] || lab_source_exact_objective(m, j, x) ||
        lab_source_exact_objective(m, h, y) || !mpq_equal(x, y))
        return 0;
    size_t p = m->col_start[j], q = m->col_start[h];
    for (;;) {
        while (p < m->col_start[j + 1] && g->cover[m->row_index[p]])
            ++p;
        while (q < m->col_start[h + 1] && g->cover[m->row_index[q]])
            ++q;
        if (p == m->col_start[j + 1] || q == m->col_start[h + 1])
            return p == m->col_start[j + 1] && q == m->col_start[h + 1];
        if (m->row_index[p] != m->row_index[q] || m->value[p] != m->value[q] ||
            lab_source_exact_coefficient(m, p, x) || lab_source_exact_coefficient(m, q, y) ||
            !mpq_equal(x, y))
            return 0;
        ++p;
        ++q;
    }
}

static int group_vertices(const Assignment *g, int right, VertexClasses *classes) {
    mpq_t x, y;
    mpq_inits(x, y, NULL);
    int complete = 1;
    for (int i = 0; i < g->n && complete; ++i) {
        int group;
        for (group = 0; group < classes->count; ++group) {
            if (lab_now() >= g->deadline) {
                complete = 0;
                break;
            }
            int rep = classes->representative[group], k;
            for (k = 0; k < g->n; ++k) {
                int a = g->edge[right ? k * g->n + i : i * g->n + k];
                int b = g->edge[right ? k * g->n + rep : rep * g->n + k];
                if (!equivalent_edges(g, a, b, x, y))
                    break;
            }
            if (k == g->n)
                break;
        }
        if (group == classes->count)
            classes->representative[classes->count++] = i;
        classes->of[i] = group;
        ++classes->size[group];
    }
    mpq_clears(x, y, NULL);
    return complete;
}

static void quotient_entry(Quotient *q, size_t row, size_t col, const mpq_t value) {
    if (mpq_sgn(value)) {
        QuotientEntry *e = &q->entries[q->nz++];
        e->row = row;
        e->col = col;
        mpq_init(e->value);
        mpq_set(e->value, value);
    }
}

static void quotient_quota(Quotient *q, size_t row, unsigned quota) {
    q->row[row].equality = 1;
    mpq_init(q->row[row].side);
    mpq_set_ui(q->row[row].side, quota, 1);
}

/* For mutually exclusive unit choices, at most one discrete term contributes.
 * For every such state the continuous activity z has the same valid lower
 * bound b. Replace z >= t by z >= max(b,t) for each state, and omit the row
 * altogether when b satisfies every state. */
static int tighten_choice_row(const Quotient *q, size_t n, const size_t *columns, mpq_t *coef,
                              mpq_t side, const double *lower, const double *upper) {
    int left = -1, right = -1, choices = 0, redundant = 0;
    mpq_t bound, term, value, base;
    mpq_inits(bound, term, value, base, NULL);
    for (size_t k = 0; k < n; ++k) {
        QuotientColumn c = q->column[columns[k]];
        if (c.left >= 0) {
            if (q->left.size[c.left] != 1 && q->right.size[c.right] != 1)
                goto done;
            if (!choices) {
                left = c.left;
                right = c.right;
            } else {
                if (left != c.left)
                    left = -1;
                if (right != c.right)
                    right = -1;
            }
            ++choices;
        } else {
            int high = mpq_sgn(coef[k]) < 0;
            double v = high ? upper[c.original] : lower[c.original];
            if (lab_source_exact_column_value(q->graph->m, c.original, high, v, value))
                goto done;
            mpq_mul(term, coef[k], value);
            mpq_add(bound, bound, term);
        }
    }
    if (choices > 1 && !(left >= 0 && q->left.size[left] == 1) &&
        !(right >= 0 && q->right.size[right] == 1))
        goto done;
    mpq_set(base, mpq_cmp(bound, side) > 0 ? bound : side);
    int restricts_state = mpq_cmp(side, bound) > 0;
    for (size_t k = 0; k < n; ++k)
        if (q->column[columns[k]].left >= 0) {
            mpq_sub(term, side, coef[k]);
            if (mpq_cmp(term, bound) < 0)
                mpq_set(term, bound);
            restricts_state |= mpq_cmp(term, bound) > 0;
            mpq_sub(coef[k], base, term);
        }
    mpq_set(side, base);
    redundant = !restricts_state;
done:
    mpq_clears(bound, term, value, base, NULL);
    return redundant;
}

static int append_side(Quotient *q, size_t i, int high, const double *lower, const double *upper) {
    const Assignment *g = q->graph;
    size_t columns[QUOTIENT_ROW_LIMIT], n = 0;
    mpq_t coef[QUOTIENT_ROW_LIMIT], side;
    mpq_init(side);
    int ok = 0;
    if (lab_source_exact_row_bound(g->m, i, high, side))
        goto done;
    if (high)
        mpq_neg(side, side);
    for (size_t p = g->offset[i]; p < g->offset[i + 1]; ++p) {
        size_t entry = g->entry[p], col = q->index[g->column[entry]];
        if (!col)
            continue;
        if (n == QUOTIENT_ROW_LIMIT)
            goto done;
        columns[n] = col - 1;
        mpq_init(coef[n]);
        ++n;
        if (lab_source_exact_coefficient(g->m, entry, coef[n - 1]))
            goto done;
        if (high)
            mpq_neg(coef[n - 1], coef[n - 1]);
    }
    if (!tighten_choice_row(q, n, columns, coef, side, lower, upper)) {
        size_t row = q->rows++;
        mpq_init(q->row[row].side);
        mpq_set(q->row[row].side, side);
        for (size_t k = 0; k < n; ++k)
            quotient_entry(q, row, columns[k], coef[k]);
    }
    ok = 1;
done:
    while (n)
        mpq_clear(coef[--n]);
    mpq_clear(side);
    return ok;
}

static int build_quotient_rows(Quotient *q, const double *lower, const double *upper) {
    const Assignment *g = q->graph;
    const LabModel *m = g->m;
    if (!group_vertices(g, 0, &q->left) || !group_vertices(g, 1, &q->right) ||
        (q->left.count == g->n && q->right.count == g->n))
        return 0;
    q->column = lab_alloc(m->cols, sizeof *q->column);
    q->index = lab_alloc(m->cols, sizeof *q->index);
    q->row = lab_alloc(2 * m->rows + 2 * (size_t)g->n, sizeof *q->row);
    q->entries = lab_alloc(2 * m->nnz, sizeof *q->entries);
    q->rows = (size_t)(q->left.count + q->right.count);
    for (int a = 0; a < q->left.count; ++a)
        quotient_quota(q, (size_t)a, (unsigned)q->left.size[a]);
    for (int b = 0; b < q->right.count; ++b)
        quotient_quota(q, (size_t)(q->left.count + b), (unsigned)q->right.size[b]);
    mpq_t one;
    mpq_init(one);
    mpq_set_ui(one, 1, 1);
    for (int a = 0; a < q->left.count; ++a)
        for (int b = 0; b < q->right.count; ++b) {
            int edge = g->edge[q->left.representative[a] * g->n + q->right.representative[b]];
            if (!edge)
                continue;
            size_t j = (size_t)edge - 1, col = q->cols++;
            q->column[col] = (QuotientColumn){a, b, j};
            q->index[j] = col + 1;
            quotient_entry(q, (size_t)a, col, one);
            quotient_entry(q, (size_t)(q->left.count + b), col, one);
        }
    mpq_clear(one);
    for (size_t j = 0; j < m->cols; ++j)
        if (!m->integer[j]) {
            q->column[q->cols] = (QuotientColumn){-1, -1, j};
            q->index[j] = ++q->cols;
        }
    for (size_t i = 0; i < m->rows; ++i) {
        if (!(i % 64) && lab_now() >= g->deadline)
            return 0;
        if (g->cover[i])
            continue;
        if ((isfinite(m->row_lower[i]) && !append_side(q, i, 0, lower, upper)) ||
            (isfinite(m->row_upper[i]) && !append_side(q, i, 1, lower, upper)))
            return 0;
    }
    return 1;
}

static int materialize_quotient(Quotient *q, const double *lower, const double *upper) {
    LabModel *m = &q->reduced;
    const LabModel *src = q->graph->m;
    m->rows = q->rows;
    m->cols = q->cols;
    m->nnz = q->nz;
    m->parsed_decimal = 1;
    m->maximize = src->maximize;
    m->row_lower = lab_alloc(m->rows, sizeof(double));
    m->row_upper = lab_alloc(m->rows, sizeof(double));
    m->col_lower = lab_alloc(m->cols, sizeof(double));
    m->col_upper = lab_alloc(m->cols, sizeof(double));
    m->integer = lab_alloc(m->cols, sizeof *m->integer);
    m->objective = lab_alloc(m->cols, sizeof(double));
    m->col_start = lab_alloc(m->cols + 1, sizeof *m->col_start);
    m->row_index = lab_alloc(m->nnz, sizeof *m->row_index);
    m->value = lab_alloc(m->nnz, sizeof(double));
    LabExactModel *s = m->exact = lab_exact_model_new(m->rows, m->cols, m->nnz);
    if (lab_source_exact_offset(src, s->offset) || lab_exact_nearest(s->offset, &m->offset))
        return 0;
    for (size_t i = 0; i < m->rows; ++i) {
        if (!(i % 256) && lab_now() >= q->graph->deadline)
            return 0;
        mpq_set(s->row_lower[i], q->row[i].side);
        if (lab_exact_nearest(s->row_lower[i], &m->row_lower[i]))
            return 0;
        m->row_upper[i] = INFINITY;
        if (q->row[i].equality) {
            mpq_set(s->row_upper[i], s->row_lower[i]);
            m->row_upper[i] = m->row_lower[i];
        }
    }
    for (size_t j = 0; j < m->cols; ++j) {
        if (!(j % 256) && lab_now() >= q->graph->deadline)
            return 0;
        QuotientColumn c = q->column[j];
        if (lab_source_exact_objective(src, c.original, s->objective[j]) ||
            lab_exact_nearest(s->objective[j], &m->objective[j]))
            return 0;
        if (c.left >= 0) {
            m->integer[j] = 1;
            int n = q->left.size[c.left] < q->right.size[c.right] ? q->left.size[c.left]
                                                                  : q->right.size[c.right];
            mpq_set_ui(s->col_upper[j], (unsigned)n, 1);
            m->col_upper[j] = n;
        } else {
            m->col_lower[j] = lower[c.original];
            m->col_upper[j] = upper[c.original];
            if ((isfinite(m->col_lower[j]) &&
                 lab_source_exact_column_value(src, c.original, 0, m->col_lower[j],
                                               s->col_lower[j])) ||
                (isfinite(m->col_upper[j]) &&
                 lab_source_exact_column_value(src, c.original, 1, m->col_upper[j],
                                               s->col_upper[j])))
                return 0;
        }
        s->col_lower_source[j] = m->col_lower[j];
        s->col_upper_source[j] = m->col_upper[j];
    }
    for (size_t k = 0; k < q->nz; ++k)
        ++m->col_start[q->entries[k].col + 1];
    for (size_t j = 0; j < m->cols; ++j)
        m->col_start[j + 1] += m->col_start[j];
    size_t *cursor = lab_alloc(m->cols, sizeof *cursor);
    memcpy(cursor, m->col_start, m->cols * sizeof *cursor);
    int ok = 1;
    for (size_t k = 0; k < q->nz && ok; ++k) {
        if (!(k % 256) && lab_now() >= q->graph->deadline) {
            ok = 0;
            break;
        }
        QuotientEntry *e = &q->entries[k];
        size_t p = cursor[e->col]++;
        m->row_index[p] = e->row;
        mpq_set(s->value[p], e->value);
        ok = !lab_exact_nearest(s->value[p], &m->value[p]) && m->value[p] != 0;
    }
    free(cursor);
    return ok;
}

static int lift_quotient(const Quotient *q, const double *x, const double *lower,
                         const double *upper, double *point, char ***decimal) {
    const Assignment *g = q->graph;
    int left[MAX_SIDE], right[MAX_SIDE], used = 0;
    memcpy(left, q->left.of, (size_t)g->n * sizeof *left);
    memcpy(right, q->right.of, (size_t)g->n * sizeof *right);
    memset(point, 0, g->m->cols * sizeof *point);
    for (size_t j = 0; j < q->cols; ++j) {
        QuotientColumn c = q->column[j];
        if (c.left < 0)
            continue;
        if (!isfinite(x[j]) || x[j] != floor(x[j]) || x[j] < 0 || x[j] > q->reduced.col_upper[j])
            return 0;
        for (int h = 0; h < (int)x[j]; ++h) {
            int a = 0, b = 0;
            while (a < g->n && left[a] != c.left)
                ++a;
            while (b < g->n && right[b] != c.right)
                ++b;
            if (a == g->n || b == g->n || !g->edge[a * g->n + b])
                return 0;
            left[a] = right[b] = -1;
            point[g->edge[a * g->n + b] - 1] = 1;
            ++used;
        }
    }
    if (used != g->n)
        return 0;
    for (size_t j = 0; j < g->m->cols; ++j)
        if (g->m->integer[j] && (point[j] < lower[j] || point[j] > upper[j]))
            return 0;
    int n = 0;
    BoundArc *arcs = lab_alloc(2 * (g->m->rows + g->m->cols), sizeof *arcs);
    int ok = build_arcs(g, point, arcs, &n) && recover(g, arcs, n, point, decimal);
    free(arcs);
    return ok;
}

int lab_assignment_quotient(const LabModel *m, const double *lower, const double *upper,
                            double deadline, size_t node_limit, size_t iteration_limit,
                            double *point, char ***decimal, LabScipMipStats *stats) {
    *stats = (LabScipMipStats){0};
    *decimal = NULL;
    if (!node_limit || !iteration_limit || !m->integer || !m->rows || !m->cols || m->rows > 8192 ||
        m->cols > 16384 || m->nnz > 100000 || lab_now() >= deadline)
        return 0;
    Assignment g = assignment_graph(m, deadline);
    Quotient q = {.graph = &g};
    int ok = 0;
    if (recognize(&g, upper) && g.nc && build_quotient_rows(&q, lower, upper) &&
        materialize_quotient(&q, lower, upper) && lab_now() < deadline) {
        double *x = lab_alloc(q.cols, sizeof *x);
        LabSearchResult result = lab_scip_mip(&q.reduced, q.reduced.col_lower, q.reduced.col_upper, NULL, NULL,
                                  NULL, deadline, node_limit, iteration_limit, x, stats);
        if (result.status == LAB_MIP_OPTIMAL)
            ok = lift_quotient(&q, x, lower, upper, point, decimal);
        free(x);
    }
    for (size_t i = 0; i < q.rows; ++i)
        mpq_clear(q.row[i].side);
    for (size_t k = 0; k < q.nz; ++k)
        mpq_clear(q.entries[k].value);
    free(q.row);
    free(q.entries);
    free(q.column);
    free(q.index);
    lab_model_free(&q.reduced);
    free_assignment_graph(&g);
    return ok;
}

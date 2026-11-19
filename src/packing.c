#include "packing.h"
#include "common.h"
#include "exact.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { MAX_ITEMS = 512, MAX_SIZE = 1000000 };
typedef struct {
    int size, column;
} Option;
typedef struct {
    int items, bins, target;
    int weight[MAX_ITEMS], capacity[MAX_ITEMS], gate[MAX_ITEMS], degree[MAX_ITEMS];
    Option *options;
    uint32_t random;
} Packing;

static uint32_t random_word(Packing *p) {
    uint32_t x = p->random;
    x ^= x << 13;
    x ^= x >> 17;
    return p->random = x ^ (x << 5);
}
static Option *option(const Packing *p, int item, int bin) {
    return &p->options[item * p->bins + bin];
}
static int coefficient(const LabModel *m, size_t k, mpq_t q, int *value) {
    if (lab_source_exact_coefficient(m, k, q) || mpz_cmp_ui(mpq_denref(q), 1) ||
        mpz_cmpabs_ui(mpq_numref(q), MAX_SIZE) > 0)
        return 0;
    *value = (int)mpz_get_si(mpq_numref(q));
    return m->value[k] == (double)*value;
}

/* Item equations equal one. Each resource row has a single paid gate. All
 * other columns consume capacity and may cover multiple item equations.
 * Singleton columns supply a restricted set of options for the heuristic. */
static int read_columns(Packing *p, const LabModel *m, const int *id, double deadline) {
    mpq_t q, charge;
    mpq_inits(q, charge, NULL);
    int ok = 0;
    for (size_t j = 0; j < m->cols; ++j) {
        if (m->col_lower[j] != 0 || m->col_upper[j] != 1 || !m->integer[j] ||
            ((j & 127U) == 0 && lab_now() >= deadline) ||
            lab_source_exact_column_bound(m, j, 0, q) || mpq_sgn(q) ||
            lab_source_exact_column_bound(m, j, 1, q) || mpq_cmp_ui(q, 1, 1) ||
            lab_source_exact_objective(m, j, q))
            goto done;
        size_t begin = m->col_start[j], end = m->col_start[j + 1];
        if (mpq_sgn(q)) {
            if (m->maximize)
                mpq_neg(q, q);
            if (mpq_sgn(q) <= 0 || (mpq_sgn(charge) && !mpq_equal(charge, q)) ||
                end - begin != 1 || id[m->row_index[begin]] >= 0)
                goto done;
            mpq_set(charge, q);
            int b = -1 - id[m->row_index[begin]], c;
            if (p->capacity[b] || !coefficient(m, begin, q, &c) || c >= 0)
                goto done;
            p->capacity[b] = -c;
            p->gate[b] = (int)j;
            continue;
        }
        int ne = 0, nr = 0, item = -1, bin = -1, size = 0;
        for (size_t k = begin; k < end; ++k) {
            int r = id[m->row_index[k]], a;
            if (!coefficient(m, k, q, &a) || a <= 0)
                goto done;
            if (r >= 0) {
                if (a != 1)
                    goto done;
                item = r;
                ++ne;
            } else {
                bin = -1 - r;
                size = a;
                ++nr;
            }
        }
        if (!ne || !nr)
            goto done;
        if (ne == 1 && nr == 1) {
            Option *o = option(p, item, bin);
            if (!o->size)
                ++p->degree[item];
            if (!o->size || size < o->size)
                *o = (Option){size, (int)j};
            if (!p->weight[item] || size < p->weight[item])
                p->weight[item] = size;
        }
    }
    for (int i = 0; i < p->items; ++i)
        if (!p->degree[i])
            goto done;
    for (int b = 0; b < p->bins; ++b)
        if (!p->capacity[b])
            goto done;
    /* Sum the item equations weighted by their minimum singleton sizes.
     * Check the inequality column by column, including every bundled option.
     * All matrix integers here were recovered exactly above. */
    for (size_t j = 0; j < m->cols; ++j) {
        if ((j & 255U) == 0 && lab_now() >= deadline)
            goto done;
        size_t begin = m->col_start[j], end = m->col_start[j + 1];
        if (end - begin == 1)
            continue; /* all non-gate columns have at least two entries */
        int64_t slack = 0;
        for (size_t k = begin; k < end; ++k) {
            int r = id[m->row_index[k]];
            slack += r >= 0 ? -p->weight[r] : (int)m->value[k];
        }
        if (slack < 0)
            goto done;
    }
    int total = 0, largest = 0;
    for (int i = 0; i < p->items; ++i)
        total += p->weight[i];
    for (int b = 0; b < p->bins; ++b)
        if (p->capacity[b] > largest)
            largest = p->capacity[b];
    p->target = (total + largest - 1) / largest;
    ok = p->target <= p->bins;
done:
    mpq_clears(q, charge, NULL);
    return ok;
}

typedef struct {
    int selected[MAX_ITEMS], load[MAX_ITEMS], bin[MAX_ITEMS], tabu[MAX_ITEMS];
} Placement;
static int place(Packing *p, Placement *s, int restart) {
    int order[MAX_ITEMS], key[MAX_ITEMS], rank[MAX_ITEMS];
    double score[MAX_ITEMS] = {0};
    for (int i = 0; i < p->items; ++i) {
        order[i] = i;
        key[i] = p->weight[i] * 100 + (int)(random_word(p) % (restart % 2 ? 1100U : 100U));
        for (int b = 0; b < p->bins; ++b)
            if (option(p, i, b)->size)
                score[b] += 1.0 / p->degree[i];
    }
    for (int b = 0; b < p->bins; ++b) {
        rank[b] = b;
        score[b] += .0001 * (random_word(p) % 1000);
    }
    for (int a = 0; a < p->target; ++a) {
        int k = a;
        for (int z = a + 1; z < p->bins; ++z)
            if (score[rank[z]] > score[rank[k]])
                k = z;
        s->selected[a] = rank[k];
        rank[k] = rank[a];
        rank[a] = s->selected[a];
    }
    for (int a = 0; a < p->items; ++a) {
        int k = a;
        for (int z = a + 1; z < p->items; ++z)
            if (key[order[z]] > key[order[k]])
                k = z;
        int i = order[k];
        order[k] = order[a];
        order[a] = i;
        int best = -1, score_best = 0;
        for (int t = 0; t < p->target; ++t) {
            int b = s->selected[t], w = option(p, i, b)->size;
            if (!w)
                continue;
            int left = p->capacity[b] - s->load[b] - w;
            int quality = left >= 0 ? left : MAX_SIZE - left;
            if (best < 0 || quality < score_best) {
                best = b;
                score_best = quality;
            }
        }
        if (best < 0)
            return 0;
        s->bin[i] = best;
        s->load[best] += option(p, i, best)->size;
    }
    return 1;
}
static int overflow(const Packing *p, const Placement *s, int bin, int change) {
    int extra = s->load[bin] + change - p->capacity[bin];
    return extra > 0 ? extra : 0;
}

typedef struct {
    int item, other, to, delta;
    unsigned ties;
} Move;
static void consider(Packing *p, Move *best, int i, int j, int b, int delta) {
    if (best->item < 0 || delta < best->delta) {
        best->delta = delta;
        best->ties = 0;
    }
    if (delta == best->delta && random_word(p) % ++best->ties == 0) {
        best->item = i;
        best->other = j;
        best->to = b;
    }
}
static Move choose(Packing *p, const Placement *s, int iter, int total, int record) {
    Move best = {.item = -1};
    for (int i = 0; i < p->items; ++i) {
        int a = s->bin[i], excess = overflow(p, s, a, 0);
        if (!excess)
            continue;
        int old = option(p, i, a)->size;
        for (int t = 0; t < p->target; ++t) {
            int b = s->selected[t], w = option(p, i, b)->size;
            if (a == b || !w)
                continue;
            int delta = overflow(p, s, a, -old) + overflow(p, s, b, w) - excess - overflow(p, s, b, 0);
            if (s->tabu[i] < iter || total + delta < record)
                consider(p, &best, i, -1, b, delta);
        }
        for (int j = 0; j < p->items; ++j) {
            int b = s->bin[j], w = option(p, i, b)->size, v = option(p, j, a)->size;
            if (a == b || !w || !v)
                continue;
            int delta = overflow(p, s, a, v - old) + overflow(p, s, b, w - option(p,j,b)->size) -
                        excess - overflow(p, s, b, 0);
            if ((s->tabu[i] < iter && s->tabu[j] < iter) || total + delta < record)
                consider(p, &best, i, j, b, delta);
        }
    }
    return best;
}
static void move_item(Packing *p, Placement *s, int i, int to, int iter) {
    s->load[s->bin[i]] -= option(p, i, s->bin[i])->size;
    s->load[to] += option(p, i, to)->size;
    s->bin[i] = to;
    s->tabu[i] = iter + 2 + (int)(random_word(p) % 5);
}

int lab_resource_packing(const LabModel *m, double deadline, size_t limit, double *point,
                         size_t *attempts) {
    *attempts = 0;
    if (!limit || m->rows < 2 || m->rows > 2 * MAX_ITEMS || m->cols > 200000 ||
        m->nnz > 1000000 || lab_now() >= deadline)
        return 0;
    Packing p = {.random = 19581};
    int id[2 * MAX_ITEMS];
    mpq_t q;
    mpq_init(q);
    for (size_t i = 0; i < m->rows; ++i) {
        if (m->row_lower[i] == 1 && m->row_upper[i] == 1 && p.items < MAX_ITEMS &&
            !lab_source_exact_row_bound(m,i,0,q) && !mpq_cmp_ui(q,1,1) &&
            !lab_source_exact_row_bound(m,i,1,q) && !mpq_cmp_ui(q,1,1))
            id[i] = p.items++;
        else if (m->row_lower[i] == -INFINITY && m->row_upper[i] == 0 && p.bins < MAX_ITEMS &&
                 !lab_source_exact_row_bound(m,i,1,q) && !mpq_sgn(q))
            id[i] = -1 - p.bins++;
        else {
            mpq_clear(q);
            return 0;
        }
    }
    mpq_clear(q);
    if (!p.items || !p.bins)
        return 0;
    p.options = lab_alloc((size_t)p.items * (size_t)p.bins, sizeof *p.options);
    int solved = 0;
    if (read_columns(&p, m, id, deadline)) {
        for (int restart = 0; restart < 64 && *attempts < limit && lab_now() < deadline; ++restart) {
            Placement s = {0};
            ++*attempts;
            if (!place(&p, &s, restart))
                continue;
            int total = 0;
            for (int b = 0; b < p.bins; ++b)
                total += overflow(&p, &s, b, 0);
            int record = total, stale = 0;
            for (int iter = 1; total && iter <= 6000 && stale < 600 && *attempts < limit &&
                               lab_now() < deadline; ++iter) {
                ++*attempts;
                Move step = choose(&p, &s, iter, total, record);
                if (step.item < 0)
                    break;
                int a = s.bin[step.item];
                move_item(&p, &s, step.item, step.to, iter);
                if (step.other >= 0)
                    move_item(&p, &s, step.other, a, iter);
                total += step.delta;
                if (total < record) {
                    record = total;
                    stale = 0;
                } else
                    ++stale;
            }
            if (total || lab_now() >= deadline)
                continue;
            memset(point, 0, m->cols * sizeof *point);
            for (int i = 0; i < p.items; ++i) {
                int b = s.bin[i];
                point[p.gate[b]] = point[option(&p, i, b)->column] = 1;
            }
            solved = 1;
            break;
        }
    }
    free(p.options);
    return solved;
}

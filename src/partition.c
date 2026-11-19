#include "partition.h"
#include "common.h"
#include "exact.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { SCALE = 1024, MEMO_BITS = 17, MASK_WORDS = 4, MAX_ROWS = MASK_WORDS * 64 };
typedef struct {
    uint64_t w[MASK_WORDS];
} Mask;
typedef struct {
    Mask mask;
    int64_t cost, reduced;
    size_t column;
} Choice;
typedef struct {
    Mask mask;
    int64_t spent;
} Visit;
typedef struct {
    Choice *choices;
    size_t offsets[MAX_ROWS + 1], *by_row;
    unsigned max_cardinality, words;
    size_t path[MAX_ROWS], witness[MAX_ROWS], length, work, limit;
    int64_t base, ceiling;
    double deadline;
    Visit *memo;
    LabPartitionResult *result;
    int stopped;
} Search;

static unsigned first_bit64(uint64_t mask) {
#if defined(__GNUC__) || defined(__clang__)
    return (unsigned)__builtin_ctzll(mask);
#else
    unsigned index = 0;
    while (!(mask & 1)) {
        mask >>= 1;
        index++;
    }
    return index;
#endif
}
static int mask_is_zero(Mask m) {
    return (m.w[0] | m.w[1] | m.w[2] | m.w[3]) == 0;
}
static int mask_eq(Mask a, Mask b) {
    return a.w[0] == b.w[0] && a.w[1] == b.w[1] && a.w[2] == b.w[2] && a.w[3] == b.w[3];
}
static int mask_cmp(Mask a, Mask b) {
    for (int k = MASK_WORDS - 1; k >= 0; --k)
        if (a.w[k] != b.w[k])
            return a.w[k] < b.w[k] ? -1 : 1;
    return 0;
}
static int mask_subset(Mask sub, Mask sup) {
    return ((sub.w[0] & sup.w[0]) == sub.w[0]) && ((sub.w[1] & sup.w[1]) == sub.w[1]) &&
           ((sub.w[2] & sup.w[2]) == sub.w[2]) && ((sub.w[3] & sup.w[3]) == sub.w[3]);
}
static Mask mask_xor(Mask a, Mask b) {
    return (Mask){{a.w[0] ^ b.w[0], a.w[1] ^ b.w[1], a.w[2] ^ b.w[2], a.w[3] ^ b.w[3]}};
}
static int mask_test(Mask m, size_t row) {
    return (m.w[row >> 6] & (UINT64_C(1) << (row & 63))) != 0;
}
static void mask_set(Mask *m, size_t row) {
    m->w[row >> 6] |= UINT64_C(1) << (row & 63);
}
static unsigned mask_first_bit(Mask m) {
    for (unsigned k = 0; k < MASK_WORDS; ++k)
        if (m.w[k])
            return k * 64U + first_bit64(m.w[k]);
    return 0;
}
static uint64_t mask_hash(Mask m) {
    uint64_t h = m.w[0] ^ (m.w[1] * UINT64_C(0x9e3779b97f4a7c15)) ^
                 (m.w[2] * UINT64_C(0xbf58476d1ce4e5b9)) ^
                 (m.w[3] * UINT64_C(0x94d049bb133111eb));
    return (h ^ (h >> 31)) * UINT64_C(11400714819323198485);
}

static int by_mask_cost(const void *l, const void *r) {
    const Choice *a = l, *b = r;
    int cmp = mask_cmp(a->mask, b->mask);
    if (cmp)
        return cmp;
    if (a->cost != b->cost)
        return a->cost < b->cost ? -1 : 1;
    return (a->column > b->column) - (a->column < b->column);
}
static int by_reduced(const void *l, const void *r) {
    const Choice *a = l, *b = r;
    if (a->reduced != b->reduced)
        return a->reduced < b->reduced ? -1 : 1;
    return by_mask_cost(l, r);
}

/* Choose an integer objective scale from the original decimal rationals.
 * Leave the bounded 64-bit search exactly the same after converting units. */
static long decimal_scale(const LabModel *m, double deadline) {
    int fractional = 0;
    for (size_t j = 0; j < m->cols; ++j)
        if (trunc(m->objective[j]) != m->objective[j]) {
            fractional = 1;
            break;
        }
    if (!fractional)
        return 1;
    mpz_t common;
    mpz_init_set_ui(common, 1);
    mpq_t q;
    mpq_init(q);
    long scale = 0;
    for (size_t j = 0; j < m->cols; ++j) {
        if (((j & 255U) == 0 && lab_now() >= deadline) || lab_source_exact_objective(m, j, q))
            goto done;
        mpz_lcm(common, common, mpq_denref(q));
        if (mpz_cmp_ui(common, 1000000000) > 0)
            goto done;
    }
    scale = mpz_get_si(common);
done:
    mpq_clear(q);
    mpz_clear(common);
    return scale;
}

static int collect(const LabModel *m, const double *row_dual, Choice *choices,
                   int64_t *price, double deadline) {
    mpq_t q;
    mpq_init(q);
    int valid = 0;
    long units = decimal_scale(m, deadline);
    if (!units)
        goto done;
    int64_t cost_cap = m->rows > 64 ? (INT64_C(1) << 36) : (INT64_C(1) << 38);
    double dual_cap = m->rows > 64 ? 0x1p36 : 0x1p38;
    for (size_t i = 0; i < m->rows; i++) {
        if (m->row_lower[i] != 1 || m->row_upper[i] != 1 ||
            lab_source_exact_row_bound(m, i, 0, q) || mpq_cmp_ui(q, 1, 1) ||
            lab_source_exact_row_bound(m, i, 1, q) || mpq_cmp_ui(q, 1, 1))
            goto done;
        double y = row_dual ? row_dual[i] * (double)units : 0;
        price[i] = isfinite(y) && fabs(y) <= dual_cap ? (int64_t)floor(y * SCALE) : 0;
    }
    for (size_t j = 0; j < m->cols; j++) {
        if (m->col_start[j] > m->col_start[j + 1] || m->col_start[j + 1] > m->nnz ||
            !m->integer[j] || m->col_lower[j] != 0 || m->col_upper[j] != 1 ||
            lab_source_exact_column_bound(m, j, 0, q) || mpq_sgn(q) ||
            lab_source_exact_column_bound(m, j, 1, q) || mpq_cmp_ui(q, 1, 1) ||
            !isfinite(m->objective[j]) || fabs(m->objective[j]) > 0x1p30 ||
            lab_source_exact_objective(m, j, q) ||
            !mpz_fits_slong_p(mpq_numref(q)) || mpz_cmp_ui(mpq_denref(q), (unsigned long)units) > 0 ||
            ((j & 255U) == 0 && lab_now() >= deadline))
            goto done;
        long den = mpz_get_si(mpq_denref(q)), num = mpz_get_si(mpq_numref(q));
        if (units % den)
            goto done;
        int64_t scale = units / den, bound = cost_cap / scale;
        if (num < -bound || num > bound)
            goto done;
        int64_t cost = (int64_t)num * scale;
        Mask mask = {{0}};
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            size_t row = m->row_index[k];
            if (row >= m->rows)
                goto done;
            if (mask_test(mask, row) || m->value[k] != 1 || lab_source_exact_coefficient(m, k, q) ||
                mpq_cmp_ui(q, 1, 1))
                goto done;
            mask_set(&mask, row);
        }
        if (mask_is_zero(mask))
            goto done;
        choices[j] = (Choice){.mask = mask, .column = j,
                              .cost = cost * (m->maximize ? -SCALE : SCALE)};
    }
    valid = 1;
done:
    mpq_clear(q);
    return valid;
}
static int64_t support(const Choice *c, const int64_t *price, unsigned words) {
    int64_t sum = 0;
    for (unsigned k = 0; k < words; ++k)
        for (uint64_t bits = c->mask.w[k]; bits; bits &= bits - 1)
            sum += price[k * 64U + first_bit64(bits)];
    return sum;
}
static int expired(Search *s) {
    if ((++s->work & 4095U) == 0 && lab_now() >= s->deadline)
        s->stopped = 1;
    return s->stopped;
}
static void explore(Search *s, Mask remaining, int64_t spent, size_t depth) {
    if (s->result->nodes == s->limit || lab_now() >= s->deadline) {
        s->stopped = 1;
        return;
    }
    s->result->nodes++;
    if (spent > s->ceiling)
        return;
    if (mask_is_zero(remaining)) {
        s->result->has_point = 1;
        s->length = depth;
        memcpy(s->witness, s->path, depth * sizeof *s->path);
        s->ceiling = spent - SCALE; /* every complete cost is integral */
        return;
    }
    /* For identical remaining rows, a smaller spent cost dominates. The
     * ceiling can only decrease within a traversal; collisions evict entries. */
    uint64_t hash = mask_hash(remaining);
    Visit *slot = &s->memo[hash >> (64 - MEMO_BITS)];
    if (mask_eq(slot->mask, remaining) && slot->spent <= spent)
        return;
    *slot = (Visit){remaining, spent};
    size_t branch = 0, best = SIZE_MAX, expensive_row = MAX_ROWS;
    int64_t budget = s->ceiling - spent;
    int64_t lower = 0, shares[MAX_ROWS], largest = -1;
    for (unsigned k = 0; k < s->words; ++k) {
        for (uint64_t rows = remaining.w[k]; rows; rows &= rows - 1) {
            size_t row = k * 64U + first_bit64(rows), count = 0;
            for (size_t p = s->offsets[row]; p < s->offsets[row + 1] && count < best; p++) {
                const Choice *c = &s->choices[s->by_row[p]];
                if (expired(s))
                    return;
                if (c->reduced > s->ceiling - spent)
                    break;
                if (mask_subset(c->mask, remaining)) {
                    /* Every completion pays for a column covering this row. Its
                     * cost can be charged to at most max_cardinality rows. */
                    if (!count) {
                        int64_t share = c->reduced / s->max_cardinality;
                        if (share > budget - lower)
                            return;
                        shares[row] = share;
                        lower += share;
                        if (c->reduced > largest) {
                            largest = c->reduced;
                            expensive_row = row;
                        }
                    }
                    ++count;
                }
            }
            if (!count)
                return;
            if (count < best) {
                best = count;
                branch = row;
            }
        }
    }
    for (size_t p = s->offsets[branch]; p < s->offsets[branch + 1]; p++) {
        const Choice *c = &s->choices[s->by_row[p]];
        if (expired(s) || c->reduced > s->ceiling - spent)
            break;
        if (mask_subset(c->mask, remaining)) {
            /* The selected column pays for its own rows. The uncovered
             * rows must still pay their shares from the parent lower bound. */
            int64_t left = s->ceiling - spent - c->reduced;
            if (expensive_row < MAX_ROWS && !mask_test(c->mask, expensive_row) && largest > left)
                continue;
            if (lower > left) {
                int64_t unpaid = lower;
                for (unsigned k = 0; k < s->words; ++k)
                    for (uint64_t bits = c->mask.w[k]; bits; bits &= bits - 1)
                        unpaid -= shares[k * 64U + first_bit64(bits)];
                if (unpaid > left)
                    continue;
            }
            s->path[depth] = c->column;
            explore(s, mask_xor(remaining, c->mask), spent + c->reduced, depth + 1);
        }
    }
}
void lab_partition_search(const LabModel *m, const double *row_dual, double deadline,
                          size_t node_limit, double *point, LabPartitionResult *r) {
    if (!r)
        return;
    *r = (LabPartitionResult){0};
    if (!m || !m->integer || !point || !node_limit || !m->rows || m->rows > MAX_ROWS ||
        !m->cols || m->cols > 200000 || m->nnz > 2000000 || !isfinite(deadline) || lab_now() >= deadline ||
        !m->row_lower || !m->row_upper || !m->col_start || !m->row_index || !m->value ||
        !m->col_lower || !m->col_upper || !m->objective || m->col_start[0] != 0 ||
        m->col_start[m->cols] != m->nnz)
        return;
    for (size_t i = 0; i < m->rows; i++)
        if (m->row_lower[i] != 1 || m->row_upper[i] != 1)
            return;
    Search s = {.deadline = deadline, .limit = node_limit, .result = r,
                .words = (unsigned)((m->rows + 63) / 64)};
    int64_t price[MAX_ROWS];
    s.choices = lab_alloc(m->cols, sizeof *s.choices);
    if (!collect(m, row_dual, s.choices, price, deadline))
        goto done;
    qsort(s.choices, m->cols, sizeof *s.choices, by_mask_cost);
    size_t n = 0;
    for (size_t j = 0; j < m->cols; j++)
        if (!n || !mask_eq(s.choices[j].mask, s.choices[n - 1].mask))
            s.choices[n++] = s.choices[j];
    /* Restore dual feasibility exactly. Decreasing a row price can never
     * increase any support sum. A path covers at most MAX_ROWS rows. */
    for (size_t j = 0; j < n; j++) {
        if ((j & 1023U) == 0 && lab_now() >= deadline)
            goto done;
        Choice *c = &s.choices[j];
        int64_t excess = support(c, price, s.words) - c->cost;
        if (excess > 0)
            price[mask_first_bit(c->mask)] -= excess;
    }
    for (size_t i = 0; i < m->rows; i++)
        s.base += price[i];
    for (size_t j = 0; j < n; j++) {
        if ((j & 1023U) == 0 && lab_now() >= deadline)
            goto done;
        Choice *c = &s.choices[j];
        c->reduced = c->cost - support(c, price, s.words);
        unsigned count = 0;
        for (unsigned k = 0; k < s.words; ++k)
            for (uint64_t bits = c->mask.w[k]; bits; bits &= bits - 1) {
                ++s.offsets[k * 64U + first_bit64(bits) + 1];
                ++count;
            }
        if (count > s.max_cardinality)
            s.max_cardinality = count;
    }
    if (lab_now() >= deadline)
        goto done;
    qsort(s.choices, n, sizeof *s.choices, by_reduced);
    for (size_t i = 0; i < m->rows; i++)
        s.offsets[i + 1] += s.offsets[i];
    s.by_row = lab_alloc(s.offsets[m->rows], sizeof *s.by_row);
    size_t next[MAX_ROWS];
    memcpy(next, s.offsets, m->rows * sizeof *next);
    for (size_t j = 0; j < n; j++)
        for (unsigned k = 0; k < s.words; ++k)
            for (uint64_t bits = s.choices[j].mask.w[k]; bits; bits &= bits - 1)
                s.by_row[next[k * 64U + first_bit64(bits)]++] = j;
    s.memo = lab_alloc(1U << MEMO_BITS, sizeof *s.memo);
    Mask all = {{0}};
    for (size_t i = 0; i < m->rows; ++i)
        mask_set(&all, i);
    /* Every completion has scaled cost base + sum(reduced). Increasing
     * ceilings find a first cover cheaply. Once a cover is found, finish
     * exhausting every strictly cheaper completion under that same ceiling. */
    for (int64_t radius = 0; !s.stopped && radius < (INT64_C(1) << 50);
         radius = radius ? 2 * radius : SCALE) {
        s.ceiling = radius + ((SCALE - s.base % SCALE) % SCALE);
        memset(s.memo, 0, (1U << MEMO_BITS) * sizeof *s.memo);
        explore(&s, all, 0, 0);
        if (r->has_point) {
            r->optimal = !s.stopped;
            memset(point, 0, m->cols * sizeof *point);
            for (size_t k = 0; k < s.length; k++)
                point[s.witness[k]] = 1;
            break;
        }
    }
done:
    free(s.memo);
    free(s.by_row);
    free(s.choices);
}

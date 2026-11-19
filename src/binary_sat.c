#define _POSIX_C_SOURCE 200809L
#include "binary_sat.h"
#include "common.h"
#include "exact.h"
#include <errno.h>
#include <kissat.h>
#include <limits.h>
#include <math.h>
#include <poll.h>
#include <signal.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

/* Emit clauses for minimal overweight subsets. Two common long cardinality
 * cases have linear encodings. Every arithmetic operation stays within the
 * exact integer range (model size and individual weights are capped). */
enum { SHORT_ROW = 16, MAX_OBJ_TERMS = 512, MAX_TERMS = 2000000, MAX_WEIGHT = 1000000000 };
typedef struct {
    int lit;
    int64_t weight;
    size_t entry;
} Term;
typedef struct {
    kissat *solver;
    int variables, base_vars;
    int *map;
    int mapped;
    size_t clauses;
    double deadline;
} Encoder;

static int stopped(const Encoder *e) {
    return lab_now() >= e->deadline;
}

/* 4.0.4 stores the termination callback without polling it. Run a decision
 * call in a child so even preprocessing is subject to the monotonic deadline.
 * A partial pipe message is never accepted as a witness. */
static int bounded_solve(Encoder *e, size_t cols, double *point) {
    int fd[2];
    if (stopped(e) || pipe(fd) != 0)
        return 0;
    pid_t pid = fork();
    if (pid == 0) {
        close(fd[0]);
        unsigned char status = (unsigned char)kissat_solve(e->solver);
        if (write(fd[1], &status, 1) != 1)
            _exit(1);
        if (status == 10 && point) {
            unsigned char chunk[2048];
            for (size_t j = 0; j < cols;) {
                size_t n = cols - j < sizeof chunk ? cols - j : sizeof chunk;
                for (size_t k = 0; k < n; ++k)
                    chunk[k] = kissat_value(e->solver, (int)(j + k) + 1) > 0;
                if (write(fd[1], chunk, n) != (ssize_t)n)
                    _exit(1);
                j += n;
            }
        }
        _exit(0);
    }
    close(fd[1]);
    if (pid < 0) {
        close(fd[0]);
        return 0;
    }
    unsigned char status = 0, chunk[2048];
    int header = 0, complete = 0;
    size_t j = 0;
    while (!complete) {
        double ms = ceil(1000 * (e->deadline - lab_now()));
        if (ms <= 0)
            break;
        struct pollfd wait = {.fd = fd[0], .events = POLLIN};
        int ready = poll(&wait, 1, (int)fmin(INT_MAX, ms));
        if (ready < 0 && errno == EINTR)
            continue;
        if (ready <= 0)
            break;
        ssize_t got = read(fd[0], header ? chunk : &status, header ? sizeof chunk : 1);
        if (got < 0 && errno == EINTR)
            continue;
        if (got <= 0)
            break;
        if (!header) {
            header = 1;
            complete = status != 10 || !point;
        } else {
            if ((size_t)got > cols - j)
                break;
            for (ssize_t k = 0; k < got; ++k)
                point[j++] = chunk[k];
            complete = j == cols;
        }
    }
    close(fd[0]);
    if (!complete)
        kill(pid, SIGKILL);
    int ws = 0;
    while (waitpid(pid, &ws, 0) < 0) {
        if (errno != EINTR)
            return 0;
    }
    return complete && WIFEXITED(ws) && WEXITSTATUS(ws) == 0 ? status : 0;
}

static int fresh_var(Encoder *e) {
    if (e->map)
        return e->base_vars + (++e->mapped);
    return ++e->variables;
}

static void literal(Encoder *e, int lit) {
    if (e->map && lit) {
        int v = abs(lit);
        int mapped;
        if (v <= e->base_vars) {
            int *slot = &e->map[v - 1];
            if (!*slot)
                *slot = ++e->mapped;
            mapped = *slot;
        } else {
            mapped = v - e->base_vars;
        }
        lit = lit < 0 ? -mapped : mapped;
    }
    kissat_add(e->solver, lit);
}

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

static int include_row_den(int64_t *scale, mpq_srcptr q) {
    enum { MAX_ROW_DEN = 1000000 };
    if (!*scale || !mpz_fits_ulong_p(mpq_denref(q)) || !mpz_fits_slong_p(mpq_numref(q))) {
        *scale = 0;
        return 0;
    }
    unsigned long den = mpz_get_ui(mpq_denref(q));
    if (!den || den > MAX_ROW_DEN) {
        *scale = 0;
        return 0;
    }
    int64_t d = (int64_t)den;
    int64_t g = gcd64(*scale, d);
    int64_t mult = d / g;
    if (*scale > MAX_ROW_DEN / mult) {
        *scale = 0;
        return 0;
    }
    *scale *= mult;
    return 1;
}

static int convert_scaled(mpq_srcptr q, double value, int64_t scale, int64_t *out) {
    if (!scale || !isfinite(value) || !mpz_fits_ulong_p(mpq_denref(q)) ||
        !mpz_fits_slong_p(mpq_numref(q)))
        return 0;
    unsigned long den = mpz_get_ui(mpq_denref(q));
    if (!den || scale % (int64_t)den)
        return 0;
    long num = mpz_get_si(mpq_numref(q));
    int64_t factor = scale / (int64_t)den;
    if (num && llabs((int64_t)num) > MAX_WEIGHT / factor)
        return 0;
    *out = (int64_t)num * factor;
    return 1;
}

static int convert(mpq_srcptr q, double value, int64_t *out) {
    return convert_scaled(q, value, 1024, out);
}

static int heavier(const void *a, const void *b) {
    int64_t x = ((const Term *)a)->weight, y = ((const Term *)b)->weight;
    return (x < y) - (x > y);
}
static int covers(Encoder *e, const Term *t, int count, int index, int64_t capacity,
                  int64_t remaining, int *clause, int used) {
    if (capacity < 0) {
        if (++e->clauses > 2000000 || (!(e->clauses % 1024) && stopped(e)))
            return 0;
        for (int k = 0; k < used; ++k)
            literal(e, clause[k]);
        literal(e, 0);
        return 1;
    }
    if (remaining <= capacity || index == count)
        return 1;
    int64_t w = t[index].weight;
    clause[used] = -t[index].lit;
    return covers(e, t, count, index + 1, capacity - w, remaining - w, clause, used + 1) &&
           covers(e, t, count, index + 1, capacity, remaining - w, clause, used);
}

static void binary_clause(Encoder *e, int a, int b) {
    literal(e, a);
    literal(e, b);
    literal(e, 0);
    ++e->clauses;
}

static void ternary_clause(Encoder *e, int a, int b, int c) {
    literal(e, a);
    literal(e, b);
    literal(e, c);
    literal(e, 0);
    ++e->clauses;
}

/* Exact weighted sequential counter for sum_{i=0..n-1} w_i * l_i <= cap
 * with w_i >= 1, 1 <= cap <= 512. Auxiliary variable s[u] encodes prefix sum >= u. */
static int sequential_counter(Encoder *e, const Term *t, int n, int64_t cap) {
    int C = (int)cap;
    int *prev = lab_alloc((size_t)C + 1, sizeof *prev);
    int *curr = lab_alloc((size_t)C + 1, sizeof *curr);
    int prev_max = 0, ok = 1;
    for (int i = 0; i < n && ok; ++i) {
        if ((i & 15) == 0 && (e->clauses > 2000000 || stopped(e))) {
            ok = 0;
            break;
        }
        int lit = t[i].lit;
        int w = (int)t[i].weight;
        if (w > C) {
            literal(e, -lit);
            literal(e, 0);
            ++e->clauses;
            continue;
        }
        int curr_max = prev_max + w;
        if (curr_max > C)
            curr_max = C;
        int last = (i == n - 1);
        if (!last) {
            for (int v = 1; v <= curr_max; ++v)
                curr[v] = fresh_var(e);
            int w_lim = w < curr_max ? w : curr_max;
            for (int v = 1; v <= w_lim; ++v)
                binary_clause(e, -lit, curr[v]);
            for (int u = 1; u <= prev_max; ++u) {
                binary_clause(e, -prev[u], curr[u]);
                if (u + w > C)
                    binary_clause(e, -lit, -prev[u]);
                else
                    ternary_clause(e, -lit, -prev[u], curr[u + w]);
            }
            memcpy(prev + 1, curr + 1, (size_t)curr_max * sizeof *prev);
            prev_max = curr_max;
        } else {
            for (int u = 1; u <= prev_max; ++u)
                if (u + w > C)
                    binary_clause(e, -lit, -prev[u]);
        }
        if (e->clauses > 2000000)
            ok = 0;
    }
    free(curr);
    free(prev);
    return ok;
}

static int at_most(Encoder *e, Term *t, int n, int64_t cap) {
    int64_t total = 0, minimum = INT64_MAX;
    for (int i = 0; i < n; ++i) {
        if (t[i].weight < 0) {
            cap -= t[i].weight;
            t[i].weight = -t[i].weight;
            t[i].lit = -t[i].lit;
        }
        total += t[i].weight;
        if (t[i].weight < minimum)
            minimum = t[i].weight;
    }
    if (cap >= total)
        return 1;
    if (cap < 0) {
        literal(e, 0);
        ++e->clauses;
        return 1;
    }
    if (total - minimum <= cap) {
        for (int i = 0; i < n; ++i)
            literal(e, -t[i].lit);
        literal(e, 0);
        ++e->clauses;
        return 1;
    }
    if (cap < minimum) {
        for (int i = 0; i < n; ++i) {
            literal(e, -t[i].lit);
            literal(e, 0);
            ++e->clauses;
        }
        return 1;
    }
    if (n > SHORT_ROW) {
        /* Strip terms heavier than cap immediately, then divide by GCD. */
        int m = 0;
        int64_t g = 0, rem_total = 0;
        for (int i = 0; i < n; ++i) {
            if (t[i].weight > cap) {
                literal(e, -t[i].lit);
                literal(e, 0);
                ++e->clauses;
                continue;
            }
            g = gcd64(g, t[i].weight);
            rem_total += t[i].weight;
            t[m++] = t[i];
        }
        if (!m || cap >= rem_total)
            return 1;
        int64_t scaled_cap = cap / g;
        for (int i = 0; i < m; ++i)
            t[i].weight /= g;
        if (scaled_cap == 1) {
            int prev = 0;
            for (int i = 0; i < m; ++i) {
                int lit = t[i].lit;
                if (prev)
                    binary_clause(e, -prev, -lit);
                if (i + 1 < m) {
                    int next = fresh_var(e);
                    if (prev)
                        binary_clause(e, -prev, next);
                    binary_clause(e, -lit, next);
                    prev = next;
                }
            }
            return e->clauses <= 2000000;
        }
        if (scaled_cap <= 512 && (int64_t)m * scaled_cap <= 120000)
            return sequential_counter(e, t, m, scaled_cap);
        return 0;
    }
    int clause[SHORT_ROW];
    qsort(t, (size_t)n, sizeof *t, heavier);
    return covers(e, t, n, 0, cap, total, clause, 0);
}

/* Exact scalar fallback for assignments with overlapping source intervals. */
static int exact_assignment(const LabModel *m, const Term *t, size_t n, size_t row, unsigned bits) {
    mpq_t sum, q;
    mpq_inits(sum, q, NULL);
    int valid = 1;
    for (size_t k = 0; k < n; ++k)
        if (bits & (1u << k)) {
            if (lab_source_exact_coefficient(m, t[k].entry, q)) {
                valid = -1;
                break;
            }
            mpq_add(sum, sum, q);
        }
    for (int upper = 0; valid == 1 && upper < 2; ++upper) {
        if (!isfinite(upper ? m->row_upper[row] : m->row_lower[row]))
            continue;
        if (lab_source_exact_row_bound(m, row, upper, q))
            valid = -1;
        else if (upper ? mpq_cmp(sum, q) > 0 : mpq_cmp(sum, q) < 0)
            valid = 0;
    }
    mpq_clears(sum, q, NULL);
    return valid;
}

/* Resolve short non-dyadic rows by interval arithmetic on every binary
 * assignment. Emit only proven violations and abandon the entire encoding
 * if any assignment straddles a source boundary. */
static int truth_table(const LabModel *m, Encoder *e, const Term *t, size_t n, size_t i) {
    const LabSourceEnvelope *v = m->envelope;
    if (n > 10 || !v || !v->complete || v->nnz != m->nnz || v->rows != m->rows)
        return 0;
    unsigned char forbidden[1024] = {0}, covered[1024] = {0};
    unsigned assignments = 1u << n;
    for (unsigned bits = 0; bits < (1u << n); ++bits) {
        LabInterval a = {0, 0};
        for (size_t k = 0; k < n; ++k)
            if (bits & (1u << k))
                a = lab_interval_add(a, v->value[t[k].entry]);
        LabInterval lo = v->row_lower[i], hi = v->row_upper[i];
        if (a.lower >= lo.upper && a.upper <= hi.lower)
            continue;
        if (!(a.upper < lo.lower || a.lower > hi.upper)) {
            int valid = exact_assignment(m, t, n, i, bits);
            if (valid < 0)
                return 0;
            if (valid)
                continue;
        }
        forbidden[bits] = 1;
    }
    /* Shorten each clause while every completion of its remaining literals
     * still violates the original row. This exposes implications to SAT unit
     * propagation even when source decimals required a truth table. */
    for (unsigned bits = 0; bits < assignments; ++bits) {
        if (!forbidden[bits] || covered[bits])
            continue;
        unsigned free_mask = 0;
        for (size_t k = 0; k < n; ++k) {
            unsigned expanded = free_mask | (1u << k);
            unsigned fixed = bits & ~expanded;
            int all_bad = 1;
            for (unsigned sub = expanded;; sub = (sub - 1) & expanded) {
                if (!forbidden[fixed | sub]) {
                    all_bad = 0;
                    break;
                }
                if (!sub)
                    break;
            }
            if (all_bad)
                free_mask = expanded;
        }
        for (unsigned sub = free_mask;; sub = (sub - 1) & free_mask) {
            covered[(bits & ~free_mask) | sub] = 1;
            if (!sub)
                break;
        }
        if (++e->clauses > 2000000 || stopped(e))
            return 0;
        for (size_t k = 0; k < n; ++k)
            if (!(free_mask & (1u << k)))
                literal(e, (bits & (1u << k)) ? -t[k].lit : t[k].lit);
        literal(e, 0);
    }
    return 1;
}

typedef struct {
    int n, lit[MAX_OBJ_TERMS];
    mpq_t weight[MAX_OBJ_TERMS], cutoff;
    int64_t scaled[MAX_OBJ_TERMS], scaled_gcd, scaled_total;
    int pb_exact;
} Objective;

/* A cheaper assignment cannot contain any subset with cost >= cutoff.
 * Costs represent positive penalties on signed literals. */
static int objective_cut(Encoder *e, const Objective *cost, int from, int *selected, int used,
                         mpq_t sum) {
    if (mpq_cmp(sum, cost->cutoff) >= 0) {
        for (int i = 0; i < used; ++i)
            literal(e, -selected[i]);
        literal(e, 0);
        ++e->clauses;
        return 1;
    }
    if (stopped(e))
        return 0;
    for (int k = from; k < cost->n; ++k) {
        selected[used] = cost->lit[k];
        mpq_add(sum, sum, cost->weight[k]);
        int ok = objective_cut(e, cost, k + 1, selected, used + 1, sum);
        mpq_sub(sum, sum, cost->weight[k]);
        if (!ok)
            return 0;
    }
    return 1;
}

static int objective_pb_cut(Encoder *e, const Objective *cost, int64_t strict_bound) {
    if (strict_bound <= 0) {
        literal(e, 0);
        ++e->clauses;
        return 1;
    }
    Term *terms = lab_alloc((size_t)cost->n, sizeof *terms);
    for (int k = 0; k < cost->n; ++k)
        terms[k] = (Term){cost->lit[k], cost->scaled[k], 0};
    int ok = at_most(e, terms, cost->n, strict_bound - 1);
    free(terms);
    return ok;
}

/* Build both sides separately from the immutable original CSC data. */
static int encode(const LabModel *m, Encoder *e, const size_t *start, const Term *terms,
                  const int64_t *row_scale, Term *row, const Objective *cost, int improving,
                  int64_t pb_cap, const size_t *indices, size_t count, int *skipped) {
    mpq_t q;
    mpq_init(q);
    int ok = 1;
    if (skipped)
        *skipped = 0;
    for (size_t index = 0; ok && index < count; ++index) {
        size_t i = indices ? indices[index] : index;
        if (!(index % 128) && stopped(e)) {
            ok = 0;
            break;
        }
        size_t n = start[i + 1] - start[i];
        int interval = !row_scale[i];
        for (size_t k = start[i]; !interval && k < start[i + 1]; ++k)
            if (!terms[k].weight)
                interval = 1;
        if (interval) {
            int row_ok = truth_table(m, e, terms + start[i], n, i);
            if (!row_ok) {
                if (skipped && !stopped(e) && e->clauses < 1900000)
                    ++(*skipped);
                else
                    ok = 0;
            }
            continue;
        }
        for (int upper = 0; upper < 2 && ok; ++upper) {
            double rhs = upper ? m->row_upper[i] : m->row_lower[i];
            if (!isfinite(rhs))
                continue;
            int64_t cap;
            if (lab_source_exact_row_bound(m, i, upper, q) ||
                !convert_scaled(q, rhs, row_scale[i], &cap)) {
                int row_ok = truth_table(m, e, terms + start[i], n, i);
                if (!row_ok) {
                    if (skipped && !stopped(e) && e->clauses < 1900000)
                        ++(*skipped);
                    else
                        ok = 0;
                }
                break;
            }
            memcpy(row, terms + start[i], n * sizeof *row);
            if (!upper) {
                cap = -cap;
                for (size_t k = 0; k < n; ++k)
                    row[k].weight = -row[k].weight;
            }
            int row_ok = at_most(e, row, (int)n, cap);
            if (!row_ok) {
                if (skipped && !stopped(e) && e->clauses < 1900000)
                    ++(*skipped);
                else
                    ok = 0;
            }
        }
    }
    if (ok && improving) {
        if (cost->n <= SHORT_ROW) {
            mpq_set_ui(q, 0, 1);
            int selected[SHORT_ROW];
            ok = objective_cut(e, cost, 0, selected, 0, q);
        } else if (cost->pb_exact && pb_cap >= 0) {
            ok = objective_pb_cut(e, cost, pb_cap + 1);
        } else {
            ok = 0;
        }
    }
    mpq_clear(q);
    return ok;
}
/* Collect a small connected row relaxation around one ordinary column.
 * Objective columns remain in its rows but do not connect distant rows.
 * UNSAT with the improving cutoff suffices to close the full model. */
static size_t neighborhood(const LabModel *m, const size_t *start, const Term *terms, size_t seed,
                           unsigned char *interior, size_t *indices, size_t maxrows,
                           uint64_t *shape) {
    enum { WIDTH = 384 };
    size_t queue[WIDTH], head = 0, tail = 0;
    unsigned char *seen = lab_alloc(m->cols, 1), *selected = lab_alloc(m->rows, 1);
    size_t count = 0;
    while (seed < m->cols && m->objective[seed] != 0)
        ++seed;
    if (seed < m->cols) {
        queue[tail++] = seed;
        seen[seed] = 1;
    }
    size_t layer_end = tail;
    *shape = 0;
    while (head < tail && tail < WIDTH) {
        size_t j = queue[head++];
        if (head <= WIDTH / 4)
            interior[j] = 1;
        for (size_t p = m->col_start[j]; p < m->col_start[j + 1]; ++p) {
            size_t i = m->row_index[p];
            if (count >= maxrows || selected[i] || start[i + 1] - start[i] > SHORT_ROW)
                continue;
            selected[i] = 1;
            indices[count++] = i;
            for (size_t k = start[i]; k < start[i + 1] && tail < WIDTH; ++k) {
                size_t v = (size_t)(terms[k].lit - 1);
                if (seen[v] || m->objective[v] != 0)
                    continue;
                queue[tail++] = v;
                seen[v] = 1;
            }
        }
        if (head == layer_end && tail < WIDTH) {
            *shape = (*shape * UINT64_C(65599) + tail) * UINT64_C(65599) + count;
            layer_end = tail;
        }
    }
    free(selected);
    free(seen);
    return count;
}

typedef struct {
    size_t seed, frequency;
    uint64_t shape;
} Region;

static int by_shape(const void *x, const void *y) {
    const Region *a = x, *b = y;
    return (a->shape > b->shape) - (a->shape < b->shape);
}
static int by_frequency(const void *x, const void *y) {
    const Region *a = x, *b = y;
    if (a->frequency != b->frequency)
        return (a->frequency > b->frequency) - (a->frequency < b->frequency);
    return (a->seed > b->seed) - (a->seed < b->seed);
}

/* Try rare local constraint structures first. Fingerprints only change
 * scheduling. Every retained lower bound comes from an actual UNSAT solve. */
static size_t rank_regions(const LabModel *m, const size_t *start, const Term *terms,
                           unsigned char *interior, Region *regions, size_t limit,
                           double deadline) {
    size_t indices[3000], count = 0;
    memset(interior, 0, m->cols);
    for (size_t j = 0; j < m->cols && count < limit; ++j) {
        if (interior[j] || m->objective[j] != 0)
            continue;
        if (lab_now() >= deadline)
            break;
        Region *r = &regions[count++];
        r->seed = j;
        neighborhood(m, start, terms, j, interior, indices, 3000, &r->shape);
    }
    qsort(regions, count, sizeof *regions, by_shape);
    for (size_t a = 0, b; a < count; a = b) {
        for (b = a + 1; b < count && regions[b].shape == regions[a].shape; ++b) {
        }
        for (size_t k = a; k < b; ++k)
            regions[k].frequency = b - a;
    }
    qsort(regions, count, sizeof *regions, by_frequency);
    return count;
}

/* Fix each priced binary at its cheapest permitted endpoint. A feasible
 * point on this face attains the objective bound from the box alone.
 * During SAT search we may omit hard inequalities. An exact second pass
 * checks both sides of every row before accepting the point. */
static int cheapest_face(const LabModel *m, const size_t *start, const Term *terms,
                         const int64_t *row_scale, Term *row, const signed char *fixed,
                         double deadline, size_t budget, double *point, size_t *calls) {
    Encoder e = {.solver = kissat_init(), .variables = (int)m->cols, .base_vars = (int)m->cols,
                 .deadline = deadline};
    kissat_reserve(e.solver, e.variables);
    mpq_t q;
    mpq_init(q);
    int valid = 0;
    for (int check = 0; check <= 1; ++check) {
        for (size_t i = 0; i < m->rows; ++i) {
            if (i % 128 == 0 && stopped(&e))
                goto done;
            for (int upper = 0; upper < 2; ++upper) {
                double rhs = upper ? m->row_upper[i] : m->row_lower[i];
                if (!isfinite(rhs))
                    continue;
                int64_t cap;
                if (lab_source_exact_row_bound(m, i, upper, q) ||
                    !convert_scaled(q, rhs, row_scale[i], &cap))
                    goto done;
                if (!upper)
                    cap = -cap;
                size_t count = 0;
                int64_t minimum = 0, maximum = 0;
                for (size_t k = start[i]; k < start[i + 1]; ++k) {
                    Term t = terms[k];
                    if (!t.weight)
                        goto done;
                    if (!upper)
                        t.weight = -t.weight;
                    size_t j = (size_t)t.lit - 1;
                    int value = check ? (int)point[j] : fixed[j];
                    if (value >= 0) {
                        cap -= value * t.weight;
                    } else {
                        row[count++] = t;
                        if (t.weight < 0)
                            minimum += t.weight;
                        else
                            maximum += t.weight;
                    }
                }
                /* MAX_TERMS * MAX_WEIGHT fits in int64_t. In the checking
                 * pass every column has a value and both extrema are zero. */
                if (minimum > cap)
                    goto done;
                if (maximum > cap && e.clauses < 1800000)
                    (void)at_most(&e, row, (int)count, cap);
            }
        }
        if (!check) {
            kissat_set_decision_limit(e.solver, (unsigned)fmin(UINT_MAX, (double)budget));
            ++*calls;
            if (bounded_solve(&e, m->cols, point) != 10)
                goto done;
            for (size_t j = 0; j < m->cols; ++j)
                if (fixed[j] >= 0)
                    point[j] = fixed[j];
        }
    }
    valid = !stopped(&e);
done:
    mpq_clear(q);
    kissat_release(e.solver);
    return valid;
}

/* Check in exact source arithmetic whether sum_{i in E} a_{ij} == k * c_j for
 * every column j and a single nonzero rational k, where E is the set of exact
 * equality rows. When true, c^T x is constant over all feasible solutions. */
static int constant_objective_via_equalities(const LabModel *m, double deadline) {
    unsigned char *is_eq = lab_alloc(m->rows ? m->rows : 1, 1);
    mpq_t q, q2, sum, ratio, cur;
    mpq_inits(q, q2, sum, ratio, cur, NULL);
    size_t eq_count = 0;
    int valid = 1;
    for (size_t i = 0; i < m->rows; ++i) {
        if ((i & 255U) == 0 && lab_now() >= deadline) {
            valid = 0;
            break;
        }
        if (isfinite(m->row_lower[i]) && m->row_lower[i] == m->row_upper[i]) {
            if (lab_source_exact_row_bound(m, i, 0, q) ||
                lab_source_exact_row_bound(m, i, 1, q2) || mpq_cmp(q, q2)) {
                valid = 0;
                break;
            }
            is_eq[i] = 1;
            ++eq_count;
        }
    }
    if (!eq_count)
        valid = 0;
    int have_ratio = 0;
    for (size_t j = 0; valid && j < m->cols; ++j) {
        if ((j & 255U) == 0 && lab_now() >= deadline) {
            valid = 0;
            break;
        }
        if (lab_source_exact_objective(m, j, q)) {
            valid = 0;
            break;
        }
        mpq_set_ui(sum, 0, 1);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            if (!is_eq[m->row_index[k]])
                continue;
            if (lab_source_exact_coefficient(m, k, cur)) {
                valid = 0;
                break;
            }
            mpq_add(sum, sum, cur);
        }
        if (!valid)
            break;
        int s_obj = mpq_sgn(q), s_sum = mpq_sgn(sum);
        if (!s_obj && !s_sum)
            continue;
        if (!s_obj || !s_sum) {
            valid = 0;
            break;
        }
        mpq_div(cur, sum, q);
        if (!have_ratio) {
            mpq_set(ratio, cur);
            have_ratio = 1;
        } else if (mpq_cmp(cur, ratio)) {
            valid = 0;
            break;
        }
    }
    mpq_clears(q, q2, sum, ratio, cur, NULL);
    free(is_eq);
    return valid && have_ratio;
}

static int order_lit(int base, int width, int neg, int v) {
    return neg ? -(base + width - v) : (base + v - 1);
}

/* Order-encoded SAT/PB translation for pure-integer models whose non-binary
 * variables all have small finite integer domains (width <= 32). Handles exact
 * feasibility, infeasibility, and constant-objective optimality. */
static int small_integer_sat(const LabModel *m, const double *lo, const double *hi, double deadline,
                             size_t budget, double *point, size_t *calls) {
    enum { MAX_DOMAIN_WIDTH = 32, MAX_ORDER_BITS = 524288 };
    mpq_t q;
    mpq_init(q);
    int64_t *col_lo = lab_alloc(m->cols, sizeof *col_lo);
    int *width = lab_alloc(m->cols, sizeof *width);
    int *base = lab_alloc(m->cols, sizeof *base);
    size_t total_bits = 0;
    int ok = 1, has_obj = 0;
    for (size_t j = 0; j < m->cols && ok; ++j) {
        if ((j & 255U) == 0 && lab_now() >= deadline) {
            ok = 0;
            break;
        }
        double l = m->col_lower[j], u = m->col_upper[j];
        if (!m->integer[j] || !isfinite(l) || !isfinite(u) || !isfinite(lo[j]) || !isfinite(hi[j]) ||
            trunc(l) != l || trunc(u) != u || trunc(lo[j]) != lo[j] || trunc(hi[j]) != hi[j] ||
            l > lo[j] || u < hi[j] || lo[j] > hi[j] ||
            lab_source_exact_column_bound(m, j, 0, q) || mpz_cmp_ui(mpq_denref(q), 1) ||
            !mpz_fits_slong_p(mpq_numref(q)) || (double)mpz_get_si(mpq_numref(q)) != l ||
            lab_source_exact_column_bound(m, j, 1, q) || mpz_cmp_ui(mpq_denref(q), 1) ||
            !mpz_fits_slong_p(mpq_numref(q)) || (double)mpz_get_si(mpq_numref(q)) != u ||
            lab_source_exact_objective(m, j, q)) {
            ok = 0;
            break;
        }
        if (mpq_sgn(q) && lo[j] != hi[j])
            has_obj = 1;
        int64_t ilo = (int64_t)lo[j], ihi = (int64_t)hi[j];
        if (llabs(ilo) > MAX_WEIGHT || llabs(ihi) > MAX_WEIGHT || ihi - ilo > MAX_DOMAIN_WIDTH) {
            ok = 0;
            break;
        }
        col_lo[j] = ilo;
        width[j] = (int)(ihi - ilo);
        if (width[j] > 0) {
            if (total_bits + (size_t)width[j] > MAX_ORDER_BITS) {
                ok = 0;
                break;
            }
            base[j] = (int)total_bits + 1;
            total_bits += (size_t)width[j];
        }
    }
    if (!ok || !total_bits) {
        mpq_clear(q);
        free(col_lo);
        free(width);
        free(base);
        return 0;
    }
    int const_obj = !has_obj || constant_objective_via_equalities(m, fmin(deadline, lab_now() + .5));
    /* Build CSR of exact converted entries. */
    size_t *start = lab_alloc(m->rows + 1, sizeof *start);
    for (size_t k = 0; k < m->nnz; ++k)
        ++start[m->row_index[k] + 1];
    size_t max_row = 0;
    for (size_t i = 0; i < m->rows; ++i) {
        if (start[i + 1] > max_row)
            max_row = start[i + 1];
        start[i + 1] += start[i];
    }
    size_t *next = lab_alloc(m->rows, sizeof *next);
    memcpy(next, start, m->rows * sizeof *next);
    Term *terms = lab_alloc(m->nnz + 1, sizeof *terms);
    for (size_t j = 0; j < m->cols && ok; ++j) {
        if ((j & 255U) == 0 && lab_now() >= deadline) {
            ok = 0;
            break;
        }
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            int64_t w = 0;
            if (lab_source_exact_coefficient(m, k, q) || !convert(q, m->value[k], &w)) {
                ok = 0;
                break;
            }
            terms[next[m->row_index[k]]++] = (Term){(int)j, w, k};
        }
    }
    Encoder e = {.solver = kissat_init(), .variables = (int)total_bits,
                 .base_vars = (int)total_bits, .deadline = deadline};
    Term *expanded = NULL;
    double *bits = NULL;
    int result = 0;
    if (!ok)
        goto done;
    kissat_reserve(e.solver, e.variables);
    /* Symmetry-breaking order clauses: b_{t} => b_{t-1}. */
    for (size_t j = 0; j < m->cols; ++j)
        for (int t = 1; t < width[j]; ++t)
            binary_clause(&e, -(base[j] + t), base[j] + t - 1);
    size_t exp_cap = max_row * MAX_DOMAIN_WIDTH + 1;
    expanded = lab_alloc(exp_cap, sizeof *expanded);
    for (size_t i = 0; i < m->rows && ok; ++i) {
        if ((i & 127U) == 0 && (e.clauses > 3000000 || stopped(&e))) {
            ok = 0;
            break;
        }
        for (int upper = 0; upper < 2 && ok; ++upper) {
            double rhs = upper ? m->row_upper[i] : m->row_lower[i];
            if (!isfinite(rhs))
                continue;
            int64_t cap = 0;
            if (lab_source_exact_row_bound(m, i, upper, q) || !convert(q, rhs, &cap)) {
                ok = 0;
                break;
            }
            if (!upper)
                cap = -cap;
            /* Collect active (non-fixed) normalized terms. */
            int k_act = 0;
            int act_col[4], act_w[4], act_neg[4];
            int64_t act_coef[4];
            int has_multi = 0;
            int64_t adj_cap = cap;
            for (size_t p = start[i]; p < start[i + 1]; ++p) {
                size_t j = (size_t)terms[p].lit;
                int64_t a = upper ? terms[p].weight : -terms[p].weight;
                if (!a)
                    continue;
                if (col_lo[j]) {
                    if (llabs(col_lo[j]) > INT64_MAX / (llabs(a) + 1)) {
                        ok = 0;
                        break;
                    }
                    adj_cap -= a * col_lo[j];
                }
                if (!width[j])
                    continue;
                if (width[j] > 1)
                    has_multi = 1;
                if (k_act < 4) {
                    act_col[k_act] = (int)j;
                    act_w[k_act] = width[j];
                    if (a < 0) {
                        act_neg[k_act] = 1;
                        act_coef[k_act] = -a;
                        adj_cap += (-a) * width[j];
                    } else {
                        act_neg[k_act] = 0;
                        act_coef[k_act] = a;
                    }
                }
                ++k_act;
            }
            if (!ok)
                break;
            if (k_act == 0) {
                if (adj_cap < 0) {
                    literal(&e, 0);
                    ++e.clauses;
                }
                continue;
            }
            if (k_act == 1) {
                int64_t v = adj_cap < 0 ? 0 : (adj_cap / act_coef[0] + 1);
                if (v == 0) {
                    literal(&e, 0);
                    ++e.clauses;
                } else if (v <= act_w[0]) {
                    literal(&e, -order_lit(base[act_col[0]], act_w[0], act_neg[0], (int)v));
                    literal(&e, 0);
                    ++e.clauses;
                }
                continue;
            }
            if (k_act == 2 && has_multi) {
                for (int v0 = 0; v0 <= act_w[0]; ++v0) {
                    int64_t rem = adj_cap - act_coef[0] * v0;
                    int64_t v1 = rem < 0 ? 0 : (rem / act_coef[1] + 1);
                    if (v1 > act_w[1])
                        continue;
                    if (v0 > 0)
                        literal(&e, -order_lit(base[act_col[0]], act_w[0], act_neg[0], v0));
                    if (v1 > 0)
                        literal(&e, -order_lit(base[act_col[1]], act_w[1], act_neg[1], (int)v1));
                    literal(&e, 0);
                    ++e.clauses;
                }
                continue;
            }
            if (k_act == 3 && has_multi) {
                /* Order by smallest domain first for minimal clause count. */
                for (int a = 0; a < 3; ++a)
                    for (int b = a + 1; b < 3; ++b)
                        if (act_w[b] < act_w[a]) {
                            int tc = act_col[a], tw = act_w[a], tn = act_neg[a];
                            int64_t tf = act_coef[a];
                            act_col[a] = act_col[b]; act_w[a] = act_w[b];
                            act_neg[a] = act_neg[b]; act_coef[a] = act_coef[b];
                            act_col[b] = tc; act_w[b] = tw;
                            act_neg[b] = tn; act_coef[b] = tf;
                        }
                if ((act_w[0] + 1) * (act_w[1] + 1) <= 1024) {
                    for (int v0 = 0; v0 <= act_w[0]; ++v0) {
                        for (int v1 = 0; v1 <= act_w[1]; ++v1) {
                            int64_t rem = adj_cap - act_coef[0] * v0 - act_coef[1] * v1;
                            int64_t v2 = rem < 0 ? 0 : (rem / act_coef[2] + 1);
                            if (v2 > act_w[2])
                                continue;
                            if (v0 > 0)
                                literal(&e, -order_lit(base[act_col[0]], act_w[0], act_neg[0], v0));
                            if (v1 > 0)
                                literal(&e, -order_lit(base[act_col[1]], act_w[1], act_neg[1], v1));
                            if (v2 > 0)
                                literal(&e, -order_lit(base[act_col[2]], act_w[2], act_neg[2], (int)v2));
                            literal(&e, 0);
                            ++e.clauses;
                        }
                    }
                    continue;
                }
            }
            /* General expansion into order-encoded binary terms. */
            int exp_n = 0;
            int64_t exp_cap_rhs = cap;
            for (size_t p = start[i]; p < start[i + 1]; ++p) {
                size_t j = (size_t)terms[p].lit;
                int64_t a = upper ? terms[p].weight : -terms[p].weight;
                if (!a)
                    continue;
                exp_cap_rhs -= a * col_lo[j];
                for (int t = 0; t < width[j]; ++t)
                    expanded[exp_n++] = (Term){base[j] + t, a, terms[p].entry};
            }
            ok = at_most(&e, expanded, exp_n, exp_cap_rhs);
        }
    }
    if (!ok || stopped(&e))
        goto done;
    unsigned slice = (unsigned)fmin(UINT_MAX, (double)budget);
    kissat_set_decision_limit(e.solver, slice);
    ++*calls;
    bits = lab_alloc(total_bits, sizeof *bits);
    int status = bounded_solve(&e, total_bits, bits);
    if (status == 20) {
        result = -1;
    } else if (status == 10) {
        for (size_t j = 0; j < m->cols; ++j) {
            int64_t val = col_lo[j];
            for (int t = 0; t < width[j]; ++t)
                if (bits[base[j] - 1 + t] != 0)
                    ++val;
            point[j] = (double)val;
        }
        result = const_obj ? 1 : 2;
    }
done:
    free(bits);
    free(expanded);
    kissat_release(e.solver);
    free(terms);
    free(next);
    free(start);
    free(base);
    free(width);
    free(col_lo);
    mpq_clear(q);
    return lab_now() < deadline ? result : 0;
}

/* Fix all priced integer variables to their cheapest finite box endpoint,
 * propagate tight row faces (forcing nonnegative penalty slacks to zero when
 * their sum is fixed to zero), and solve the resulting binary face in Kissat.
 * A verified point on this face attains the box lower bound and is optimal. */
static int cheapest_face_general_int(const LabModel *m, const double *lo, const double *hi,
                                     double deadline, size_t budget, double *point, size_t *calls) {
    if (!m->cols || m->cols > MAX_TERMS || m->rows > MAX_TERMS || m->nnz > MAX_TERMS ||
        lab_now() >= deadline)
        return 0;
    double face_dl = fmin(deadline, lab_now() + fmin(8.0, .35 * (deadline - lab_now())));
    mpq_t q;
    mpq_init(q);
    double *clo = lab_alloc(m->cols, sizeof *clo);
    double *chi = lab_alloc(m->cols, sizeof *chi);
    int ok = 1, any_priced = 0;
    for (size_t j = 0; j < m->cols && ok; ++j) {
        if ((j & 1023U) == 0 && lab_now() >= face_dl) {
            ok = 0;
            break;
        }
        if (!m->integer[j] || lo[j] > hi[j] || lab_source_exact_objective(m, j, q)) {
            ok = 0;
            break;
        }
        clo[j] = lo[j];
        chi[j] = hi[j];
        int sign = mpq_sgn(q);
        if (m->maximize)
            sign = -sign;
        if (sign > 0) {
            if (!isfinite(clo[j]) || trunc(clo[j]) != clo[j] || fabs(clo[j]) > MAX_WEIGHT) {
                ok = 0;
                break;
            }
            chi[j] = clo[j];
            any_priced = 1;
        } else if (sign < 0) {
            if (!isfinite(chi[j]) || trunc(chi[j]) != chi[j] || fabs(chi[j]) > MAX_WEIGHT) {
                ok = 0;
                break;
            }
            clo[j] = chi[j];
            any_priced = 1;
        }
    }
    if (!ok || !any_priced) {
        mpq_clear(q);
        free(clo);
        free(chi);
        return 0;
    }

    size_t *start = lab_alloc(m->rows + 1, sizeof *start);
    for (size_t k = 0; k < m->nnz; ++k)
        ++start[m->row_index[k] + 1];
    size_t max_row = 0;
    for (size_t i = 0; i < m->rows; ++i) {
        if (start[i + 1] > max_row)
            max_row = start[i + 1];
        start[i + 1] += start[i];
    }
    size_t *next = lab_alloc(m->rows, sizeof *next);
    memcpy(next, start, m->rows * sizeof *next);
    size_t *rcol = lab_alloc(m->nnz + 1, sizeof *rcol);
    double *rval = lab_alloc(m->nnz + 1, sizeof *rval);
    int64_t *rcoef = lab_alloc(m->nnz + 1, sizeof *rcoef);
    int64_t *row_scale = lab_alloc(m->rows ? m->rows : 1, sizeof *row_scale);
    for (size_t i = 0; i < m->rows; ++i) {
        row_scale[i] = 1;
        for (int upper = 0; upper < 2; ++upper) {
            double rhs = upper ? m->row_upper[i] : m->row_lower[i];
            if (isfinite(rhs)) {
                if (lab_source_exact_row_bound(m, i, upper, q))
                    row_scale[i] = 0;
                else
                    (void)include_row_den(&row_scale[i], q);
            }
        }
    }
    for (size_t j = 0; j < m->cols && ok; ++j) {
        if ((j & 1023U) == 0 && lab_now() >= face_dl) {
            ok = 0;
            break;
        }
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            if (row_scale[i]) {
                if (lab_source_exact_coefficient(m, k, q))
                    row_scale[i] = 0;
                else
                    (void)include_row_den(&row_scale[i], q);
            }
        }
    }
    for (size_t j = 0; j < m->cols && ok; ++j) {
        if ((j & 1023U) == 0 && lab_now() >= face_dl) {
            ok = 0;
            break;
        }
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            size_t pos = next[i]++;
            rcol[pos] = j;
            rval[pos] = m->value[k];
            int64_t w = 0;
            if (row_scale[i] &&
                (lab_source_exact_coefficient(m, k, q) ||
                 !convert_scaled(q, m->value[k], row_scale[i], &w)))
                row_scale[i] = 0;
            rcoef[pos] = w;
        }
    }

    /* Propagate tight row faces to fix implied penalty slacks. */
    int infeasible = 0;
    for (int pass = 0; pass < 8 && ok && !infeasible; ++pass) {
        int changed = 0;
        for (size_t i = 0; i < m->rows; ++i) {
            if ((i & 1023U) == 0 && lab_now() >= face_dl) {
                ok = 0;
                break;
            }
            double min_a = 0.0, max_a = 0.0;
            int inf_min = 0, inf_max = 0;
            for (size_t p = start[i]; p < start[i + 1]; ++p) {
                size_t j = rcol[p];
                double a = rval[p];
                if (a > 0) {
                    if (isfinite(clo[j]))
                        min_a += a * clo[j];
                    else
                        inf_min = 1;
                    if (isfinite(chi[j]))
                        max_a += a * chi[j];
                    else
                        inf_max = 1;
                } else if (a < 0) {
                    if (isfinite(chi[j]))
                        min_a += a * chi[j];
                    else
                        inf_min = 1;
                    if (isfinite(clo[j]))
                        max_a += a * clo[j];
                    else
                        inf_max = 1;
                }
            }
            if (!inf_min && isfinite(m->row_upper[i])) {
                if (min_a > m->row_upper[i] + 1e-6) {
                    infeasible = 1;
                    break;
                }
                if (fabs(min_a - m->row_upper[i]) <= 1e-9) {
                    for (size_t p = start[i]; p < start[i + 1]; ++p) {
                        size_t j = rcol[p];
                        if (clo[j] == chi[j])
                            continue;
                        if (rval[p] > 0 && isfinite(clo[j])) {
                            chi[j] = clo[j];
                            changed = 1;
                        } else if (rval[p] < 0 && isfinite(chi[j])) {
                            clo[j] = chi[j];
                            changed = 1;
                        }
                    }
                }
            }
            if (!inf_max && isfinite(m->row_lower[i])) {
                if (max_a < m->row_lower[i] - 1e-6) {
                    infeasible = 1;
                    break;
                }
                if (fabs(max_a - m->row_lower[i]) <= 1e-9) {
                    for (size_t p = start[i]; p < start[i + 1]; ++p) {
                        size_t j = rcol[p];
                        if (clo[j] == chi[j])
                            continue;
                        if (rval[p] > 0 && isfinite(chi[j])) {
                            clo[j] = chi[j];
                            changed = 1;
                        } else if (rval[p] < 0 && isfinite(clo[j])) {
                            chi[j] = clo[j];
                            changed = 1;
                        }
                    }
                }
            }
        }
        if (!changed)
            break;
    }

    /* Map remaining unfixed binary variables to compact SAT IDs 1..nfree. */
    int *var_id = lab_alloc(m->cols, sizeof *var_id);
    size_t nfree = 0;
    for (size_t j = 0; j < m->cols && ok && !infeasible; ++j) {
        if (clo[j] == chi[j]) {
            if (!isfinite(clo[j]) || trunc(clo[j]) != clo[j] || fabs(clo[j]) > MAX_WEIGHT)
                ok = 0;
        } else if (clo[j] == 0.0 && chi[j] == 1.0) {
            var_id[j] = (int)(++nfree);
        } else {
            ok = 0;
        }
    }

    int solved = 0;
    if (ok && !infeasible && nfree > 0 && nfree <= 524288) {
        Encoder e = {.solver = kissat_init(), .variables = (int)nfree, .base_vars = (int)nfree,
                     .deadline = face_dl};
        kissat_reserve(e.solver, e.variables);
        Term *row = lab_alloc(max_row + 1, sizeof *row);
        int enc_ok = 1;
        for (size_t i = 0; i < m->rows && enc_ok; ++i) {
            if ((i & 255U) == 0 && (e.clauses > 2000000 || stopped(&e))) {
                enc_ok = 0;
                break;
            }
            if (!row_scale[i])
                continue;
            for (int upper = 0; upper < 2 && enc_ok; ++upper) {
                double rhs = upper ? m->row_upper[i] : m->row_lower[i];
                if (!isfinite(rhs))
                    continue;
                int64_t cap = 0;
                if (lab_source_exact_row_bound(m, i, upper, q) ||
                    !convert_scaled(q, rhs, row_scale[i], &cap))
                    continue;
                if (!upper)
                    cap = -cap;
                size_t cnt = 0;
                int64_t min_s = 0, max_s = 0;
                int row_valid = 1;
                for (size_t p = start[i]; p < start[i + 1]; ++p) {
                    size_t j = rcol[p];
                    int64_t w = upper ? rcoef[p] : -rcoef[p];
                    if (!var_id[j]) {
                        int64_t fx = (int64_t)clo[j];
                        if (fx && w && llabs(w) > (INT64_MAX / 4) / llabs(fx)) {
                            row_valid = 0;
                            break;
                        }
                        cap -= fx * w;
                    } else if (w) {
                        row[cnt++] = (Term){var_id[j], w, 0};
                        if (w < 0)
                            min_s += w;
                        else
                            max_s += w;
                    }
                }
                if (!row_valid)
                    continue;
                if (min_s > cap) {
                    enc_ok = 0;
                    break;
                }
                if (max_s > cap && e.clauses < 1900000)
                    (void)at_most(&e, row, (int)cnt, cap);
            }
        }
        free(row);
        if (enc_ok) {
            double *bits = lab_alloc(nfree, sizeof *bits);
            kissat_set_decision_limit(e.solver, (unsigned)fmin(UINT_MAX, (double)budget));
            ++*calls;
            if (bounded_solve(&e, nfree, bits) == 10) {
                for (size_t j = 0; j < m->cols; ++j)
                    point[j] = var_id[j] ? bits[var_id[j] - 1] : clo[j];
                LabModel view = *m;
                view.col_lower = (double *)lo;
                view.col_upper = (double *)hi;
                double eval_obj = 0, eval_viol = 0;
                if (!lab_exact_evaluate(&view, point, &eval_obj, &eval_viol) && eval_viol <= 1e-7)
                    solved = 1;
            }
            free(bits);
        }
        kissat_release(e.solver);
    }

    free(var_id);
    free(row_scale);
    free(rcoef);
    free(rval);
    free(rcol);
    free(next);
    free(start);
    free(clo);
    free(chi);
    mpq_clear(q);
    return solved && lab_now() < deadline;
}

int lab_binary_sat(const LabModel *m, const double *lo, const double *hi, double deadline,
                   size_t budget, double *point, size_t *calls) {
    *calls = 0;
    if (!m->cols || m->cols > MAX_TERMS || m->rows > MAX_TERMS || m->nnz > MAX_TERMS || !budget ||
        lab_now() >= deadline)
        return 0;
    int all_int = 1, has_nonbinary = 0;
    for (size_t j = 0; j < m->cols; ++j) {
        if (!m->integer[j]) {
            all_int = 0;
            break;
        }
        double l = m->col_lower[j], u = m->col_upper[j];
        if ((l != 0 && l != 1) || (u != 0 && u != 1) || (lo[j] != 0 && lo[j] != 1) ||
            (hi[j] != 0 && hi[j] != 1))
            has_nonbinary = 1;
    }
    if (!all_int)
        return 0;
    if (has_nonbinary) {
        if (cheapest_face_general_int(m, lo, hi, deadline, budget, point, calls))
            return 1;
        return small_integer_sat(m, lo, hi, deadline, budget, point, calls);
    }
    Objective cost = {0};
    cost.pb_exact = 1;
    mpq_init(cost.cutoff);
    mpq_t q;
    mpq_init(q);
    size_t *start = NULL, *next = NULL;
    Term *terms = NULL, *row = NULL;
    int64_t *row_scale = NULL;
    unsigned char *interior = NULL;
    double *trial = NULL;
    Region *regions = NULL;
    int result = 0, incumbent = 0, face_only = 0, try_face = 0;
    signed char *fixed = lab_alloc(m->cols, sizeof *fixed);
    int const_obj = constant_objective_via_equalities(m, fmin(deadline, lab_now() + .5));
    for (size_t j = 0; j < m->cols; ++j) {
        double l = m->col_lower[j], u = m->col_upper[j];
        if (!m->integer[j] || (l != 0 && l != 1) || (u != 0 && u != 1) ||
            (lo[j] != 0 && lo[j] != 1) || (hi[j] != 0 && hi[j] != 1) ||
            l > lo[j] || u < hi[j] || lo[j] > hi[j] ||
            lab_source_exact_column_bound(m, j, 0, q) || mpq_cmp_ui(q, (unsigned long)l, 1) ||
            lab_source_exact_column_bound(m, j, 1, q) || mpq_cmp_ui(q, (unsigned long)u, 1) ||
            lab_source_exact_objective(m, j, q) || (!(j % 256) && lab_now() >= deadline))
            goto done;
        int sign = const_obj ? 0 : mpq_sgn(q);
        if (m->maximize)
            sign = -sign;
        fixed[j] = lo[j] == hi[j] ? (signed char)lo[j] : sign > 0 ? 0 : sign < 0 ? 1 : -1;
        if (!const_obj && (l == u || (sign && cost.n == SHORT_ROW)))
            try_face = 1;
        if (sign && cost.n == MAX_OBJ_TERMS)
            face_only = 1;
        if (!sign || face_only || lo[j] == hi[j])
            continue;
        int lit = (int)j + 1;
        if ((mpq_sgn(q) < 0) != (m->maximize != 0))
            lit = -lit;
        cost.lit[cost.n] = lit;
        mpq_init(cost.weight[cost.n]);
        mpq_abs(cost.weight[cost.n], q);
        int64_t sw = 0;
        if (cost.pb_exact && convert(cost.weight[cost.n], fabs(m->objective[j]), &sw) && sw > 0) {
            cost.scaled[cost.n] = sw;
            cost.scaled_gcd = gcd64(cost.scaled_gcd, sw);
            if (cost.scaled_total > INT64_MAX - sw)
                cost.pb_exact = 0;
            else
                cost.scaled_total += sw;
        } else {
            cost.pb_exact = 0;
        }
        ++cost.n;
    }
    if (cost.n > SHORT_ROW &&
        (!cost.pb_exact || !cost.scaled_gcd || cost.scaled_total / cost.scaled_gcd > 512))
        face_only = 1;
    double orig_deadline = deadline;
    if (face_only)
        deadline = fmin(deadline, lab_now() + fmin(3.0, .25 * (deadline - lab_now())));
    start = lab_alloc(m->rows + 1, sizeof *start);
    for (size_t k = 0; k < m->nnz; ++k)
        ++start[m->row_index[k] + 1];
    size_t width = 0;
    for (size_t i = 0; i < m->rows; ++i) {
        if (start[i + 1] > width)
            width = start[i + 1];
        start[i + 1] += start[i];
    }
    next = lab_alloc(m->rows, sizeof *next);
    memcpy(next, start, m->rows * sizeof *next);
    terms = lab_alloc(m->nnz + 1, sizeof *terms);
    row = lab_alloc(width + 1, sizeof *row);
    row_scale = lab_alloc(m->rows ? m->rows : 1, sizeof *row_scale);
    for (size_t i = 0; i < m->rows; ++i) {
        row_scale[i] = 1;
        for (int upper = 0; upper < 2; ++upper) {
            double rhs = upper ? m->row_upper[i] : m->row_lower[i];
            if (isfinite(rhs)) {
                if (lab_source_exact_row_bound(m, i, upper, q))
                    row_scale[i] = 0;
                else
                    (void)include_row_den(&row_scale[i], q);
            }
        }
    }
    for (size_t j = 0; j < m->cols; ++j) {
        if (!(j % 256) && lab_now() >= deadline)
            goto done;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            if (row_scale[i]) {
                if (lab_source_exact_coefficient(m, k, q))
                    row_scale[i] = 0;
                else
                    (void)include_row_den(&row_scale[i], q);
            }
        }
    }
    for (size_t j = 0; j < m->cols; ++j) {
        if (!(j % 256) && lab_now() >= deadline)
            goto done;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            int64_t w = 0;
            if (row_scale[i] &&
                (lab_source_exact_coefficient(m, k, q) ||
                 !convert_scaled(q, m->value[k], row_scale[i], &w)))
                w = 0;
            terms[next[i]++] = (Term){(int)j + 1, w, k};
        }
    }
    if (try_face || face_only) {
        double face_dl = face_only
                             ? fmin(orig_deadline, lab_now() + fmin(8.0, .35 * (orig_deadline - lab_now())))
                             : fmin(orig_deadline, lab_now() + fmin(3.0, .25 * (orig_deadline - lab_now())));
        int faced = cheapest_face(m, start, terms, row_scale, row, fixed, face_dl, budget, point, calls);
        if (faced) {
            result = 1;
            goto done;
        }
        if (face_only) {
            /* Before giving up on a large-objective pure-binary model, run one
             * fast feasibility check; UNSAT on a relaxation proves infeasibility. */
            double feas_dl = fmin(orig_deadline, lab_now() + fmin(2.0, .20 * (orig_deadline - lab_now())));
            if (lab_now() < feas_dl && budget) {
                Encoder ef = {.solver = kissat_init(), .variables = (int)m->cols,
                              .base_vars = (int)m->cols, .deadline = feas_dl};
                kissat_reserve(ef.solver, ef.variables);
                for (size_t j = 0; j < m->cols; ++j)
                    if (lo[j] == hi[j]) {
                        literal(&ef, lo[j] == 1 ? (int)j + 1 : -((int)j + 1));
                        literal(&ef, 0);
                    }
                int skipped_rows = 0;
                if (encode(m, &ef, start, terms, row_scale, row, &cost, 0, -1, NULL, m->rows,
                           &skipped_rows)) {
                    kissat_set_decision_limit(ef.solver, (unsigned)fmin(UINT_MAX, (double)budget));
                    ++*calls;
                    int st = bounded_solve(&ef, m->cols, NULL);
                    if (st == 20)
                        result = -1;
                }
                kissat_release(ef.solver);
            }
            goto done;
        }
    }
    size_t seed = 0, regions_count = 0;
    double local_deadline = 0;
    interior = lab_alloc(m->cols, 1);
    trial = lab_alloc(m->cols, sizeof *trial);
    int pb_mode = cost.n > SHORT_ROW && cost.pb_exact && cost.scaled_gcd > 0;
    int64_t pb_lo = 0, pb_hi = pb_mode ? cost.scaled_total / cost.scaled_gcd : 0;
    while (lab_now() < deadline && budget) {
        Encoder e = {.solver = kissat_init(), .variables = (int)m->cols, .base_vars = (int)m->cols,
                     .deadline = deadline};
        size_t indices[3000], count = m->rows;
        int local = !pb_mode && incumbent && cost.n > 1 && m->cols > 512 && seed < regions_count &&
                    lab_now() < local_deadline && budget > 1000;
        int64_t pb_probe = -1;
        if (pb_mode && incumbent) {
            if (pb_lo >= pb_hi) {
                kissat_release(e.solver);
                result = 1;
                break;
            }
            pb_probe = pb_hi - pb_lo <= 2 ? pb_hi - 1 : pb_lo + (pb_hi - pb_lo) / 2;
        }
        if (local) {
            uint64_t shape;
            count = neighborhood(m, start, terms, regions[seed++].seed, interior, indices, 3000,
                                 &shape);
            e.deadline = fmin(local_deadline, lab_now() + .15);
            e.map = lab_alloc(m->cols, sizeof *e.map);
        } else {
            kissat_reserve(e.solver, e.variables);
            for (size_t j = 0; j < m->cols; ++j) {
                if (lo[j] == hi[j]) {
                    literal(&e, lo[j] == 1 ? (int)j + 1 : -((int)j + 1));
                    literal(&e, 0);
                }
            }
        }
        int skipped_rows = 0;
        int allow_skip = !incumbent && !local;
        int encoded = encode(m, &e, start, terms, row_scale, row, &cost, incumbent,
                             pb_probe >= 0 ? pb_probe * cost.scaled_gcd : -1,
                             local ? indices : NULL, count, allow_skip ? &skipped_rows : NULL);
        free(e.map);
        if (!encoded) {
            kissat_release(e.solver);
            if (local)
                continue;
            break;
        }
        /* Reserve a disjoint part of the decision budget for each call. */
        size_t want = pb_mode ? (budget > 8 ? budget / 4 : budget) : (budget / 2 + budget % 2);
        unsigned slice = (unsigned)fmin(UINT_MAX, (double)want);
        if (local && slice > 1024)
            slice = 1024;
        budget -= slice;
        kissat_set_decision_limit(e.solver, slice);
        ++*calls;
        int status = bounded_solve(&e, m->cols, local ? NULL : trial);
        if (local && status != 20) {
            kissat_release(e.solver);
            continue;
        }
        if (status == 10) {
            if (skipped_rows > 0) {
                /* SAT on a relaxed encoding does not verify skipped rows. */
                kissat_release(e.solver);
                break;
            }
            memcpy(point, trial, m->cols * sizeof *point);
            mpq_set_ui(cost.cutoff, 0, 1);
            int64_t cur_scaled = 0;
            for (int k = 0; k < cost.n; ++k)
                if ((point[abs(cost.lit[k]) - 1] != 0) == (cost.lit[k] > 0)) {
                    mpq_add(cost.cutoff, cost.cutoff, cost.weight[k]);
                    if (pb_mode)
                        cur_scaled += cost.scaled[k];
                }
            incumbent = 1;
            if (pb_mode) {
                pb_hi = cur_scaled / cost.scaled_gcd;
                if (pb_hi <= pb_lo)
                    result = 1;
            } else {
                seed = 0;
                local_deadline = lab_now() + .6 * (deadline - lab_now());
                if (!mpq_sgn(cost.cutoff))
                    result = 1;
                else if (cost.n > 1 && m->cols > 512 && !regions) {
                    regions = lab_alloc(4096, sizeof *regions);
                    regions_count =
                        rank_regions(m, start, terms, interior, regions, 4096, local_deadline);
                }
            }
        } else if (status == 20) {
            if (pb_mode && incumbent) {
                pb_lo = pb_probe + 1;
                if (pb_lo >= pb_hi)
                    result = 1;
            } else {
                result = incumbent ? 1 : -1;
            }
        }
        kissat_release(e.solver);
        if (result || (status != 10 && !(pb_mode && incumbent && status == 20)))
            break;
    }
done:
    for (int k = 0; k < cost.n; ++k)
        mpq_clear(cost.weight[k]);
    mpq_clear(cost.cutoff);
    mpq_clear(q);
    free(fixed);
    free(regions);
    free(interior);
    free(trial);
    free(start);
    free(next);
    free(terms);
    free(row_scale);
    free(row);
    if (result && lab_now() < deadline)
        return result;
    return incumbent ? 2 : 0;
}

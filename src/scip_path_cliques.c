#include "scip_path_cliques.h"
#include "common.h"
#include "exact.h"
#include "scip/cons_linear.h"
#include <math.h>
#include <stdlib.h>

enum { MAX_ARCS = 4096, MAX_DEGREE = 96, MAX_CUTS = 1024, MAX_CUT_TERMS = 32000 };

typedef struct {
    size_t first, second;
} Membership;

typedef struct {
    size_t choice, to;
    int next_in, next_out;
    mpq_t earliest, latest;
} ConditionalArc;

static int shares(Membership a, Membership b) {
    return (a.first && (a.first == b.first || a.first == b.second)) ||
           (a.second && (a.second == b.first || a.second == b.second));
}

/* The selected columns share one source-checked at-most-one inequality. */
static int common_cover(const Membership *members, const size_t *columns, int count) {
    if (count <= 1)
        return 1;
    Membership possible = members[columns[0]];
    for (int k = 1; k < count; ++k) {
        Membership m = members[columns[k]];
        if (possible.first != m.first && possible.first != m.second)
            possible.first = 0;
        if (possible.second != m.first && possible.second != m.second)
            possible.second = 0;
    }
    return possible.first || possible.second;
}

static int directed_row(const LabModel *m, const int *col, const int *entry, size_t n) {
    int positive = 0, negative = 0;
    for (size_t k = 0; k < n; ++k) {
        int e = entry[k];
        if (!m->integer[col[e]]) {
            positive += m->value[e] == 1;
            negative += m->value[e] == -1;
            if (fabs(m->value[e]) != 1)
                return 0;
        }
    }
    return positive == 1 && negative == 1;
}

static void remember_covers(const LabModel *m, const size_t *start, const int *entry,
                            const int *col, const double *lower, const double *upper,
                            const double *ru, Membership *members, double deadline) {
    mpq_t q;
    mpq_init(q);
    for (size_t i = 0; i < m->rows; ++i) {
        if (!(i % 128) && lab_now() >= deadline)
            break;
        if (ru[i] != 1 || m->row_upper[i] != 1 || start[i + 1] - start[i] < 2 ||
            lab_source_exact_row_bound(m, i, 1, q) || mpq_cmp_ui(q, 1, 1))
            continue;
        int valid = 1;
        for (size_t p = start[i]; p < start[i + 1]; ++p) {
            int e = entry[p], j = col[e];
            if ((!(p % 256) && lab_now() >= deadline) || !m->integer[j] || lower[j] < 0 ||
                m->value[e] != 1 || lab_source_exact_column_bound(m, (size_t)j, 0, q) ||
                mpq_sgn(q) < 0 || lab_source_exact_coefficient(m, (size_t)e, q) ||
                mpq_cmp_ui(q, 1, 1)) {
                valid = 0;
                break;
            }
        }
        if (valid)
            for (size_t p = start[i]; p < start[i + 1]; ++p) {
                size_t j = (size_t)col[entry[p]];
                if (upper[j] < 1)
                    continue;
                if (!members[j].first)
                    members[j].first = i + 1;
                else if (!members[j].second)
                    members[j].second = i + 1;
            }
    }
    mpq_clear(q);
}

/* With the trigger set to one, upper-bound the contribution of every
 * remaining discrete term. An edge in the same cover must be zero. */
static int read_arc(const LabModel *m, size_t i, int hi, size_t trigger, const int *entries,
                    size_t n, const int *col, const Membership *members, const double *lower,
                    const double *upper, ConditionalArc *arc, size_t *origin) {
    mpq_t lag, q, b;
    mpq_inits(lag, q, b, NULL);
    size_t from = m->cols, to = m->cols;
    int ok = 0;
    if (lab_source_exact_row_bound(m, i, hi, lag))
        goto done;
    if (hi)
        mpq_neg(lag, lag);
    for (size_t k = 0; k < n; ++k) {
        size_t e = (size_t)entries[k], j = (size_t)col[e];
        if (lab_source_exact_coefficient(m, e, q))
            goto done;
        if (hi)
            mpq_neg(q, q);
        if (j == trigger) {
            if (mpq_sgn(q) >= 0)
                goto done;
            mpq_sub(lag, lag, q);
        } else if (m->integer[j]) {
            if (shares(members[j], members[trigger]))
                continue;
            int high = mpq_sgn(q) > 0;
            double bound = high ? upper[j] : lower[j];
            if (lab_source_exact_column_value(m, j, high, bound, b))
                goto done;
            mpq_mul(q, q, b);
            mpq_sub(lag, lag, q);
        } else if (mpq_cmp_si(q, -1, 1) == 0)
            from = j;
        else if (mpq_cmp_ui(q, 1, 1) == 0)
            to = j;
        else
            goto done;
    }
    if (from == m->cols || to == m->cols ||
        lab_source_exact_column_value(m, from, 0, lower[from], q) ||
        lab_source_exact_column_value(m, to, 1, upper[to], b))
        goto done;
    *origin = from;
    arc->choice = trigger;
    arc->to = to;
    mpq_inits(arc->earliest, arc->latest, NULL);
    mpq_add(arc->earliest, q, lag);
    mpq_sub(arc->latest, b, lag);
    ok = 1;
done:
    mpq_clears(lag, q, b, NULL);
    return ok;
}

static SCIP_RETCODE emit_cliques(SCIP *scip, SCIP_VAR **vars, const Membership *members,
                                 const ConditionalArc *arc, int count, const int *first_in,
                                 const int *first_out, unsigned char *seen, double deadline,
                                 size_t *added) {
    size_t columns[2 * MAX_DEGREE], total = 0;
    SCIP_VAR *selected[2 * MAX_DEGREE];
    double ones[2 * MAX_DEGREE];
    for (int k = 0; k < 2 * MAX_DEGREE; ++k)
        ones[k] = 1;
    for (int e = 1; e <= count && *added < MAX_CUTS && total < MAX_CUT_TERMS; ++e) {
        if (lab_now() >= deadline)
            break;
        size_t v = arc[e].to;
        int ni = 0, no = 0, duplicate = 0, degree = 0;
        for (int h = first_in[v]; h && degree <= MAX_DEGREE; h = arc[h].next_in) {
            ++degree;
            int cmp = mpq_cmp(arc[h].earliest, arc[e].earliest);
            if (h < e && !cmp)
                duplicate = 1;
            if (cmp >= 0 && ni < MAX_DEGREE)
                columns[ni++] = arc[h].choice;
        }
        if (degree > MAX_DEGREE || duplicate || !common_cover(members, columns, ni))
            continue;
        degree = 0;
        for (int h = first_out[v]; h && degree <= MAX_DEGREE; h = arc[h].next_out) {
            ++degree;
            if (mpq_cmp(arc[h].latest, arc[e].earliest) < 0 && no < MAX_DEGREE)
                columns[ni + no++] = arc[h].choice;
        }
        if (degree > MAX_DEGREE || !no || !common_cover(members, columns + ni, no))
            continue;
        for (int k = 0; k < ni + no; ++k) {
            duplicate |= seen[columns[k]];
            seen[columns[k]] = 1;
            selected[k] = vars[columns[k]];
        }
        for (int k = 0; k < ni + no; ++k)
            seen[columns[k]] = 0;
        if (duplicate)
            continue;
        SCIP_CONS *cons = NULL;
        char name[48];
        snprintf(name, sizeof name, "conditional_path_%zu", *added);
        SCIP_RETCODE rc = SCIPcreateConsBasicLinear(scip, &cons, name, ni + no, selected, ones,
                                                    -SCIPinfinity(scip), 1);
        if (rc != SCIP_OKAY)
            return rc;
        rc = SCIPaddCons(scip, cons);
        SCIP_RETCODE released = SCIPreleaseCons(scip, &cons);
        if (rc != SCIP_OKAY || released != SCIP_OKAY)
            return rc != SCIP_OKAY ? rc : released;
        ++*added;
        total += (size_t)(ni + no);
    }
    return SCIP_OKAY;
}

SCIP_RETCODE lab_scip_path_cliques(SCIP *scip, const LabModel *m, SCIP_VAR **vars,
                                   const size_t *start, const int *entry, const double *lower,
                                   const double *upper, const double *rl, const double *ru,
                                   double deadline, size_t *added) {
    *added = 0;
    if (m->nnz > 1000000 || m->cols > 150000 || lab_now() >= deadline)
        return SCIP_OKAY;
    int potentials = 0;
    for (size_t j = 0; j < m->cols && potentials < 2; ++j)
        potentials += !m->integer[j];
    if (potentials < 2)
        return SCIP_OKAY;
    size_t count_rows = 0;
    for (size_t i = 0; i < m->rows; ++i) {
        size_t n = start[i + 1] - start[i];
        if (n == 3 || n == 4)
            ++count_rows;
    }
    if (count_rows < 2)
        return SCIP_OKAY;
    int *col = lab_alloc(m->nnz, sizeof *col);
    size_t eligible = 0;
    /* A skipped zero in the imported matrix may hide a source-decimal
     * residual. Without that term the derived implication would be unsound. */
    for (size_t j = 0; j < m->cols; ++j) {
        if (!(j % 256) && lab_now() >= deadline)
            goto finish;
        for (size_t e = m->col_start[j]; e < m->col_start[j + 1]; ++e) {
            if (m->value[e] == 0)
                goto finish;
            col[e] = (int)j;
        }
    }
    for (size_t i = 0; i < m->rows && eligible < 2; ++i) {
        size_t n = start[i + 1] - start[i];
        if ((n == 3 || n == 4) && directed_row(m, col, entry + start[i], n))
            ++eligible;
    }
    if (eligible < 2 || lab_now() >= deadline)
        goto finish;
    Membership *members = lab_alloc(m->cols, sizeof *members);
    ConditionalArc *arcs = lab_alloc(MAX_ARCS + 1, sizeof *arcs);
    int *first_in = lab_alloc(m->cols, sizeof *first_in),
        *first_out = lab_alloc(m->cols, sizeof *first_out);
    int count = 0;
    remember_covers(m, start, entry, col, lower, upper, ru, members, deadline);
    for (size_t i = 0; i < m->rows && count < MAX_ARCS; ++i) {
        if (!(i % 32) && lab_now() >= deadline)
            break;
        size_t n = start[i + 1] - start[i];
        const int *entries = entry + start[i];
        if ((n != 3 && n != 4) || !directed_row(m, col, entries, n))
            continue;
        for (int hi = 0; hi < 2; ++hi) {
            double side = hi ? ru[i] : rl[i];
            if (!isfinite(side) || side != (hi ? m->row_upper[i] : m->row_lower[i]))
                continue;
            for (size_t k = 0; k < n && count < MAX_ARCS; ++k) {
                size_t j = (size_t)col[entries[k]], origin = 0;
                if (!m->integer[j] || lower[j] != 0 || upper[j] != 1 ||
                    (hi ? m->value[entries[k]] <= 0 : m->value[entries[k]] >= 0))
                    continue;
                if (read_arc(m, i, hi, j, entries, n, col, members, lower, upper, &arcs[count + 1],
                             &origin)) {
                    ++count;
                    size_t to = arcs[count].to;
                    arcs[count].next_in = first_in[to];
                    first_in[to] = count;
                    arcs[count].next_out = first_out[origin];
                    first_out[origin] = count;
                }
            }
        }
    }
    unsigned char *seen = lab_alloc(m->cols, sizeof *seen);
    SCIP_RETCODE rc =
        emit_cliques(scip, vars, members, arcs, count, first_in, first_out, seen, deadline, added);
    for (int e = 1; e <= count; ++e)
        mpq_clears(arcs[e].earliest, arcs[e].latest, NULL);
    free(seen);
    free(arcs);
    free(members);
    free(first_in);
    free(first_out);
    free(col);
    return rc;
finish:
    free(col);
    return SCIP_OKAY;
}

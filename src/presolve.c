#include "presolve.h"
#include "common.h"
#include "exact.h"
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static LabInterval point(double x) {
    return (LabInterval){x, x};
}
static int valid_interval(LabInterval x) {
    return !isnan(x.lower) && !isnan(x.upper) && x.lower <= x.upper;
}
static LabInterval rational_interval(const mpq_t q) {
    double x = mpq_get_d(q);
    if (isinf(x))
        return x > 0 ? (LabInterval){DBL_MAX, INFINITY}
                     : (LabInterval){-INFINITY, -DBL_MAX};
    mpq_t binary;
    mpq_init(binary);
    mpq_set_d(binary, x);
    int sign = mpq_cmp(binary, q);
    mpq_clear(binary);
    return (LabInterval){sign > 0 ? nextafter(x, -INFINITY) : x,
                         sign < 0 ? nextafter(x, INFINITY) : x};
}
static void diagnostic(LabPresolveResult *r, const char *status, const char *message) {
    r->status = status;
    snprintf(r->message, sizeof r->message, "%s", message);
}
static int source_data(const LabModel *m, LabInterval *a, LabInterval *rl,
                        LabInterval *ru, LabInterval *box) {
    const LabExactModel *e = m->exact;
    const LabSourceEnvelope *s = m->envelope;
    if (e && (e->rows != m->rows || e->cols != m->cols || e->nnz != m->nnz))
        return -1;
    int envelope = s && s->complete && s->rows == m->rows &&
                   s->cols == m->cols && s->nnz == m->nnz;
    if (!e && s && !envelope)
        return s->complete ? -1 : 0;
    if ((e && ((m->nnz && !e->value) || (m->rows && (!e->row_lower || !e->row_upper)) ||
                (m->cols && (!e->col_lower || !e->col_upper || !e->col_lower_source ||
                             !e->col_upper_source)))) ||
        (envelope && ((m->nnz && !s->value) || (m->rows && (!s->row_lower || !s->row_upper)) ||
                      (m->cols && (!s->col_lower || !s->col_upper || !s->col_lower_source ||
                                   !s->col_upper_source)))))
        return -1;
    if (!e && !envelope && m->parsed_decimal)
        return 0;
    for (size_t k = 0; k < m->nnz; k++)
        a[k] = e ? rational_interval(e->value[k])
                 : envelope ? s->value[k] : point(m->value[k]);
    for (size_t i = 0; i < m->rows; i++) {
        rl[i] = e && isfinite(m->row_lower[i]) ? rational_interval(e->row_lower[i])
                : envelope ? s->row_lower[i] : point(m->row_lower[i]);
        ru[i] = e && isfinite(m->row_upper[i]) ? rational_interval(e->row_upper[i])
                : envelope ? s->row_upper[i] : point(m->row_upper[i]);
    }
    for (size_t j = 0; j < m->cols; j++) {
        double l = m->col_lower[j], u = m->col_upper[j];
        LabInterval li = e ? (isfinite(e->col_lower_source[j])
                              ? rational_interval(e->col_lower[j])
                              : point(e->col_lower_source[j]))
                         : envelope ? s->col_lower[j] : point(l);
        LabInterval ui = e ? (isfinite(e->col_upper_source[j])
                              ? rational_interval(e->col_upper[j])
                              : point(e->col_upper_source[j]))
                         : envelope ? s->col_upper[j] : point(u);
        box[j] = (LabInterval){li.lower, ui.upper};
    }
    for (size_t k = 0; k < m->nnz; k++)
        if (!valid_interval(a[k]) || a[k].lower == INFINITY || a[k].upper == -INFINITY)
            return -1;
    for (size_t i = 0; i < m->rows; i++)
        if (!valid_interval(rl[i]) || !valid_interval(ru[i]))
            return -1;
    return 1;
}
static void overlay_source_box(const LabModel *m, const double *lower, const double *upper,
                                LabInterval *box, LabPresolveResult *r, int relax) {
    for (size_t j = 0; j < m->cols; j++) {
        double source_l = m->col_lower[j], source_u = m->col_upper[j];
        if (m->exact) {
            source_l = m->exact->col_lower_source[j];
            source_u = m->exact->col_upper_source[j];
        } else if (m->envelope) {
            source_l = m->envelope->col_lower_source[j];
            source_u = m->envelope->col_upper_source[j];
        }
        if (lower[j] != source_l)
            box[j].lower = fmax(box[j].lower, lower[j]);
        if (upper[j] != source_u)
            box[j].upper = fmin(box[j].upper, upper[j]);
        if (m->integer && m->integer[j] && !relax) {
            double old_lower = box[j].lower, old_upper = box[j].upper;
            if (isfinite(box[j].lower))
                box[j].lower = ceil(box[j].lower);
            if (isfinite(box[j].upper))
                box[j].upper = floor(box[j].upper);
            r->tightenings += (size_t)(old_lower < box[j].lower) +
                              (size_t)(old_upper > box[j].upper);
        }
    }
}
typedef struct {
    LabInterval finite_min, finite_max;
    size_t infinite_min, infinite_max;
} Activity;
static int row_activity(size_t begin, size_t end, const size_t *entries,
                         const size_t *column, const LabInterval *a,
                         const LabInterval *box, LabInterval *term, Activity *activity,
                         double deadline) {
    memset(activity, 0, sizeof *activity);
    for (size_t p = begin; p < end; p++) {
        size_t k = entries[p];
        LabInterval t = lab_interval_multiply(a[k], box[column[k]]);
        term[p - begin] = t;
        if (isfinite(t.lower))
            activity->finite_min = lab_interval_add(activity->finite_min, point(t.lower));
        else
            activity->infinite_min++;
        if (isfinite(t.upper))
            activity->finite_max = lab_interval_add(activity->finite_max, point(t.upper));
        else
            activity->infinite_max++;
        if ((p & 1023U) == 0 && lab_now() >= deadline)
            return -1;
    }
    return 0;
}
static int tighten(LabInterval *box, double candidate, int upper, int integer,
                     LabPresolveResult *r) {
    if (isnan(candidate))
        return 0;
    if (integer && isfinite(candidate))
        candidate = upper ? floor(candidate) : ceil(candidate);
    if (upper) {
        if (candidate >= box->upper || candidate == -INFINITY)
            return 0;
        box->upper = candidate;
    } else {
        if (candidate <= box->lower || candidate == INFINITY)
            return 0;
        box->lower = candidate;
    }
    r->tightenings++;
    return 1;
}
#include "clique_probing.inc"
/* Reduce exact integer equalities modulo two. Only point source enclosures
 * are eligible: a decimal that merely rounds to an integer has no parity.
 * Limit packed storage and row-reduction work on large sparse LPs. */
static int integer_point(LabInterval a) {
    return a.lower == a.upper && isfinite(a.lower) && trunc(a.lower) == a.lower;
}
static int odd(double integer) {
    return fmod(fabs(integer), 2.0) == 1.0;
}
static int propagate_parity(const LabModel *m, const LabInterval *a,
                            const LabInterval *rl, const LabInterval *ru,
                            const size_t *start, const size_t *entries, const size_t *column,
                            LabInterval *box, LabPresolveResult *r, double deadline) {
    size_t n = m->cols;
    if (!m->integer || !n || n > 1024 || !m->rows || m->rows > 4096 || m->nnz > 262144)
        return 0;
    size_t words = (n + 63) / 64;
    uint64_t *basis = lab_alloc(n * words, sizeof *basis);
    uint64_t *row = lab_alloc(words, sizeof *row);
    unsigned char *rhs = lab_alloc(n, 1);
    int conflict = 0;
    for (size_t i = 0; i < m->rows; i++) {
        if (lab_now() >= deadline)
            goto done;
        if (!integer_point(rl[i]) || !integer_point(ru[i]) || rl[i].lower != ru[i].lower)
            continue;
        unsigned char bit = (unsigned char)odd(rl[i].lower);
        int usable = 1;
        memset(row, 0, words * sizeof *row);
        for (size_t p = start[i]; p < start[i+1]; p++) {
            size_t k = entries[p], j = column[k];
            if (!integer_point(a[k]) || (a[k].lower != 0 && !m->integer[j])) {
                usable = 0;
                break;
            }
            if (odd(a[k].lower)) {
                if (integer_point(box[j]))
                    bit ^= (unsigned char)odd(box[j].lower);
                else
                    row[j/64] ^= UINT64_C(1) << (j%64);
            }
        }
        if (!usable)
            continue;
        size_t j;
        for (j = 0; j < n; j++) {
            uint64_t mask = UINT64_C(1) << (j%64);
            if (!(row[j/64] & mask))
                continue;
            uint64_t *pivot = basis + j * words;
            if (!(pivot[j/64] & mask)) {
                memcpy(pivot, row, words * sizeof *pivot);
                rhs[j] = bit;
                break;
            }
            for (size_t w = j/64; w < words; w++)
                row[w] ^= pivot[w];
            bit ^= rhs[j];
        }
        if (j == n && bit) {
            conflict = 1;
            goto done;
        }
    }
    /* Substitute determined parities into earlier equations. A row still
     * depending on another variable cannot justify an endpoint change. */
    for (size_t j = n; j-- > 0;) {
        if (lab_now() >= deadline)
            goto done;
        uint64_t mask = UINT64_C(1) << (j%64);
        const uint64_t *pivot = basis + j * words;
        if (pivot[j/64] != mask)
            continue;
        size_t w = j/64 + 1;
        while (w < words && !pivot[w])
            w++;
        if (w != words)
            continue;
        double lo = box[j].lower, hi = box[j].upper;
        if (fabs(lo) < 0x1p53 && odd(lo) != rhs[j]) {
            box[j].lower = lo + 1;
            r->tightenings++;
        }
        if (fabs(hi) < 0x1p53 && odd(hi) != rhs[j]) {
            box[j].upper = hi - 1;
            r->tightenings++;
        }
        if (box[j].lower > box[j].upper) {
            conflict = 1;
            goto done;
        }
        for (size_t i = 0; i < j; i++) {
            uint64_t *word = &basis[i * words + j/64];
            if (*word & mask) {
                *word ^= mask;
                rhs[i] ^= rhs[j];
            }
        }
    }
done:
    free(row);
    free(basis);
    free(rhs);
    return conflict;
}

void lab_presolve(const LabModel *m, const double *lower, const double *upper,
                  const LabPresolveOptions *options, LabPresolveResult *r) {
    memset(r, 0, sizeof *r);
    double started = lab_now();
    LabPresolveOptions opt = options ? *options : (LabPresolveOptions){1, 16, 0};
    diagnostic(r, "input_error", "Invalid presolve model, options or bound overlay");
    if (!m || !isfinite(opt.time_limit) || opt.time_limit <= 0 || !opt.max_passes ||
        fegetround() != FE_TONEAREST || (m->cols && (!m->col_start || !m->col_lower ||
        !m->col_upper)) || (m->rows && (!m->row_lower || !m->row_upper)) ||
        (m->nnz && (!m->row_index || !m->value)))
        goto done;
    if (!m->cols && m->nnz)
        goto done;
    if (!lower)
        lower = m->col_lower;
    if (!upper)
        upper = m->col_upper;
    if (m->cols && (m->col_start[0] || m->col_start[m->cols] != m->nnz))
        goto done;
    for (size_t j = 0; j < m->cols; j++)
        if (isnan(lower[j]) || isnan(upper[j]) || lower[j] == INFINITY ||
            upper[j] == -INFINITY || lower[j] < m->col_lower[j] ||
            upper[j] > m->col_upper[j] || isnan(m->col_lower[j]) ||
            isnan(m->col_upper[j]) || m->col_lower[j] == INFINITY ||
            m->col_upper[j] == -INFINITY || m->col_start[j] > m->col_start[j + 1])
            goto done;
    for (size_t i = 0; i < m->rows; i++)
        if (isnan(m->row_lower[i]) || isnan(m->row_upper[i]) ||
            m->row_lower[i] == INFINITY || m->row_upper[i] == -INFINITY)
            goto done;
    for (size_t k = 0; k < m->nnz; k++)
        if (m->row_index[k] >= m->rows || !isfinite(m->value[k]))
            goto done;
    r->lower = lab_alloc(m->cols, sizeof(double));
    r->upper = lab_alloc(m->cols, sizeof(double));
    r->redundant = lab_alloc(m->rows, 1);
    if (m->cols) {
        memcpy(r->lower, lower, m->cols * sizeof(double));
        memcpy(r->upper, upper, m->cols * sizeof(double));
    }
    diagnostic(r, "unresolved", "Source provenance unavailable; propagation abstained");
    LabInterval *a = lab_alloc(m->nnz, sizeof(LabInterval));
    LabInterval *rl = lab_alloc(m->rows, sizeof(LabInterval));
    LabInterval *ru = lab_alloc(m->rows, sizeof(LabInterval));
    LabInterval *box = lab_alloc(m->cols, sizeof(LabInterval));
    LabInterval *initial = lab_alloc(m->cols, sizeof(LabInterval));
    size_t *start = lab_alloc(m->rows + 1, sizeof(size_t));
    size_t *column = lab_alloc(m->nnz, sizeof(size_t));
    size_t *entries = lab_alloc(m->nnz, sizeof(size_t));
    size_t *at = lab_alloc(m->rows, sizeof(size_t));
    int source = source_data(m, a, rl, ru, box);
    if (source <= 0) {
        if (source < 0)
            diagnostic(r, "input_error", "Invalid source-envelope metadata");
        goto cleanup;
    }
    r->source_available = 1;
    overlay_source_box(m, lower, upper, box, r, opt.relax);
    if (m->cols)
        memcpy(initial, box, m->cols * sizeof(LabInterval));
    for (size_t j = 0; j < m->cols; j++)
        if (box[j].lower > box[j].upper) {
            diagnostic(r, "infeasible", "Source bounds and integrality yield a contradictory box");
            goto cleanup;
        }
    for (size_t j = 0; j < m->cols; j++)
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            column[k] = j;
            start[m->row_index[k] + 1]++;
        }
    size_t longest = 0;
    for (size_t i = 0; i < m->rows; i++) {
        longest = start[i + 1] > longest ? start[i + 1] : longest;
        start[i + 1] += start[i];
        at[i] = start[i];
    }
    for (size_t k = 0; k < m->nnz; k++)
        entries[at[m->row_index[k]]++] = k;
    LabInterval *term = lab_alloc(longest, sizeof(LabInterval));
    double deadline = started + opt.time_limit;
    int changed = 0;
    int infeasible = !opt.relax && propagate_parity(m, a, rl, ru, start, entries, column,
                                                    box, r, deadline);
    for (size_t pass = 0; pass < opt.max_passes && !infeasible; pass++) {
        changed = 0;
        int completed = 1;
        for (size_t i = 0; i < m->rows; i++) {
            Activity activity;
            if (lab_now() >= deadline || row_activity(start[i], start[i + 1], entries,
                    column, a, box, term, &activity, deadline)) {
                r->limited = 1;
                completed = 0;
                break;
            }
            double minimum = activity.infinite_min ? -INFINITY : activity.finite_min.lower;
            double maximum = activity.infinite_max ? INFINITY : activity.finite_max.upper;
            if (maximum < rl[i].lower || minimum > ru[i].upper) {
                infeasible = 1;
                break;
            }
            for (size_t p = start[i]; p < start[i + 1]; p++) {
                size_t k = entries[p], j = column[k];
                LabInterval coefficient = a[k], own = term[p - start[i]];
                if (!(coefficient.lower > 0 || coefficient.upper < 0))
                    continue;
                double other_min = -INFINITY, other_max = INFINITY;
                if (!activity.infinite_min)
                    other_min = lab_interval_subtract(point(activity.finite_min.lower),
                                                       point(own.lower)).lower;
                else if (activity.infinite_min == 1 && !isfinite(own.lower))
                    other_min = activity.finite_min.lower;
                if (!activity.infinite_max)
                    other_max = lab_interval_subtract(point(activity.finite_max.upper),
                                                       point(own.upper)).upper;
                else if (activity.infinite_max == 1 && !isfinite(own.upper))
                    other_max = activity.finite_max.upper;
                if (isfinite(ru[i].upper) && isfinite(other_min)) {
                    LabInterval numerator = lab_interval_subtract(point(ru[i].upper), point(other_min));
                    LabInterval bound = lab_interval_divide(numerator, coefficient);
                    int up = coefficient.lower > 0;
                    changed |= tighten(&box[j], up ? bound.upper : bound.lower,
                                        up, m->integer && m->integer[j] && !opt.relax, r);
                }
                if (isfinite(rl[i].lower) && isfinite(other_max)) {
                    LabInterval numerator = lab_interval_subtract(point(rl[i].lower), point(other_max));
                    LabInterval bound = lab_interval_divide(numerator, coefficient);
                    int up = coefficient.upper < 0;
                    changed |= tighten(&box[j], up ? bound.upper : bound.lower,
                                        up, m->integer && m->integer[j] && !opt.relax, r);
                }
                if (box[j].lower > box[j].upper) {
                    infeasible = 1;
                    break;
                }
            }
            if (infeasible)
                break;
        }
        if (completed)
            r->passes++;
        if (r->limited || !changed)
            break;
    }
    if (!infeasible && !opt.relax)
        probe_cliques(m, a, rl, ru, start, entries, column, box, r, deadline);
    if (changed && r->passes == opt.max_passes)
        r->limited = 1;
    for (size_t j = 0; j < m->cols; j++) {
        r->lower[j] = fmax(lower[j], box[j].lower);
        r->upper[j] = fmin(upper[j], box[j].upper);
        if (isfinite(box[j].lower) && box[j].lower == box[j].upper &&
            r->lower[j] == box[j].lower && r->upper[j] == box[j].upper)
            r->fixed++;
        /* A row may be removed only under an outer box containing the actual
         * published numerical box; unpublishable internal tightening is omitted. */
        if (r->lower[j] > lower[j])
            initial[j].lower = fmax(initial[j].lower, r->lower[j]);
        if (r->upper[j] < upper[j])
            initial[j].upper = fmin(initial[j].upper, r->upper[j]);
    }
    if (infeasible)
        diagnostic(r, "infeasible", "Source constraints and variable bounds imply infeasibility");
    else {
        diagnostic(r, "unresolved", r->limited ? "Safe propagation stopped at its resource limit"
                                                : "Safe propagation reached a fixed point");
        for (size_t i = 0; i < m->rows && lab_now() < deadline; i++) {
            Activity activity;
            if (row_activity(start[i], start[i + 1], entries, column, a, initial, term,
                              &activity, deadline))
                break;
            double minimum = activity.infinite_min ? -INFINITY : activity.finite_min.lower;
            double maximum = activity.infinite_max ? INFINITY : activity.finite_max.upper;
            if (minimum >= rl[i].upper && maximum <= ru[i].lower) {
                r->redundant[i] = 1;
                r->redundant_rows++;
            }
        }
    }
    free(term);
cleanup:
    free(a); free(rl); free(ru); free(box); free(initial);
    free(start); free(column); free(entries); free(at);
done:
    r->elapsed = lab_now() - started;
}
void lab_presolve_result_free(LabPresolveResult *r) {
    if (!r)
        return;
    free(r->lower);
    free(r->upper);
    free(r->redundant);
    memset(r, 0, sizeof *r);
}

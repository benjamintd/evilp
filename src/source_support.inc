/* Rigorous enclosure-only support-gap fallback. This is an optimizer stopping
 * safeguard, not an independent feasibility, bound or optimality certificate.
 * Only available==1 and a zero return may accept this fallback. */
#include <fenv.h>
/* The enclosure scan uses one interval per row and linear work. It does
 * not allocate a rational basis or a dense tableau. Limit row scratch to
 * 16 MB and bound both the sparse scan and its elapsed time independently. */
enum { SOURCE_SUPPORT_ROWS = 1000000, SOURCE_SUPPORT_COLS = 2000000, SOURCE_SUPPORT_NNZ = 4000000 };
static int source_support_finite(LabInterval q) {
    return isfinite(q.lower) && isfinite(q.upper) && q.lower <= q.upper;
}
static int source_support_absolute(LabInterval q, LabInterval *out) {
    if (!source_support_finite(q))
        return -1;
    *out = q.lower >= 0   ? q
           : q.upper <= 0 ? (LabInterval){-q.upper, -q.lower}
                          : (LabInterval){0, fmax(-q.lower, q.upper)};
    return 0;
}
static int source_support_add_absolute(LabInterval *gap, LabInterval term) {
    LabInterval magnitude;
    if (source_support_absolute(term, &magnitude))
        return -1;
    *gap = lab_interval_add(*gap, magnitude);
    return source_support_finite(*gap) ? 0 : -1;
}
/* Matching source identities retain the original decimal enclosure. A changed
 * bound is an explicit binary overlay, including a finite side added to an
 * originally infinite side. No uncertain source side is replaced by a point. */
static int source_support_column_bound(const LabModel *m, size_t j, int upper, LabInterval *out) {
    const LabSourceEnvelope *e = m->envelope;
    double value = upper ? m->col_upper[j] : m->col_lower[j];
    double source = upper ? e->col_upper_source[j] : e->col_lower_source[j];
    if (!isfinite(value) || isnan(source) || (upper ? source == -INFINITY : source == INFINITY))
        return -1;
    *out =
        value == source ? (upper ? e->col_upper[j] : e->col_lower[j]) : (LabInterval){value, value};
    return source_support_finite(*out) ? 0 : -1;
}
static int source_support_interval_gap(Work *w, const double *point, const double *row_dual,
                                       double deadline, int *available) {
    if (!available)
        return 0;
    *available = 0;
    if (!w || !w->model || !isfinite(deadline) || fegetround() != FE_TONEAREST)
        return 0;
    deadline = fmin(deadline, lab_now() + 1);
    const LabModel *m = w->model;
    const LabSourceEnvelope *e = m->envelope;
    if (m->rows > SOURCE_SUPPORT_ROWS || m->cols > SOURCE_SUPPORT_COLS ||
        m->nnz > SOURCE_SUPPORT_NNZ || !e || !e->complete || !e->objective_complete ||
        e->rows != m->rows || e->cols != m->cols || e->nnz != m->nnz || !m->col_start ||
        m->col_start[0] || m->col_start[m->cols] != m->nnz ||
        (m->rows &&
         (!m->row_lower || !m->row_upper || !row_dual || !e->row_lower || !e->row_upper)) ||
        (m->cols &&
         (!point || !m->objective || !m->col_lower || !m->col_upper || !e->objective ||
          !e->col_lower || !e->col_upper || !e->col_lower_source || !e->col_upper_source)) ||
        (m->nnz && (!m->row_index || !e->value)) || lab_now() >= deadline)
        return 0;
    for (size_t i = 0; i < m->rows; i++) {
        if (!isfinite(row_dual[i]) || isnan(m->row_lower[i]) || isnan(m->row_upper[i]) ||
            m->row_lower[i] == INFINITY || m->row_upper[i] == -INFINITY ||
            m->row_lower[i] > m->row_upper[i] || ((i & 255U) == 0 && lab_now() >= deadline))
            return 0;
        /* Finite row overlays have no source identity API. Only original
         * finite sides (or a mask to infinity) are admitted here. */
        if (isfinite(m->row_lower[i]) &&
            (!source_support_finite(e->row_lower[i]) || m->row_lower[i] < e->row_lower[i].lower ||
             m->row_lower[i] > e->row_lower[i].upper))
            return 0;
        if (isfinite(m->row_upper[i]) &&
            (!source_support_finite(e->row_upper[i]) || m->row_upper[i] < e->row_upper[i].lower ||
             m->row_upper[i] > e->row_upper[i].upper))
            return 0;
    }
    for (size_t j = 0; j < m->cols; j++) {
        if (!isfinite(point[j]) || !isfinite(m->objective[j]) ||
            !source_support_finite(e->objective[j]) || isnan(m->col_lower[j]) ||
            isnan(m->col_upper[j]) || m->col_lower[j] == INFINITY || m->col_upper[j] == -INFINITY ||
            m->col_lower[j] > m->col_upper[j] || isnan(e->col_lower_source[j]) ||
            isnan(e->col_upper_source[j]) || e->col_lower_source[j] == INFINITY ||
            e->col_upper_source[j] == -INFINITY || m->col_start[j] > m->col_start[j + 1] ||
            m->col_start[j + 1] > m->nnz || ((j & 255U) == 0 && lab_now() >= deadline))
            return 0;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++)
            if (m->row_index[k] >= m->rows || !source_support_finite(e->value[k]) ||
                ((k & 255U) == 0 && lab_now() >= deadline))
                return 0;
    }
    LabInterval *activity = lab_alloc(m->rows, sizeof *activity);
    LabInterval gap = {0, 0};
    int valid = 1;
    for (size_t j = 0; j < m->cols && valid; j++) {
        if (point[j] == 0)
            continue; /* finite source coefficient times zero is exactly zero */
        LabInterval x = {point[j], point[j]};
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            size_t i = m->row_index[k];
            if (row_dual[i] == 0)
                continue; /* this row contributes zero to the support gap */
            LabInterval term = lab_interval_multiply(e->value[k], x);
            activity[i] = lab_interval_add(activity[i], term);
            if (!source_support_finite(term) || !source_support_finite(activity[i]) ||
                ((k & 255U) == 0 && lab_now() >= deadline)) {
                valid = 0;
                break;
            }
        }
    }
    for (size_t i = 0; i < m->rows && valid; i++) {
        double y = row_dual[i];
        if (y == 0)
            continue;
        double current = y > 0 ? m->row_lower[i] : m->row_upper[i];
        LabInterval side = y > 0 ? e->row_lower[i] : e->row_upper[i];
        if (!isfinite(current) || !source_support_finite(side)) {
            valid = 0;
            break;
        }
        LabInterval slack = lab_interval_subtract(activity[i], side);
        LabInterval term = lab_interval_multiply((LabInterval){y, y}, slack);
        if (source_support_add_absolute(&gap, term) || lab_now() >= deadline)
            valid = 0;
    }
    for (size_t j = 0; j < m->cols && valid; j++) {
        LabInterval reduced = e->objective[j];
        if (m->maximize)
            reduced = (LabInterval){-reduced.upper, -reduced.lower};
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            double y = row_dual[m->row_index[k]];
            if (y == 0)
                continue;
            LabInterval term = lab_interval_multiply(e->value[k], (LabInterval){y, y});
            reduced = lab_interval_subtract(reduced, term);
            if (!source_support_finite(term) || !source_support_finite(reduced) ||
                ((k & 255U) == 0 && lab_now() >= deadline)) {
                valid = 0;
                break;
            }
        }
        if (!valid)
            break;
        if (reduced.lower == 0 && reduced.upper == 0)
            continue;
        /* At zero, either support choice gives a zero contribution. A sign
         * interval crossing zero requires both finite sides; bound the union
         * of the positive/lower and negative/upper support contributions. */
        LabInterval contribution = {0, 0};
        for (int upper = 0; upper < 2; upper++) {
            if (upper ? reduced.lower >= 0 : reduced.upper <= 0)
                continue;
            LabInterval bound;
            if (source_support_column_bound(m, j, upper, &bound)) {
                valid = 0;
                break;
            }
            LabInterval sign_part = upper ? (LabInterval){reduced.lower, fmin(0, reduced.upper)}
                                          : (LabInterval){fmax(0, reduced.lower), reduced.upper};
            LabInterval slack = lab_interval_subtract((LabInterval){point[j], point[j]}, bound);
            LabInterval term = lab_interval_multiply(sign_part, slack), magnitude;
            if (source_support_absolute(term, &magnitude)) {
                valid = 0;
                break;
            }
            contribution.upper = fmax(contribution.upper, magnitude.upper);
        }
        if (valid && (source_support_add_absolute(&gap, contribution) || lab_now() >= deadline))
            valid = 0;
    }
    free(activity);
    if (!valid || !source_support_finite(gap) || lab_now() >= deadline)
        return 0;
    *available = 1;
    return gap.upper > 1e-7;
}

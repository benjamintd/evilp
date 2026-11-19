/* Bounded source-objective decisions. This file never changes .exact or the
 * floating factor: exact signs select ordinary guarded primal pivots, and an
 * exact source dual can close a numerically rejected support-gap check. */
enum { SOURCE_PRICE_ROWS = 4096, SOURCE_PRICE_COLS = 131072, SOURCE_PRICE_NNZ = 2000000 };
static int source_price_q(const mpq_t q) {
    return mpz_sgn(mpq_denref(q)) > 0 &&
           mpz_sizeinbase(mpq_numref(q), 2) <= LAB_RATIONAL_MAX_BITS &&
           mpz_sizeinbase(mpq_denref(q), 2) <= LAB_RATIONAL_MAX_BITS;
}
static int source_price_guard(const mpq_t q, double deadline) {
    return source_price_q(q) && lab_now() < deadline;
}
/* Reuse the prior binary-gap path only where it already applied, and where
 * complete source enclosures prove every datum used by that path is binary.
 * Small models that previously used exact support checks retain those checks. */
static int source_price_binary_equivalent(const Work *w) {
    const LabModel *m = w->model;
    if (exact_source(m))
        return 0;
    if (!m->parsed_decimal && !m->native_exact && !m->objective_exact)
        return 1;
    const LabSourceEnvelope *e = m->envelope;
    if (!e || !e->complete || !e->objective_complete || e->rows != m->rows || e->cols != m->cols ||
        e->nnz != m->nnz || !e->value || !e->objective || !e->row_lower || !e->row_upper ||
        !e->col_lower || !e->col_upper || !e->col_lower_source || !e->col_upper_source)
        return 0;
    for (size_t k = 0; k < m->nnz; k++)
        if (e->value[k].lower != m->value[k] || e->value[k].upper != m->value[k])
            return 0;
    for (size_t i = 0; i < m->rows; i++) {
        if (isfinite(m->row_lower[i]) &&
            (e->row_lower[i].lower != m->row_lower[i] || e->row_lower[i].upper != m->row_lower[i]))
            return 0;
        if (isfinite(m->row_upper[i]) &&
            (e->row_upper[i].lower != m->row_upper[i] || e->row_upper[i].upper != m->row_upper[i]))
            return 0;
    }
    for (size_t j = 0; j < m->cols; j++) {
        if (e->objective[j].lower != m->objective[j] || e->objective[j].upper != m->objective[j])
            return 0;
        if (isfinite(m->col_lower[j]) && m->col_lower[j] == e->col_lower_source[j] &&
            (e->col_lower[j].lower != m->col_lower[j] || e->col_lower[j].upper != m->col_lower[j]))
            return 0;
        if (isfinite(m->col_upper[j]) && m->col_upper[j] == e->col_upper_source[j] &&
            (e->col_upper[j].lower != m->col_upper[j] || e->col_upper[j].upper != m->col_upper[j]))
            return 0;
        if (w->transform[j] != 'X' && w->transform[j] != 'F') {
            int upper = w->transform[j] == 'U';
            double source = upper ? e->col_upper_source[j] : e->col_lower_source[j];
            LabInterval bound = upper ? e->col_upper[j] : e->col_lower[j];
            if (w->shift[j] == source && (bound.lower != source || bound.upper != source))
                return 0;
        }
    }
    return 1;
}
static int source_objective_risk(const Work *w, double deadline) {
    const LabModel *m = w->model;
    if (!lab_source_objective_needs_exact(m))
        return 0;
    if (m->cols > SOURCE_PRICE_COLS)
        return 1;
    mpq_t cost, binary, error, lo, hi, width, total, tolerance;
    mpq_inits(cost, binary, error, lo, hi, width, total, tolerance, NULL);
    mpq_set_ui(tolerance, 1, 10000000);
    int risk = 0;
    for (size_t j = 0; j < m->cols; j++) {
        if (lab_now() >= deadline || lab_source_exact_objective(m, j, cost) ||
            !source_price_guard(cost, deadline)) {
            risk = 1;
            break;
        }
        mpq_set_d(binary, m->objective[j]);
        mpq_sub(error, cost, binary);
        mpq_abs(error, error);
        if (!mpq_sgn(error))
            continue;
        if (!isfinite(m->col_lower[j]) || !isfinite(m->col_upper[j]) ||
            lab_source_exact_column_bound(m, j, 0, lo) ||
            lab_source_exact_column_bound(m, j, 1, hi)) {
            risk = 1;
            break;
        }
        mpq_sub(width, hi, lo);
        if (mpq_sgn(width) < 0 || !source_price_guard(width, deadline) || !source_price_q(error)) {
            risk = 1;
            break;
        }
        mpq_mul(error, error, width);
        if (!source_price_guard(error, deadline)) {
            risk = 1;
            break;
        }
        mpq_add(total, total, error);
        if (!source_price_q(total) || mpq_cmp(total, tolerance) > 0) {
            risk = 1;
            break;
        }
    }
    mpq_clears(cost, binary, error, lo, hi, width, total, tolerance, NULL);
    return risk;
}
static int source_price_cost(const Work *w, size_t j, mpq_t q) {
    if (j >= w->structural) {
        mpq_set_ui(q, 0, 1);
        return 0;
    }
    if (lab_source_exact_objective(w->model, w->original[j], q) || !source_price_q(q))
        return -1;
    mpq_t scale;
    mpq_init(scale);
    mpq_set_d(scale, w->sign[j]);
    mpq_mul(q, q, scale);
    if (w->model->maximize)
        mpq_neg(q, q);
    mpq_clear(scale);
    return source_price_q(q) ? 0 : -1;
}
static int source_price_entry(const Work *w, size_t j, size_t k, mpq_t q) {
    const Column *c = &w->a[j];
    size_t i = c->row[k];
    if (j >= w->structural || w->source_row[i] == SIZE_MAX) {
        mpq_set_d(q, c->value[k]);
        return 0;
    }
    const LabModel *m = w->model;
    size_t original = w->original[j], source = w->source_row[i];
    size_t lo = m->col_start[original], hi = m->col_start[original + 1];
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (m->row_index[mid] < source)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo == m->col_start[original + 1] || m->row_index[lo] != source ||
        lab_source_exact_coefficient(m, lo, q) || !source_price_q(q))
        return -1;
    mpq_t scale;
    mpq_init(scale);
    mpq_set_d(scale, w->sign[j]);
    mpq_mul(q, q, scale);
    if (!source_price_q(q)) {
        mpq_clear(scale);
        return -1;
    }
    mpq_set_d(scale, w->row_scale[i]);
    mpq_mul(q, q, scale);
    mpq_clear(scale);
    return source_price_q(q) ? 0 : -1;
}
/* An exact coefficient that vanished from the floating canonical column
 * cannot silently disappear from an exact reduced-cost calculation. */
static int source_price_coverage(const Work *w, double deadline) {
    const LabModel *m = w->model;
    size_t *head = lab_alloc(m->rows, sizeof(size_t)), *next = lab_alloc(w->m, sizeof(size_t));
    size_t *seen = lab_alloc(w->m, sizeof *seen);
    mpq_t q;
    mpq_init(q);
    int valid = 1;
    for (size_t i = 0; i < m->rows; i++)
        head[i] = SIZE_MAX;
    for (size_t i = 0; i < w->m; i++)
        if (w->source_row[i] != SIZE_MAX) {
            next[i] = head[w->source_row[i]];
            head[w->source_row[i]] = i;
        }
    for (size_t j = 0; j < w->structural && valid; j++) {
        if (lab_now() >= deadline) {
            valid = 0;
            break;
        }
        const Column *c = &w->a[j];
        for (size_t k = 0; k < c->count; k++)
            seen[c->row[k]] = j + 1;
        size_t original = w->original[j];
        for (size_t k = m->col_start[original]; k < m->col_start[original + 1]; k++) {
            /* Only omitted coordinates need a source zero check here. The
             * dual solve and reduced-cost scan fetch represented coefficients. */
            for (size_t i = head[m->row_index[k]]; i != SIZE_MAX; i = next[i]) {
                if (seen[i] != j + 1 &&
                    (lab_source_exact_coefficient(m, k, q) ||
                     !source_price_guard(q, deadline) || mpq_sgn(q))) {
                    valid = 0;
                    break;
                }
            }
            if (!valid)
                break;
        }
    }
    mpq_clear(q);
    free(head);
    free(next);
    free(seen);
    return valid ? 0 : -1;
}
static void source_price_free(mpq_t *q, size_t n) {
    if (q)
        for (size_t i = 0; i < n; i++)
            mpq_clear(q[i]);
    free(q);
}
static int source_price_dual(Work *w, double deadline, mpq_t **out) {
    const LabModel *m = w->model;
    size_t n = w->m, nnz = 0;
    if (w->options->no_exact || n > LAB_RATIONAL_SPARSE_MAX_DIMENSION ||
        m->rows > SOURCE_PRICE_ROWS || m->cols > SOURCE_PRICE_COLS || m->nnz > SOURCE_PRICE_NNZ ||
        w->n > SOURCE_PRICE_NNZ)
        return -1;
    for (size_t t = 0; t < n; t++)
        nnz += w->a[w->basis[t]].count;
    if (nnz > LAB_RATIONAL_SPARSE_MAX_NONZEROS || lab_now() >= deadline ||
        source_price_coverage(w, deadline))
        return -1;
    size_t *start = lab_alloc(n + 1, sizeof(size_t)), *at = lab_alloc(n, sizeof(size_t)),
           *row = lab_alloc(nnz, sizeof(size_t));
    mpq_t *value = lab_alloc(nnz, sizeof(mpq_t)), *rhs = lab_alloc(n, sizeof(mpq_t)),
          *pi = lab_alloc(n, sizeof(mpq_t));
    for (size_t k = 0; k < nnz; k++)
        mpq_init(value[k]);
    for (size_t i = 0; i < n; i++) {
        mpq_init(rhs[i]);
        mpq_init(pi[i]);
    }
    for (size_t t = 0; t < n; t++)
        for (size_t k = 0; k < w->a[w->basis[t]].count; k++)
            start[w->a[w->basis[t]].row[k] + 1]++;
    for (size_t i = 0; i < n; i++) {
        start[i + 1] += start[i];
        at[i] = start[i];
    }
    int valid = 1;
    for (size_t t = 0; t < n && valid; t++) {
        size_t j = w->basis[t];
        if (source_price_cost(w, j, rhs[t])) {
            valid = 0;
            break;
        }
        for (size_t k = 0; k < w->a[j].count; k++) {
            size_t p = at[w->a[j].row[k]]++;
            row[p] = t;
            if (source_price_entry(w, j, k, value[p]) || !source_price_guard(value[p], deadline)) {
                valid = 0;
                break;
            }
        }
    }
    if (valid && lab_rational_sparse_solve(n, start, row, value, rhs, pi, deadline))
        valid = 0;
    for (size_t i = 0; i < n && valid; i++)
        if (!source_price_guard(pi[i], deadline))
            valid = 0;
    /* Check B^T*pi=cB afresh, independently of sparse-solve success. */
    mpq_t sum, term, entry;
    mpq_inits(sum, term, entry, NULL);
    for (size_t t = 0; t < n && valid; t++) {
        mpq_set_ui(sum, 0, 1);
        size_t j = w->basis[t];
        for (size_t k = 0; k < w->a[j].count; k++) {
            if (source_price_entry(w, j, k, entry)) {
                valid = 0;
                break;
            }
            mpq_mul(term, entry, pi[w->a[j].row[k]]);
            if (!source_price_guard(term, deadline)) {
                valid = 0;
                break;
            }
            mpq_add(sum, sum, term);
            if (!source_price_guard(sum, deadline)) {
                valid = 0;
                break;
            }
        }
        if (mpq_cmp(sum, rhs[t]) || !source_price_q(sum) || lab_now() >= deadline)
            valid = 0;
    }
    mpq_clears(sum, term, entry, NULL);
    free(start);
    free(at);
    free(row);
    source_price_free(value, nnz);
    source_price_free(rhs, n);
    if (!valid) {
        source_price_free(pi, n);
        return -1;
    }
    *out = pi;
    return 0;
}
/* 0 source dual feasible, 1 real source-improving entering column, -1 abstain. */
static int source_price_scan(Work *w, double deadline, mpq_t *pi, size_t *entering) {
    mpq_t reduced, entry, term, best;
    mpq_inits(reduced, entry, term, best, NULL);
    *entering = SIZE_MAX;
    int valid = 1;
    for (size_t j = 0; j < w->n && valid; j++) {
        if (w->artificial[j] || is_basic(w, j))
            continue;
        if (source_price_cost(w, j, reduced)) {
            valid = 0;
            break;
        }
        const Column *c = &w->a[j];
        for (size_t k = 0; k < c->count; k++) {
            if (source_price_entry(w, j, k, entry)) {
                valid = 0;
                break;
            }
            mpq_mul(term, entry, pi[c->row[k]]);
            if (!source_price_guard(term, deadline)) {
                valid = 0;
                break;
            }
            mpq_sub(reduced, reduced, term);
            if (!source_price_guard(reduced, deadline)) {
                valid = 0;
                break;
            }
        }
        if (!source_price_q(reduced) || lab_now() >= deadline) {
            valid = 0;
            break;
        }
        if (mpq_sgn(reduced) < 0 && (*entering == SIZE_MAX || mpq_cmp(reduced, best) < 0)) {
            *entering = j;
            mpq_set(best, reduced);
        }
    }
    mpq_clears(reduced, entry, term, best, NULL);
    return !valid ? -1 : *entering == SIZE_MAX ? 0 : 1;
}
static int source_objective_price(Work *w, double deadline, size_t *entering) {
    w->result->source_pricing_attempts++;
    mpq_t *pi = NULL;
    if (source_price_dual(w, deadline, &pi))
        return -1;
    int result = source_price_scan(w, deadline, pi, entering);
    source_price_free(pi, w->m);
    if (!result)
        w->result->source_pricing_completions++;
    return result;
}
/* Exact source support/complementarity gap. Input row multipliers are exact;
 * a zero gap is distinct from merely rounding a source dual into binary64. */
static int source_price_gap(Work *w, const double *point, mpq_t *dual, double deadline,
                            int *available) {
    const LabModel *m = w->model;
    *available = 0;
    if (m->rows > SOURCE_PRICE_ROWS || m->cols > SOURCE_PRICE_COLS || m->nnz > SOURCE_PRICE_NNZ)
        return 0;
    mpq_t *activity = lab_alloc(m->rows, sizeof(mpq_t));
    for (size_t i = 0; i < m->rows; i++)
        mpq_init(activity[i]);
    mpq_t x, coefficient, term, reduced, side, gap, variation, anchor, tolerance;
    mpq_inits(x, coefficient, term, reduced, side, gap, variation, anchor, tolerance, NULL);
    int valid = 1;
    for (size_t j = 0; j < m->cols && valid; j++) {
        mpq_set_d(x, point[j]);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            if (lab_source_exact_coefficient(m, k, coefficient) ||
                !source_price_guard(coefficient, deadline)) {
                valid = 0;
                break;
            }
            mpq_mul(term, coefficient, x);
            if (!source_price_guard(term, deadline)) {
                valid = 0;
                break;
            }
            mpq_add(activity[m->row_index[k]], activity[m->row_index[k]], term);
            if (!source_price_guard(activity[m->row_index[k]], deadline)) {
                valid = 0;
                break;
            }
        }
        if (!valid)
            break;
        if (w->transform[j] != 'X') {
            if (lab_source_exact_objective(m, j, coefficient) ||
                !source_price_guard(coefficient, deadline)) {
                valid = 0;
                break;
            }
            if (w->transform[j] == 'F')
                mpq_set_ui(anchor, 0, 1);
            else if (lab_source_exact_column_value(m, j, w->transform[j] == 'U', w->shift[j],
                                                   anchor)) {
                valid = 0;
                break;
            }
            mpq_sub(term, x, anchor);
            if (!source_price_guard(term, deadline)) {
                valid = 0;
                break;
            }
            mpq_mul(term, term, coefficient);
            mpq_abs(term, term);
            if (!source_price_guard(term, deadline)) {
                valid = 0;
                break;
            }
            mpq_add(variation, variation, term);
            if (!source_price_guard(variation, deadline)) {
                valid = 0;
                break;
            }
        }
        if (lab_now() >= deadline) {
            valid = 0;
            break;
        }
    }
    for (size_t i = 0; i < m->rows && valid; i++) {
        int sign = mpq_sgn(dual[i]);
        if (!sign)
            continue;
        if (!source_price_guard(dual[i], deadline)) {
            valid = 0;
            break;
        }
        if (lab_source_exact_row_bound(m, i, sign < 0, side) ||
            !source_price_guard(side, deadline)) {
            valid = 0;
            break;
        }
        mpq_sub(term, activity[i], side);
        if (!source_price_guard(term, deadline)) {
            valid = 0;
            break;
        }
        mpq_mul(term, term, dual[i]);
        mpq_abs(term, term);
        if (!source_price_guard(term, deadline)) {
            valid = 0;
            break;
        }
        mpq_add(gap, gap, term);
        if (!source_price_guard(gap, deadline)) {
            valid = 0;
            break;
        }
    }
    for (size_t j = 0; j < m->cols && valid; j++) {
        if (lab_source_exact_objective(m, j, reduced) || !source_price_guard(reduced, deadline)) {
            valid = 0;
            break;
        }
        if (m->maximize)
            mpq_neg(reduced, reduced);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            if (lab_source_exact_coefficient(m, k, coefficient) ||
                !source_price_guard(coefficient, deadline)) {
                valid = 0;
                break;
            }
            if (!source_price_guard(dual[m->row_index[k]], deadline)) {
                valid = 0;
                break;
            }
            mpq_mul(term, coefficient, dual[m->row_index[k]]);
            if (!source_price_guard(term, deadline)) {
                valid = 0;
                break;
            }
            mpq_sub(reduced, reduced, term);
            if (!source_price_guard(reduced, deadline)) {
                valid = 0;
                break;
            }
        }
        if (!valid)
            break;
        int sign = mpq_sgn(reduced);
        if (!sign)
            continue;
        if (lab_source_exact_column_bound(m, j, sign < 0, side) ||
            !source_price_guard(side, deadline)) {
            valid = 0;
            break;
        }
        mpq_set_d(x, point[j]);
        mpq_sub(term, x, side);
        if (!source_price_guard(term, deadline)) {
            valid = 0;
            break;
        }
        mpq_mul(term, term, reduced);
        mpq_abs(term, term);
        if (!source_price_guard(term, deadline)) {
            valid = 0;
            break;
        }
        mpq_add(gap, gap, term);
        if (!source_price_q(gap) || lab_now() >= deadline) {
            valid = 0;
            break;
        }
    }
    mpq_set_ui(tolerance, 1, 1);
    mpq_add(tolerance, tolerance, variation);
    mpq_set_ui(term, 1, 10000000);
    mpq_mul(tolerance, tolerance, term);
    if (!source_price_guard(tolerance, deadline))
        valid = 0;
    *available = valid ? 1 : -1;
    int bad = valid && mpq_cmp(gap, tolerance) > 0;
    source_price_free(activity, m->rows);
    mpq_clears(x, coefficient, term, reduced, side, gap, variation, anchor, tolerance, NULL);
    return bad;
}
static int source_numeric_gap(Work *w, const double *point, const double *row_dual,
                              int *available) {
    const LabModel *m = w->model;
    *available = 0;
    if (m->rows > SOURCE_PRICE_ROWS || m->cols > SOURCE_PRICE_COLS || m->nnz > SOURCE_PRICE_NNZ)
        return 0;
    w->result->source_support_scans++;
    mpq_t *dual = lab_alloc(m->rows, sizeof(mpq_t));
    for (size_t i = 0; i < m->rows; i++) {
        mpq_init(dual[i]);
        mpq_set_d(dual[i], row_dual[i]);
    }
    double deadline = fmin(w->started + w->options->time_limit, lab_now() + 1);
    int bad = source_price_gap(w, point, dual, deadline, available);
    source_price_free(dual, m->rows);
    return bad;
}
static int source_gap_closure(Work *w, const double *point) {
    double deadline = fmin(w->started + w->options->time_limit, lab_now() + 1);
    mpq_t *pi = NULL;
    if (source_price_dual(w, deadline, &pi))
        return 0;
    size_t entering;
    int scan = source_price_scan(w, deadline, pi, &entering);
    if (scan) {
        source_price_free(pi, w->m);
        return 0;
    }
    const LabModel *m = w->model;
    mpq_t *dual = lab_alloc(m->rows, sizeof(mpq_t));
    for (size_t i = 0; i < m->rows; i++)
        mpq_init(dual[i]);
    mpq_t scale, term;
    mpq_inits(scale, term, NULL);
    int valid = 1;
    for (size_t i = 0; i < w->m && valid; i++)
        if (w->source_row[i] != SIZE_MAX) {
            mpq_set_d(scale, w->row_scale[i]);
            mpq_mul(term, pi[i], scale);
            if (!source_price_guard(term, deadline)) {
                valid = 0;
                break;
            }
            mpq_add(dual[w->source_row[i]], dual[w->source_row[i]], term);
            if (!source_price_guard(dual[w->source_row[i]], deadline)) {
                valid = 0;
                break;
            }
        }
    int available = 0, bad = 1;
    if (valid)
        bad = source_price_gap(w, point, dual, deadline, &available);
    if (available == 1 && !bad) {
        w->result->dual_rational = lab_alloc(m->rows, sizeof(char *));
        w->result->dual_rows = m->rows;
        for (size_t i = 0; i < m->rows; i++)
            w->result->dual_rational[i] = mpq_get_str(NULL, 10, dual[i]);
    }
    source_price_free(pi, w->m);
    source_price_free(dual, m->rows);
    mpq_clears(scale, term, NULL);
    if (available == 1 && !bad)
        w->result->source_dual_closures++;
    return available == 1 && !bad;
}

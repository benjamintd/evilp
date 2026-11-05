#include "farkas_candidate.h"
#include "exact.h"
#include "common.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum {
    MAX_ROWS = 4000000, MAX_COLUMNS = 2000000, MAX_ENTRIES = 32000000,
    MAX_TERMS = 8000000, MAX_BITS = 8192
};
static int bounded(const mpq_t q) {
    return mpz_sgn(mpq_denref(q)) > 0 &&
           mpz_sizeinbase(mpq_numref(q), 2) <= MAX_BITS &&
           mpz_sizeinbase(mpq_denref(q), 2) <= MAX_BITS;
}
static void message(LabFarkasCandidate *r, const char *text) {
    snprintf(r->message, sizeof r->message, "%s", text);
}
static int interval_valid(LabInterval x) {
    return !isnan(x.lower) && !isnan(x.upper) && x.lower <= x.upper;
}
static int provenance(const LabModel *m) {
    const LabExactModel *e = m->exact;
    const LabSourceEnvelope *s = m->envelope;
    if (e) {
        if (e->rows != m->rows || e->cols != m->cols || e->nnz != m->nnz ||
            (m->nnz && !e->value) || (m->rows && (!e->row_lower || !e->row_upper)) ||
            (m->cols && (!e->col_lower || !e->col_upper || !e->col_lower_source ||
                         !e->col_upper_source)))
            return -1;
        return 1;
    }
    if (s) {
        if (!s->complete)
            return 0;
        if (s->rows != m->rows || s->cols != m->cols || s->nnz != m->nnz ||
            (m->nnz && !s->value) || (m->rows && (!s->row_lower || !s->row_upper)) ||
            (m->cols && (!s->col_lower || !s->col_upper || !s->col_lower_source ||
                         !s->col_upper_source)))
            return -1;
        return 1;
    }
    return m->parsed_decimal ? 0 : 1;
}
static int coefficient(const LabModel *m, size_t k, mpq_t lo, mpq_t hi) {
    if (m->exact) {
        mpq_set(lo, m->exact->value[k]);
        mpq_set(hi, lo);
    } else if (m->envelope) {
        LabInterval x = m->envelope->value[k];
        if (!interval_valid(x) || !isfinite(x.lower) || !isfinite(x.upper))
            return 0;
        mpq_set_d(lo, x.lower);
        mpq_set_d(hi, x.upper);
    } else {
        mpq_set_d(lo, m->value[k]);
        mpq_set(hi, lo);
    }
    return bounded(lo) && bounded(hi);
}
static int row_side(const LabModel *m, size_t i, int upper, mpq_t side) {
    double current = upper ? m->row_upper[i] : m->row_lower[i];
    if (!isfinite(current))
        return 0;
    if (m->exact)
        mpq_set(side, upper ? m->exact->row_upper[i] : m->exact->row_lower[i]);
    else if (m->envelope) {
        LabInterval x = upper ? m->envelope->row_upper[i] : m->envelope->row_lower[i];
        double endpoint = upper ? x.upper : x.lower;
        if (!interval_valid(x) || !isfinite(endpoint))
            return 0;
        mpq_set_d(side, endpoint);
    } else
        mpq_set_d(side, current);
    return bounded(side);
}
typedef struct {
    mpq_t value;
    int infinity; /* -1,-infinity;0,finite;+1,+infinity */
} Endpoint;
static int column_side(const LabModel *m, size_t j, int upper, Endpoint *out) {
    double current = upper ? m->col_upper[j] : m->col_lower[j];
    double identity = current;
    out->infinity = 0;
    if (m->exact) {
        const LabExactModel *e = m->exact;
        identity = upper ? e->col_upper_source[j] : e->col_lower_source[j];
        if (!isfinite(identity))
            out->infinity = identity > 0 ? 1 : -1;
        else
            mpq_set(out->value, upper ? e->col_upper[j] : e->col_lower[j]);
    } else if (m->envelope) {
        const LabSourceEnvelope *e = m->envelope;
        LabInterval x = upper ? e->col_upper[j] : e->col_lower[j];
        identity = upper ? e->col_upper_source[j] : e->col_lower_source[j];
        if (!interval_valid(x))
            return 0;
        double endpoint = upper ? x.upper : x.lower;
        if (isfinite(endpoint))
            mpq_set_d(out->value, endpoint);
        else
            out->infinity = endpoint > 0 ? 1 : -1;
    } else if (isfinite(current))
        mpq_set_d(out->value, current);
    else
        out->infinity = current > 0 ? 1 : -1;
    /* Match LP view/evaluator semantics: a changed side replaces that source
     * side, including a widening to infinity. Only an unchanged binary identity
     * retains its original decimal enclosure. */
    if (current != identity) {
        if (isfinite(current)) {
            mpq_set_d(out->value, current);
            out->infinity = 0;
        } else
            out->infinity = current > 0 ? 1 : -1;
    }
    return out->infinity || bounded(out->value);
}
static int support_product(const mpq_t q, const Endpoint *x, mpq_t product,
                             mpq_t maximum, int *has_maximum) {
    int sign = mpq_sgn(q);
    if (x->infinity && sign) {
        if (sign * x->infinity > 0)
            return 0;
        return 1; /* a negative-infinite endpoint cannot increase the maximum */
    }
    if (!sign)
        mpq_set_ui(product, 0, 1);
    else
        mpq_mul(product, q, x->value);
    if (!bounded(product))
        return 0;
    if (!*has_maximum || mpq_cmp(product, maximum) > 0) {
        mpq_set(maximum, product);
        *has_maximum = 1;
    }
    return 1;
}
static int construct_candidate(const LabModel *m, size_t terms, const size_t *source_row,
                                 const double *scale, const double *pi, double deadline,
                                 int project_rows, int simplify_terms, LabFarkasCandidate *r) {
    if (!r)
        return -1;
    memset(r, 0, sizeof *r);
    double started = lab_now();
    message(r, "Invalid Farkas candidate input");
    int outcome = -1;
    if (!m || !isfinite(deadline) || (terms && (!source_row || !scale || !pi)) ||
        (m->cols && (!m->col_start || !m->col_lower || !m->col_upper)) ||
        (m->rows && (!m->row_lower || !m->row_upper)) ||
        (m->nnz && (!m->row_index || !m->value)) || (!m->cols && m->nnz))
        goto done;
    if (m->rows > MAX_ROWS || m->cols > MAX_COLUMNS || m->nnz > MAX_ENTRIES ||
        terms > MAX_TERMS) {
        message(r, "Farkas candidate storage limit");
        outcome = 0;
        goto done;
    }
    if (lab_now() >= deadline) {
        message(r, "Farkas candidate deadline exhausted");
        outcome = 0;
        goto done;
    }
    if (m->cols && (m->col_start[0] || m->col_start[m->cols] != m->nnz))
        goto done;
    for (size_t j = 0; j < m->cols; j++) {
        if (m->col_start[j] > m->col_start[j + 1] || isnan(m->col_lower[j]) ||
            isnan(m->col_upper[j]) || m->col_lower[j] == INFINITY ||
            m->col_upper[j] == -INFINITY)
            goto done;
        if ((j & 1023U) == 0 && lab_now() >= deadline)
            goto expired;
    }
    for (size_t i = 0; i < m->rows; i++) {
        if (isnan(m->row_lower[i]) || isnan(m->row_upper[i]))
            goto done;
        if ((i & 1023U) == 0 && lab_now() >= deadline)
            goto expired;
    }
    for (size_t i = 0; i < terms; i++) {
        if (!isfinite(scale[i]) || !isfinite(pi[i]) ||
            (source_row[i] != SIZE_MAX && source_row[i] >= m->rows))
            goto done;
        if ((i & 1023U) == 0 && lab_now() >= deadline)
            goto expired;
    }
    for (size_t k = 0; k < m->nnz; k++) {
        if (m->row_index[k] >= m->rows || !isfinite(m->value[k]))
            goto done;
        if ((k & 1023U) == 0 && lab_now() >= deadline)
            goto expired;
    }
    int source = provenance(m);
    if (source <= 0) {
        outcome = source;
        message(r, source ? "Malformed original-source provenance"
                           : "Original-source provenance is unavailable");
        goto done;
    }
    if (m->exact || m->envelope)
        for (size_t j = 0; j < m->cols; j++) {
            double lo = m->exact ? m->exact->col_lower_source[j]
                                 : m->envelope->col_lower_source[j];
            double hi = m->exact ? m->exact->col_upper_source[j]
                                 : m->envelope->col_upper_source[j];
            if (isnan(lo) || isnan(hi) || lo == INFINITY || hi == -INFINITY) {
                message(r, "Invalid original-source bound identity");
                goto done;
            }
            if ((j & 1023U) == 0 && lab_now() >= deadline)
                goto expired;
        }
    outcome = 0;
    message(r, "No strict source-enclosure contradiction established");
    mpq_t *y = calloc(m->rows ? m->rows : 1, sizeof(mpq_t));
    if (!y) {
        message(r, "Cannot allocate candidate multipliers");
        goto done;
    }
    size_t initialized = 0;
    for (; initialized < m->rows; initialized++) {
        if ((initialized & 1023U) == 0 && lab_now() >= deadline)
            break;
        mpq_init(y[initialized]);
    }
    mpq_t qlo, qhi, alo, ahi, product_lo, product_hi, factor, value, maximum, margin, cutoff;
    mpq_inits(qlo, qhi, alo, ahi, product_lo, product_hi, factor, value, maximum, margin, cutoff, NULL);
    Endpoint box_lo = {0}, box_hi = {0};
    mpq_init(box_lo.value);
    mpq_init(box_hi.value);
    if (initialized != m->rows)
        goto cleanup;
    if (simplify_terms) {
        /* Candidate construction only: exact comparisons avoid underflow in
         * max(abs(pi))*2^-40 and preserve the inclusive cutoff at every scale. */
        for (size_t i = 0; i < terms; i++) {
            mpq_set_d(value, pi[i]);
            mpq_abs(value, value);
            if (mpq_cmp(value, cutoff) > 0)
                mpq_set(cutoff, value);
            if ((i & 1023U) == 0 && lab_now() >= deadline)
                goto cleanup;
        }
        mpq_div_2exp(cutoff, cutoff, 40);
    }
    for (size_t i = 0; i < terms; i++) {
        if ((i & 1023U) == 0 && lab_now() >= deadline)
            goto cleanup;
        if (simplify_terms && pi[i] != 0) {
            mpq_set_d(value, pi[i]);
            mpq_abs(factor, value);
            if (mpq_cmp(factor, cutoff) <= 0) {
                r->removed_terms++;
                continue;
            }
        }
        if (source_row[i] != SIZE_MAX) {
            mpq_set_d(value, pi[i]);
            mpq_set_d(factor, scale[i]);
            mpq_mul(value, value, factor);
            mpq_add(y[source_row[i]], y[source_row[i]], value);
            if (!bounded(y[source_row[i]]))
                goto cleanup;
        }
    }
    if (project_rows)
        for (size_t i = 0; i < m->rows; i++) {
            int sign = mpq_sgn(y[i]);
            if (sign && !isfinite(sign > 0 ? m->row_lower[i] : m->row_upper[i])) {
                mpq_set_ui(y[i], 0, 1);
                r->projected_rows++;
            }
            if ((i & 1023U) == 0 && lab_now() >= deadline)
                goto cleanup;
        }
    for (size_t i = 0; i < m->rows; i++) {
        int sign = mpq_sgn(y[i]);
        if (sign) {
            if (!row_side(m, i, sign < 0, value)) {
                message(r, "Selected source row side is unavailable or outside the arithmetic limit");
                goto cleanup;
            }
            mpq_mul(value, value, y[i]);
            mpq_add(margin, margin, value);
            if (!bounded(value) || !bounded(margin))
                goto cleanup;
        }
        if ((i & 1023U) == 0 && lab_now() >= deadline)
            goto cleanup;
    }
    for (size_t j = 0; j < m->cols; j++) {
        mpq_set_ui(qlo, 0, 1);
        mpq_set_ui(qhi, 0, 1);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            size_t i = m->row_index[k];
            if (mpq_sgn(y[i])) {
                if (!coefficient(m, k, alo, ahi)) {
                    message(r, "Source coefficient enclosure is unavailable or outside the arithmetic limit");
                    goto cleanup;
                }
                mpq_mul(product_lo, alo, y[i]);
                mpq_mul(product_hi, ahi, y[i]);
                if (mpq_sgn(y[i]) < 0)
                    mpq_swap(product_lo, product_hi);
                mpq_add(qlo, qlo, product_lo);
                mpq_add(qhi, qhi, product_hi);
                if (!bounded(product_lo) || !bounded(product_hi) || !bounded(qlo) || !bounded(qhi))
                    goto cleanup;
            }
            if ((k & 1023U) == 0 && lab_now() >= deadline)
                goto cleanup;
        }
        if (!column_side(m, j, 0, &box_lo) ||
            !column_side(m, j, 1, &box_hi) || box_lo.infinity > 0 ||
            box_hi.infinity < 0 || (!box_lo.infinity && !box_hi.infinity &&
                                    mpq_cmp(box_lo.value, box_hi.value) > 0))
            goto cleanup;
        int has_maximum = 0;
        if (!support_product(qlo, &box_lo, value, maximum, &has_maximum) ||
            !support_product(qlo, &box_hi, value, maximum, &has_maximum) ||
            !support_product(qhi, &box_lo, value, maximum, &has_maximum) ||
            !support_product(qhi, &box_hi, value, maximum, &has_maximum) || !has_maximum) {
            message(r, "Source box support is unbounded or outside the arithmetic limit");
            goto cleanup;
        }
        mpq_sub(margin, margin, maximum);
        if (!bounded(margin) || lab_now() >= deadline)
            goto cleanup;
    }
    if (mpq_sgn(margin) <= 0 || lab_now() >= deadline)
        goto cleanup;
    char **text = calloc(m->rows ? m->rows : 1, sizeof(char *));
    if (!text)
        goto cleanup;
    size_t count = 0;
    for (; count < m->rows; count++) {
        if ((count & 1023U) == 0 && lab_now() >= deadline)
            break;
        text[count] = mpq_get_str(NULL, 10, y[count]);
    }
    if (count != m->rows) {
        for (size_t i = 0; i < count; i++)
            free(text[i]);
        free(text);
        goto cleanup;
    }
    r->rows = m->rows;
    r->multiplier = text;
    r->margin_lower_bound = mpq_get_str(NULL, 10, margin);
    r->margin_kind = "source_enclosure";
    message(r, "Strict positive original-source contradiction lower bound");
    outcome = 1;
cleanup:
    mpq_clear(box_lo.value); mpq_clear(box_hi.value);
    mpq_clears(qlo, qhi, alo, ahi, product_lo, product_hi, factor, value, maximum, margin, cutoff, NULL);
    for (size_t i = 0; i < initialized; i++)
        mpq_clear(y[i]);
    free(y);
    if (!outcome && lab_now() >= deadline)
        message(r, "Farkas candidate deadline exhausted");
    goto done;
expired:
    outcome = 0;
    message(r, "Farkas candidate deadline exhausted");
done:
    r->elapsed = lab_now() - started;
    return outcome;
}
int lab_farkas_candidate(const LabModel *m, size_t terms, const size_t *source_row,
                          const double *scale, const double *pi, double deadline,
                          LabFarkasCandidate *r) {
    return construct_candidate(m, terms, source_row, scale, pi, deadline, 0, 0, r);
}
int lab_farkas_projected_candidate(const LabModel *m, size_t terms, const size_t *source_row,
                                    const double *scale, const double *pi, double deadline,
                                    LabFarkasCandidate *r) {
    return construct_candidate(m, terms, source_row, scale, pi, deadline, 1, 0, r);
}
int lab_farkas_simplified_candidate(const LabModel *m, size_t terms, const size_t *source_row,
                                     const double *scale, const double *pi, double deadline,
                                     LabFarkasCandidate *r) {
    return construct_candidate(m, terms, source_row, scale, pi, deadline, 1, 1, r);
}
void lab_farkas_candidate_free(LabFarkasCandidate *r) {
    if (!r)
        return;
    for (size_t i = 0; i < r->rows; i++)
        free(r->multiplier[i]);
    free(r->multiplier);
    free(r->margin_lower_bound);
    memset(r, 0, sizeof *r);
}

#include "exact.h"
#include "common.h"
#include <ctype.h>
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int lab_exact_decimal(mpq_t value, const char *text) {
    /* Optional provenance is deliberately bounded, including spellings like
     * 0e999999999 whose floating conversion need not overflow. */
    enum { DIGITS = 1024, EXPONENT = 4096 };
    if (!text || strlen(text) > DIGITS)
        return -1;
    const char *s = text;
    int negative = *s == '-';
    if (*s == '-' || *s == '+')
        s++;
    char digits[DIGITS + 1];
    size_t count = 0;
    long fraction = 0;
    int point = 0;
    while (isdigit((unsigned char)*s) || *s == '.') {
        if (*s == '.') {
            if (point)
                return -1;
            point = 1;
        } else {
            digits[count++] = *s;
            fraction += point;
        }
        s++;
    }
    if (!count)
        return -1;
    digits[count] = 0;
    long exponent = 0;
    if (*s == 'e' || *s == 'E' || *s == 'd' || *s == 'D') {
        s++;
        int neg_exp = *s == '-';
        if (*s == '-' || *s == '+')
            s++;
        if (!isdigit((unsigned char)*s))
            return -1;
        while (isdigit((unsigned char)*s)) {
            exponent = 10 * exponent + (*s++ - '0');
            if (exponent > EXPONENT)
                return -1;
        }
        if (neg_exp)
            exponent = -exponent;
    }
    if (*s)
        return -1;
    exponent -= fraction;
    if (exponent < -EXPONENT || exponent > EXPONENT)
        return -1;
    if (mpz_set_str(mpq_numref(value), digits, 10))
        return -1;
    if (negative)
        mpz_neg(mpq_numref(value), mpq_numref(value));
    mpz_set_ui(mpq_denref(value), 1);
    mpz_t power;
    mpz_init(power);
    mpz_ui_pow_ui(power, 10, (unsigned long)(exponent < 0 ? -exponent : exponent));
    if (exponent < 0)
        mpz_set(mpq_denref(value), power);
    else
        mpz_mul(mpq_numref(value), mpq_numref(value), power);
    mpz_clear(power);
    mpq_canonicalize(value);
    return 0;
}

static mpq_t *rationals(size_t count) {
    mpq_t *values = lab_alloc(count, sizeof *values);
    for (size_t k = 0; k < count; k++)
        mpq_init(values[k]);
    return values;
}
LabExactModel *lab_exact_model_new(size_t rows, size_t cols, size_t nnz) {
    LabExactModel *m = lab_alloc(1, sizeof *m);
    m->rows = rows;
    m->cols = cols;
    m->nnz = nnz;
    m->value = rationals(nnz);
    m->row_lower = rationals(rows);
    m->row_upper = rationals(rows);
    m->col_lower = rationals(cols);
    m->col_upper = rationals(cols);
    m->objective = rationals(cols);
    mpq_init(m->offset);
    m->col_lower_source = lab_alloc(cols, sizeof(double));
    m->col_upper_source = lab_alloc(cols, sizeof(double));
    return m;
}
static void free_rationals(mpq_t *values, size_t count) {
    for (size_t k = 0; k < count; k++)
        mpq_clear(values[k]);
    free(values);
}
void lab_exact_model_free(LabExactModel *m) {
    if (!m)
        return;
    free_rationals(m->value, m->nnz);
    free_rationals(m->row_lower, m->rows);
    free_rationals(m->row_upper, m->rows);
    free_rationals(m->col_lower, m->cols);
    free_rationals(m->col_upper, m->cols);
    free_rationals(m->objective, m->cols);
    mpq_clear(m->offset);
    free(m->col_lower_source);
    free(m->col_upper_source);
    free(m);
}

static void original_bound(const LabModel *m, size_t j, int upper, mpq_t q) {
    const LabExactModel *e = m->exact;
    double current = upper ? m->col_upper[j] : m->col_lower[j];
    if (e && current == (upper ? e->col_upper_source[j] : e->col_lower_source[j]))
        mpq_set(q, upper ? e->col_upper[j] : e->col_lower[j]);
    else
        mpq_set_d(q, current);
}
/* GMP's mpq_get_d truncates. Select the closest adjacent finite double for an
 * objective, or round upward for a nonnegative violation. */
static int rational_double(const mpq_t q, double *result, int upward) {
    mpq_t converted, error, other;
    mpq_inits(converted, error, other, NULL);
    mpq_set_d(converted, DBL_MAX);
    mpq_abs(error, q);
    int bad = mpq_cmp(error, converted) > 0;
    double value = 0;
    if (!bad) {
        value = mpq_get_d(q);
        bad = !isfinite(value);
    }
    if (!bad) {
        mpq_set_d(converted, value);
        int comparison = mpq_cmp(converted, q);
        if (upward && comparison < 0)
            value = nextafter(value, INFINITY);
        else if (!upward && comparison) {
            double adjacent = nextafter(value, comparison < 0 ? INFINITY : -INFINITY);
            if (isfinite(adjacent)) {
                mpq_sub(error, converted, q);
                mpq_abs(error, error);
                mpq_set_d(other, adjacent);
                mpq_sub(other, other, q);
                mpq_abs(other, other);
                if (mpq_cmp(other, error) < 0)
                    value = adjacent;
            }
        }
        bad = !isfinite(value);
    }
    mpq_clears(converted, error, other, NULL);
    if (bad)
        return -1;
    *result = value == 0 ? 0 : value;
    return 0;
}
int lab_exact_nearest(const mpq_t value, double *result) {
    return result ? rational_double(value, result, 0) : -1;
}
static void larger(mpq_t largest, const mpq_t candidate) {
    if (mpq_cmp(candidate, largest) > 0)
        mpq_set(largest, candidate);
}
static int evaluate_source(const LabModel *m, const double *x, double *objective,
                           double *violation) {
    const LabExactModel *e = m->exact;
    if (e->rows != m->rows || e->cols != m->cols || e->nnz != m->nnz)
        return -1;
    mpq_t *activity = rationals(m->rows);
    mpq_t value, term, side, cost, worst;
    mpq_inits(value, term, side, cost, worst, NULL);
    mpq_set(cost, e->offset);
    for (size_t j = 0; j < m->cols; j++) {
        mpq_set_d(value, x[j]);
        mpq_mul(term, e->objective[j], value);
        mpq_add(cost, cost, term);
        if (isfinite(m->col_lower[j])) {
            original_bound(m, j, 0, side);
            mpq_sub(term, side, value);
            larger(worst, term);
        }
        if (isfinite(m->col_upper[j])) {
            original_bound(m, j, 1, side);
            mpq_sub(term, value, side);
            larger(worst, term);
        }
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++) {
            mpq_mul(term, e->value[k], value);
            mpq_add(activity[m->row_index[k]], activity[m->row_index[k]], term);
        }
    }
    for (size_t i = 0; i < m->rows; i++) {
        if (isfinite(m->row_lower[i])) {
            mpq_sub(term, e->row_lower[i], activity[i]);
            larger(worst, term);
        }
        if (isfinite(m->row_upper[i])) {
            mpq_sub(term, activity[i], e->row_upper[i]);
            larger(worst, term);
        }
    }
    int bad = rational_double(cost, objective, 0) || rational_double(worst, violation, 1);
    free_rationals(activity, m->rows);
    mpq_clears(value, term, side, cost, worst, NULL);
    return bad ? -1 : 0;
}

/* A two-component expansion plus an outward-rounded bound on discarded terms.
 * No extra precision from long double is assumed. */
typedef struct {
    double high, low, error;
} Accumulator;
static int two_sum(double a, double b, double *sum, double *error) {
    double s = a + b;
    double b_virtual = s - a;
    double e = (a - (s - b_virtual)) + (b - b_virtual);
    if (!isfinite(s) || !isfinite(e))
        return -1;
    *sum = s;
    *error = e;
    return 0;
}
static int add_error(Accumulator *s, double error) {
    if (error > 0) {
        s->error = nextafter(s->error + error, INFINITY);
        if (!isfinite(s->error))
            return -1;
    }
    return 0;
}
static int accumulate(Accumulator *s, double term) {
    double high, low, first_error, discarded;
    if (two_sum(s->high, term, &high, &first_error) ||
        two_sum(s->low, first_error, &low, &discarded) ||
        two_sum(high, low, &s->high, &s->low) || add_error(s, fabs(discarded)))
        return -1;
    return 0;
}
static int product(Accumulator *s, double a, double b) {
    double rounded = a * b;
    double error = fma(a, b, -rounded);
    if (!isfinite(rounded) || !isfinite(error) || accumulate(s, rounded) ||
        accumulate(s, error))
        return -1;
    /* Away from underflow the FMA product error is exact. Near underflow its
     * lost tail has magnitude below one smallest subnormal. */
    if (a != 0 && b != 0 && fabs(rounded) < scalbn(DBL_MIN, DBL_MANT_DIG) &&
        add_error(s, DBL_TRUE_MIN))
        return -1;
    return 0;
}
static int enclosure(const Accumulator *s, double *lower, double *upper) {
    double value = s->high + s->low;
    if (!isfinite(value))
        return -1;
    *lower = s->low == 0 ? value : nextafter(value, -INFINITY);
    *upper = s->low == 0 ? value : nextafter(value, INFINITY);
    if (s->error > 0) {
        *lower = nextafter(*lower - s->error, -INFINITY);
        *upper = nextafter(*upper + s->error, INFINITY);
    }
    return isfinite(*lower) && isfinite(*upper) ? 0 : -1;
}
static int violation_difference(double a, double b, double *worst) {
    double difference, error;
    if (two_sum(a, -b, &difference, &error))
        return -1;
    if (error > 0)
        difference = nextafter(difference, INFINITY);
    if (!isfinite(difference))
        return -1;
    *worst = fmax(*worst, difference);
    return 0;
}
static int evaluate_binary(const LabModel *m, const double *x, double *objective,
                           double *violation) {
    if (fegetround() != FE_TONEAREST || !isfinite(m->offset))
        return -1;
    Accumulator *activity = lab_alloc(m->rows, sizeof *activity);
    Accumulator cost = {.high = m->offset};
    double worst = 0;
    int bad = 0;
    for (size_t j = 0; j < m->cols && !bad; j++) {
        bad = product(&cost, m->objective[j], x[j]);
        if (!bad && isfinite(m->col_lower[j]))
            bad = violation_difference(m->col_lower[j], x[j], &worst);
        if (!bad && isfinite(m->col_upper[j]))
            bad = violation_difference(x[j], m->col_upper[j], &worst);
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1] && !bad; k++)
            bad = product(&activity[m->row_index[k]], m->value[k], x[j]);
    }
    for (size_t i = 0; i < m->rows && !bad; i++) {
        double lower, upper;
        bad = enclosure(&activity[i], &lower, &upper);
        if (!bad && isfinite(m->row_lower[i]))
            bad = violation_difference(m->row_lower[i], lower, &worst);
        if (!bad && isfinite(m->row_upper[i]))
            bad = violation_difference(upper, m->row_upper[i], &worst);
    }
    double value = cost.high + cost.low;
    if (!isfinite(value))
        bad = 1;
    if (!bad && cost.error > 0) {
        /* A discarded third component can become significant after later
         * cancellation. Recover only this scalar binary objective exactly;
         * sparse row activities still use the bounded FMA enclosures. */
        mpq_t total, coefficient, point;
        mpq_inits(total, coefficient, point, NULL);
        mpq_set_d(total, m->offset);
        for (size_t j = 0; j < m->cols; j++) {
            mpq_set_d(coefficient, m->objective[j]);
            mpq_set_d(point, x[j]);
            mpq_mul(coefficient, coefficient, point);
            mpq_add(total, total, coefficient);
        }
        bad = rational_double(total, &value, 0);
        mpq_clears(total, coefficient, point, NULL);
    }
    free(activity);
    if (bad)
        return -1;
    *objective = value == 0 ? 0 : value;
    *violation = worst;
    return 0;
}
int lab_exact_evaluate(const LabModel *m, const double *x, double *objective,
                       double *max_violation) {
    if (!objective || !max_violation)
        return -1;
    *objective = NAN;
    *max_violation = INFINITY;
    if (!m || (m->cols && !x))
        return -1;
    for (size_t j = 0; j < m->cols; j++)
        if (!isfinite(x[j]) || isnan(m->col_lower[j]) || isnan(m->col_upper[j]) ||
            m->col_lower[j] == INFINITY || m->col_upper[j] == -INFINITY)
            return -1;
    for (size_t i = 0; i < m->rows; i++)
        if (isnan(m->row_lower[i]) || isnan(m->row_upper[i]) ||
            m->row_lower[i] == INFINITY || m->row_upper[i] == -INFINITY)
            return -1;
    int failed = m->exact ? evaluate_source(m, x, objective, max_violation)
                          : evaluate_binary(m, x, objective, max_violation);
    if (failed) {
        *objective = NAN;
        *max_violation = INFINITY;
    }
    return failed;
}

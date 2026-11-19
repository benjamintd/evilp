#include "exact.h"
#include "common.h"
#include <ctype.h>
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <errno.h>
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
static int evaluate_enclosed_source(const LabModel *,const double *,double *,double *);
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
                          : m->parsed_decimal ? evaluate_enclosed_source(m,x,objective,max_violation)
                          : evaluate_binary(m, x, objective, max_violation);
    if (failed) {
        *objective = NAN;
        *max_violation = INFINITY;
    }
    return failed;
}

static LabInterval whole_interval(void) {
    return (LabInterval){-INFINITY, INFINITY};
}
static LabInterval scalar_sum_interval(double a, double b) {
    if (isinf(a) || isinf(b)) {
        double sum = a + b;
        return isnan(sum) ? whole_interval() : (LabInterval){sum, sum};
    }
    double sum, error;
    if (two_sum(a, b, &sum, &error)) {
        sum = a + b;
        if (isinf(sum))
            return sum > 0 ? (LabInterval){DBL_MAX, INFINITY}
                           : (LabInterval){-INFINITY, -DBL_MAX};
        return whole_interval();
    }
    return (LabInterval){error < 0 ? nextafter(sum, -INFINITY) : sum,
                         error > 0 ? nextafter(sum, INFINITY) : sum};
}
LabInterval lab_interval_add(LabInterval a, LabInterval b) {
    if (fegetround() != FE_TONEAREST || isnan(a.lower) || isnan(a.upper) ||
        isnan(b.lower) || isnan(b.upper) || a.lower > a.upper || b.lower > b.upper)
        return whole_interval();
    LabInterval lo = scalar_sum_interval(a.lower, b.lower);
    LabInterval hi = scalar_sum_interval(a.upper, b.upper);
    return (LabInterval){lo.lower, hi.upper};
}
LabInterval lab_interval_subtract(LabInterval a, LabInterval b) {
    return lab_interval_add(a, (LabInterval){-b.upper, -b.lower});
}
static LabInterval scalar_product_interval(double a, double b) {
    if (a == 0 || b == 0)
        return (LabInterval){0, 0};
    double p = a * b;
    if (isinf(a) || isinf(b))
        return (LabInterval){p, p};
    if (isinf(p))
        return p > 0 ? (LabInterval){DBL_MAX, INFINITY}
                     : (LabInterval){-INFINITY, -DBL_MAX};
    double error = fma(a, b, -p);
    if (!isfinite(error))
        return whole_interval();
    int tiny = fabs(p) < scalbn(DBL_MIN, DBL_MANT_DIG);
    return (LabInterval){error < 0 || tiny ? nextafter(p, -INFINITY) : p,
                         error > 0 || tiny ? nextafter(p, INFINITY) : p};
}
LabInterval lab_interval_multiply(LabInterval a, LabInterval b) {
    if (fegetround() != FE_TONEAREST || isnan(a.lower) || isnan(a.upper) ||
        isnan(b.lower) || isnan(b.upper) || a.lower > a.upper || b.lower > b.upper)
        return whole_interval();
    LabInterval result = {INFINITY, -INFINITY};
    double av[] = {a.lower, a.upper}, bv[] = {b.lower, b.upper};
    for (size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < 2; j++) {
            LabInterval p = scalar_product_interval(av[i], bv[j]);
            result.lower = fmin(result.lower, p.lower);
            result.upper = fmax(result.upper, p.upper);
        }
    return result;
}
static LabInterval reciprocal_interval(double value) {
    if (isinf(value))
        return (LabInterval){0, 0};
    double r = 1 / value;
    if (isinf(r))
        return r > 0 ? (LabInterval){DBL_MAX, INFINITY}
                     : (LabInterval){-INFINITY, -DBL_MAX};
    /* Adjacent doubles enclose correctly rounded division even when an FMA
     * residual would itself underflow and lose its direction. */
    return (LabInterval){nextafter(r, -INFINITY), nextafter(r, INFINITY)};
}
LabInterval lab_interval_divide(LabInterval a, LabInterval b) {
    if (isnan(b.lower) || isnan(b.upper) || b.lower > b.upper ||
        (b.lower <= 0 && b.upper >= 0))
        return whole_interval();
    LabInterval lo = reciprocal_interval(b.upper), hi = reciprocal_interval(b.lower);
    return lab_interval_multiply(a, (LabInterval){lo.lower, hi.upper});
}
LabInterval lab_source_number_interval(const char *text, double rounded) {
    if (!text || !isfinite(rounded))
        return whole_interval();
    int nonzero = 0, integer = 1;
    const char *s = text;
    if (*s == '+' || *s == '-')
        s++;
    const char *digits = s;
    for (; *s; s++) {
        if (*s == 'e' || *s == 'E' || *s == 'd' || *s == 'D') {
            integer = 0;
            break;
        }
        if (*s == '.')
            integer = 0;
        else if (*s >= '1' && *s <= '9')
            nonzero = 1;
    }
    if (!nonzero)
        return (LabInterval){0, 0};
    if (integer && strlen(digits) <= 20) {
        errno = 0;
        char *end;
        unsigned long long magnitude = strtoull(digits, &end, 10);
        double absolute = fabs(rounded);
        if (!errno && !*end && absolute < 0x1p64 &&
            (unsigned long long)absolute == magnitude)
            return (LabInterval){rounded, rounded};
    }
    /* Only a successful bounded exact comparison may collapse this interval.
     * Allocation/parser limits retain the previous outward enclosure. */
    mpq_t source, binary;
    mpq_inits(source, binary, NULL);
    int represented = !lab_exact_decimal(source, text);
    if (represented) {
        mpq_set_d(binary, rounded);
        represented = !mpq_cmp(source, binary);
    }
    mpq_clears(source, binary, NULL);
    if (represented)
        return (LabInterval){rounded, rounded};
    return (LabInterval){nextafter(rounded, -INFINITY), nextafter(rounded, INFINITY)};
}
LabSourceEnvelope *lab_source_envelope_new(size_t rows, size_t cols, size_t nnz) {
    LabSourceEnvelope *s = lab_alloc(1, sizeof *s);
    s->rows = rows;
    s->cols = cols;
    s->nnz = nnz;
    s->value = lab_alloc(nnz, sizeof(LabInterval));
    s->row_lower = lab_alloc(rows, sizeof(LabInterval));
    s->row_upper = lab_alloc(rows, sizeof(LabInterval));
    s->col_lower = lab_alloc(cols, sizeof(LabInterval));
    s->col_upper = lab_alloc(cols, sizeof(LabInterval));
    s->col_lower_source = lab_alloc(cols, sizeof(double));
    s->col_upper_source = lab_alloc(cols, sizeof(double));
    s->complete = 1;
    s->objective=lab_alloc(cols,sizeof(LabInterval));
    return s;
}
void lab_source_envelope_free(LabSourceEnvelope *s) {
    if (!s)
        return;
    for (size_t k = 0; k < s->exception_count; ++k)
        mpq_clear(s->exceptions[k].value);
    free(s->exceptions);
    free(s->value);
    free(s->objective);
    free(s->row_lower);
    free(s->row_upper);
    free(s->col_lower);
    free(s->col_upper);
    free(s->col_lower_source);
    free(s->col_upper_source);
    free(s);
}

/* Exact source accessors never infer a rational from an uncertain enclosure.
 * A distinct column value is an explicit binary bound/shift overlay. */
static int source_envelope_valid(const LabModel *m) {
    const LabSourceEnvelope *s = m->envelope;
    return s && s->complete && s->rows == m->rows && s->cols == m->cols &&
        s->nnz == m->nnz && s->value && s->row_lower && s->row_upper &&
        s->col_lower && s->col_upper && s->col_lower_source && s->col_upper_source;
}
static const LabExactModel *native_source_snapshot(const LabModel *m) {
    return m->exact ? m->exact : m->native_exact;
}
static int source_exact_valid(const LabModel *m) {
    const LabExactModel *e=native_source_snapshot(m);
    return e && e->rows==m->rows && e->cols==m->cols && e->nnz==m->nnz;
}
static int exact_point(LabInterval side, double value) {
    return isfinite(value) && side.lower == value && side.upper == value;
}
int lab_source_exact_coefficient(const LabModel *m, size_t k, mpq_t out) {
    if (!m || k >= m->nnz || !isfinite(m->value[k]) || ((m->exact || m->native_exact) && !source_exact_valid(m))) return -1;
    if (source_exact_valid(m)) mpq_set(out,native_source_snapshot(m)->value[k]);
    else if (!m->parsed_decimal || (source_envelope_valid(m) &&
             exact_point(m->envelope->value[k], m->value[k]))) mpq_set_d(out, m->value[k]);
    else if (source_envelope_valid(m)) {
        const LabSourceEnvelope *s = m->envelope;
        size_t lo = 0, hi = s->exception_count;
        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if (s->exceptions[mid].entry < k) lo = mid + 1;
            else hi = mid;
        }
        if (lo == s->exception_count || s->exceptions[lo].entry != k ||
            s->exceptions[lo].rounded != m->value[k])
            return -1;
        mpq_set(out, s->exceptions[lo].value);
    }
    else return -1;
    return 0;
}
int lab_source_exact_row_bound(const LabModel *m, size_t i, int upper, mpq_t out) {
    if (!m || i >= m->rows || (upper != 0 && upper != 1) || ((m->exact || m->native_exact) && !source_exact_valid(m))) return -1;
    double value = upper ? m->row_upper[i] : m->row_lower[i];
    if (!isfinite(value)) return -1;
    if (source_exact_valid(m)) mpq_set(out,upper ? native_source_snapshot(m)->row_upper[i] : native_source_snapshot(m)->row_lower[i]);
    else if (!m->parsed_decimal || (source_envelope_valid(m) && exact_point(
        upper ? m->envelope->row_upper[i] : m->envelope->row_lower[i], value))) mpq_set_d(out,value);
    else return -1;
    return 0;
}
int lab_source_exact_column_value(const LabModel *m, size_t j, int upper,
                                  double value, mpq_t out) {
    if (!m || j >= m->cols || (upper != 0 && upper != 1) || !isfinite(value) ||
        ((m->exact || m->native_exact) && !source_exact_valid(m))) return -1;
    if (source_exact_valid(m)) {
        const LabExactModel *e = native_source_snapshot(m);
        if (value == (upper ? e->col_upper_source[j] : e->col_lower_source[j]))
            mpq_set(out, upper ? e->col_upper[j] : e->col_lower[j]);
        else mpq_set_d(out,value);
    } else if (!m->parsed_decimal) mpq_set_d(out,value);
    else if (source_envelope_valid(m)) {
        const LabSourceEnvelope *s = m->envelope;
        if (value == (upper ? s->col_upper_source[j] : s->col_lower_source[j]) &&
            !exact_point(upper ? s->col_upper[j] : s->col_lower[j],value)) return -1;
        mpq_set_d(out,value);
    } else return -1;
    return 0;
}
int lab_source_exact_column_bound(const LabModel *m, size_t j, int upper, mpq_t out) {
    if (!m || j >= m->cols) return -1;
    return lab_source_exact_column_value(m,j,upper,upper ? m->col_upper[j] : m->col_lower[j],out);
}

LabSourceObjective *lab_source_objective_new(size_t cols) {
    LabSourceObjective *o=lab_alloc(1,sizeof *o);o->cols=cols;
    o->value=rationals(cols);o->source=lab_alloc(cols,sizeof(double));mpq_init(o->offset);return o;
}
void lab_source_objective_free(LabSourceObjective *o) {
    if(!o)return;free_rationals(o->value,o->cols);free(o->source);mpq_clear(o->offset);free(o);
}
int lab_source_exact_objective(const LabModel *m,size_t j,mpq_t out) {
    if(!m||j>=m->cols||!isfinite(m->objective[j]))return -1;
    if(m->exact) {
        if(m->exact->cols!=m->cols)return -1;
        mpq_set(out,m->exact->objective[j]);return 0;
    }
    const LabSourceObjective *o=m->objective_exact;
    if(o) {
        if(o->cols!=m->cols||o->source[j]!=m->objective[j])return -1;
        mpq_set(out,o->value[j]);return 0;
    }
    if(!m->parsed_decimal) {mpq_set_d(out,m->objective[j]);return 0;}
    const LabSourceEnvelope *e=m->envelope;
    if(e&&e->complete&&e->objective_complete&&e->cols==m->cols&&e->objective&&
       e->objective[j].lower==m->objective[j]&&e->objective[j].upper==m->objective[j]) {
        mpq_set_d(out,m->objective[j]);return 0;
    }
    return -1;
}
int lab_source_exact_offset(const LabModel *m,mpq_t out) {
    if(!m||!isfinite(m->offset))return -1;
    if(m->exact){mpq_set(out,m->exact->offset);return 0;}
    if(m->objective_exact) {
        if(m->objective_exact->cols!=m->cols||m->objective_exact->offset_source!=m->offset)return -1;
        mpq_set(out,m->objective_exact->offset);return 0;
    }
    if(!m->parsed_decimal){mpq_set_d(out,m->offset);return 0;}
    const LabSourceEnvelope *e=m->envelope;
    if(e&&e->complete&&e->objective_complete&&e->offset.lower==m->offset&&e->offset.upper==m->offset) {
        mpq_set_d(out,m->offset);return 0;
    }
    return -1;
}
int lab_source_objective_needs_exact(const LabModel *m) {
    if(!m)return 1;
    if(!m->parsed_decimal&&!m->exact&&!m->objective_exact)return 0;
    const LabSourceEnvelope *e=m->envelope;
    if(e&&e->complete&&e->objective_complete&&e->cols==m->cols&&e->objective) {
        for(size_t j=0;j<m->cols;j++)if(e->objective[j].lower!=m->objective[j]||e->objective[j].upper!=m->objective[j])return 1;
        return 0;
    }
    if(!m->exact&&!m->objective_exact)return 1;
    mpq_t q,d;mpq_inits(q,d,NULL);int needed=0;
    for(size_t j=0;j<m->cols;j++) {
        if(lab_source_exact_objective(m,j,q)){needed=1;break;}
        mpq_set_d(d,m->objective[j]);if(mpq_cmp(q,d)){needed=1;break;}
    }
    mpq_clears(q,d,NULL);return needed;
}
static int source_scalar(const LabModel *m,const double *left,const double *right,double *out) {
    if(!m||!out||(m->cols&&!left))return -1;
    int fast=!lab_source_objective_needs_exact(m);
    if(fast&&!right) {
        mpq_t q,d;mpq_inits(q,d,NULL);
        if(lab_source_exact_offset(m,q))fast=0;
        else {mpq_set_d(d,m->offset);if(mpq_cmp(q,d))fast=0;}
        mpq_clears(q,d,NULL);
    }
    if(fast) {
        Accumulator total={0};
        if(!right) {
            mpq_t offset;mpq_init(offset);int bad=lab_source_exact_offset(m,offset);
            double d=0;if(!bad)bad=rational_double(offset,&d,0);mpq_clear(offset);
            if(bad)return -1;total.high=d;
        }
        int bad=0;
        for(size_t j=0;j<m->cols&&!bad;j++) {
            if(!isfinite(left[j])||(right&&!isfinite(right[j])))return -1;
            if(right) {
                double hi,lo;if(two_sum(left[j],-right[j],&hi,&lo))return -1;
                bad=product(&total,m->objective[j],hi)||product(&total,m->objective[j],lo);
            } else bad=product(&total,m->objective[j],left[j]);
        }
        if(bad)return -1;
        if(total.error==0){double value=total.high+total.low;if(!isfinite(value))return -1;*out=value==0?0:value;return 0;}
    }
    /* Charge the exact-work budget for nonzero products. When a complete
     * finite source enclosure is available, multiplying any value in it by
     * zero gives exactly zero, even for an unknown decimal coefficient. */
    size_t exact_products = 0;
    const LabSourceEnvelope *e = m->envelope;
    int enclosed = e && e->complete && e->objective_complete &&
                   e->cols == m->cols && e->objective;
    mpq_t total,c,x,y,term;mpq_inits(total,c,x,y,term,NULL);int bad=0;
    if(!right)bad=lab_source_exact_offset(m,total);
    for(size_t j=0;j<m->cols&&!bad;j++) {
        if(!isfinite(left[j])||(right&&!isfinite(right[j]))){bad=1;break;}
        if (enclosed) {
            LabInterval cost = e->objective[j];
            int inactive = left[j] == (right ? right[j] : 0);
            if ((inactive || (cost.lower == 0 && cost.upper == 0)) &&
                isfinite(cost.lower) && isfinite(cost.upper) &&
                cost.lower <= m->objective[j] && m->objective[j] <= cost.upper)
                continue;
        }
        if (++exact_products > 524288) {bad=1;break;}
        bad=lab_source_exact_objective(m,j,c);if(bad)break;
        mpq_set_d(x,left[j]);if(right){mpq_set_d(y,right[j]);mpq_sub(x,x,y);}
        mpq_mul(term,c,x);mpq_add(total,total,term);
        if(mpz_sizeinbase(mpq_numref(total),2)>16384||mpz_sizeinbase(mpq_denref(total),2)>16384)bad=1;
    }
    double value=0;if(!bad)bad=rational_double(total,&value,0);
    mpq_clears(total,c,x,y,term,NULL);if(bad)return -1;*out=value==0?0:value;return 0;
}
int lab_source_objective_value(const LabModel *m,const double *x,double *out) {return source_scalar(m,x,NULL,out);}
int lab_source_objective_difference(const LabModel *m,const double *x,const double *y,double *out) {
    if(!y)return -1;return source_scalar(m,x,y,out);
}

/* Evaluate all source enclosures first. Only ambiguous row scalars are rebuilt
 * exactly, within the independently bounded native source budget. This never
 * changes the global exact-basis eligibility or acceptance tolerances. */
static double row_enclosure_violation(const LabModel *m,size_t i,LabInterval activity) {
    const LabSourceEnvelope *e=m->envelope;
    double worst=0;
    if(isfinite(m->row_lower[i]))worst=fmax(worst,lab_interval_subtract(e->row_lower[i],activity).upper);
    if(isfinite(m->row_upper[i]))worst=fmax(worst,lab_interval_subtract(activity,e->row_upper[i]).upper);
    return worst;
}
static int scalar_bits_fit(const mpq_t q) {
    return mpz_sizeinbase(mpq_numref(q),2)<=16384 && mpz_sizeinbase(mpq_denref(q),2)<=16384;
}
static int evaluate_enclosed_source(const LabModel *m,const double *x,double *objective,double *violation) {
    if(fegetround()!=FE_TONEAREST||!source_envelope_valid(m))return -1;
    const LabSourceEnvelope *e=m->envelope;
    Accumulator *binary=lab_alloc(m->rows,sizeof *binary);
    LabInterval *uncertainty=lab_alloc(m->rows,sizeof *uncertainty);
    double *row_worst=lab_alloc(m->rows,sizeof *row_worst);
    size_t *selected=lab_alloc(m->rows,sizeof *selected);
    double worst=0;int bad=0;size_t count=0;
    mpq_t bound,point,term;mpq_inits(bound,point,term,NULL);
    for(size_t j=0;j<m->cols&&!bad;j++) {
        LabInterval value={x[j],x[j]};
        for(int upper=0;upper<2;upper++) {
            double side=upper?m->col_upper[j]:m->col_lower[j];if(!isfinite(side))continue;
            LabInterval source=upper?e->col_upper[j]:e->col_lower[j];
            if(side!=(upper?e->col_upper_source[j]:e->col_lower_source[j]))source=(LabInterval){side,side};
            double error=(upper?lab_interval_subtract(value,source):lab_interval_subtract(source,value)).upper;
            if(error>1e-12&&!lab_source_exact_column_bound(m,j,upper,bound)) {
                mpq_set_d(point,x[j]);if(upper)mpq_sub(term,point,bound);else mpq_sub(term,bound,point);
                if(mpq_sgn(term)<=0)error=0;else if(rational_double(term,&error,1))bad=1;
            }
            worst=fmax(worst,error);
        }
        for(size_t k=m->col_start[j];k<m->col_start[j+1]&&!bad;k++) {
            size_t i=m->row_index[k];if(i>=m->rows){bad=1;break;}
            if(!isfinite(m->row_lower[i])&&!isfinite(m->row_upper[i]))continue;
            bad=product(&binary[i],m->value[k],x[j]);
            LabInterval delta=lab_interval_subtract(e->value[k],(LabInterval){m->value[k],m->value[k]});
            uncertainty[i]=lab_interval_add(uncertainty[i],lab_interval_multiply(delta,value));
        }
    }
    for(size_t i=0;i<m->rows&&!bad;i++) {
        selected[i]=SIZE_MAX;
        if(!isfinite(m->row_lower[i])&&!isfinite(m->row_upper[i]))continue;
        double lo,hi;bad=enclosure(&binary[i],&lo,&hi);if(bad)break;
        LabInterval activity=lab_interval_add((LabInterval){lo,hi},uncertainty[i]);
        row_worst[i]=row_enclosure_violation(m,i,activity);
        if(row_worst[i]>1e-12&&count<2048)selected[i]=count++;
    }
    mpq_t *activity=rationals(count);unsigned char *known=lab_alloc(count,sizeof *known);
    memset(known,1,count);size_t terms=0;
    for(size_t j=0;j<m->cols&&!bad;j++) {
        if(x[j]==0)continue;mpq_set_d(point,x[j]);
        for(size_t k=m->col_start[j];k<m->col_start[j+1];k++) {
            size_t r=selected[m->row_index[k]];if(r==SIZE_MAX||!known[r])continue;
            if(terms++>=100000||lab_source_exact_coefficient(m,k,term)){known[r]=0;continue;}
            mpq_mul(term,term,point);mpq_add(activity[r],activity[r],term);
            if(!scalar_bits_fit(activity[r]))known[r]=0;
        }
    }
    for(size_t i=0;i<m->rows&&!bad;i++) {
        size_t r=selected[i];
        if(r!=SIZE_MAX&&known[r]) {
            mpq_set_ui(point,0,1);int available=1;
            for(int upper=0;upper<2;upper++) {
                double side=upper?m->row_upper[i]:m->row_lower[i];if(!isfinite(side))continue;
                if(lab_source_exact_row_bound(m,i,upper,bound)){available=0;break;}
                if(upper)mpq_sub(term,activity[r],bound);else mpq_sub(term,bound,activity[r]);
                larger(point,term);
            }
            if(available&&rational_double(point,&row_worst[i],1))bad=1;
        }
        worst=fmax(worst,row_worst[i]);
    }
    if(!bad)bad=!isfinite(worst)||lab_source_objective_value(m,x,objective);
    free_rationals(activity,count);free(known);free(binary);free(uncertainty);free(row_worst);free(selected);
    mpq_clears(bound,point,term,NULL);if(bad)return -1;*violation=worst;return 0;
}

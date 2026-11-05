#include "exact.h"
#include "common.h"
#include <ctype.h>
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
    free(m->col_lower_source);
    free(m->col_upper_source);
    free(m);
}

#include "solution.h"
#include "common.h"
#include "exact.h"
#include <math.h>
#include <stdlib.h>

/* Choose the shortest significant-digit decimal that roundtrips to the same
 * double. This avoids exporting .29999999999999999 for input .3. Independent
 * rational checking still decides whether that decimal is a valid witness. */
static void decimal_text(char text[64], double x) {
    for (int digits = 1; digits <= 17; digits++) {
        snprintf(text, 64, "%.*g", digits, x == 0 ? 0 : x);
        char *end;
        double parsed = strtod(text, &end);
        if (!*end && parsed == x)
            break;
    }
}
void lab_solution_decimal(FILE *f, double x) {
    char text[64];
    decimal_text(text, x);
    fputs(text, f);
}
static int solution_scalar_fits(const mpq_t q) {
    return mpz_sgn(mpq_denref(q)) > 0 && mpz_sizeinbase(mpq_numref(q), 2) <= 16384 &&
           mpz_sizeinbase(mpq_denref(q), 2) <= 16384;
}
/* The checker evaluates the exported decimal coordinates, which can differ
 * from their internal binary representatives. Accumulate that same source
 * objective before emitting its optional claim. A rounded double can lose a
 * material improvement next to a large constant. On unavailable provenance or
 * resource limits, omit the optional claim and let the checker score the point. */
static void write_source_objective(FILE *f, const LabModel *m, const double *point) {
    if (m->cols > 65536)
        return;
    mpq_t total, cost, coordinate, term;
    mpq_inits(total, cost, coordinate, term, NULL);
    int bad = lab_source_exact_offset(m, total) || !solution_scalar_fits(total);
    for (size_t j = 0; j < m->cols && !bad; j++) {
        char text[64];
        decimal_text(text, point[j]);
        bad = lab_source_exact_objective(m, j, cost) || lab_exact_decimal(coordinate, text) ||
              !solution_scalar_fits(cost);
        if (bad)
            break;
        mpq_mul(term, cost, coordinate);
        if (!solution_scalar_fits(term)) {
            bad = 1;
            break;
        }
        mpq_add(total, total, term);
        bad = !solution_scalar_fits(total);
    }
    mpz_t denominator, numerator;
    mpz_inits(denominator, numerator, NULL);
    if (!bad) {
        mpz_set(denominator, mpq_denref(total));
        mp_bitcnt_t twos = mpz_scan1(denominator, 0);
        mpz_fdiv_q_2exp(denominator, denominator, twos);
        unsigned long fives = 0;
        while (mpz_divisible_ui_p(denominator, 5)) {
            mpz_divexact_ui(denominator, denominator, 5);
            fives++;
        }
        unsigned long places = (unsigned long)twos > fives ? (unsigned long)twos : fives;
        if (mpz_cmp_ui(denominator, 1) || places > 4096) {
            bad = 1;
        } else {
            mpz_set(numerator, mpq_numref(total));
            mpz_mul_2exp(numerator, numerator, places - (unsigned long)twos);
            mpz_ui_pow_ui(denominator, 5, places - fives);
            mpz_mul(numerator, numerator, denominator);
            size_t digits = mpz_sizeinbase(numerator, 10);
            if (digits <= 4080) {
                char text[4098];
                mpz_get_str(text, 10, numerator);
                if (places)
                    fprintf(f, "=obj= %se-%lu\n", text, places);
                else
                    fprintf(f, "=obj= %s\n", text);
            }
        }
    }
    mpz_clears(denominator, numerator, NULL);
    mpq_clears(total, cost, coordinate, term, NULL);
}
int lab_solution_write_point(FILE *f, const LabModel *model, const double *point) {
    if (!f || !model || (model->cols && !point))
        return -1;
    for (size_t j = 0; j < model->cols; j++)
        if (!isfinite(point[j]))
            return -1;
    write_source_objective(f, model, point);
    for (size_t j = 0; j < model->cols; j++) {
        fprintf(f, "%s ", model->col_name[j]);
        lab_solution_decimal(f, point[j]);
        fputc('\n', f);
    }
    return ferror(f) ? -1 : 0;
}

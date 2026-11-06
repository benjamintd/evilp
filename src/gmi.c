#include "gmi.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#if defined(__FAST_MATH__)
#error "Exact-to-binary cut conversion requires IEEE finite/nextafter semantics"
#endif
_Static_assert(FLT_RADIX == 2 && DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024 &&
                   DBL_MIN_EXP == -1021,
               "This conversion prototype targets IEEE binary64 doubles");

typedef struct {
    const LabGmiBudget *budget;
    size_t bits;
} Guard;

static LabGmiStatus tick(const Guard *g) {
    double now = g->budget->now(g->budget->context);
    if (!isfinite(now))
        return LAB_GMI_INVALID;
    return now >= g->budget->deadline ? LAB_GMI_DEADLINE : LAB_GMI_OK;
}

static LabGmiStatus begin(Guard *g, size_t n, const LabGmiBudget *budget) {
    if (!budget || !budget->now || !isfinite(budget->deadline) ||
        !budget->max_bits || !budget->max_total_bits)
        return LAB_GMI_INVALID;
    if (n > budget->max_terms || n > SIZE_MAX / sizeof(mpq_t))
        return LAB_GMI_RESOURCE;
    g->budget = budget;
    g->bits = 0;
    return tick(g);
}

static LabGmiStatus check_q(Guard *g, const mpq_t q, int account) {
    if (mpz_sgn(mpq_denref(q)) <= 0)
        return LAB_GMI_INVALID;
    size_t numerator = mpz_sizeinbase(mpq_numref(q), 2);
    size_t denominator = mpz_sizeinbase(mpq_denref(q), 2);
    if (numerator > g->budget->max_bits || denominator > g->budget->max_bits)
        return LAB_GMI_RESOURCE;
    if (account) {
        if (numerator > SIZE_MAX - denominator ||
            numerator + denominator > g->budget->max_total_bits - g->bits)
            return LAB_GMI_RESOURCE;
        g->bits += numerator + denominator;
    }
    mpz_t gcd;
    mpz_init(gcd);
    mpz_gcd(gcd, mpq_numref(q), mpq_denref(q));
    int canonical = mpz_cmp_ui(gcd, 1) == 0;
    mpz_clear(gcd);
    return canonical ? tick(g) : LAB_GMI_INVALID;
}

void lab_gmi_cut_init(LabGmiCut *cut) {
    cut->n = 0;
    cut->coefficient = NULL;
    mpq_init(cut->rhs);
}

void lab_gmi_cut_clear(LabGmiCut *cut) {
    for (size_t j = 0; j < cut->n; ++j)
        mpq_clear(cut->coefficient[j]);
    free(cut->coefficient);
    mpq_clear(cut->rhs);
    cut->n = 0;
    cut->coefficient = NULL;
}

void lab_gmi_binary_init(LabGmiBinaryCut *cut) {
    cut->n = 0;
    cut->coefficient = NULL;
    cut->rhs = 0;
    mpq_inits(cut->error_lower, cut->adjusted_rhs, NULL);
}

void lab_gmi_binary_clear(LabGmiBinaryCut *cut) {
    free(cut->coefficient);
    mpq_clears(cut->error_lower, cut->adjusted_rhs, NULL);
    cut->n = 0;
    cut->coefficient = NULL;
}

void lab_gmi_box_init(LabGmiBox *box) {
    box->has_lower = box->has_upper = 0;
    mpq_inits(box->lower, box->upper, NULL);
}

void lab_gmi_box_clear(LabGmiBox *box) {
    mpq_clears(box->lower, box->upper, NULL);
}

const char *lab_gmi_status_name(LabGmiStatus status) {
    switch (status) {
    case LAB_GMI_OK: return "ok";
    case LAB_GMI_INTEGRAL: return "integral_beta";
    case LAB_GMI_RESOURCE: return "resource_limit";
    case LAB_GMI_DEADLINE: return "deadline";
    case LAB_GMI_INVALID: return "invalid_input";
    case LAB_GMI_UNBOUNDED_ROUNDING: return "unbounded_rounding_error";
    case LAB_GMI_NONFINITE_ROUNDING: return "nonfinite_rounding";
    }
    return "invalid_status";
}

static void fraction(mpq_t out, const mpq_t q, mpz_t floor_q, mpq_t temporary) {
    mpz_fdiv_q(floor_q, mpq_numref(q), mpq_denref(q));
    mpq_set_z(temporary, floor_q);
    mpq_sub(out, q, temporary);
}

LabGmiStatus lab_gmi_build(LabGmiCut *out, size_t n, const mpq_t *a,
                         const unsigned char *is_integer, const mpq_t beta,
                         const LabGmiBudget *budget) {
    if (!out || !beta || (n && (!a || !is_integer)))
        return LAB_GMI_INVALID;
    Guard guard;
    LabGmiStatus status = begin(&guard, n, budget);
    if (status != LAB_GMI_OK || (status = check_q(&guard, beta, 1)) != LAB_GMI_OK)
        return status;
    for (size_t j = 0; j < n; ++j) {
        if (is_integer[j] > 1)
            return LAB_GMI_INVALID;
        if ((status = check_q(&guard, a[j], 1)) != LAB_GMI_OK)
            return status;
    }

    LabGmiCut result;
    lab_gmi_cut_init(&result);
    mpq_t one_minus_f, fj, temporary;
    mpq_inits(one_minus_f, fj, temporary, NULL);
    mpz_t floor_q;
    mpz_init(floor_q);
    fraction(result.rhs, beta, floor_q, temporary);
    if (!mpq_sgn(result.rhs)) {
        status = LAB_GMI_INTEGRAL;
        goto finish;
    }
    mpq_set_ui(one_minus_f, 1, 1);
    mpq_sub(one_minus_f, one_minus_f, result.rhs);
    if ((status = check_q(&guard, result.rhs, 1)) != LAB_GMI_OK ||
        (status = check_q(&guard, one_minus_f, 0)) != LAB_GMI_OK)
        goto finish;
    result.coefficient = n ? calloc(n, sizeof *result.coefficient) : NULL;
    if (n && !result.coefficient) {
        status = LAB_GMI_RESOURCE;
        goto finish;
    }
    for (size_t j = 0; j < n; ++j) {
        mpq_init(result.coefficient[j]);
        ++result.n;
        if (is_integer[j]) {
            fraction(fj, a[j], floor_q, temporary);
            if (mpq_cmp(fj, result.rhs) <= 0) {
                mpq_set(result.coefficient[j], fj);
            } else {
                mpq_set_ui(temporary, 1, 1);
                mpq_sub(temporary, temporary, fj);
                mpq_mul(result.coefficient[j], result.rhs, temporary);
                mpq_div(result.coefficient[j], result.coefficient[j], one_minus_f);
            }
        } else if (mpq_sgn(a[j]) >= 0) {
            mpq_set(result.coefficient[j], a[j]);
        } else {
            mpq_mul(result.coefficient[j], result.rhs, a[j]);
            mpq_neg(result.coefficient[j], result.coefficient[j]);
            mpq_div(result.coefficient[j], result.coefficient[j], one_minus_f);
        }
        if ((status = check_q(&guard, result.coefficient[j], 1)) != LAB_GMI_OK)
            goto finish;
    }
    status = tick(&guard);
    if (status == LAB_GMI_OK) {
        LabGmiCut old = *out;
        *out = result;
        result = old;
    }
finish:
    mpz_clear(floor_q);
    mpq_clears(one_minus_f, fj, temporary, NULL);
    lab_gmi_cut_clear(&result);
    return status;
}

/* Obtain adjacent finite binary64 candidates bracketing q where possible.
 * Every claim about direction is checked exactly after mpq_set_d. */
static size_t adjacent(const mpq_t q, double candidate[2], mpq_t temporary) {
    double d = mpq_get_d(q);
    if (!isfinite(d))
        d = mpq_sgn(q) < 0 ? -DBL_MAX : DBL_MAX;
    candidate[0] = d;
    mpq_set_d(temporary, d);
    int comparison = mpq_cmp(temporary, q);
    if (!comparison)
        return 1;
    double other = nextafter(d, comparison < 0 ? INFINITY : -INFINITY);
    if (!isfinite(other))
        return 1;
    candidate[1] = other;
    return 2;
}

LabGmiStatus lab_gmi_round_box(LabGmiBinaryCut *out, size_t n, const mpq_t *a,
                             const mpq_t rhs, const LabGmiBox *box,
                             const LabGmiBudget *budget) {
    if (!out || !rhs || (n && (!a || !box)))
        return LAB_GMI_INVALID;
    Guard guard;
    LabGmiStatus status = begin(&guard, n, budget);
    if (status != LAB_GMI_OK || (status = check_q(&guard, rhs, 1)) != LAB_GMI_OK)
        return status;
    for (size_t j = 0; j < n; ++j) {
        if ((box[j].has_lower != 0 && box[j].has_lower != 1) ||
            (box[j].has_upper != 0 && box[j].has_upper != 1))
            return LAB_GMI_INVALID;
        if ((status = check_q(&guard, a[j], 1)) != LAB_GMI_OK)
            return status;
        if (box[j].has_lower &&
            (status = check_q(&guard, box[j].lower, 1)) != LAB_GMI_OK)
            return status;
        if (box[j].has_upper &&
            (status = check_q(&guard, box[j].upper, 1)) != LAB_GMI_OK)
            return status;
        if (box[j].has_lower && box[j].has_upper &&
            mpq_cmp(box[j].lower, box[j].upper) > 0)
            return LAB_GMI_INVALID;
    }
    LabGmiBinaryCut result;
    lab_gmi_binary_init(&result);
    result.coefficient = n ? calloc(n, sizeof *result.coefficient) : NULL;
    if (n && !result.coefficient) {
        lab_gmi_binary_clear(&result);
        return LAB_GMI_RESOURCE;
    }
    result.n = n;
    mpq_t binary, delta, magnitude, best_magnitude, term, selected_term;
    mpq_inits(binary, delta, magnitude, best_magnitude, term, selected_term, NULL);
    for (size_t j = 0; j < n; ++j) {
        double candidate[2];
        size_t count = adjacent(a[j], candidate, binary);
        int selected = 0;
        for (size_t k = 0; k < count; ++k) {
            double d = candidate[k];
            /* No coefficient deletion, sign reversal, or infinite arithmetic. */
            if (mpq_sgn(a[j]) && (d == 0 || (d < 0) != (mpq_sgn(a[j]) < 0)))
                continue;
            mpq_set_d(binary, d);
            mpq_sub(delta, binary, a[j]);
            int sign = mpq_sgn(delta);
            if ((sign > 0 && !box[j].has_lower) || (sign < 0 && !box[j].has_upper))
                continue;
            if ((status = check_q(&guard, delta, 0)) != LAB_GMI_OK)
                goto finish;
            if (sign)
                mpq_mul(term, delta, sign > 0 ? box[j].lower : box[j].upper);
            else
                mpq_set_ui(term, 0, 1);
            if ((status = check_q(&guard, term, 0)) != LAB_GMI_OK)
                goto finish;
            mpq_abs(magnitude, delta);
            if (!selected || mpq_cmp(magnitude, best_magnitude) < 0) {
                selected = 1;
                result.coefficient[j] = d;
                mpq_set(best_magnitude, magnitude);
                mpq_set(selected_term, term);
            }
        }
        if (!selected) {
            status = LAB_GMI_UNBOUNDED_ROUNDING;
            goto finish;
        }
        mpq_add(result.error_lower, result.error_lower, selected_term);
        if ((status = check_q(&guard, result.error_lower, 0)) != LAB_GMI_OK)
            goto finish;
    }
    mpq_add(result.adjusted_rhs, rhs, result.error_lower);
    if ((status = check_q(&guard, result.error_lower, 1)) != LAB_GMI_OK ||
        (status = check_q(&guard, result.adjusted_rhs, 1)) != LAB_GMI_OK)
        goto finish;
    double candidate[2];
    size_t count = adjacent(result.adjusted_rhs, candidate, binary);
    int selected = 0;
    for (size_t k = 0; k < count; ++k) {
        mpq_set_d(binary, candidate[k]);
        if (mpq_cmp(binary, result.adjusted_rhs) <= 0 &&
            (!selected || candidate[k] > result.rhs)) {
            selected = 1;
            result.rhs = candidate[k];
        }
    }
    if (!selected) {
        status = LAB_GMI_NONFINITE_ROUNDING;
        goto finish;
    }
    status = tick(&guard);
    if (status == LAB_GMI_OK) {
        LabGmiBinaryCut old = *out;
        *out = result;
        result = old;
    }
finish:
    mpq_clears(binary, delta, magnitude, best_magnitude, term, selected_term, NULL);
    lab_gmi_binary_clear(&result);
    return status;
}

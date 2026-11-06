#include "gmi.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <fenv.h>

typedef struct { size_t ticks, expires; } Clock;
static double now(void *context) {
    Clock *clock = context;
    return ++clock->ticks >= clock->expires ? 100 : 0;
}
static double bad_now(void *context) {
    (void)context;
    return NAN;
}

static void set(mpq_t q, const char *s) {
    assert(!mpq_set_str(q, s, 10));
    mpq_canonicalize(q);
}

static void eq(const mpq_t q, const char *s) {
    mpq_t expected;
    mpq_init(expected);
    set(expected, s);
    assert(!mpq_cmp(q, expected));
    mpq_clear(expected);
}

int main(void) {
    Clock clock = {0, 1000000};
    LabGmiBudget budget = {100, 16384, 1U << 20, 100, now, &clock};
    mpq_t beta, a[4], saved[4];
    mpq_init(beta);
    for (size_t j = 0; j < 4; ++j)
        mpq_inits(a[j], saved[j], NULL);
    const unsigned char integer[4] = {1, 1, 0, 0};
    set(beta, "7/3");
    set(a[0], "-5/4"); set(a[1], "-2/3");
    set(a[2], "-7/5"); set(a[3], "2/7");
    for (size_t j = 0; j < 4; ++j)
        mpq_set(saved[j], a[j]);
    LabGmiCut cut;
    lab_gmi_cut_init(&cut);
    assert(lab_gmi_build(&cut, 4, (const mpq_t *)a, integer, beta, &budget) == LAB_GMI_OK);
    eq(cut.rhs, "1/3");
    eq(cut.coefficient[0], "1/8"); eq(cut.coefficient[1], "1/3");
    eq(cut.coefficient[2], "7/10"); eq(cut.coefficient[3], "2/7");
    eq(beta, "7/3");
    for (size_t j = 0; j < 4; ++j)
        assert(!mpq_cmp(a[j], saved[j]));

    /* Output/input aliasing is supported, including its old rhs. */
    assert(lab_gmi_build(&cut, 4, (const mpq_t *)cut.coefficient, integer, cut.rhs,
                         &budget) == LAB_GMI_OK);
    eq(cut.coefficient[0], "1/8"); eq(cut.rhs, "1/3");
    mpq_t old_rhs;
    mpq_init(old_rhs);
    mpq_set(old_rhs, cut.rhs);
    mpq_t *old_coefficients = cut.coefficient;
    set(beta, "-2");
    assert(lab_gmi_build(&cut, 4, (const mpq_t *)a, integer, beta, &budget) == LAB_GMI_INTEGRAL);
    assert(cut.coefficient == old_coefficients && !mpq_cmp(old_rhs, cut.rhs));
    set(beta, "-7/3");
    assert(lab_gmi_build(&cut, 0, NULL, NULL, beta, &budget) == LAB_GMI_OK);
    eq(cut.rhs, "2/3"); assert(!cut.n);
    old_coefficients = cut.coefficient;
    mpq_set(old_rhs, cut.rhs);
    budget.max_terms = 3;
    assert(lab_gmi_build(&cut, 4, (const mpq_t *)a, integer, beta, &budget) == LAB_GMI_RESOURCE);
    budget.max_terms = 100;
    budget.max_bits = 1;
    assert(lab_gmi_build(&cut, 4, (const mpq_t *)a, integer, beta, &budget) == LAB_GMI_RESOURCE);
    budget.max_bits = 16384;
    budget.max_total_bits = 2;
    assert(lab_gmi_build(&cut, 4, (const mpq_t *)a, integer, beta, &budget) == LAB_GMI_RESOURCE);
    budget.max_total_bits = 1U << 20;
    clock.expires = clock.ticks + 3;
    assert(lab_gmi_build(&cut, 4, (const mpq_t *)a, integer, beta, &budget) == LAB_GMI_DEADLINE);
    clock.expires = 1000000;
    assert(cut.coefficient == old_coefficients && !mpq_cmp(old_rhs, cut.rhs));
    LabGmiBudget invalid_budget = budget;
    invalid_budget.now = bad_now;
    assert(lab_gmi_build(&cut, 4, (const mpq_t *)a, integer, beta, &invalid_budget) == LAB_GMI_INVALID);
    invalid_budget = budget; invalid_budget.max_bits = 0;
    assert(lab_gmi_build(&cut, 4, (const mpq_t *)a, integer, beta, &invalid_budget) == LAB_GMI_INVALID);
    unsigned char bad_flags[4] = {2, 0, 0, 0};
    assert(lab_gmi_build(&cut, 4, (const mpq_t *)a, bad_flags, beta, &budget) == LAB_GMI_INVALID);
    mpz_set_ui(mpq_numref(a[0]), 2); mpz_set_ui(mpq_denref(a[0]), 4);
    assert(lab_gmi_build(&cut, 4, (const mpq_t *)a, integer, beta, &budget) == LAB_GMI_INVALID);
    mpq_canonicalize(a[0]);

    LabGmiBox box[4];
    for (size_t j = 0; j < 4; ++j) {
        lab_gmi_box_init(&box[j]);
        box[j].has_lower = 1;
    }
    LabGmiBinaryCut binary;
    lab_gmi_binary_init(&binary);
    set(a[0], "1/3"); set(beta, "1/3");
    assert(lab_gmi_round_box(&binary, 1, (const mpq_t *)a, beta, box, &budget) == LAB_GMI_OK);
    mpq_set_d(saved[0], binary.coefficient[0]);
    assert(mpq_cmp(saved[0], a[0]) > 0);
    eq(binary.error_lower, "0");
    mpq_set_d(saved[1], binary.rhs);
    assert(mpq_cmp(saved[1], beta) <= 0);
    double *old_binary = binary.coefficient;
    double old_binary_rhs = binary.rhs;
    box[0].has_lower = 0;
    assert(lab_gmi_round_box(&binary, 1, (const mpq_t *)a, beta, box, &budget) ==
           LAB_GMI_UNBOUNDED_ROUNDING);
    assert(binary.coefficient == old_binary && binary.rhs == old_binary_rhs);
    set(a[0], "1/2");
    assert(lab_gmi_round_box(&binary, 1, (const mpq_t *)a, beta, box, &budget) == LAB_GMI_OK);
    assert(binary.coefficient[0] == .5);

    /* Direction is established by exact comparisons, independent of fenv. */
    int original_round = fegetround();
    const int modes[] = {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST};
    set(a[0], "1/3"); box[0].has_lower = 1;
    for (size_t k = 0; k < sizeof modes / sizeof *modes; ++k) {
        assert(!fesetround(modes[k]));
        assert(lab_gmi_round_box(&binary, 1, (const mpq_t *)a, beta, box, &budget) == LAB_GMI_OK);
        mpq_set_d(saved[0], binary.coefficient[0]);
        assert(mpq_cmp(saved[0], a[0]) >= 0);
        mpq_set_d(saved[1], binary.rhs);
        assert(mpq_cmp(saved[1], beta) <= 0);
    }
    assert(!fesetround(original_round));

    /* Subnormal exact coefficient is retained rather than silently deleted. */
    mpq_set_ui(a[0], 1, 1); mpz_mul_2exp(mpq_denref(a[0]), mpq_denref(a[0]), 2000);
    box[0].has_lower = 1;
    assert(lab_gmi_round_box(&binary, 1, (const mpq_t *)a, beta, box, &budget) == LAB_GMI_OK);
    assert(binary.coefficient[0] == DBL_TRUE_MIN);
    old_binary = binary.coefficient;
    box[0].has_lower = 0; box[0].has_upper = 1;
    assert(lab_gmi_round_box(&binary, 1, (const mpq_t *)a, beta, box, &budget) ==
           LAB_GMI_UNBOUNDED_ROUNDING);
    assert(binary.coefficient == old_binary);

    /* Large coefficient is representable conservatively on a fixed zero box. */
    mpq_set_ui(a[0], 1, 1); mpz_mul_2exp(mpq_numref(a[0]), mpq_numref(a[0]), 2000);
    box[0].has_lower = 1;
    assert(lab_gmi_round_box(&binary, 1, (const mpq_t *)a, beta, box, &budget) == LAB_GMI_OK);
    assert(binary.coefficient[0] == DBL_MAX);
    box[0].has_upper = 0;
    assert(lab_gmi_round_box(&binary, 1, (const mpq_t *)a, beta, box, &budget) ==
           LAB_GMI_UNBOUNDED_ROUNDING);

    set(a[0], "1/10"); set(beta, "-10000000000000000000");
    box[0].has_lower = box[0].has_upper = 1;
    set(box[0].lower, "-100000000000000000000");
    mpq_set(box[0].upper, box[0].lower);
    assert(lab_gmi_round_box(&binary, 1, (const mpq_t *)a, beta, box, &budget) == LAB_GMI_OK);
    assert(mpq_sgn(binary.error_lower) < 0);
    mpq_set_d(saved[0], binary.coefficient[0]); mpq_mul(saved[0], saved[0], box[0].lower);
    assert(mpq_cmp(saved[0], beta) < 0); /* Naive conversion cuts this exact witness. */
    mpq_set_d(saved[1], binary.rhs);
    assert(mpq_cmp(saved[0], saved[1]) >= 0);
    old_binary = binary.coefficient;
    clock.expires = clock.ticks + 8;
    assert(lab_gmi_round_box(&binary, 1, (const mpq_t *)a, beta, box, &budget) == LAB_GMI_DEADLINE);
    assert(binary.coefficient == old_binary);
    clock.expires = 1000000;
    set(box[0].upper, "-100000000000000000001");
    assert(lab_gmi_round_box(&binary, 1, (const mpq_t *)a, beta, box, &budget) == LAB_GMI_INVALID);
    assert(binary.coefficient == old_binary);
    mpq_set_ui(beta, 1, 1); mpz_mul_2exp(mpq_numref(beta), mpq_numref(beta), 2000);
    mpq_neg(beta, beta);
    assert(lab_gmi_round_box(&binary, 0, NULL, beta, NULL, &budget) == LAB_GMI_NONFINITE_ROUNDING);

    lab_gmi_binary_clear(&binary);
    for (size_t j = 0; j < 4; ++j) {
        lab_gmi_box_clear(&box[j]);
        mpq_clears(a[j], saved[j], NULL);
    }
    mpq_clears(beta, old_rhs, NULL);
    lab_gmi_cut_clear(&cut);
    puts("Exact GMI and safe binary conversion contracts passed.");
    return 0;
}

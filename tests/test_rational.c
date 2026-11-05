#include "rational.h"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>

static void require(int condition, const char *message) {
    if (!condition) {
        fprintf(stderr, "exact reconstruction contract: %s\n", message);
        exit(1);
    }
}
int main(void) {
    mpq_t a[9], b[3], sum, product;
    for (size_t i = 0; i < 9; i++)
        mpq_init(a[i]);
    for (size_t i = 0; i < 3; i++)
        mpq_init(b[i]);
    mpq_inits(sum, product, NULL);
    /* The first pivot requires a row exchange. Solutions are non-dyadic. */
    const long matrix[] = {0, 2, 1, 3, -1, 2, 1, 1, -4};
    for (size_t i = 0; i < 9; i++)
        mpq_set_si(a[i], matrix[i], 1);
    for (size_t i = 0; i < 3; i++)
        mpq_set_si(b[i], (long)i + 1, 1);
    require(!lab_rational_solve(3, a, b, lab_now() + 5), "exact nonsingular solve");
    for (size_t i = 0; i < 3; i++) {
        mpq_set_ui(sum, 0, 1);
        for (size_t j = 0; j < 3; j++) {
            mpq_set_si(product, matrix[i * 3 + j], 1);
            mpq_mul(product, product, b[j]);
            mpq_add(sum, sum, product);
        }
        require(mpq_cmp_si(sum, (long)i + 1, 1) == 0, "exact original equations");
    }
    for (size_t i = 0; i < 9; i++)
        mpq_set_ui(a[i], 0, 1);
    require(lab_rational_solve(3, a, b, lab_now() + 5), "singularity rejects proof");
    mpq_set_ui(a[0], 1, 1);
    mpz_mul_2exp(mpq_numref(a[0]), mpq_numref(a[0]), LAB_RATIONAL_MAX_BITS + 1);
    require(lab_rational_solve(1, a, b, lab_now() + 5), "bit budget is bounded");
    require(lab_rational_solve(0, NULL, NULL, lab_now() - 1), "expired budget");
    require(lab_rational_solve(LAB_RATIONAL_MAX_DIMENSION + 1, NULL, NULL, lab_now() + 1),
            "dimension cap");
    require(!lab_rational_solve(0, NULL, NULL, lab_now() + 1), "empty exact solve");
    for (size_t i = 0; i < 9; i++)
        mpq_clear(a[i]);
    for (size_t i = 0; i < 3; i++)
        mpq_clear(b[i]);
    mpq_clears(sum, product, NULL);
    puts("Bounded exact reconstruction contracts passed.");
    return 0;
}

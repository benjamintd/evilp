#ifndef LAB_GMI_H
#define LAB_GMI_H

#include <gmp.h>
#include <stddef.h>

/* All rationals must be initialized and canonical. */
typedef enum {
    LAB_GMI_OK = 0,
    LAB_GMI_INTEGRAL,
    LAB_GMI_RESOURCE,
    LAB_GMI_DEADLINE,
    LAB_GMI_INVALID,
    LAB_GMI_UNBOUNDED_ROUNDING,
    LAB_GMI_NONFINITE_ROUNDING
} LabGmiStatus;

typedef struct {
    size_t max_terms, max_bits, max_total_bits;
    double deadline;
    double (*now)(void *context);
    void *context;
} LabGmiBudget;

typedef struct {
    size_t n;
    mpq_t *coefficient;
    mpq_t rhs;
} LabGmiCut;

typedef struct {
    int has_lower, has_upper;
    mpq_t lower, upper;
} LabGmiBox;

typedef struct {
    size_t n;
    double *coefficient;
    double rhs;
    mpq_t error_lower; /* exact inf_box sum((double coefficient - exact)*x) */
    mpq_t adjusted_rhs; /* exact original rhs + error_lower */
} LabGmiBinaryCut;

void lab_gmi_cut_init(LabGmiCut *cut);
void lab_gmi_cut_clear(LabGmiCut *cut);
void lab_gmi_binary_init(LabGmiBinaryCut *cut);
void lab_gmi_binary_clear(LabGmiBinaryCut *cut);
void lab_gmi_box_init(LabGmiBox *box);
void lab_gmi_box_clear(LabGmiBox *box);
const char *lab_gmi_status_name(LabGmiStatus status);

/* x_B + sum(a_j*z_j) = beta; x_B integer; all z_j >= 0.
 * is_integer[j] is exactly 0 or 1. No implicit integrality inference.
 * Failure leaves initialized out unchanged. Inputs may alias old out.
 * An integral beta is unavailable, including the zero-dimensional case.
 */
LabGmiStatus lab_gmi_build(LabGmiCut *out, size_t n, const mpq_t *a,
                         const unsigned char *is_integer, const mpq_t beta,
                         const LabGmiBudget *budget);

/* General exact cut sum(a_j*x_j) >= rhs, valid throughout the explicit box.
 * Box bounds need not be nonnegative. This is a validity conversion, not a
 * floating dot-product error bound. Finite, same-sign, nonzero doubles are
 * required for every nonzero exact coefficient. Missing required box sides
 * cause abstention. Exact zero coefficients remain zero. Failure is atomic.
 */
LabGmiStatus lab_gmi_round_box(LabGmiBinaryCut *out, size_t n, const mpq_t *a,
                             const mpq_t rhs, const LabGmiBox *box,
                             const LabGmiBudget *budget);

#endif

#ifndef LAB_FACTOR_H
#define LAB_FACTOR_H
#include <stddef.h>
/* Reference backend. No numerical policy or LP model semantics live here. */
typedef struct {
    size_t n;
    double *lu;
    size_t *pivot;
} LabFactor;
int lab_factorize(LabFactor *f, const double *matrix, size_t n);
int lab_factor_solve(const LabFactor *f, const double *rhs, double *x, int transpose);
void lab_factor_free(LabFactor *f);
#endif

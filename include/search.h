#ifndef LAB_SEARCH_H
#define LAB_SEARCH_H
#include "model.h"
#include <stddef.h>

/* These statistics rank branches only. They establish no bound or proof and
 * must never be used to discard a node. Initialize a zeroed object once, then
 * free it before reinitialization. Original variable indices are preserved. */
typedef struct {
    long double down_mean, up_mean;
    size_t down_count, up_count;
} LabSearchCost;
typedef struct {
    size_t cols;
    LabSearchCost *cost;
} LabSearch;
typedef enum { LAB_SEARCH_DOWN = -1, LAB_SEARCH_UP = 1 } LabSearchDirection;
typedef enum {
    LAB_SEARCH_CHOSEN = 0,
    LAB_SEARCH_INTEGRAL = 1, /* no fraction beyond tolerance; not feasibility */
    LAB_SEARCH_BAD = -1, /* invalid point/bounds/API; retain and resolve node */
    LAB_SEARCH_UNSUPPORTED = -2 /* integer value outside exact double spacing */
} LabSearchChoice;
#define LAB_SEARCH_RELIABLE 2

/* 0 success; -1 allocation, overflow or invalid initialization. */
int lab_search_init(LabSearch *search, size_t cols);
void lab_search_free(LabSearch *search);
/* 1 recorded, 0 ignored nonfinite/negative/unrepresentable observation,
 * -1 invalid API arguments. Zero degradation is a valid observation. Pass an
 * offset-free normalized child-minus-parent bound difference, never raw totals.
 * Distance is the positive movement to floor/ceil of the parent value (<=1).
 * Infeasible children and unsolved children supply no observation. */
int lab_search_observe(LabSearch *search, size_t variable, LabSearchDirection direction,
                       double distance, long double degradation);
/* Continuous variables are excluded. Bounds describe the current node in
 * original units. An absent/incomplete history uses most fractionality;
 * learned ranking maximizes balanced predicted bound improvements. Ties use
 * fractionality then original index. History never changes eligibility.
 * A NULL or dimension-mismatched history falls back safely. */
LabSearchChoice lab_search_choose(const LabModel *model, const double *primal,
                                  const double *lower, const double *upper,
                                  const LabSearch *search, double integrality_tolerance,
                                  size_t *variable);
#endif

#ifndef LAB_PRICING_H
#define LAB_PRICING_H
#include <stddef.h>
#include <stdint.h>

#define LAB_PRICING_SHORT_STREAK 8

/* Phase-local ranking state. It never decides eligibility, feasibility,
 * optimality, unboundedness, pruning, or whether a pivot is numerically safe. */
typedef struct {
    int bland_requested, bland_active;
    int adaptive;
    uint64_t history[64];
    size_t history_count, history_next;
    unsigned consecutive_short;
    size_t accepted_pivots, zero_steps, short_steps, fallback_transitions;
} LabPricing;
typedef struct {
    size_t column;
    long double reduced;
} LabPricingChoice;
typedef struct {
    int valid;
    size_t position, entering_column;
} LabPricingPending;
enum {
    LAB_PRICING_ZERO = 1,
    LAB_PRICING_SHORT = 2, /* positive short steps; zeros counted separately */
    LAB_PRICING_FALLBACK = 4
};

/* Begin each primal simplex phase afresh, including forced Bland mode. */
void lab_pricing_begin(LabPricing *state, int bland_requested);
/* Adaptive mode tolerates degeneracy until a repeated basis or 512
 * consecutive nonprogress steps. Fingerprint collisions only select the
 * conservative rule; fingerprints never establish a mathematical claim. */
void lab_pricing_begin_adaptive(LabPricing *state, int bland_requested);
int lab_pricing_watch_basis(LabPricing *state, size_t n, const size_t *basis);
void lab_pricing_scan_begin(LabPricingChoice *choice);
/* Caller supplies EXISTING eligibility after all unchanged finite/tolerance
 * checks. Entire pricing scan must still run. Return 1 changed best, 0 unchanged,
 * -1 invalid API/nonfinite data. Ties use the lower canonical column index. */
int lab_pricing_consider(const LabPricing *state, LabPricingChoice *choice,
                          size_t column, long double reduced, int eligible);
/* Record only a validated accepted canonical step. Zero triggers conservative
 * permanent fallback immediately; eight consecutive POSITIVE short steps do
 * likewise. Pass the existing LAB_PRIMAL_TOLERANCE as short_threshold; it is
 * used solely to choose a ranking policy, never as a new feasibility cutoff.
 * Return event bits above, or -1 without state changes on invalid input. */
int lab_pricing_record_step(LabPricing *state, double accepted_step, double short_threshold);

/* Queue the FINAL selected pivot position/entering index. Rejected trials must
 * not survive into this event. In the deferred no-update/Bland-ratio path,
 * keep it pending until the next refactor has validated the new basis. */
int lab_pricing_queue(LabPricingPending *pending, size_t position, size_t entering_column);
void lab_pricing_discard(LabPricingPending *pending);
/* Return 1 recorded, 0 absent/deferred, -1 incoherent event. A deferred event
 * changes no state. A committed event reads xb at its accepted position and
 * verifies the entering column actually occupies that position. It is cleared
 * only on successful recording, so double finish calls cannot double count. */
int lab_pricing_finish(LabPricing *state, LabPricingPending *pending, int basis_validated,
                       size_t basis_size, const size_t *basis, const double *xb,
                       double short_threshold, unsigned *events);
#endif

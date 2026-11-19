#include "pricing.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "pricing contract line %d: %s\n", __LINE__, #condition); \
        exit(1); \
    } \
} while (0)

static void ranking(void) {
    LabPricing state;
    LabPricingChoice best;
    lab_pricing_begin(&state, 0);
    lab_pricing_scan_begin(&best);
    CHECK(best.column == SIZE_MAX);
    CHECK(lab_pricing_consider(&state, &best, 9, -4, 1) == 1);
    CHECK(lab_pricing_consider(&state, &best, 1, -2, 1) == 0);
    CHECK(lab_pricing_consider(&state, &best, 5, -4, 1) == 1);
    CHECK(lab_pricing_consider(&state, &best, 0, -100, 0) == 0);
    CHECK(best.column == 5 && best.reduced == -4);
    CHECK(lab_pricing_consider(&state, &best, 3, -5, 1) == 1);
    CHECK(best.column == 3);
    CHECK(lab_pricing_consider(&state, &best, 2, NAN, 0) == -1);
    CHECK(best.column == 3 && best.reduced == -5);
    CHECK(lab_pricing_consider(&state, &best, SIZE_MAX, -6, 1) == -1);
    lab_pricing_begin(&state, 1);
    lab_pricing_scan_begin(&best);
    CHECK(lab_pricing_consider(&state, &best, 9, -4, 1) == 1);
    CHECK(lab_pricing_consider(&state, &best, 1, -2, 1) == 1);
    CHECK(lab_pricing_consider(&state, &best, 5, -100, 1) == 0);
    CHECK(best.column == 1 && best.reduced == -2);
    /* The helper adds no absolute reduced-cost cutoff of its own. */
    lab_pricing_scan_begin(&best);
    CHECK(lab_pricing_consider(&state, &best, 8, -1e-300L, 1) == 1);
    CHECK(best.column == 8 && best.reduced < 0);
}
static void accepted_steps(void) {
    LabPricing state;
    lab_pricing_begin(&state, 0);
    const double threshold = 1e-8;
    CHECK(lab_pricing_record_step(&state, .5, threshold) == 0);
    CHECK(state.accepted_pivots == 1 && !state.bland_active);
    CHECK(lab_pricing_record_step(&state, 0, threshold) ==
          (LAB_PRICING_ZERO | LAB_PRICING_FALLBACK));
    CHECK(state.accepted_pivots == 2 && state.zero_steps == 1 && state.short_steps == 0);
    CHECK(state.bland_active && state.fallback_transitions == 1);
    CHECK(lab_pricing_record_step(&state, 1, threshold) == 0);
    CHECK(state.bland_active && state.fallback_transitions == 1);
    CHECK(lab_pricing_record_step(&state, -1, threshold) == -1);
    CHECK(lab_pricing_record_step(&state, NAN, threshold) == -1);
    CHECK(lab_pricing_record_step(&state, INFINITY, threshold) == -1);
    CHECK(lab_pricing_record_step(&state, 1, -1) == -1);
    CHECK(lab_pricing_record_step(&state, 1, NAN) == -1);
    CHECK(state.accepted_pivots == 3 && state.zero_steps == 1 && state.short_steps == 0);
    CHECK(state.fallback_transitions == 1);
    /* New phase/cold solve starts with fresh policy and fresh local counters. */
    lab_pricing_begin(&state, 0);
    CHECK(!state.bland_active && !state.accepted_pivots && !state.fallback_transitions);
    for (unsigned i = 0; i < LAB_PRICING_SHORT_STREAK - 1; i++)
        CHECK(lab_pricing_record_step(&state, threshold, threshold) == LAB_PRICING_SHORT);
    CHECK(!state.bland_active && state.consecutive_short == 7);
    CHECK(lab_pricing_record_step(&state, 1, threshold) == 0);
    CHECK(state.consecutive_short == 0);
    for (unsigned i = 0; i < LAB_PRICING_SHORT_STREAK; i++) {
        int events = lab_pricing_record_step(&state, DBL_TRUE_MIN, threshold);
        CHECK((events & LAB_PRICING_SHORT) != 0);
        CHECK(((events & LAB_PRICING_FALLBACK) != 0) == (i + 1 == LAB_PRICING_SHORT_STREAK));
    }
    CHECK(state.bland_active && state.fallback_transitions == 1 && state.zero_steps == 0);
    CHECK(state.short_steps == 15 && state.accepted_pivots == 16);
    lab_pricing_begin(&state, 1);
    CHECK(lab_pricing_record_step(&state, 0, threshold) == LAB_PRICING_ZERO);
    CHECK(state.bland_requested && state.bland_active && !state.fallback_transitions);
    state.accepted_pivots = state.zero_steps = SIZE_MAX;
    CHECK(lab_pricing_record_step(&state, 0, threshold) == LAB_PRICING_ZERO);
    CHECK(state.accepted_pivots == SIZE_MAX && state.zero_steps == SIZE_MAX);
}
static void final_validated_event(void) {
    LabPricing state;
    LabPricingPending pending = {0};
    lab_pricing_begin(&state, 0);
    size_t basis[] = {3, 7, 9};
    double xb[] = {0, .5, 1};
    unsigned events = 99;
    /* A zero-step ordinary trial is rejected. It must not survive as an event. */
    CHECK(lab_pricing_queue(&pending, 0, 7) == 0);
    CHECK(lab_pricing_finish(&state, &pending, 1, 3, basis, xb, 1e-8, &events) == -1);
    CHECK(pending.valid && !state.accepted_pivots && !state.bland_active && events == 0);
    lab_pricing_discard(&pending);
    /* Harris alternative actually accepts position1 at step.5. A stale
     * original leaving position would have incorrectly triggered fallback. */
    CHECK(lab_pricing_queue(&pending, 1, 7) == 0);
    CHECK(lab_pricing_queue(&pending, 2, 9) == -1);
    CHECK(lab_pricing_finish(&state, &pending, 0, 0, NULL, NULL, 1e-8, &events) == 0);
    CHECK(pending.valid && !state.accepted_pivots && events == 0);
    CHECK(lab_pricing_finish(&state, &pending, 1, 3, basis, xb, 1e-8, &events) == 1);
    CHECK(!pending.valid && state.accepted_pivots == 1 && events == 0 && !state.bland_active);
    CHECK(lab_pricing_finish(&state, &pending, 1, 3, basis, xb, 1e-8, &events) == 0);
    CHECK(state.accepted_pivots == 1);
    /* Deferred no-update/Bland-ratio path: only post-refactor xb is read. */
    CHECK(lab_pricing_queue(&pending, 2, 9) == 0);
    CHECK(lab_pricing_finish(&state, &pending, 0, 3, basis, xb, 1e-8, &events) == 0);
    CHECK(pending.valid && state.accepted_pivots == 1);
    xb[2] = 0;
    CHECK(lab_pricing_finish(&state, &pending, 1, 3, basis, xb, 1e-8, &events) == 1);
    CHECK(state.accepted_pivots == 2 && state.zero_steps == 1 && state.bland_active);
    CHECK(events == (LAB_PRICING_ZERO | LAB_PRICING_FALLBACK));
    CHECK(lab_pricing_queue(&pending, 1, 7) == 0);
    xb[1] = NAN;
    CHECK(lab_pricing_finish(&state, &pending, 1, 3, basis, xb, 1e-8, &events) == -1);
    CHECK(pending.valid && state.accepted_pivots == 2 && events == 0);
    lab_pricing_discard(&pending);
    CHECK(lab_pricing_queue(&pending, SIZE_MAX, 7) == -1);
    CHECK(lab_pricing_consider(NULL, NULL, 0, 0, 0) == -1);
    CHECK(lab_pricing_record_step(NULL, 0, 0) == -1);
    CHECK(lab_pricing_finish(NULL, &pending, 0, 0, NULL, NULL, 0, &events) == -1);
}
static void adaptive_guards(void) {
    LabPricing state;
    lab_pricing_begin_adaptive(&state, 0);
    size_t basis[] = {1, 2, 3};
    CHECK(lab_pricing_watch_basis(&state, 3, basis) == 0);
    CHECK(lab_pricing_record_step(&state, 0, 1e-8) == LAB_PRICING_ZERO);
    CHECK(!state.bland_active && state.zero_steps == 1);
    basis[1] = 4;
    CHECK(lab_pricing_watch_basis(&state, 3, basis) == 0);
    basis[1] = 2;
    CHECK(lab_pricing_watch_basis(&state, 3, basis) == 1);
    CHECK(state.bland_active && state.fallback_transitions == 1);
    CHECK(lab_pricing_watch_basis(&state, 3, basis) == 0);
    lab_pricing_begin_adaptive(&state, 0);
    for (size_t i = 0; i < 511; i++)
        CHECK(lab_pricing_record_step(&state, 0, 1e-8) == LAB_PRICING_ZERO);
    CHECK(!state.bland_active);
    CHECK(lab_pricing_record_step(&state, 1, 1e-8) == 0);
    CHECK(state.consecutive_short == 0);
    for (size_t i = 0; i < 512; i++) {
        int events = lab_pricing_record_step(&state, 1e-10, 1e-8);
        CHECK(((events & LAB_PRICING_FALLBACK) != 0) == (i == 511));
    }
    CHECK(state.bland_active);
    lab_pricing_begin_adaptive(&state, 1);
    CHECK(state.bland_requested && state.bland_active && !state.history_count);
    lab_pricing_begin_adaptive(&state, 0);
    for (size_t i = 0; i < 130; i++) {
        basis[0] = i;
        CHECK(lab_pricing_watch_basis(&state, 3, basis) == 0);
    }
    CHECK(state.history_count == 64 && !state.bland_active);
    CHECK(lab_pricing_watch_basis(&state, 3, basis) == 1);
    CHECK(lab_pricing_watch_basis(NULL, 0, NULL) == -1);
    CHECK(lab_pricing_watch_basis(&state, 1, NULL) == -1);
}
int main(void) {
    adaptive_guards();
    ranking();
    accepted_steps();
    final_validated_event();
    puts("Primal pricing helper contracts passed.");
    return 0;
}

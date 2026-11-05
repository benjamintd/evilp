#include "pricing.h"
#include <math.h>
#include <stdint.h>
#include <string.h>

static void increment(size_t *value) {
    if (*value != SIZE_MAX)
        (*value)++;
}
void lab_pricing_begin(LabPricing *state, int bland_requested) {
    if (!state)
        return;
    memset(state, 0, sizeof *state);
    state->bland_requested = state->bland_active = bland_requested != 0;
}
void lab_pricing_scan_begin(LabPricingChoice *choice) {
    if (!choice)
        return;
    choice->column = SIZE_MAX;
    choice->reduced = 0;
}
int lab_pricing_consider(const LabPricing *state, LabPricingChoice *choice,
                          size_t column, long double reduced, int eligible) {
    if (!state || !choice || column == SIZE_MAX || !isfinite(reduced) ||
        (choice->column != SIZE_MAX && !isfinite(choice->reduced)))
        return -1;
    if (!eligible)
        return 0;
    int selected = choice->column == SIZE_MAX;
    if (!selected && state->bland_active)
        selected = column < choice->column;
    if (!selected && !state->bland_active)
        selected = reduced < choice->reduced ||
                   (reduced == choice->reduced && column < choice->column);
    if (!selected)
        return 0;
    choice->column = column;
    choice->reduced = reduced;
    return 1;
}
int lab_pricing_record_step(LabPricing *state, double accepted_step, double short_threshold) {
    if (!state || !isfinite(accepted_step) || accepted_step < 0 ||
        !isfinite(short_threshold) || short_threshold < 0)
        return -1;
    increment(&state->accepted_pivots);
    int events = 0, trigger = 0;
    if (accepted_step == 0) {
        increment(&state->zero_steps);
        state->consecutive_short = 0;
        trigger = 1;
        events |= LAB_PRICING_ZERO;
    } else if (accepted_step <= short_threshold) {
        increment(&state->short_steps);
        if (state->consecutive_short < LAB_PRICING_SHORT_STREAK)
            state->consecutive_short++;
        trigger = state->consecutive_short == LAB_PRICING_SHORT_STREAK;
        events |= LAB_PRICING_SHORT;
    } else {
        state->consecutive_short = 0;
    }
    if (trigger && !state->bland_active) {
        state->bland_active = 1;
        increment(&state->fallback_transitions);
        events |= LAB_PRICING_FALLBACK;
    }
    return events;
}
int lab_pricing_queue(LabPricingPending *pending, size_t position, size_t entering_column) {
    if (!pending || pending->valid || position == SIZE_MAX || entering_column == SIZE_MAX)
        return -1;
    pending->valid = 1;
    pending->position = position;
    pending->entering_column = entering_column;
    return 0;
}
void lab_pricing_discard(LabPricingPending *pending) {
    if (pending)
        memset(pending, 0, sizeof *pending);
}
int lab_pricing_finish(LabPricing *state, LabPricingPending *pending, int basis_validated,
                       size_t basis_size, const size_t *basis, const double *xb,
                       double short_threshold, unsigned *events) {
    if (events)
        *events = 0;
    if (!state || !pending || !events)
        return -1;
    if (!pending->valid || !basis_validated)
        return 0;
    if (!basis || !xb || pending->position >= basis_size ||
        basis[pending->position] != pending->entering_column)
        return -1;
    int recorded = lab_pricing_record_step(state, xb[pending->position], short_threshold);
    if (recorded < 0)
        return -1;
    *events = (unsigned)recorded;
    lab_pricing_discard(pending);
    return 1;
}

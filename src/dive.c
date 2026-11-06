#include "dive.h"
#include "common.h"
#include "exact.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    size_t column;
    double lower, upper, sibling;
    int has_sibling, sibling_tried;
} Decision;
typedef struct {
    const LabModel *original;
    LabModel view;
    const LabLpOptions *lp_options;
    const LabDiveOptions *options;
    LabDiveResult *result;
    LabLpSession *session;
    double started, deadline;
    double *domain_lower, *domain_upper, *lower, *upper, *point, *rounded;
    Decision *decision;
    size_t depth;
} Dive;

LabDiveOptions lab_dive_options(double seconds) {
    LabDiveOptions o = {.time_limit = seconds, .max_steps = 128, .max_backtracks = 8,
        .max_rows = 20000, .max_columns = 20000, .max_nonzeros = 1000000};
    return o;
}
const char *lab_dive_status_name(LabDiveStatus status) {
    switch (status) {
    case LAB_DIVE_CANDIDATE: return "candidate";
    case LAB_DIVE_NOT_FOUND: return "not_found";
    case LAB_DIVE_TIME_LIMIT: return "time_limit";
    case LAB_DIVE_STEP_LIMIT: return "step_limit";
    case LAB_DIVE_LP_ABSTAIN: return "lp_abstain";
    case LAB_DIVE_INVALID: return "invalid_input";
    case LAB_DIVE_RESOURCE: return "resource_limit";
    }
    return "invalid_status";
}
static void status(Dive *d, LabDiveStatus s, const char *message) {
    d->result->status = s;
    snprintf(d->result->message, sizeof d->result->message, "%s", message);
}
static int budget(Dive *d) {
    if (lab_now() >= d->deadline) {
        status(d, LAB_DIVE_TIME_LIMIT, "Primal diving time budget exhausted");
        return -1;
    }
    return 0;
}
static int set_bounds(Dive *d, size_t j, double lo, double hi) {
    if (isnan(lo) || isnan(hi) || lo > hi || lo < d->domain_lower[j] ||
        hi > d->domain_upper[j] ||
        (d->session && lab_lp_session_set_bounds(d->session, j, lo, hi))) {
        status(d, LAB_DIVE_INVALID, "A dive bound leaves its immutable root domain");
        return -1;
    }
    return 0;
}
/* Return 0 optimal, 1 infeasible on this tentative path, -1 abstention.
 * No heuristic LP status is propagated into the caller's search status. */
static int solve(Dive *d) {
    if (budget(d))
        return -1;
    LabLpResult lp = {0};
    if (d->session)
        lab_lp_session_solve(d->session, &lp);
    else {
        LabLpOptions options = *d->lp_options;
        options.basis_in = options.basis_out = NULL;
        options.time_limit = fmin(options.time_limit, d->deadline - lab_now());
        LabModel view = d->view;
        view.col_lower = d->lower;
        view.col_upper = d->upper;
        lab_lp_solve(&view, &options, &lp);
    }
    d->result->lp_solves++;
    d->result->lp_iterations += lp.iterations;
    d->result->warm_starts += lp.warm_start != 0;
    d->result->last_lp_status = lp.status;
    int result = -1;
    if (!strcmp(lp.status, "infeasible"))
        result = 1;
    else if (!strcmp(lp.status, "optimal") && lp.has_primal) {
        memcpy(d->point, lp.primal, d->view.cols * sizeof(double));
        result = 0;
    } else
        status(d, LAB_DIVE_LP_ABSTAIN, "An unresolved heuristic LP was discarded");
    lab_lp_result_free(&lp);
    if (budget(d))
        return -1;
    return result;
}
/* A rounded vector needs no extra LP if its full original-model evaluation
 * already succeeds. Otherwise, sequential fixing can move continuous values
 * and other integer coordinates before another checked rounding attempt. */
static int rounded_candidate(Dive *d) {
    d->result->rounded_checks++;
    for (size_t j = 0; j < d->view.cols; ++j) {
        double x = d->point[j];
        if (!isfinite(x))
            return 0;
        if (d->original->integer[j] && fabs(x) > 9007199254740991.0) {
            status(d, LAB_DIVE_LP_ABSTAIN, "Integer dive value exceeds the exactly spaced range");
            return -1;
        }
        d->rounded[j] = d->original->integer[j] ? round(x) : x;
        if (d->rounded[j] < d->lower[j] || d->rounded[j] > d->upper[j])
            return 0;
    }
    LabModel view = *d->original;
    view.col_lower = d->domain_lower;
    view.col_upper = d->domain_upper;
    double objective, violation;
    if (lab_exact_evaluate(&view, d->rounded, &objective, &violation) || violation > 1e-7)
        return 0;
    if (budget(d))
        return -1;
    d->result->primal = malloc(d->view.cols * sizeof(double));
    if (d->view.cols && !d->result->primal) {
        status(d, LAB_DIVE_RESOURCE, "Cannot allocate a heuristic candidate");
        return -1;
    }
    memcpy(d->result->primal, d->rounded, d->view.cols * sizeof(double));
    d->result->objective = objective;
    status(d, LAB_DIVE_CANDIDATE, "Original-model checks accepted a diving candidate");
    return 1;
}
static size_t choose(Dive *d, double *target, double *sibling, int *has_sibling) {
    size_t chosen = SIZE_MAX;
    double best_distance = INFINITY;
    long double best_change = INFINITY;
    for (size_t j = 0; j < d->view.cols; ++j) {
        if (!d->original->integer[j] || d->lower[j] == d->upper[j])
            continue;
        double x = d->point[j];
        if (!isfinite(x) || fabs(x) > 9007199254740991.0)
            continue;
        double down = floor(x), up = ceil(x), distance = fabs(x - round(x));
        if (distance <= 1e-7)
            continue;
        int down_ok = down >= d->lower[j] && down <= d->upper[j];
        int up_ok = up >= d->lower[j] && up <= d->upper[j];
        if (!down_ok && !up_ok)
            continue;
        long double cost = (d->original->maximize ? -1.L : 1.L) *
                           d->original->objective[j];
        double preferred = x - down < up - x ? down :
                           x - down > up - x ? up : cost < 0 ? up : down;
        if ((preferred == down && !down_ok) || (preferred == up && !up_ok))
            preferred = preferred == down ? up : down;
        double travel = fabs(preferred - x);
        long double change = cost * ((long double)preferred - x);
        if (chosen == SIZE_MAX || travel < best_distance ||
            (travel == best_distance && (change < best_change ||
             (change == best_change && j < chosen)))) {
            chosen = j;
            *target = preferred;
            *sibling = preferred == down ? up : down;
            *has_sibling = preferred == down ? up_ok : down_ok;
            best_distance = travel;
            best_change = change;
        }
    }
    return chosen;
}
static int backtrack(Dive *d) {
    while (d->depth) {
        if (budget(d))
            return -1;
        Decision *f = &d->decision[d->depth - 1];
        size_t j = f->column;
        if (set_bounds(d, j, f->lower, f->upper))
            return -1;
        d->lower[j] = f->lower;
        d->upper[j] = f->upper;
        if (f->has_sibling && !f->sibling_tried &&
            d->result->backtracks < d->options->max_backtracks &&
            d->result->steps < d->options->max_steps) {
            f->sibling_tried = 1;
            if (set_bounds(d, j, f->sibling, f->sibling))
                return -1;
            d->lower[j] = d->upper[j] = f->sibling;
            d->result->backtracks++;
            d->result->steps++;
            return 0;
        }
        d->depth--;
    }
    status(d, LAB_DIVE_NOT_FOUND, "Bounded diving paths did not produce a candidate");
    return 1;
}
/* If only tiny fractional distances remain, solving with all integers fixed
 * can move coupled continuous coordinates. Keep the preceding dive domain so
 * an infeasible rounding can still backtrack, without changing the root box. */
static int final_resolve(Dive *d) {
    size_t n = d->view.cols;
    double *lower = malloc(n * sizeof(double)), *upper = malloc(n * sizeof(double));
    if (n && (!lower || !upper)) {
        free(lower);
        free(upper);
        status(d, LAB_DIVE_RESOURCE, "Cannot preserve the final rounding domain");
        return -1;
    }
    memcpy(lower, d->lower, n * sizeof(double));
    memcpy(upper, d->upper, n * sizeof(double));
    int outcome = -1;
    for (size_t j = 0; j < n; ++j)
        if (d->original->integer[j]) {
            double value = round(d->point[j]);
            if (value < lower[j] || value > upper[j] || set_bounds(d, j, value, value))
                goto done;
            d->lower[j] = d->upper[j] = value;
        }
    d->result->final_resolves++;
    int solved = solve(d);
    if (!solved) {
        int candidate = rounded_candidate(d);
        if (candidate > 0)
            outcome = 0;
        else if (!candidate)
            status(d, LAB_DIVE_LP_ABSTAIN, "Fixed-integer heuristic witness failed source checks");
    } else if (solved > 0)
        outcome = 1;
done:
    if (outcome == 1)
        for (size_t j = 0; j < n; ++j)
            if (set_bounds(d, j, lower[j], upper[j])) {
                outcome = -1;
                break;
            }
    memcpy(d->lower, lower, n * sizeof(double));
    memcpy(d->upper, upper, n * sizeof(double));
    free(lower);
    free(upper);
    return outcome;
}
void lab_primal_dive(const LabModel *original, const LabModel *lp_model,
                     const double *root_lower, const double *root_upper,
                     const double *initial, const LabLpOptions *lp_options,
                     const LabDiveOptions *options, LabDiveResult *out) {
    if (!out)
        return;
    memset(out, 0, sizeof *out);
    Dive d = {.original = original, .lp_options = lp_options, .options = options,
              .result = out, .started = lab_now()};
    status(&d, LAB_DIVE_INVALID, "Invalid primal diving input");
    if (!original || !lp_model || !lp_options || !options ||
        !isfinite(options->time_limit) || options->time_limit <= 0 ||
        !options->max_steps || options->max_steps > 4096 || options->max_backtracks > 256 ||
        original->cols != lp_model->cols || !isfinite(lp_options->time_limit) ||
        lp_options->time_limit <= 0 || !lp_options->iteration_limit ||
        (original->cols && (!root_lower || !root_upper || !initial)))
        goto cleanup;
    d.deadline = d.started + options->time_limit;
    size_t n = original->cols;
    if (lp_model->rows > options->max_rows || n > options->max_columns ||
        lp_model->nnz > options->max_nonzeros || n > SIZE_MAX / sizeof(double)) {
        status(&d, LAB_DIVE_RESOURCE, "Primal diving model-size cap reached");
        goto cleanup;
    }
    d.view = *lp_model;
    d.domain_lower = malloc(n * sizeof(double));
    d.domain_upper = malloc(n * sizeof(double));
    d.lower = malloc(n * sizeof(double));
    d.upper = malloc(n * sizeof(double));
    d.point = malloc(n * sizeof(double));
    d.rounded = malloc(n * sizeof(double));
    d.decision = calloc(options->max_steps, sizeof *d.decision);
    if ((n && (!d.domain_lower || !d.domain_upper || !d.lower || !d.upper ||
               !d.point || !d.rounded)) || !d.decision) {
        status(&d, LAB_DIVE_RESOURCE, "Cannot allocate separate primal diving state");
        goto cleanup;
    }
    int integers = 0;
    for (size_t j = 0; j < n; ++j) {
        double lo = root_lower[j], hi = root_upper[j], x = initial[j];
        if (isnan(lo) || isnan(hi) || lo > hi || lo == INFINITY || hi == -INFINITY ||
            lo < lp_model->col_lower[j] || hi > lp_model->col_upper[j] ||
            !isfinite(x) || x < lo - 1e-7 || x > hi + 1e-7 ||
            original->integer[j] != lp_model->integer[j] ||
            (original->integer[j] && fabs(x) > 9007199254740991.0))
            goto cleanup;
        integers |= original->integer[j] != 0;
        d.domain_lower[j] = d.lower[j] = lo;
        d.domain_upper[j] = d.upper[j] = hi;
        d.point[j] = x;
    }
    if (!integers) {
        status(&d, LAB_DIVE_NOT_FOUND, "No integer coordinates to dive on");
        goto cleanup;
    }
    d.view.col_lower = d.domain_lower;
    d.view.col_upper = d.domain_upper;
    if (budget(&d))
        goto cleanup;
    if (!options->cold) {
        LabLpOptions lp = *lp_options;
        lp.basis_in = lp.basis_out = NULL;
        lp.time_limit = fmin(lp.time_limit, d.deadline - lab_now());
        d.session = lab_lp_session_create(&d.view, &lp);
        if (!d.session) {
            status(&d, LAB_DIVE_LP_ABSTAIN, "Cannot create an independent diving session");
            goto cleanup;
        }
    }
    for (;;) {
        if (budget(&d))
            break;
        int candidate = rounded_candidate(&d);
        if (candidate)
            break;
        if (out->steps >= options->max_steps) {
            status(&d, LAB_DIVE_STEP_LIMIT, "Primal diving step cap reached");
            break;
        }
        double target = 0, sibling = 0;
        int has_sibling = 0;
        size_t j = choose(&d, &target, &sibling, &has_sibling);
        if (j == SIZE_MAX) {
            int final = final_resolve(&d);
            if (final <= 0 || backtrack(&d))
                break;
            int solved;
            while ((solved = solve(&d)) > 0)
                if (backtrack(&d))
                    goto cleanup;
            if (solved < 0)
                break;
            continue;
        }
        Decision *f = &d.decision[d.depth++];
        *f = (Decision){.column = j, .lower = d.lower[j], .upper = d.upper[j],
            .sibling = sibling, .has_sibling = has_sibling};
        if (set_bounds(&d, j, target, target))
            break;
        d.lower[j] = d.upper[j] = target;
        out->steps++;
        int solved;
        while ((solved = solve(&d)) > 0)
            if (backtrack(&d))
                goto cleanup;
        if (solved < 0)
            break;
    }
cleanup:
    lab_lp_session_free(d.session);
    free(d.domain_lower);
    free(d.domain_upper);
    free(d.lower);
    free(d.upper);
    free(d.point);
    free(d.rounded);
    free(d.decision);
    out->elapsed = lab_now() - d.started;
}
void lab_dive_result_free(LabDiveResult *result) {
    if (!result)
        return;
    free(result->primal);
    memset(result, 0, sizeof *result);
}

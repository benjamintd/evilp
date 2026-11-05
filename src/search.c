#include "search.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int lab_search_init(LabSearch *search, size_t cols) {
    if (!search || search->cols || search->cost || cols > SIZE_MAX / sizeof(LabSearchCost))
        return -1;
    LabSearchCost *cost = cols ? calloc(cols, sizeof *cost) : NULL;
    if (cols && !cost)
        return -1;
    search->cols = cols;
    search->cost = cost;
    return 0;
}
void lab_search_free(LabSearch *search) {
    if (!search)
        return;
    free(search->cost);
    memset(search, 0, sizeof *search);
}
static int valid_mean(long double mean) {
    return isfinite(mean) && mean >= 0;
}
int lab_search_observe(LabSearch *search, size_t variable, LabSearchDirection direction,
                       double distance, long double degradation) {
    if (!search || !search->cost || variable >= search->cols ||
        (direction != LAB_SEARCH_DOWN && direction != LAB_SEARCH_UP))
        return -1;
    if (!isfinite(distance) || distance <= 0 || distance > 1 ||
        !isfinite(degradation) || degradation < 0)
        return 0;
    long double sample = degradation / distance;
    if (!valid_mean(sample))
        return 0;
    LabSearchCost *cost = &search->cost[variable];
    long double *mean = direction == LAB_SEARCH_DOWN ? &cost->down_mean : &cost->up_mean;
    size_t *count = direction == LAB_SEARCH_DOWN ? &cost->down_count : &cost->up_count;
    if (*count && !valid_mean(*mean))
        return 0;
    if (!*count) {
        *mean = sample;
        *count = 1;
    } else {
        size_t next = *count == SIZE_MAX ? SIZE_MAX : *count + 1;
        /* Convex online mean: no growing sum/product and no counter wrap. */
        long double updated = *mean + (sample - *mean) / (long double)next;
        if (!valid_mean(updated))
            return 0;
        *mean = updated;
        *count = next;
    }
    return 1;
}
/* Validate eligibility independently of every learned statistic. The caller
 * still owns candidate checking and construction of the actual child boxes. */
static LabSearchChoice eligible(const LabModel *model, const double *primal,
                                const double *lower, const double *upper,
                                size_t j, double tolerance, double *down, double *up) {
    double x = primal[j], lo = lower[j], hi = upper[j];
    if (!isfinite(x) || isnan(lo) || isnan(hi) || lo == INFINITY || hi == -INFINITY || lo > hi)
        return LAB_SEARCH_BAD;
    /* Match the engine's existing original-unit feasibility allowance. This
     * only rejects incoherent inputs; it never establishes their feasibility. */
    if ((isfinite(lo) && (long double)lo - x > 1e-7L) ||
        (isfinite(hi) && (long double)x - hi > 1e-7L))
        return LAB_SEARCH_BAD;
    if (!model->integer[j])
        return LAB_SEARCH_INTEGRAL;
    if (fabs(x) > 9007199254740991.0)
        return LAB_SEARCH_UNSUPPORTED;
    double left = floor(x), right = ceil(x);
    *down = x - left;
    *up = right - x;
    if (fmin(*down, *up) <= tolerance)
        return LAB_SEARCH_INTEGRAL;
    int has_left = left >= lo && left < hi;
    int has_right = right <= hi && right > lo;
    return has_left || has_right ? LAB_SEARCH_CHOSEN : LAB_SEARCH_BAD;
}
static int reliable(const LabSearchCost *cost) {
    return cost->down_count >= LAB_SEARCH_RELIABLE && cost->up_count >= LAB_SEARCH_RELIABLE &&
           valid_mean(cost->down_mean) && valid_mean(cost->up_mean);
}
static long double prediction(long double mean, double distance) {
    /* A log compares the same positive products without overflowing or
     * underflowing. A genuine zero remains the least improvement. */
    return mean == 0 ? -INFINITY : logl(mean) + logl(distance);
}
LabSearchChoice lab_search_choose(const LabModel *model, const double *primal,
                                  const double *lower, const double *upper,
                                  const LabSearch *search, double tolerance,
                                  size_t *variable) {
    if (!variable)
        return LAB_SEARCH_BAD;
    *variable = SIZE_MAX;
    if (!model || !isfinite(tolerance) || tolerance < 0 || tolerance >= .5 ||
        (model->cols && (!model->integer || !primal || !lower || !upper)))
        return LAB_SEARCH_BAD;
    double most_fractional = -1;
    for (size_t j = 0; j < model->cols; j++) {
        double down = 0, up = 0;
        LabSearchChoice status = eligible(model, primal, lower, upper, j, tolerance, &down, &up);
        if (status < 0) {
            *variable = SIZE_MAX;
            return status;
        }
        if (status == LAB_SEARCH_CHOSEN && fmin(down, up) > most_fractional) {
            most_fractional = fmin(down, up);
            *variable = j;
        }
    }
    if (*variable == SIZE_MAX)
        return LAB_SEARCH_INTEGRAL;
    if (!search || search->cols != model->cols || !search->cost)
        return LAB_SEARCH_CHOSEN;
    long double global_down = 0, global_up = 0;
    size_t count = 0;
    for (size_t j = 0; j < model->cols; j++)
        if (model->integer[j] && reliable(&search->cost[j])) {
            count++;
            global_down += (search->cost[j].down_mean - global_down) / (long double)count;
            global_up += (search->cost[j].up_mean - global_up) / (long double)count;
        }
    if (!count || !valid_mean(global_down) || !valid_mean(global_up) ||
        (global_down == 0 && global_up == 0))
        return LAB_SEARCH_CHOSEN;
    size_t best = SIZE_MAX;
    long double best_small = -INFINITY, best_large = -INFINITY;
    double best_fraction = -1;
    for (size_t j = 0; j < model->cols; j++) {
        double down = 0, up = 0;
        if (eligible(model, primal, lower, upper, j, tolerance, &down, &up) != LAB_SEARCH_CHOSEN)
            continue;
        const LabSearchCost *cost = &search->cost[j];
        long double down_mean = cost->down_count >= LAB_SEARCH_RELIABLE && valid_mean(cost->down_mean)
                                    ? cost->down_mean : global_down;
        long double up_mean = cost->up_count >= LAB_SEARCH_RELIABLE && valid_mean(cost->up_mean)
                                  ? cost->up_mean : global_up;
        long double d = prediction(down_mean, down), u = prediction(up_mean, up);
        long double small = fminl(d, u), large = fmaxl(d, u);
        double fraction = fmin(down, up);
        if (best == SIZE_MAX || small > best_small ||
            (small == best_small && (large > best_large ||
             (large == best_large && fraction > best_fraction)))) {
            best = j;
            best_small = small;
            best_large = large;
            best_fraction = fraction;
        }
    }
    if (best != SIZE_MAX)
        *variable = best;
    return LAB_SEARCH_CHOSEN;
}

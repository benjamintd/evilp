#include "search.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static void require(int condition, const char *message) {
    if (!condition) {
        fprintf(stderr, "search contract: %s\n", message);
        exit(1);
    }
}
static void train(LabSearch *search, size_t j, long double down, long double up) {
    for (int k = 0; k < LAB_SEARCH_RELIABLE; k++) {
        require(lab_search_observe(search, j, LAB_SEARCH_DOWN, .25, down*.25L) == 1,
                "down history recorded");
        require(lab_search_observe(search, j, LAB_SEARCH_UP, .25, up*.25L) == 1,
                "up history recorded");
    }
}
int main(void) {
    LabSearch search = {0};
    require(lab_search_init(&search, SIZE_MAX) == -1 && !search.cost && !search.cols,
            "allocation size overflow rejected without mutation");
    require(!lab_search_init(&search, 0), "empty statistics initialize");
    lab_search_free(&search);
    lab_search_free(&search);
    lab_search_free(NULL);
    require(!lab_search_init(&search, 5), "statistics initialize");
    require(lab_search_init(&search, 5) == -1, "owned allocation cannot be overwritten");
    unsigned char integers[] = {1, 0, 1, 1, 1};
    double x[] = {.2, .5, .49, .5, .5}, lo[] = {0, 0, 0, 0, 0}, hi[] = {1, 1, 1, 1, 1};
    LabModel model = {.cols = 5, .integer = integers};
    size_t choice;
    require(lab_search_choose(&model, x, lo, hi, &search, 1e-7, &choice) == LAB_SEARCH_CHOSEN && choice == 3,
            "most fractional fallback excludes continuous columns and breaks ties by index");
    x[3] = .49;
    x[4] = .49;
    require(lab_search_choose(&model, x, lo, hi, NULL, 1e-7, &choice) == LAB_SEARCH_CHOSEN && choice == 2,
            "fallback tie retains first original index");
    lab_search_free(&search);

    unsigned char two_integer[] = {1, 1};
    double p[] = {.5, .4}, lower[] = {0, 0}, upper[] = {1, 1};
    model.cols = 2;
    model.integer = two_integer;
    require(!lab_search_init(&search, 2), "two-column history initializes");
    train(&search, 0, 100, 1);
    train(&search, 1, 10, 10);
    require(lab_search_choose(&model, p, lower, upper, &search, 1e-7, &choice) == LAB_SEARCH_CHOSEN && choice == 1,
            "balanced expected improvement beats one spectacular weak-sided branch");
    model.offset = 1e300;
    model.maximize = 1;
    double objective[] = {1e300, -1e-300};
    model.objective = objective;
    require(lab_search_choose(&model, p, lower, upper, &search, 1e-7, &choice) == LAB_SEARCH_CHOSEN && choice == 1,
            "objective offsets, constants and sense cannot alter normalized history ranking");
    model.offset = -1e300;
    model.maximize = 0;
    require(lab_search_choose(&model, p, lower, upper, &search, 1e-7, &choice) == LAB_SEARCH_CHOSEN && choice == 1,
            "large negative offset also irrelevant");
    lab_search_free(&search);
    require(!lab_search_init(&search, 2), "partial-history fixture initializes");
    require(lab_search_observe(&search, 0, LAB_SEARCH_DOWN, .5, 1000) == 1,
            "single partial observation accepted");
    require(lab_search_choose(&model, p, lower, upper, &search, 1e-7, &choice) == LAB_SEARCH_CHOSEN && choice == 0,
            "insufficient history preserves most-fractional fallback");
    train(&search, 1, 10, 10);
    p[0] = .1;
    p[1] = .5;
    require(lab_search_choose(&model, p, lower, upper, &search, 1e-7, &choice) == LAB_SEARCH_CHOSEN && choice == 1,
            "global reliable means support an unseen side without starvation");
    lab_search_free(&search);

    require(!lab_search_init(&search, 2), "observation guard fixture initializes");
    require(lab_search_observe(&search, 0, LAB_SEARCH_DOWN, .25, 0) == 1 && search.cost[0].down_count == 1,
            "zero degradation is learned safely");
    require(lab_search_observe(&search, 0, LAB_SEARCH_DOWN, 0, 1) == 0 &&
            lab_search_observe(&search, 0, LAB_SEARCH_DOWN, -1, 1) == 0 &&
            lab_search_observe(&search, 0, LAB_SEARCH_DOWN, 2, 1) == 0 &&
            lab_search_observe(&search, 0, LAB_SEARCH_DOWN, NAN, 1) == 0 &&
            lab_search_observe(&search, 0, LAB_SEARCH_DOWN, .5, -1) == 0 &&
            lab_search_observe(&search, 0, LAB_SEARCH_DOWN, .5, NAN) == 0 &&
            lab_search_observe(&search, 0, LAB_SEARCH_DOWN, .5, INFINITY) == 0,
            "invalid observations ignored without manufacturing evidence");
    require(search.cost[0].down_count == 1 && search.cost[0].down_mean == 0,
            "ignored observations leave history unchanged");
    require(lab_search_observe(&search, 2, LAB_SEARCH_UP, .5, 1) == -1 &&
            lab_search_observe(&search, 0, (LabSearchDirection)0, .5, 1) == -1,
            "invalid indices/directions rejected");
    require(lab_search_observe(&search, 0, LAB_SEARCH_DOWN, DBL_TRUE_MIN, LDBL_MAX) == 0,
            "overflowing per-unit estimate ignored");
    for (size_t k = 0; k < 20000; k++)
        require(lab_search_observe(&search, 1, LAB_SEARCH_UP, .25, (long double)(k % 7)) == 1,
                "long observation stream retained");
    require(search.cost[1].up_count == 20000 && isfinite(search.cost[1].up_mean) &&
            search.cost[1].up_mean >= 0 && search.cost[1].up_mean <= 24,
            "online mean has no growing-sum overflow");
    search.cost[0].down_mean = 1;
    search.cost[0].down_count = SIZE_MAX;
    require(lab_search_observe(&search, 0, LAB_SEARCH_DOWN, 1, LDBL_MAX) == 1 &&
            search.cost[0].down_count == SIZE_MAX && isfinite(search.cost[0].down_mean),
            "saturated counts do not wrap or overflow means");
    for (size_t j = 0; j < 2; j++)
        search.cost[j] = (LabSearchCost){.down_mean = LDBL_TRUE_MIN, .up_mean = LDBL_TRUE_MIN,
                                       .down_count = 2, .up_count = 2};
    p[0] = .25;
    p[1] = .5;
    require(lab_search_choose(&model, p, lower, upper, &search, 0, &choice) == LAB_SEARCH_CHOSEN && choice == 1,
            "log scores preserve subnormal improvement ordering");
    for (size_t j = 0; j < 2; j++) {
        search.cost[j].down_mean = LDBL_MAX;
        search.cost[j].up_mean = LDBL_MAX;
    }
    require(lab_search_choose(&model, p, lower, upper, &search, 0, &choice) == LAB_SEARCH_CHOSEN && choice == 1,
            "very large learned costs remain finite comparisons");
    p[0] = p[1] = .5;
    require(lab_search_choose(&model, p, lower, upper, &search, 0, &choice) == LAB_SEARCH_CHOSEN && choice == 0,
            "learned score ties retain lowest original index");
    lab_search_free(&search);

    double tiny = ldexp(1, -30);
    p[0] = tiny;
    p[1] = 1-tiny;
    require(lab_search_choose(&model, p, lower, upper, NULL, 2*tiny, &choice) == LAB_SEARCH_INTEGRAL && choice == SIZE_MAX,
            "near integer tolerance returns no branch, not a feasibility claim");
    require(lab_search_choose(&model, p, lower, upper, NULL, 0, &choice) == LAB_SEARCH_CHOSEN && choice == 0,
            "zero-tolerance fallback detects actual fractions deterministically");
    p[0] = .5;
    p[1] = 0;
    lower[0] = .2;
    upper[0] = 1.8;
    require(lab_search_choose(&model, p, lower, upper, NULL, 0, &choice) == LAB_SEARCH_CHOSEN && choice == 0,
            "one feasible progressing child is enough to branch");
    lower[0] = 1.2;
    upper[0] = 1.8;
    p[0] = 1.5;
    require(lab_search_choose(&model, p, lower, upper, NULL, 0, &choice) == LAB_SEARCH_BAD,
            "no progressing child is explicit bad input, never an infeasibility proof");
    lower[0] = 0;
    upper[0] = INFINITY;
    p[0] = 9007199254740992.0;
    require(lab_search_choose(&model, p, lower, upper, NULL, 0, &choice) == LAB_SEARCH_UNSUPPORTED,
            "integer exact-spacing boundary is explicit");
    two_integer[0] = 0;
    p[0] = 0x1p60;
    p[1] = .5;
    require(lab_search_choose(&model, p, lower, upper, NULL, 0, &choice) == LAB_SEARCH_CHOSEN && choice == 1,
            "large continuous coordinate is excluded from integer spacing guard");
    p[0] = NAN;
    require(lab_search_choose(&model, p, lower, upper, NULL, 0, &choice) == LAB_SEARCH_BAD,
            "nonfinite points cannot silently look integral");
    p[0] = 0;
    lower[0] = NAN;
    require(lab_search_choose(&model, p, lower, upper, NULL, 0, &choice) == LAB_SEARCH_BAD,
            "nonfinite bounds rejected");
    lower[0] = 0;
    require(lab_search_choose(&model, p, lower, upper, NULL, NAN, &choice) == LAB_SEARCH_BAD &&
            lab_search_choose(&model, p, lower, upper, NULL, -.1, &choice) == LAB_SEARCH_BAD &&
            lab_search_choose(&model, p, lower, upper, NULL, .5, &choice) == LAB_SEARCH_BAD,
            "invalid integrality tolerances rejected");
    require(lab_search_choose(&model, p, lower, upper, NULL, 0, NULL) == LAB_SEARCH_BAD,
            "missing result pointer rejected");
    LabModel empty = {0};
    require(lab_search_choose(&empty, NULL, NULL, NULL, NULL, 0, &choice) == LAB_SEARCH_INTEGRAL,
            "empty model has no branch candidates");
    puts("Branch selection and offset-free pseudocost contracts passed.");
    return 0;
}

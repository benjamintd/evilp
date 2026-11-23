#include "lpi/lpi.h"
#include <math.h>
#include <stdio.h>

#define CHECK(e) do { if (!(e)) { fprintf(stderr, "line %d: %s\n", __LINE__, #e); return 1; } } while (0)
#define OK(e) CHECK((e) == SCIP_OKAY)

/* The original wide interval remains finite when the column box, objective,
 * row scale or exposed row sense changes. No optimization library is used
 * as the oracle: each one-variable optimum is obtained from its two bounds. */
static int verify(SCIP_LPI *lp, double lo, double hi, double x, double coefficient,
                  double expected_objective) {
    double left, right, objective, primal, activity;
    OK(SCIPlpiGetSides(lp, 0, 0, &left, &right));
    CHECK(left == lo && right == hi);
    OK(SCIPlpiSolveDual(lp));
    CHECK(SCIPlpiIsOptimal(lp));
    OK(SCIPlpiGetSol(lp, &objective, &primal, NULL, &activity, NULL));
    CHECK(primal == x && objective == expected_objective);
    CHECK(activity == coefficient * x);
    return 0;
}

int main(void) {
    SCIP_LPI *lp = NULL;
    int column = 0, row = 0;
    double a = 1, cost = 1, lower = -100, upper = 100, lhs = -1e20, rhs = 3.25;
    OK(SCIPlpiCreate(&lp, NULL, "wide finite row", SCIP_OBJSEN_MAXIMIZE));
    OK(SCIPlpiLoadColLP(lp, SCIP_OBJSEN_MAXIMIZE, 1, &cost, &lower, &upper, NULL,
                        1, &lhs, &rhs, NULL, 1, &column, &row, &a));
    CHECK(!verify(lp, lhs, rhs, 3.25, 1, 3.25));
    rhs = -0.125;
    OK(SCIPlpiChgSides(lp, 1, &row, &lhs, &rhs));
    CHECK(!verify(lp, lhs, rhs, -0.125, 1, -0.125));
    OK(SCIPlpiScaleRow(lp, row, -2));
    CHECK(!verify(lp, 0.25, 2e20, -0.125, -2, -0.125));

    lhs = -1;
    rhs = 4;
    OK(SCIPlpiChgSides(lp, 1, &row, &lhs, &rhs));
    CHECK(!verify(lp, lhs, rhs, 0.5, -2, 0.5));
    lhs = -2e20;
    rhs = 0;
    OK(SCIPlpiChgSides(lp, 1, &row, &lhs, &rhs));
    CHECK(!verify(lp, lhs, rhs, 100, -2, 100));

    /* Expanding the column domain makes the far side of the finite row bind. */
    upper = 2e20;
    OK(SCIPlpiChgBounds(lp, 1, &column, &lower, &upper));
    CHECK(!verify(lp, lhs, rhs, 1e20, -2, 1e20));
    cost = -1;
    OK(SCIPlpiChgObj(lp, 1, &column, &cost));
    CHECK(!verify(lp, lhs, rhs, 0, -2, 0));
    /* Distinct finite sides separated by less than the backend epsilon still
     * describe an interval in the caller's LP. */
    lhs = 1.0;
    rhs = nextafter(lhs, 2.0);
    OK(SCIPlpiChgSides(lp, 1, &row, &lhs, &rhs));
    double read_lo, read_hi;
    OK(SCIPlpiGetSides(lp, 0, 0, &read_lo, &read_hi));
    CHECK(read_lo == lhs && read_hi == rhs);
    OK(SCIPlpiFree(&lp));
    return 0;
}

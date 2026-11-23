#include "lpi/lpi.h"
#include <stdio.h>

#define CHECK(x) do { if (!(x)) { fprintf(stderr, "line %d: %s\n", __LINE__, #x); return 1; } } while (0)
#define OK(x) CHECK((x) == SCIP_OKAY)

enum { ROWS = 136, COLS = 7 };

/* x0 is fixed by an existing equality. Add structural columns after the first
 * logical column, then append two sparse batches separated by a solve. */
int main(void) {
    SCIP_LPI *lp = NULL;
    OK(SCIPlpiCreate(&lp, NULL, "batch rows", SCIP_OBJSEN_MINIMIZE));
    const double inf = SCIPlpiInfinity(lp);
    double obj[COLS] = {1}, lo[COLS] = {0}, hi[COLS] = {3, 3, 3, 3, 3, 3, 3};
    int start = 0, column = 0;
    double one = 1;
    OK(SCIPlpiAddCols(lp, 1, obj, lo, hi, NULL, 0, NULL, NULL, NULL));
    OK(SCIPlpiAddRows(lp, 1, &one, &one, NULL, 1, &start, &column, &one));
    OK(SCIPlpiAddCols(lp, COLS-1, obj+1, lo+1, hi+1, NULL, 0, NULL, NULL, NULL));

    double expected[ROWS][COLS] = {{0}}, lhs[ROWS], rhs[ROWS];
    for (int batch = 0; batch < 2; ++batch) {
        int beg[ROWS/2], idx[ROWS/2 * COLS], nz = 0;
        double val[ROWS/2 * COLS];
        int first = batch * (ROWS/2);
        for (int k = 0; k < ROWS/2; ++k) {
            int i = first + k;
            beg[k] = nz;
            /* Last column stays empty until after both batches. Vary the
             * first occupied slot, sign and order, including empty rows. */
            for (int j = COLS-2; j >= 0 && i % 9 != 0; --j) {
                if ((i + j) % 3 == 0) {
                    idx[nz] = j;
                    expected[i][j] = 0.25 * (j + 1) * (i % 2 ? -1 : 1);
                    val[nz++] = expected[i][j];
                }
            }
            double activity = expected[i][0];
            lhs[i] = i % 4 == 1 ? -inf : activity;
            rhs[i] = i % 4 == 2 ? inf : activity + (i % 4 == 3 ? 1 : 0);
        }
        OK(SCIPlpiAddRows(lp, ROWS/2, lhs+first, rhs+first, NULL, nz, beg, idx, val));
        OK(SCIPlpiSolveDual(lp));
        CHECK(SCIPlpiIsOptimal(lp));
        double objective;
        OK(SCIPlpiGetSol(lp, &objective, NULL, NULL, NULL, NULL));
        CHECK(objective == 1);
    }

    int beg[ROWS+1], idx[ROWS * COLS], count;
    double val[ROWS * COLS], actual_lo[ROWS], actual_hi[ROWS];
    OK(SCIPlpiGetRows(lp, 1, ROWS, actual_lo, actual_hi, &count, beg, idx, val));
    beg[ROWS] = count;
    for (int i = 0; i < ROWS; ++i) {
        CHECK(actual_lo[i] == lhs[i] && actual_hi[i] == rhs[i]);
        for (int p = beg[i]; p < beg[i+1]; ++p) {
            CHECK(idx[p] >= 0 && idx[p] < COLS && val[p] != 0);
            CHECK(expected[i][idx[p]] == val[p]);
            expected[i][idx[p]] = 0;
        }
        for (int j = 0; j < COLS; ++j)
            CHECK(expected[i][j] == 0);
    }

    /* The getters legitimately return NULL sparse arrays for an empty slice. */
    OK(SCIPlpiGetRows(lp, 1, 1, actual_lo, actual_hi, &count, beg, idx, val));
    CHECK(count == 0 && actual_lo[0] == 0 && actual_hi[0] == 0);
    OK(SCIPlpiGetCols(lp, COLS-1, COLS-1, actual_lo, actual_hi, &count, beg, idx, val));
    CHECK(count == 0 && actual_lo[0] == 0 && actual_hi[0] == 3);
    OK(SCIPlpiGetCols(lp, 0, 0, actual_lo, actual_hi, &count, beg, idx, val));
    CHECK(count > 0 && actual_lo[0] == 0 && actual_hi[0] == 3);
    OK(SCIPlpiScaleRow(lp, 0, 2));
    OK(SCIPlpiGetSides(lp, 0, 0, actual_lo, actual_hi));
    CHECK(actual_lo[0] == 2 && actual_hi[0] == 2);
    OK(SCIPlpiScaleCol(lp, 0, 2));
    double coefficient;
    OK(SCIPlpiGetCoef(lp, 0, 0, &coefficient));
    CHECK(coefficient == 4);

    /* Reuse the formerly empty column, then remove a contradictory row. */
    column = COLS-1;
    double four = 4;
    OK(SCIPlpiAddRows(lp, 1, &four, &four, NULL, 1, &start, &column, &one));
    OK(SCIPlpiSolveDual(lp));
    CHECK(SCIPlpiIsPrimalInfeasible(lp));
    OK(SCIPlpiDelRows(lp, ROWS+1, ROWS+1));
    OK(SCIPlpiSolveDual(lp));
    CHECK(SCIPlpiIsOptimal(lp));
    OK(SCIPlpiFree(&lp));
    return 0;
}

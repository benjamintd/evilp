#include "lpi/lpi.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "failed at line %d: %s\n", __LINE__, #x);                              \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)
#define OK(x) CHECK((x) == SCIP_OKAY)

enum { N = 223, WIDTH = 17 };

int main(void) {
    SCIP_LPI *lp = NULL;
    int beg[N], ind[N * WIDTH], rows[N], cstat[N], rstat[N];
    double value[N * WIDTH], obj[N], lo[N] = {0}, hi[N], rhs[N] = {0};
    double target[N], x[N];
    OK(SCIPlpiCreate(&lp, NULL, "refactor", SCIP_OBJSEN_MINIMIZE));

    /* Every row has diagonal 32 and 16 entries of magnitude 1, so Ax=b
       determines x uniquely. Use binary fractions for the independent oracle. */
    for (int j = 0; j < N; ++j) {
        beg[j] = j * WIDTH;
        rows[j] = j;
        cstat[j] = SCIP_BASESTAT_LOWER;
        rstat[j] = SCIP_BASESTAT_BASIC;
        obj[j] = 1 + j % 7;
        hi[j] = 4;
        for (int k = 0; k < WIDTH; ++k) {
            ind[beg[j] + k] = (j + 19 * k) % N;
            value[beg[j] + k] = k == 0 ? 32 : (j + k) % 2 ? 1 : -1;
        }
    }
    OK(SCIPlpiLoadColLP(lp, SCIP_OBJSEN_MINIMIZE, N, obj, lo, hi, NULL, N, rhs, rhs, NULL,
                        N * WIDTH, beg, ind, value));

    for (int round = 0; round < 5; ++round) {
        memset(rhs, 0, sizeof rhs);
        double expected = 0;
        for (int j = 0; j < N; ++j) {
            target[j] = 0.25 * (1 + (j + round) % 13);
            if (round == 2 && j == N / 2)
                target[j] = 5;
            expected += obj[j] * target[j];
            for (int k = beg[j]; k < beg[j] + WIDTH; ++k)
                rhs[ind[k]] += value[k] * target[j];
        }
        OK(SCIPlpiChgSides(lp, N, rows, rhs, rhs));
        OK(SCIPlpiSetBase(lp, cstat, rstat));
        OK(SCIPlpiSolveDual(lp));
        if (round == 2) {
            CHECK(SCIPlpiIsPrimalInfeasible(lp));
            continue;
        }
        CHECK(SCIPlpiIsOptimal(lp));
        double observed;
        OK(SCIPlpiGetSol(lp, &observed, x, NULL, NULL, NULL));
        CHECK(fabs(observed - expected) < 1e-7);
        for (int j = 0; j < N; ++j)
            CHECK(fabs(x[j] - target[j]) < 1e-8);
    }
    /* The zero-nonzero API permits NULL sparse arrays. Leave a nonzero
       count in the first scratch slot, then reuse that slot for an empty row. */
    double zero = 0, one = 1, eight = 8;
    int start = 0, indices[] = {0, 1};
    double vals[] = {1, 1};
    OK(SCIPlpiAddRows(lp, 1, &zero, &eight, NULL, 2, &start, indices, vals));
    OK(SCIPlpiAddRows(lp, 1, &zero, &zero, NULL, 0, NULL, NULL, NULL));
    OK(SCIPlpiAddCols(lp, 1, &zero, &zero, &zero, NULL, 0, NULL, NULL, NULL));
    OK(SCIPlpiSolveDual(lp));
    CHECK(SCIPlpiIsOptimal(lp));
    OK(SCIPlpiAddRows(lp, 1, &one, &one, NULL, 0, NULL, NULL, NULL));
    OK(SCIPlpiSolveDual(lp));
    CHECK(SCIPlpiIsPrimalInfeasible(lp));
    OK(SCIPlpiFree(&lp));
    return 0;
}

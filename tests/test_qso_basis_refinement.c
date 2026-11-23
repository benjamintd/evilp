#include "lpi/lpi.h"
#include <stdio.h>
#include <string.h>

/* With both fixed variables equal to one, the equality
   small*x + 2^35*z + y = 2^35 forces y = -small exactly. */
static int check_cancellation(double small, double y_lower)
{
    const double budget = 0x1p35;
    SCIP_LPI *lp = NULL;
    double objective[] = {0, 0, 1}, lower[] = {1, 1, y_lower}, upper[] = {1, 1, 1};
    double a[] = {small, budget, 1}, solution[3], value;
    int offsets[] = {0, 1, 2}, indices[] = {0, 0, 0};
    int columns[] = {SCIP_BASESTAT_LOWER, SCIP_BASESTAT_LOWER, SCIP_BASESTAT_BASIC};
    int rows[] = {SCIP_BASESTAT_LOWER};
    int ok = 0;
    if (SCIPlpiCreate(&lp, NULL, "cancellation", SCIP_OBJSEN_MINIMIZE) != SCIP_OKAY)
        return 0;
    if (SCIPlpiLoadColLP(lp, SCIP_OBJSEN_MINIMIZE, 3, objective, lower, upper, NULL,
                         1, &budget, &budget, NULL, 3, offsets, indices, a) != SCIP_OKAY ||
        SCIPlpiSetBase(lp, columns, rows) != SCIP_OKAY ||
        SCIPlpiSolveDual(lp) != SCIP_OKAY)
        goto cleanup;
    if (-small < y_lower) {
        ok = !SCIPlpiIsOptimal(lp) && SCIPlpiIsPrimalInfeasible(lp);
        if (!ok)
            fprintf(stderr, "Expected infeasibility, got optimal=%d, primal-infeasible=%d\n",
                    SCIPlpiIsOptimal(lp), SCIPlpiIsPrimalInfeasible(lp));
        goto cleanup;
    }
    if (!SCIPlpiIsOptimal(lp) || SCIPlpiGetSol(lp, &value, solution, NULL, NULL, NULL) != SCIP_OKAY)
        goto cleanup;
    ok = solution[0] == 1 && solution[1] == 1 && solution[2] == -small && value == -small;
    if (!ok)
        fprintf(stderr, "Expected y and objective %a, got %a and %a\n", -small, solution[2], value);
cleanup:
    return SCIPlpiFree(&lp) == SCIP_OKAY && ok;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--nonnegative") == 0)
        return !check_cancellation(0x1.8p-20, 0);
    if (argc != 1) return 2;
    return !(check_cancellation(0x1.8p-20, -1) && check_cancellation(-0x1.8p-20, -1)
             && check_cancellation(0, -1));
}

#include "scip/scip.h"
#include <stdlib.h>
#include <string.h>

/* Linked only into the contract executable. SCIP solves the actual model. */
SCIP_Real test_reported_dual(SCIP *scip) {
    const char *bound = getenv("TEST_SCIP_BOUND");
    if (bound && strcmp(bound, "missing") == 0)
        return SCIPinfinity(scip);
    return SCIPgetDualbound(scip) - 1.0 / 1024.0;
}

SCIP_STATUS test_reported_status(SCIP *scip) {
    const char *status = getenv("TEST_SCIP_STATUS");
    if (status && strcmp(status, "time_limit") == 0)
        return SCIP_STATUS_TIMELIMIT;
    return SCIPgetStatus(scip);
}

#define _POSIX_C_SOURCE 200809L
#ifdef __APPLE__
#define _DARWIN_C_SOURCE
#endif
#include "common.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <stdio.h>
#include <time.h>

/* Supplied by the locked QSopt_ex C library. */
extern double ILLutil_deadline_zeit(void);

#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "failed at line %d: %s\n", __LINE__, #x);                              \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)
#define OK(x) CHECK((x) == SCIP_OKAY)

static double platform_now(void) {
    struct timespec ts;
#ifdef __APPLE__
    clockid_t id = CLOCK_UPTIME_RAW;
#else
    clockid_t id = CLOCK_MONOTONIC;
#endif
    if (clock_gettime(id, &ts))
        return -1;
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

int main(void) {
    SCIP *parent = NULL, *child = NULL;
    SCIP_CLOCK *clock = NULL;
    const char *param = "heuristics/subsciptimefrac";
    OK(SCIPcreate(&parent));
    OK(SCIPcreate(&child));
    OK(SCIPincludeDefaultPlugins(parent));
    OK(SCIPincludeDefaultPlugins(child));
    OK(SCIPcreateProbBasic(parent, "time"));
    OK(SCIPcreateProbBasic(child, "heuristic"));
    OK(SCIPcreateWallClock(parent, &clock));
    double before = platform_now();
    double lab = lab_now();
    double qs = ILLutil_deadline_zeit();
    OK(SCIPstartClock(parent, clock));
    double after = platform_now();
    CHECK(before >= 0 && before <= lab && lab <= qs && qs <= after);

    /* A descheduled solver must consume its wall budget. Bracket both clock reads
       so the test also works when the process is preempted between them. */
    const struct timespec pause = {0, 40000000};
    CHECK(nanosleep(&pause, NULL) == 0);
    double end_before = platform_now();
    OK(SCIPstopClock(parent, clock));
    double elapsed = SCIPgetClockTime(parent, clock);
    double end_after = platform_now();
    CHECK(elapsed >= end_before - after - 0.000002 && elapsed <= end_after - before + 0.000002);
    CHECK(ILLutil_deadline_zeit() - qs >= 0.039);
    OK(SCIPfreeClock(parent, &clock));

    SCIP_Real limit;
    SCIP_Longint nodes;
    OK(SCIPsetRealParam(parent, "limits/time", 80.0));
    OK(SCIPcopyLimitsHeur(parent, child));
    OK(SCIPgetRealParam(child, "limits/time", &limit));
    CHECK(limit == 80.0);
    OK(SCIPaddRealParam(parent, param, "test fraction", NULL, FALSE, 1.0, 0.0, 1.0, NULL, NULL));
    OK(SCIPsetRealParam(parent, param, 0.25));
    OK(SCIPcopyLimitsHeur(parent, child));
    OK(SCIPgetRealParam(child, "limits/time", &limit));
    CHECK(limit == 20.0);

    /* A heuristic with a stricter limit and custom node budget retains both. */
    OK(SCIPsetRealParam(child, "limits/time", 2.0));
    OK(SCIPsetLongintParam(child, "limits/nodes", 7));
    OK(SCIPcapHeurTime(parent, child));
    OK(SCIPgetRealParam(child, "limits/time", &limit));
    OK(SCIPgetLongintParam(child, "limits/nodes", &nodes));
    CHECK(limit == 2.0 && nodes == 7);
    OK(SCIPsetRealParam(parent, "limits/time", SCIPinfinity(parent)));
    OK(SCIPcapHeurTime(parent, child));
    OK(SCIPgetRealParam(child, "limits/time", &limit));
    CHECK(limit == 2.0);
    OK(SCIPsetRealParam(parent, param, 0.0));
    OK(SCIPcapHeurTime(parent, child));
    OK(SCIPgetRealParam(child, "limits/time", &limit));
    CHECK(limit == 0.0);
    OK(SCIPfree(&child));
    OK(SCIPfree(&parent));
    return 0;
}

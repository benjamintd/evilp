#define _POSIX_C_SOURCE 200809L
#include "presolve.h"
#include "exact.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>

static void require(int condition, const char *message) {
    if (!condition) {
        fprintf(stderr, "presolve contract: %s\n", message);
        exit(1);
    }
}
static LabPresolveOptions options(void) {
    return (LabPresolveOptions){5, 32, 0};
}
static uint32_t rng = UINT32_C(9917223);
static long integer(long low, long high) {
    rng = UINT32_C(1664525) * rng + UINT32_C(1013904223);
    return low + (long)(rng % (uint32_t)(high - low + 1));
}
static void generated_integer_oracles(void) {
    size_t start[] = {0, 3, 6, 9}, row[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    unsigned char integrality[] = {1, 1, 1};
    double lo[] = {-2, -2, -2}, hi[] = {3, 3, 3};
    double a[9], rl[3], ru[3];
    LabModel m = {.rows = 3, .cols = 3, .nnz = 9, .col_start = start, .row_index = row,
                  .integer = integrality, .col_lower = lo, .col_upper = hi,
                  .row_lower = rl, .row_upper = ru, .value = a};
    size_t improvements = 0;
    for (size_t sample = 0; sample < 1000; sample++) {
        long exact[9], lower[3], upper[3];
        for (size_t k = 0; k < 9; k++)
            a[k] = (double)(exact[k] = integer(-3, 3));
        for (size_t i = 0; i < 3; i++) {
            lower[i] = integer(-9, 0);
            upper[i] = integer(0, 9);
            rl[i] = sample % 3 == 0 ? -INFINITY : (double)lower[i];
            ru[i] = sample % 3 == 1 ? INFINITY : (double)upper[i];
        }
        LabPresolveResult r;
        LabPresolveOptions opt = options();
        lab_presolve(&m, NULL, NULL, &opt, &r);
        require(strcmp(r.status, "input_error") != 0 && r.source_available,
                "valid binary models propagate");
        improvements += r.tightenings;
        size_t feasible = 0;
        for (long x = -2; x <= 3; x++)
            for (long y = -2; y <= 3; y++)
                for (long z = -2; z <= 3; z++) {
                    long point[] = {x, y, z}, activity[3] = {0};
                    int original = 1, retained = 1;
                    for (size_t i = 0; i < 3; i++) {
                        for (size_t j = 0; j < 3; j++)
                            activity[i] += exact[3 * j + i] * point[j];
                        if ((isfinite(rl[i]) && activity[i] < lower[i]) ||
                            (isfinite(ru[i]) && activity[i] > upper[i]))
                            original = 0;
                    }
                    for (size_t j = 0; j < 3; j++)
                        if ((double)point[j] < r.lower[j] || (double)point[j] > r.upper[j])
                            retained = 0;
                    if (original) {
                        feasible++;
                        require(retained, "no exact feasible integer point is cut off");
                    }
                    if (retained && strcmp(r.status, "infeasible"))
                        for (size_t i = 0; i < 3; i++)
                            if (r.redundant[i])
                                require((!isfinite(rl[i]) || activity[i] >= lower[i]) &&
                                        (!isfinite(ru[i]) || activity[i] <= upper[i]),
                                        "removed row holds across the returned box");
                }
        if (!strcmp(r.status, "infeasible"))
            require(!feasible, "claimed empty box has no independently enumerated point");
        require(lo[0] == -2 && lo[1] == -2 && lo[2] == -2 &&
                hi[0] == 3 && hi[1] == 3 && hi[2] == 3,
                "original column bounds stay immutable");
        for (size_t k = 0; k < 9; k++)
            require(a[k] == (double)exact[k], "original CSC stays immutable");
        lab_presolve_result_free(&r);
    }
    require(improvements > 100, "oracle campaign exercises real tightening");
}
static void one_sided_bounds_and_constants(void) {
    size_t start[] = {0, 2, 4}, row[] = {0, 1, 0, 1};
    double a[] = {1, 1, 1, -1}, lo[] = {0, 1}, hi[] = {INFINITY, 1};
    double rl[] = {-INFINITY, 1}, ru[] = {3, INFINITY};
    unsigned char integrality[] = {1, 0};
    LabModel m = {.rows = 2, .cols = 2, .nnz = 4, .col_start = start, .row_index = row,
                  .value = a, .col_lower = lo, .col_upper = hi, .row_lower = rl,
                  .row_upper = ru, .integer = integrality};
    LabPresolveOptions opt = options();
    LabPresolveResult r;
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(!strcmp(r.status, "unresolved") && r.lower[0] == 2 && r.upper[0] == 2 &&
                r.fixed == 2 && r.redundant_rows == 2,
            "one unbounded activity can be excluded to fix an integer column");
    lab_presolve_result_free(&r);
    double changed_lo[] = {3, 1};
    lab_presolve(&m, changed_lo, NULL, &opt, &r);
    require(!strcmp(r.status, "infeasible"), "tightening can expose an empty row");
    lab_presolve_result_free(&r);
    double cl[] = {1}, cu[] = {INFINITY};
    LabModel constant = {.rows = 1, .row_lower = cl, .row_upper = cu};
    lab_presolve(&constant, NULL, NULL, &opt, &r);
    require(!strcmp(r.status, "infeasible"), "contradictory exact constant row");
    lab_presolve_result_free(&r);
    constant.parsed_decimal = 1;
    lab_presolve(&constant, NULL, NULL, &opt, &r);
    require(!strcmp(r.status, "unresolved") && !r.source_available,
            "unknown decimal provenance cannot justify even constant-row pruning");
    lab_presolve_result_free(&r);
}
static void large_sparse_chain(void) {
    size_t n = 10000, rows = n - 1, nnz = 2 * rows;
    size_t *start = calloc(n + 1, sizeof(size_t)), *row = calloc(nnz, sizeof(size_t));
    double *a = calloc(nnz, sizeof(double)), *lo = calloc(n, sizeof(double));
    double *hi = calloc(n, sizeof(double)), *rl = calloc(rows, sizeof(double));
    double *ru = calloc(rows, sizeof(double));
    unsigned char *integrality = calloc(n, 1);
    require(start && row && a && lo && hi && rl && ru && integrality, "large test allocation");
    size_t at = 0;
    for (size_t j = 0; j < n; j++) {
        start[j] = at;
        hi[j] = (double)j;
        integrality[j] = 1;
        if (j) {
            row[at] = j - 1;
            a[at++] = 1;
        }
        if (j + 1 < n) {
            row[at] = j;
            a[at++] = -1;
        }
    }
    start[n] = at;
    for (size_t i = 0; i < rows; i++) {
        rl[i] = 1;
        ru[i] = INFINITY;
    }
    LabModel m = {.rows = rows, .cols = n, .nnz = nnz, .col_start = start,
                  .row_index = row, .value = a, .col_lower = lo, .col_upper = hi,
                  .row_lower = rl, .row_upper = ru, .integer = integrality};
    LabPresolveOptions opt = options();
    LabPresolveResult r;
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(!strcmp(r.status, "unresolved") && r.fixed == n && r.redundant_rows == rows,
            "large sparse chain collapses without a dense matrix");
    for (size_t j = 0; j < n; j++)
        require(r.lower[j] == (double)j && r.upper[j] == (double)j && lo[j] == 0,
                "large sparse bound propagation is exact and immutable");
    lab_presolve_result_free(&r);
    free(start); free(row); free(a); free(lo); free(hi); free(rl); free(ru); free(integrality);
}
static void read_text(LabModel *m, const char *text, char path[64]) {
    strcpy(path, "/tmp/milp-presolve-XXXXXX");
    int fd = mkstemp(path);
    require(fd >= 0, "temporary MPS");
    FILE *f = fdopen(fd, "w");
    require(f && fputs(text, f) >= 0 && !fclose(f), "write temporary MPS");
    if (lab_model_read(m, path)) {
        fprintf(stderr, "%s\n", m->error);
        require(0, "read source model");
    }
}
static void generated_decimal_oracles(void) {
    size_t changes = 0;
    for (size_t sample = 0; sample < 100; sample++) {
        long a[9], rhs[3];
        char text[16384], path[64];
        size_t at = (size_t)snprintf(text, sizeof text, "NAME FRACTIONAL\nROWS\n N O\n");
        for (size_t i = 0; i < 129; i++)
            at += (size_t)snprintf(text + at, sizeof text - at, " L R%zu\n", i);
        at += (size_t)snprintf(text + at, sizeof text - at, "COLUMNS\n");
        for (size_t j = 0; j < 3; j++) {
            at += (size_t)snprintf(text + at, sizeof text - at, " x%zu O 0\n", j);
            for (size_t i = 0; i < 3; i++) {
                a[3 * j + i] = integer(-9, 9);
                at += (size_t)snprintf(text + at, sizeof text - at,
                                       " x%zu R%zu %.1f\n", j, i, (double)a[3 * j + i] / 10);
            }
        }
        at += (size_t)snprintf(text + at, sizeof text - at, "RHS\n");
        for (size_t i = 0; i < 3; i++) {
            rhs[i] = integer(-12, 12);
            at += (size_t)snprintf(text + at, sizeof text - at, " b R%zu %.1f\n",
                                   i, (double)rhs[i] / 10);
        }
        at += (size_t)snprintf(text + at, sizeof text - at, "BOUNDS\n");
        for (size_t j = 0; j < 3; j++)
            at += (size_t)snprintf(text + at, sizeof text - at, " BV b x%zu\n", j);
        snprintf(text + at, sizeof text - at, "ENDATA\n");
        LabModel m = {0};
        read_text(&m, text, path);
        require(!m.exact && m.envelope && m.envelope->complete,
                "large-row decimal models exercise envelope-only source arithmetic");
        LabPresolveOptions opt = options();
        LabPresolveResult r;
        lab_presolve(&m, NULL, NULL, &opt, &r);
        require(r.source_available && strcmp(r.status, "input_error"),
                "decimal envelope model is supported");
        changes += r.tightenings;
        size_t feasible = 0;
        for (long x = 0; x <= 1; x++)
            for (long y = 0; y <= 1; y++)
                for (long z = 0; z <= 1; z++) {
                    long p[] = {x, y, z};
                    int valid = 1;
                    for (size_t i = 0; i < 3; i++) {
                        long activity = 0;
                        for (size_t j = 0; j < 3; j++)
                            activity += a[3 * j + i] * p[j];
                        if (activity > rhs[i])
                            valid = 0;
                    }
                    if (valid) {
                        feasible++;
                        require(strcmp(r.status, "infeasible"),
                                "source-decimal oracle forbids a false infeasibility claim");
                        for (size_t j = 0; j < 3; j++)
                            require((double)p[j] >= r.lower[j] && (double)p[j] <= r.upper[j],
                                    "source-decimal oracle point survives tightening");
                    }
                }
        if (!strcmp(r.status, "infeasible"))
            require(!feasible, "source contradiction agrees with exact decimal enumeration");
        unlink(path);
        lab_model_free(&m);
        lab_presolve_result_free(&r);
    }
    require(changes > 20, "decimal oracle campaign exercises propagation");
}
static void independent_point_check(const char *checker, const char *model,
                                     const char *solution) {
    if (!checker)
        return;
    char path[] = "/tmp/milp-presolve-solution-XXXXXX";
    int fd = mkstemp(path);
    require(fd >= 0, "checker solution file");
    FILE *f = fdopen(fd, "w");
    require(f && fputs(solution, f) >= 0 && !fclose(f), "write checker point");
    pid_t pid = fork();
    require(pid >= 0, "checker process");
    if (!pid) {
        int quiet = open("/dev/null", O_WRONLY);
        if (quiet >= 0) {
            dup2(quiet, STDOUT_FILENO);
            close(quiet);
        }
        execl(checker, checker, "--model", model, "--solution", path, "--exact", (char *)NULL);
        _exit(127);
    }
    int status;
    require(waitpid(pid, &status, 0) == pid && WIFEXITED(status) && WEXITSTATUS(status) == 0,
            "independent checker validates the retained original-model point");
    unlink(path);
}
static void cancellation_and_source_boxes(const char *checker) {
    LabPresolveOptions opt = options();
    LabPresolveResult r;
    LabModel m = {0};
    char path[64];
    read_text(&m,
        "NAME HIDDENCOEF\nROWS\n N O\n G R\nCOLUMNS\n"
        " x O 0 R 10000000000000001\n x R -10000000000000000\n"
        "RHS\n b R 1\nBOUNDS\n FX b x 1\nENDATA\n", path);
    require(m.exact && m.envelope && m.value[0] == 1,
            "exact coalescing restores the lost source coefficient");
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(strcmp(r.status, "infeasible") && r.lower[0] <= 1 && r.upper[0] >= 1,
            "zero-coalescing cannot create a false constant contradiction");
    independent_point_check(checker, path, "x 1\n");
    unlink(path); lab_model_free(&m); lab_presolve_result_free(&r);
    read_text(&m,
        "NAME SOURCEBOX\nROWS\n N O\n L R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R 1.00000000000000001\nBOUNDS\n"
        " LO b x 0.99999999999999999\n UP b x 1.00000000000000001\nENDATA\n", path);
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(strcmp(r.status, "infeasible") && !r.fixed && !r.redundant_rows,
            "rounded fixed-looking continuous box cannot narrow the outer source domain");
    independent_point_check(checker, path, "x 1\n");
    unlink(path); lab_model_free(&m); lab_presolve_result_free(&r);
    read_text(&m,
        "NAME FREEOVERLAY\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R -5\nBOUNDS\n FR b x\nENDATA\n", path);
    double lower[] = {-4};
    lab_presolve(&m, lower, NULL, &opt, &r);
    require(strcmp(r.status, "infeasible") && r.lower[0] == -4,
            "new finite overlay does not mistake an infinite-source placeholder for zero");
    unlink(path); lab_model_free(&m); lab_presolve_result_free(&r);
}
static void parity_equalities(void) {
    /* Free integer columns absorb even residuals. Neither equality alone
     * fixes y, and replacing v by a continuous variable removes the parity
     * condition on x entirely. */
    size_t start[] = {0,2,3,4,5}, row[] = {0,1,0,0,1};
    double a[] = {1,1,1,2,-2}, lo[] = {0,0,-INFINITY,-INFINITY};
    double hi[] = {1,1,INFINITY,INFINITY}, b[] = {-1,0};
    unsigned char integer[] = {1,1,1,1};
    LabModel m = {.rows=2, .cols=4, .nnz=5, .col_start=start, .row_index=row,
                  .value=a, .integer=integer, .col_lower=lo, .col_upper=hi,
                  .row_lower=b, .row_upper=b};
    LabPresolveOptions opt = options();
    LabPresolveResult r;
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(!strcmp(r.status, "unresolved") && r.upper[0] == 0 && r.lower[1] == 1,
            "combined integer equalities determine both binary variables");
    lab_presolve_result_free(&r);
    opt.relax = 1;
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(r.lower[1] == 0 && r.upper[0] == 1, "parity does not constrain the LP");
    lab_presolve_result_free(&r);
    integer[3] = 0;
    opt.relax = 0;
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(r.lower[1] == 0 && r.upper[0] == 1, "continuous even coefficient prevents inference");
    lab_presolve_result_free(&r);
    integer[3] = 1;
    a[1] = 2;
    b[1] = -1;
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(!strcmp(r.status, "infeasible"), "odd RHS contradicts an all-even integer row");
    lab_presolve_result_free(&r);
    /* Rounding this source coefficient to 2 would incorrectly reject the
     * exact feasible point x=10^17, v=-10^17. */
    char path[64];
    LabModel decimal = {0};
    read_text(&decimal, "NAME DECIMAL\nROWS\n N O\n E R\nCOLUMNS\n"
        " x R 2.00000000000000001\n v R 2\nRHS\n b R 1\nBOUNDS\n"
        " LI b x 0\n LI b v -1e18\nENDATA\n", path);
    lab_presolve(&decimal, NULL, NULL, &opt, &r);
    require(strcmp(r.status, "infeasible"), "rounded decimal coefficient has no parity");
    unlink(path);
    lab_model_free(&decimal);
    lab_presolve_result_free(&r);
}

static void ranges_relaxation_and_limits(void) {
    const char *kinds[] = {"L", "G", "E", "E"};
    const long ranges[] = {2, -2, 2, -2}, expected_lo[] = {3, 5, 5, 3};
    const long expected_hi[] = {5, 7, 7, 5};
    for (size_t case_id = 0; case_id < 4; case_id++) {
        char text[512], path[64];
        snprintf(text, sizeof text,
                 "NAME RANGE\nROWS\n N O\n %s R\nCOLUMNS\n x O 0 R 1\n"
                 "RHS\n b R 5\nRANGES\n g R %ld\nBOUNDS\n LI b x 0\n"
                 " UI b x 10\nENDATA\n", kinds[case_id], ranges[case_id]);
        LabModel m = {0};
        read_text(&m, text, path);
        LabPresolveOptions opt = options();
        LabPresolveResult r;
        lab_presolve(&m, NULL, NULL, &opt, &r);
        require(!strcmp(r.status, "unresolved") &&
                r.lower[0] == (double)expected_lo[case_id] &&
                r.upper[0] == (double)expected_hi[case_id] && r.redundant[0],
                "all signed range conventions produce source-safe integer bounds");
        unlink(path); lab_model_free(&m); lab_presolve_result_free(&r);
    }
    size_t start[] = {0, 0};
    double lo[] = {.2}, hi[] = {.8};
    unsigned char integrality[] = {1};
    LabModel m = {.cols = 1, .col_start = start, .col_lower = lo, .col_upper = hi,
                  .integer = integrality};
    LabPresolveOptions opt = options();
    LabPresolveResult r;
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(!strcmp(r.status, "infeasible"), "empty integer interval is detected");
    lab_presolve_result_free(&r);
    opt.relax = 1;
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(!strcmp(r.status, "unresolved") && r.lower[0] == .2 && r.upper[0] == .8,
            "LP relaxation mode must not round integer metadata");
    lab_presolve_result_free(&r);
    size_t rows[] = {0};
    double a[] = {1}, rl[] = {0}, ru[] = {1};
    start[1] = 1;
    lo[0] = 0; hi[0] = 2;
    m.rows = m.nnz = 1; m.row_index = rows; m.value = a;
    m.row_lower = rl; m.row_upper = ru;
    opt.time_limit = 1e-12;
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(!strcmp(r.status, "unresolved") && r.limited &&
                r.lower[0] == 0 && r.upper[0] == 2,
            "expired budget preserves a safe unresolved box");
    lab_presolve_result_free(&r);
    opt = options();
    rows[0] = 1;
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(!strcmp(r.status, "input_error"), "out-of-range CSC row rejected");
    lab_presolve_result_free(&r);
    rows[0] = 0;
    double widened[] = {-1};
    lab_presolve(&m, widened, NULL, &opt, &r);
    require(!strcmp(r.status, "input_error"), "numerical overlay cannot widen original domain");
    lab_presolve_result_free(&r);
}
static void duplicate_clique_column(void) {
    /* x + y <= 1 and x + x + y >= 2 admit exactly (x,y) = (1,0). */
    size_t start[] = {0, 3, 5}, row[] = {0, 1, 1, 0, 1};
    double a[] = {1, 1, 1, 1, 1}, lo[] = {0, 0}, hi[] = {1, 1};
    double rl[] = {-INFINITY, 2}, ru[] = {1, INFINITY};
    unsigned char integer[] = {1, 1};
    LabModel m = {.rows = 2, .cols = 2, .nnz = 5, .col_start = start,
                  .row_index = row, .value = a, .col_lower = lo, .col_upper = hi,
                  .row_lower = rl, .row_upper = ru, .integer = integer};
    LabPresolveOptions opt = options();
    LabPresolveResult r;
    lab_presolve(&m, NULL, NULL, &opt, &r);
    require(!strcmp(r.status, "unresolved") && r.upper[0] == 1 && r.upper[1] == 0,
            "clique probing keeps every coefficient of the column fixed to one");
    lab_presolve_result_free(&r);
}

int main(int argc, char **argv) {
    duplicate_clique_column();
    generated_integer_oracles();
    parity_equalities();
    one_sided_bounds_and_constants();
    large_sparse_chain();
    generated_decimal_oracles();
    cancellation_and_source_boxes(argc > 1 ? argv[1] : NULL);
    ranges_relaxation_and_limits();
    puts("Source-safe presolve:1100 exact enumerated MILPs,10000-column sparse chain, overlays and cancellation pass.");
    return 0;
}

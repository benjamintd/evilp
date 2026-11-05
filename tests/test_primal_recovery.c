#define _POSIX_C_SOURCE 200809L
#include "lp.h"
#include <stdio.h>
#include <gmp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void require(int condition, const char *message) {
    if (!condition) {
        fprintf(stderr, "primal recovery contract failed: %s\n", message);
        exit(1);
    }
}
static void read_text(LabModel *model, const char *text) {
    char path[] = "/tmp/milp-primal-recovery-XXXXXX";
    int fd = mkstemp(path);
    require(fd >= 0, "create temporary original-decimal model");
    FILE *stream = fdopen(fd, "w");
    require(stream != NULL, "open temporary model stream");
    require(fputs(text, stream) >= 0 && !fclose(stream), "write original model");
    int failed = lab_model_read(model, path);
    unlink(path);
    if (failed)
        fprintf(stderr, "%s\n", model->error);
    require(!failed, "read original model");
}
static LabLpOptions options(int dense, int no_scaling) {
    LabLpOptions result = {0};
    result.time_limit = 10;
    result.iteration_limit = 10000;
    result.max_dimension = 128;
    result.model_sha256 = "0000000000000000000000000000000000000000000000000000000000000000";
    result.dense = dense;
    result.no_scaling = no_scaling;
    return result;
}
/* This oracle does not read the engine's exact provenance or canonical matrix.
 * Source arithmetic is x + (1/10)*2^53 >= 900719925474099.25, hence x>=1/20.
 * Both the submitted row scale and its normalized original units are checked.
 * mpq_set_d evaluates the returned binary point exactly, without relying on
 * long double having more precision on a particular platform. */
static void point(const LabLpResult *result, double lo, double hi, const char *scale,
                  double expected) {
    require(result->has_primal && result->primal != NULL, "optimal result exports an original point");
    require(isfinite(result->primal[0]) && isfinite(result->primal[1]), "finite original point");
    require(result->primal[1] == 9007199254740992.0, "fixed column survives reconstruction");
    mpq_t x, z, value, residual, tolerance, bound, factor;
    mpq_inits(x, z, value, residual, tolerance, bound, factor, NULL);
    mpq_set_d(x, result->primal[0]);
    mpq_set_d(z, result->primal[1]);
    mpq_set_ui(tolerance, 1, 10000000);
    mpq_set_ui(value, 1, 10);
    mpq_mul(value, value, z);
    mpq_add(value, value, x);
    require(!mpq_set_str(bound, "3602879701896397/4", 10), "oracle row bound");
    mpq_sub(residual, bound, value);
    if (mpq_cmp(residual, tolerance) > 0)
        gmp_fprintf(stderr, "normalized exact row violation=%Qd, x=%Qd, bounds=[%g,%g], "
                           "row_scale=%s, warm=%d, scaled_rows=%zu, scaled_columns=%zu\n",
                    residual, x, lo, hi, scale, result->warm_start, result->scaled_rows,
                    result->scaled_columns);
    require(mpq_cmp(residual, tolerance) <= 0, "normalized original row tolerance remains 1e-7");
    require(!mpq_set_str(factor, scale, 10), "oracle row scale");
    mpq_mul(residual, residual, factor);
    require(mpq_cmp(residual, tolerance) <= 0, "submitted original row tolerance remains 1e-7");
    mpq_set_d(bound, lo);
    mpq_sub(residual, bound, x);
    require(mpq_cmp(residual, tolerance) <= 0, "active lower bound checked exactly");
    mpq_set_d(bound, hi);
    mpq_sub(residual, x, bound);
    require(mpq_cmp(residual, tolerance) <= 0, "active upper bound checked exactly");
    mpq_set_d(bound, expected);
    mpq_sub(residual, x, bound);
    mpq_abs(residual, residual);
    require(mpq_cmp(residual, tolerance) <= 0, "point cost agrees with independently known optimum");
    require(isfinite(result->objective) && fabs(result->objective - result->primal[0]) <= 1e-7,
            "reported objective agrees with original point");
    require(isfinite(result->primal_residual) && result->primal_residual <= 1e-7,
            "reported original residual keeps the existing guard");
    mpq_clears(x, z, value, residual, tolerance, bound, factor, NULL);
}
static void outcome(const LabLpResult *result, double lo, double hi, const char *scale,
                    int feasible, int mandatory) {
    require(result->status != NULL, "solve returns an explicit status");
    if (!strcmp(result->status, "optimal")) {
        require(feasible, "an exactly contradictory overlay cannot export an optimum");
        point(result, lo, hi, scale, fmax(.05, lo));
        return;
    }
    if (mandatory)
        fprintf(stderr, "expected recoverable optimum; got %s (%s)\n", result->status,
                result->message);
    require(!mandatory, "well-conditioned source gap must be recovered");
    require(!result->has_primal, "uncertain or infeasible result exports no primal");
    require(!strcmp(result->status, "numerical_failure") ||
                (!feasible && !strcmp(result->status, "infeasible")),
            "uncertainty remains explicit rather than a false model claim");
}
static void compare(LabLpSession *session, const LabModel *original, const LabLpOptions *opt,
                    double lo, double hi, const char *scale, int feasible, int mandatory) {
    double lower[] = {lo, 9007199254740992.0}, upper[] = {hi, 9007199254740992.0};
    LabModel overlay = *original;
    overlay.col_lower = lower;
    overlay.col_upper = upper;
    LabLpResult warm = {0}, cold = {0};
    lab_lp_session_solve(session, &warm);
    lab_lp_solve(&overlay, opt, &cold);
    outcome(&warm, lo, hi, scale, feasible, mandatory);
    outcome(&cold, lo, hi, scale, feasible, mandatory);
    lab_lp_result_free(&warm);
    lab_lp_result_free(&cold);
    require(original->col_lower[0] == 0 && original->col_upper[0] == 1,
            "session overlays leave original bounds immutable");
}
static void sequence(int dense, int no_scaling, int units) {
    const char *coefficient[] = {"1", "1e-20", "1e20"};
    const char *fixed_coefficient[] = {".1", "1e-21", "1e19"};
    const char *rhs[] = {"900719925474099.25", "900719925474099.25e-20",
                         "900719925474099.25e20"};
    const char *scale[] = {"1", "1/100000000000000000000", "100000000000000000000"};
    char text[1024];
    /* Separate record fields explicitly; test source text is independent of
     * the optimizer's model writer and exact-arithmetic helper routines. */
    snprintf(text, sizeof text,
             "NAME RECOVERY\nROWS\n N O\n G R\nCOLUMNS\n x O 1 R %s\n z O 0 R %s\n"
             "RHS\n b R %s\nBOUNDS\n LO b x 0\n UP b x 1\n FX b z 9007199254740992\nENDATA\n",
             coefficient[units], fixed_coefficient[units], rhs[units]);
    LabModel model = {0};
    read_text(&model, text);
    require(model.rows == 1 && model.cols == 2, "oracle source dimensions");
    LabLpOptions opt = options(dense, no_scaling);
    LabLpSession *session = lab_lp_session_create(&model, &opt);
    require(session != NULL, "create source model session");
    int mandatory = units == 0;
    compare(session, &model, &opt, 0, 1, scale[units], 1, mandatory);
    const double lower[] = {.25, 0, .25, 0, 0};
    const double upper[] = {1, .5, .25, .015625, 1};
    for (size_t i = 0; i < sizeof lower / sizeof *lower; i++) {
        require(!lab_lp_session_set_bounds(session, 0, lower[i], upper[i]), "valid finite overlay");
        int feasible = upper[i] >= .05;
        compare(session, &model, &opt, lower[i], upper[i], scale[units], feasible,
                feasible && mandatory);
    }
    lab_lp_session_free(session);
    lab_model_free(&model);
}
static void canonical_negative_bound_slack(int dense, int no_scaling) {
    /* The original width 2^39 induces a canonical column scale 2^39. Tightening
     * the lower bound by 1 creates a basic slack -2^-39, smaller than a 1e-8
     * canonical tolerance but a full unit violation in the original variable.
     * The implementation may repair the basis or rebuild; either must return
     * the analytic interval optimum rather than publish the clamped point 0. */
    LabModel model = {0};
    read_text(&model, "NAME NEGATIVE\nROWS\n N O\nCOLUMNS\n x O 1\n"
                      "BOUNDS\n LO b x 0\n UP b x 549755813888\nENDATA\n");
    LabLpOptions opt = options(dense, no_scaling);
    LabLpSession *session = lab_lp_session_create(&model, &opt);
    require(session != NULL, "create wide-column session");
    const double lower[] = {0, 1, 4096, 0};
    for (size_t i = 0; i < sizeof lower / sizeof *lower; i++) {
        require(!lab_lp_session_set_bounds(session, 0, lower[i], 549755813888.0),
                "wide-column overlay");
        double lo[] = {lower[i]}, hi[] = {549755813888.0};
        LabModel overlay = model;
        overlay.col_lower = lo;
        overlay.col_upper = hi;
        LabLpResult warm = {0}, cold = {0};
        lab_lp_session_solve(session, &warm);
        lab_lp_solve(&overlay, &opt, &cold);
        const LabLpResult *results[] = {&warm, &cold};
        for (size_t k = 0; k < sizeof results / sizeof *results; k++) {
            const LabLpResult *r = results[k];
            require(r->status != NULL && !strcmp(r->status, "optimal") && r->has_primal && r->primal,
                    "wide-column original bound must be restored");
            require(isfinite(r->primal[0]) && fabs(r->primal[0] - lower[i]) <= 1e-7,
                    "wide-column point satisfies original units and exact known optimum");
            require(fabs(r->objective - lower[i]) <= 1e-7,
                    "wide-column objective matches original interval");
        }
        lab_lp_result_free(&warm);
        lab_lp_result_free(&cold);
    }
    require(model.col_lower[0] == 0 && model.col_upper[0] == 549755813888.0,
            "wide-column source bounds remain immutable");
    lab_lp_session_free(session);
    lab_model_free(&model);
}
int main(void) {
    for (int dense = 0; dense < 2; dense++)
        for (int no_scaling = 0; no_scaling < 2; no_scaling++) {
            canonical_negative_bound_slack(dense, no_scaling);
            for (int units = 0; units < 3; units++)
                sequence(dense, no_scaling, units);
        }
    puts("primal recovery session contracts passed");
    return 0;
}

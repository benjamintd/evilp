#include "bb.h"
#include "common.h"
#include "lp.h"
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int error(const char *status, const char *message, int code) {
    printf("{\"schema\":1,\"status\":");
    lab_json_string(stdout, status);
    printf(",\"message\":");
    lab_json_string(stdout, message);
    puts(",\"certifies_optimality\":false,\"has_primal\":false,\"objective\":null}");
    return code;
}
static int positive_size(const char *s, size_t *n, size_t cap) {
    if (*s < '0' || *s > '9')
        return -1;
    char *end;
    errno = 0;
    unsigned long long v = strtoull(s, &end, 10);
    if (errno || *end || !v || v > cap)
        return -1;
    *n = (size_t)v;
    return 0;
}
/* Choose the shortest significant-digit decimal that roundtrips to the same
 * double. This avoids exporting .29999999999999999 for input .3. Independent
 * rational checking still decides whether that decimal is a valid witness. */
static void decimal(FILE *f, double x) {
    char text[64];
    for (int digits = 1; digits <= 17; digits++) {
        snprintf(text, sizeof text, "%.*g", digits, x == 0 ? 0 : x);
        char *end;
        double parsed = strtod(text, &end);
        if (!*end && parsed == x)
            break;
    }
    fputs(text, f);
}
static int write_vector(const char *path, const LabModel *m, const LabLpResult *r, int dual) {
    if (!path)
        return 0;
    FILE *f = fopen(path, "wx");
    if (!f)
        return -1;
    if (dual) {
        for (size_t i = 0; i < m->rows; i++) {
            fprintf(f, "%s ", m->row_name[i]);
            decimal(f, r->dual[i]);
            fputc('\n', f);
        }
    } else {
        fprintf(f, "=obj= ");
        decimal(f, r->objective);
        fputc('\n', f);
        for (size_t j = 0; j < m->cols; j++) {
            fprintf(f, "%s ", m->col_name[j]);
            decimal(f, r->primal[j]);
            fputc('\n', f);
        }
    }
    int bad = ferror(f);
    if (fclose(f))
        bad = 1;
    return bad ? -1 : 0;
}
static int status_code(const char *status) {
    if (!strcmp(status, "input_error"))
        return 2;
    if (!strcmp(status, "unsupported"))
        return 3;
    if (!strcmp(status, "optimal") || !strcmp(status, "infeasible") || !strcmp(status, "unbounded"))
        return 0;
    return 4;
}
static int write_farkas(const char *path, const LabModel *m, const LabLpResult *r) {
    if (!path || !r->farkas)
        return 0;
    FILE *f = fopen(path, "wx");
    if (!f)
        return -1;
    fputs("MILPLAB_FARKAS_1\n", f);
    for (size_t i = 0; i < r->farkas_rows; i++)
        fprintf(f, "%s %s\n", m->row_name[i], r->farkas[i]);
    int bad = ferror(f);
    if (fclose(f))
        bad = 1;
    return bad ? -1 : 0;
}
int main(int argc, char **argv) {
    if (argc == 2 && !strcmp(argv[1], "--help")) {
        puts("milpsolve MODEL (--relax | --mip) [--solution NEW_FILE]\n"
             "          [--time-limit SECONDS] [--iteration-limit N] [--max-dimension N]\n"
             "          [--backend sparse|dense]\n"
             "          [--no-scaling] [--no-refinement] [--no-exact] [--bland-ratio]\n"
             "          LP: [--dual NEW_FILE] [--farkas NEW_FILE] [--basis NEW_FILE] [--read-basis FILE]\n"
             "          MILP: [--node-limit N] [--cold] [--trace NEW_FILE]\n"
             "C17 LP and numerical branch-and-bound prototype. Select problem mode explicitly.\n"
             "JSON statuses are numerical solver claims; use milpcheck for validation.\n"
             "Default standard-form dimension limit4096; sparse hard cap32768,dense hard cap4096.");
        return 0;
    }
    const char *path = NULL, *solution = NULL, *dual = NULL, *trace_path = NULL, *farkas = NULL;
    int relax = 0, mip = 0, cold = 0, node_limit_set = 0;
    size_t node_limit = 100000;
    LabLpOptions options = {.time_limit = 60, .iteration_limit = 100000, .max_dimension = 4096};
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--relax")) {
            relax = 1;
        } else if (!strcmp(a, "--mip")) {
            mip = 1;
        } else if (!strcmp(a, "--cold")) {
            cold = 1;
        } else if (!strcmp(a, "--no-scaling")) {
            options.no_scaling = 1;
        } else if (!strcmp(a, "--no-refinement")) {
            options.no_refinement = 1;
        } else if (!strcmp(a, "--no-exact")) {
            options.no_exact = 1;
        } else if (!strcmp(a, "--bland-ratio")) {
            options.bland_ratio = 1;
        } else if (!strcmp(a, "--solution") || !strcmp(a, "--dual") || !strcmp(a, "--basis") ||
                   !strcmp(a, "--read-basis") || !strcmp(a, "--trace") || !strcmp(a, "--farkas")) {
            if (++i == argc)
                return error("input_error", "Missing option value", 2);
            if (!strcmp(a, "--solution"))
                solution = argv[i];
            else if (!strcmp(a, "--dual"))
                dual = argv[i];
            else if (!strcmp(a, "--basis"))
                options.basis_out = argv[i];
            else if (!strcmp(a, "--trace"))
                trace_path = argv[i];
            else if (!strcmp(a, "--farkas"))
                farkas = argv[i];
            else
                options.basis_in = argv[i];
        } else if (!strcmp(a, "--backend")) {
            if (++i == argc)
                return error("input_error", "Missing backend name", 2);
            if (!strcmp(argv[i], "dense"))
                options.dense = 1;
            else if (!strcmp(argv[i], "sparse"))
                options.dense = 0;
            else
                return error("input_error", "Backend must be sparse or dense", 2);
        } else if (!strcmp(a, "--time-limit")) {
            if (++i == argc)
                return error("input_error", "Missing time limit", 2);
            char *end;
            errno = 0;
            options.time_limit = strtod(argv[i], &end);
            if (errno || end == argv[i] || *end || !isfinite(options.time_limit) ||
                options.time_limit <= 0 || options.time_limit > 1e9)
                return error("input_error", "Invalid time limit", 2);
        } else if (!strcmp(a, "--iteration-limit") || !strcmp(a, "--max-dimension") ||
                   !strcmp(a, "--node-limit")) {
            if (++i == argc)
                return error("input_error", "Missing integer limit", 2);
            size_t *v = !strcmp(a, "--iteration-limit") ? &options.iteration_limit
                        : !strcmp(a, "--node-limit")    ? &node_limit
                                                        : &options.max_dimension;
            if (!strcmp(a, "--node-limit"))
                node_limit_set = 1;
            if (positive_size(argv[i], v, !strcmp(a, "--max-dimension") ? 32768 : SIZE_MAX))
                return error("input_error", "Invalid positive integer limit", 2);
        } else if (a[0] == '-' || path) {
            return error("input_error", "Unknown option or extra model path", 2);
        } else {
            path = a;
        }
    }
    if (!path || relax == mip)
        return error("input_error", "Provide a model and exactly one of --relax or --mip", 2);
    if (mip && (dual || farkas || options.basis_in || options.basis_out))
        return error("input_error",
                     "Dual/basis files are LP-only; MILP tree proofs are not implemented", 2);
    if (!mip && (cold || node_limit_set || trace_path))
        return error("input_error", "Cold search, node limits and traces require --mip", 2);
    if (options.dense && options.max_dimension > 4096)
        return error("input_error", "Dense backend hard dimension limit is4096", 2);
    double started = lab_now();
    char hash[65];
    if (lab_sha256(path, hash))
        return error("input_error", "Cannot hash model", 2);
    options.model_sha256 = hash;
    LabModel model = {0};
    if (lab_model_read(&model, path)) {
        int code = error(model.unsupported ? "unsupported" : "input_error", model.error,
                         model.unsupported ? 3 : 2);
        lab_model_free(&model);
        return code;
    }
    /* Reading counts toward the CLI budget; each simplex iteration checks it. */
    options.time_limit -= lab_now() - started;
    if (mip) {
        LabBbOptions bb_options = {.lp = options,
                                   .time_limit = options.time_limit,
                                   .node_limit = node_limit,
                                   .cold = cold,
                                   .trace_path = trace_path};
        LabBbResult r;
        lab_bb_solve(&model, &bb_options, &r);
        if (r.has_primal) {
            LabLpResult exported = {.objective = r.objective, .primal = r.primal};
            if (write_vector(solution, &model, &exported, 0)) {
                r.status = "input_error";
                snprintf(r.message, sizeof r.message, "Cannot create fresh solution file");
            }
        }
        printf("{\"schema\":1,\"solver\":\"milp-lab\",\"status\":");
        lab_json_string(stdout, r.status);
        printf(",\"message\":");
        lab_json_string(stdout, r.message);
        printf(",\"relaxation\":false,\"certifies_optimality\":false,\"bound_certified\":false,"
               "\"solver_status_is_claim\":true,\"model_sha256\":");
        lab_json_string(stdout, hash);
        printf(",\"objective_sense\":\"%s\",\"has_primal\":%s,\"objective\":",
               model.maximize ? "max" : "min", r.has_primal ? "true" : "false");
        if (r.has_primal)
            printf("%.17g", r.objective);
        else
            fputs("null", stdout);
        printf(",\"best_bound\":");
        if (r.has_bound && isfinite(r.best_bound))
            printf("%.17g", r.best_bound);
        else
            fputs("null", stdout);
        printf(",\"rows\":%zu,\"columns\":%zu,\"nonzeros\":%zu,\"nodes\":%zu,"
               "\"lp_solves\":%zu,\"lp_iterations\":%zu,\"refactorizations\":%zu,"
               "\"warm_starts\":%zu,\"cold_starts\":%zu,\"pruned_nodes\":%zu,"
               "\"dual_iterations\":%zu,\"cold_fallbacks\":%zu,\"factor_nnz_peak\":%zu,"
               "\"standard_nnz_peak\":%zu,"
               "\"solve_wall_s\":%.9g,\"total_wall_s\":%.9g,\"backend\":\"%s\",\"cold\":%s,"
               "\"scaling_enabled\":%s,\"refinement_enabled\":%s,\"refinement_steps\":%zu,"
               "\"scaled_rows\":%zu,\"scaled_columns\":%zu,\"initial_backward_error\":%.17g,"
               "\"max_backward_error\":%.17g,\"exact_enabled\":%s,\"exact_attempts\":%zu,"
               "\"exact_infeasibilities\":%zu,\"ratio_test\":\"%s\",\"exact_model_available\":%s,"
               "\"rejected_pivots\":%zu,\"harris_rescues\":%zu}\n",
               model.rows, model.cols, model.nnz, r.nodes, r.lp_solves, r.lp_iterations,
               r.refactorizations, r.warm_starts, r.cold_starts, r.pruned_nodes, r.dual_iterations,
               r.cold_fallbacks, r.factor_nnz_peak, r.standard_nnz_peak, r.elapsed,
               lab_now() - started, options.dense ? "dense_lu_reference" : "sparse_lu",
               cold ? "true" : "false", options.no_scaling ? "false" : "true",
               options.no_refinement ? "false" : "true", r.refinement_steps, r.scaled_rows,
               r.scaled_columns, r.initial_backward_error, r.max_backward_error,
               options.no_exact ? "false" : "true", r.exact_attempts, r.exact_infeasibilities,
               options.bland_ratio ? "bland" : "guarded_bland_harris", model.exact ? "true" : "false",
               r.rejected_pivots, r.harris_rescues);
        int code = status_code(r.status);
        lab_bb_result_free(&r);
        lab_model_free(&model);
        return code;
    }
    LabLpResult r;
    lab_lp_solve(&model, &options, &r);
    if (write_farkas(farkas, &model, &r)) {
        r.status = "input_error";
        snprintf(r.message, sizeof r.message, "Cannot create fresh Farkas file");
    }
    if (r.has_primal &&
        (write_vector(solution, &model, &r, 0) || write_vector(dual, &model, &r, 1))) {
        r.status = "input_error";
        snprintf(r.message, sizeof r.message, "Cannot create fresh solution/dual file");
    }
    printf("{\"schema\":1,\"solver\":\"milp-lab\",\"status\":");
    lab_json_string(stdout, r.status);
    printf(",\"message\":");
    lab_json_string(stdout, r.message);
    printf(",\"relaxation\":true,\"certifies_optimality\":false,\"solver_status_is_claim\":true,"
           "\"model_sha256\":");
    lab_json_string(stdout, hash);
    printf(",\"objective_sense\":\"%s\",\"has_primal\":%s,\"objective\":",
           model.maximize ? "max" : "min", r.has_primal ? "true" : "false");
    if (r.has_primal)
        printf("%.17g", r.objective);
    else
        printf("null");
    printf(
        ",\"rows\":%zu,\"columns\":%zu,\"nonzeros\":%zu,\"standard_rows\":%zu,"
        "\"standard_columns\":%zu,\"redundant_rows\":%zu,\"iterations\":%zu,"
        "\"refactorizations\":%zu,\"primal_residual\":%.17g,\"dual_residual\":%.17g,"
        "\"solve_wall_s\":%.9g,\"total_wall_s\":%.9g,\"backend\":\"%s\","
        "\"warm_start\":%s,\"dual_iterations\":%zu,\"cold_fallbacks\":%zu,"
        "\"factor_nnz_peak\":%zu,\"standard_nnz\":%zu,\"scaling_enabled\":%s,"
        "\"refinement_enabled\":%s,\"refinement_steps\":%zu,\"scaled_rows\":%zu,"
        "\"scaled_columns\":%zu,\"initial_backward_error\":%.17g,\"max_backward_error\":%.17g,"
        "\"exact_enabled\":%s,\"exact_attempts\":%zu,\"exact_infeasibilities\":%zu,"
        "\"has_farkas_candidate\":%s,\"certifies_infeasibility\":false,"
        "\"ratio_test\":\"%s\",\"exact_model_available\":%s,"
        "\"rejected_pivots\":%zu,\"harris_rescues\":%zu}\n",
        model.rows, model.cols, model.nnz, r.standard_rows, r.standard_columns, r.redundant_rows,
        r.iterations, r.refactorizations, r.primal_residual, r.dual_residual, r.elapsed,
        lab_now() - started, options.dense ? "dense_lu_reference" : "sparse_lu",
        r.warm_start ? "true" : "false", r.dual_iterations, r.cold_fallbacks, r.factor_nnz_peak,
        r.standard_nnz, options.no_scaling ? "false" : "true",
        options.no_refinement ? "false" : "true", r.refinement_steps, r.scaled_rows,
        r.scaled_columns, r.initial_backward_error, r.max_backward_error,
        options.no_exact ? "false" : "true", r.exact_attempts, r.exact_infeasibilities,
        r.farkas ? "true" : "false", options.bland_ratio ? "bland" : "guarded_bland_harris",
        model.exact ? "true" : "false", r.rejected_pivots, r.harris_rescues);
    int code = status_code(r.status);
    lab_lp_result_free(&r);
    lab_model_free(&model);
    return code;
}

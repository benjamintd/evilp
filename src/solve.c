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
int main(int argc, char **argv) {
    if (argc == 2 && !strcmp(argv[1], "--help")) {
        puts("milpsolve MODEL --relax [--solution NEW_FILE] [--dual NEW_FILE]\n"
             "          [--time-limit SECONDS] [--iteration-limit N] [--max-dimension N]\n"
             "          [--basis NEW_FILE] [--read-basis FILE]\n"
             "C17 revised-simplex LP reference engine. Integrality is relaxed explicitly.\n"
             "JSON statuses are numerical solver claims; use milpcheck for validation.\n"
             "Dense backend: default maximum standard-form dimension 512, hard cap 4096.");
        return 0;
    }
    const char *path = NULL, *solution = NULL, *dual = NULL;
    int relax = 0;
    LabLpOptions options = {60, 100000, 512, NULL, NULL, NULL};
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--relax")) {
            relax = 1;
        } else if (!strcmp(a, "--solution") || !strcmp(a, "--dual") || !strcmp(a, "--basis") ||
                   !strcmp(a, "--read-basis")) {
            if (++i == argc)
                return error("input_error", "Missing option value", 2);
            if (!strcmp(a, "--solution"))
                solution = argv[i];
            else if (!strcmp(a, "--dual"))
                dual = argv[i];
            else if (!strcmp(a, "--basis"))
                options.basis_out = argv[i];
            else
                options.basis_in = argv[i];
        } else if (!strcmp(a, "--time-limit")) {
            if (++i == argc)
                return error("input_error", "Missing time limit", 2);
            char *end;
            errno = 0;
            options.time_limit = strtod(argv[i], &end);
            if (errno || end == argv[i] || *end || !isfinite(options.time_limit) ||
                options.time_limit <= 0 || options.time_limit > 1e9)
                return error("input_error", "Invalid time limit", 2);
        } else if (!strcmp(a, "--iteration-limit") || !strcmp(a, "--max-dimension")) {
            if (++i == argc)
                return error("input_error", "Missing integer limit", 2);
            size_t *v =
                !strcmp(a, "--iteration-limit") ? &options.iteration_limit : &options.max_dimension;
            if (positive_size(argv[i], v, !strcmp(a, "--max-dimension") ? 4096 : SIZE_MAX))
                return error("input_error", "Invalid positive integer limit", 2);
        } else if (a[0] == '-' || path) {
            return error("input_error", "Unknown option or extra model path", 2);
        } else {
            path = a;
        }
    }
    if (!path || !relax)
        return error("input_error", "Provide a model and --relax; MILP search is not implemented",
                     2);
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
    LabLpResult r;
    lab_lp_solve(&model, &options, &r);
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
    printf(",\"rows\":%zu,\"columns\":%zu,\"nonzeros\":%zu,\"standard_rows\":%zu,"
           "\"standard_columns\":%zu,\"redundant_rows\":%zu,\"iterations\":%zu,"
           "\"refactorizations\":%zu,\"primal_residual\":%.17g,\"dual_residual\":%.17g,"
           "\"solve_wall_s\":%.9g,\"total_wall_s\":%.9g,\"backend\":\"dense_lu_reference\"}\n",
           model.rows, model.cols, model.nnz, r.standard_rows, r.standard_columns, r.redundant_rows,
           r.iterations, r.refactorizations, r.primal_residual, r.dual_residual, r.elapsed,
           lab_now() - started);
    int code = !strcmp(r.status, "input_error")   ? 2
               : !strcmp(r.status, "unsupported") ? 3
               : !strcmp(r.status, "optimal") || !strcmp(r.status, "infeasible") ||
                       !strcmp(r.status, "unbounded")
                   ? 0
                   : 4;
    lab_lp_result_free(&r);
    lab_model_free(&model);
    return code;
}

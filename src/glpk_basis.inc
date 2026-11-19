/* GLPK proposes a starting basis. A standalone caller may also accept a
 * point and dual after checking both against the original source model. */
#ifdef LAB_GLPK_BASIS
#include <glpk.h>
#include <limits.h>

typedef struct {
    double lower, upper;
    size_t lower_row, upper_row;
} SuggestionBounds;

/* In a native variable-bound equation a*z + s*t = b, t >= 0.
 * A positive s gives a*z <= b, a negative s gives a*z >= b. */
static int collect_suggestion_bound(const Work *w, size_t row, double a,
                                    SuggestionBounds *bounds) {
    const Column *slack = &w->a[w->structural + row];
    if (slack->count != 1 || a == 0)
        return -1;
    double bound = w->b[row] / a;
    if (!isfinite(bound))
        return -1;
    if ((a > 0) == (slack->value[0] > 0)) {
        if (bound < bounds->upper) {
            bounds->upper = bound;
            bounds->upper_row = row;
        }
    } else if (bound > bounds->lower) {
        bounds->lower = bound;
        bounds->lower_row = row;
    }
    return bounds->lower > bounds->upper ? -1 : 0;
}

/* Lift GLPK's row-bounded, variable-bounded basis to native equality form.
 * Omitted variable-bound rows start with their slack basic. A variable at a
 * nonzero nonbasic bound replaces the slack of its active bound row. */
static void lift_suggested_basis(Work *w, glp_prob *problem, const size_t *source_row,
                                 size_t rows, const SuggestionBounds *bounds) {
    size_t *proposal = lab_alloc(w->m, sizeof *proposal);
    for (size_t i = 0; i < w->m; i++)
        proposal[i] = w->structural + i;
    for (size_t k = 0; k < rows; k++) {
        int id = glp_get_bhead(problem, (int)k + 1);
        if (id <= (int)rows) {
            size_t slack = w->structural + source_row[id - 1];
            /* A basic equality logical has zero displacement. Use its
             * Phase I artificial until native cleanup establishes rank. */
            proposal[source_row[k]] = w->a[slack].count ? slack : slack + w->m;
        } else {
            proposal[source_row[k]] = (size_t)(id - (int)rows - 1);
        }
    }
    for (size_t j = 0; j < w->structural; j++) {
        int status = glp_get_col_stat(problem, (int)j + 1);
        size_t side = SIZE_MAX;
        if (status == GLP_NU)
            side = bounds[j].upper_row;
        else if ((status == GLP_NL || status == GLP_NS) && bounds[j].lower != 0)
            side = bounds[j].lower_row;
        if (side != SIZE_MAX)
            proposal[side] = j;
    }
    size_t *original = w->basis;
    const char *saved_status = w->result->status;
    char saved_message[sizeof w->result->message];
    memcpy(saved_message, w->result->message, sizeof saved_message);
    w->basis = proposal;
    if (refactor(w, 1)) {
        w->result->glpk_basis_status = "native_basis_rejected";
        w->basis = original;
        w->result->status = saved_status;
        memcpy(w->result->message, saved_message, sizeof saved_message);
        free(proposal);
    } else {
        w->result->glpk_basis_status = "accepted";
        w->result->glpk_bases_accepted++;
        free(original);
    }
}

/* At an optimal vertex a fixed equality row can still be basic. Exchange
 * it at zero primal step for a true structural variable or inequality
 * slack, using the dual ratio test to retain the proposed dual solution. */
static int exchange_fixed_rows(Work *w, glp_prob *problem, const size_t *source_row,
                               size_t rows, int *index, double *value, double deadline) {
    unsigned exchanges = 0;
    for (size_t k = 0; k < rows && exchanges < 128; k++) {
        if (lab_now() >= deadline || w->result->iterations >= w->options->iteration_limit ||
            glp_get_dual_stat(problem) != GLP_FEAS)
            break;
        int row = (int)k + 1;
        if (w->a[w->structural + source_row[k]].count ||
            glp_get_row_stat(problem, row) != GLP_BS)
            continue;
        int length = glp_eval_tab_row(problem, row, index, value);
        int pivot = glp_dual_rtest(problem, length, index, value, 1, 1e-9);
        if (!pivot)
            pivot = glp_dual_rtest(problem, length, index, value, -1, 1e-9);
        if (!pivot)
            continue;
        int entering = index[pivot];
        glp_set_row_stat(problem, row, GLP_NS);
        if (entering <= (int)rows)
            glp_set_row_stat(problem, entering, GLP_BS);
        else
            glp_set_col_stat(problem, entering - (int)rows, GLP_BS);
        w->result->glpk_basis_iterations++;
        w->result->iterations++;
        exchanges++;
        if (glp_warm_up(problem) || glp_get_prim_stat(problem) != GLP_FEAS)
            return -1;
    }
    return 0;
}

/* A standalone caller needs a checked point and bound, without needing a
 * saved basis for a future bound change. Map the proposal back to the
 * original rows and columns before checking feasibility and dual support. */
static int accept_simplex_proposal(Work *w, glp_prob *p, const size_t *source_row,
                                   size_t rows, double deadline) {
    const LabModel *m = w->model;
    if (!w->standalone || w->options->no_exact || m->exact || !m->envelope ||
        source_objective_risk(w, deadline))
        return 0;
    double *x = lab_alloc(m->cols, sizeof *x), *y = lab_alloc(m->rows, sizeof *y);
    memcpy(x, w->shift, m->cols * sizeof *x);
    for (size_t j = 0; j < w->structural; ++j)
        x[w->original[j]] += w->sign[j] * glp_get_col_prim(p, (int)j + 1);
    for (size_t k = 0; k < rows; ++k) {
        size_t i = source_row[k];
        y[w->source_row[i]] += w->row_scale[i] * glp_get_row_dual(p, (int)k + 1);
    }
    double obj, residual;
    int available = 0;
    if (lab_exact_evaluate(m, x, &obj, &residual) || residual > 1e-7 ||
        source_support_interval_gap(w, x, y, deadline, &available) || available != 1 ||
        lab_now() >= deadline) {
        free(x);
        free(y);
        return 0;
    }
    LabLpResult *r = w->result;
    r->primal = x;
    r->dual = y;
    r->objective = obj;
    r->primal_residual = residual;
    r->has_primal = 1;
    r->source_support_scans++;
    r->glpk_basis_status = "source_verified_point_and_bound";
    message(w, "optimal", "Original-model primal and absolute dual support gap checked");
    return 1;
}

static void suggest_glpk_basis(Work *w) {
    size_t m = w->m, n = w->structural;
    if (!w->options->glpk_basis || w->options->basis_in || !m || !n ||
        m > 250000 || n > 250000 || w->result->standard_nnz > 2000000 ||
        w->result->iterations >= w->options->iteration_limit)
        return;
    w->result->glpk_basis_status = "setup_incomplete";
    double deadline = w->started + w->options->time_limit;
    if (lab_now() >= deadline)
        return;
    size_t *source_row = lab_alloc(m, sizeof *source_row);
    int *row_id = lab_alloc(m, sizeof *row_id);
    size_t rows = 0;
    for (size_t i = 0; i < m; i++)
        if (w->source_row[i] != SIZE_MAX) {
            source_row[rows++] = i;
            row_id[i] = (int)rows;
        }
    glp_prob *problem = glp_create_prob();
    if (rows)
        glp_add_rows(problem, (int)rows);
    glp_add_cols(problem, (int)n);
    for (size_t k = 0; k < rows; k++) {
        size_t i = source_row[k];
        const Column *slack = &w->a[n + i];
        int type = !slack->count ? GLP_FX : slack->value[0] > 0 ? GLP_UP : GLP_LO;
        glp_set_row_bnds(problem, (int)k + 1, type, w->b[i], w->b[i]);
    }
    size_t scratch = rows > n ? rows : n;
    int *index = lab_alloc(scratch + 1, sizeof *index);
    double *value = lab_alloc(scratch + 1, sizeof *value);
    size_t *seen = lab_alloc(m, sizeof *seen);
    SuggestionBounds *bounds = lab_alloc(n, sizeof *bounds);
    for (size_t j = 0; j < n; j++) {
        if ((j & 255) == 0 && lab_now() >= deadline)
            goto done;
        bounds[j] = (SuggestionBounds){0, INFINITY, SIZE_MAX, SIZE_MAX};
        const Column *c = &w->a[j];
        int length = 0;
        for (size_t k = 0; k < c->count; k++) {
            size_t i = c->row[k];
            if (seen[i] == j + 1)
                goto done; /* GLPK requires unique row indices */
            seen[i] = j + 1;
            if (row_id[i]) {
                index[++length] = row_id[i];
                value[length] = c->value[k];
            } else if (collect_suggestion_bound(w, i, c->value[k], &bounds[j])) {
                goto done;
            }
        }
        glp_set_mat_col(problem, (int)j + 1, length, index, value);
        double lo = bounds[j].lower, hi = bounds[j].upper;
        int type = lo == hi ? GLP_FX : isfinite(hi) ? GLP_DB : GLP_LO;
        glp_set_col_bnds(problem, (int)j + 1, type, lo, isfinite(hi) ? hi : 0);
        glp_set_obj_coef(problem, (int)j + 1, w->cost[j]);
    }
    glp_smcp parameters;
    glp_init_smcp(&parameters);
    parameters.msg_lev = GLP_MSG_OFF;
    parameters.meth = GLP_DUALP;
    parameters.tol_bnd = 1e-9;
    parameters.tol_dj = 1e-9;
    double ms = (deadline - lab_now()) * 800; /* reserve time for native validation */
    if (ms < 1)
        goto done;
    parameters.tm_lim = (int)fmin(INT_MAX, ms);
    size_t remaining = w->options->iteration_limit - w->result->iterations;
    parameters.it_lim = remaining < INT_MAX ? (int)remaining : INT_MAX;
    int status = glp_simplex(problem, &parameters);
    size_t pivots = (size_t)glp_get_it_cnt(problem);
    w->result->glpk_basis_iterations += pivots;
    w->result->iterations += pivots;
    if (status || glp_get_status(problem) != GLP_OPT) {
        w->result->glpk_basis_status = "no_optimal_basis";
    } else if (lab_now() >= deadline) {
        w->result->glpk_basis_status = "deadline";
    } else if (accept_simplex_proposal(w, problem, source_row, rows, deadline)) {
        goto done;
    } else if (rows && !glp_bf_exists(problem) && glp_factorize(problem)) {
        w->result->glpk_basis_status = "factorization_failed";
    } else if (exchange_fixed_rows(w, problem, source_row, rows, index, value, deadline)) {
        w->result->glpk_basis_status = "exchange_failed";
    } else if (lab_now() >= deadline || w->result->iterations >= w->options->iteration_limit) {
        w->result->glpk_basis_status = "budget_exhausted";
    } else {
        lift_suggested_basis(w, problem, source_row, rows, bounds);
    }
done:
    free(bounds);
    free(seen);
    free(value);
    free(index);
    glp_delete_prob(problem);
    free(row_id);
    free(source_row);
}
#endif

#include "bb.h"
#include "common.h"
#include "exact.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Bounds are a persistent path, not one full model copy per pending node.
 * Nodes retain their parent until both child paths have been released. */
typedef struct Node {
    struct Node *parent, *next;
    size_t variable, id, depth, references;
    double lower, upper;
    long double bound;
    int has_bound;
} Node;
typedef struct {
    const LabModel *model;
    const LabBbOptions *options;
    LabBbResult *result;
    double started, sense;
    long double incumbent, constant;
    double *base_lower, *base_upper, *reference, *lower, *upper, *repair_lower, *repair_upper;
    Node *pending;
    size_t next_id;
    FILE *trace;
#ifdef LAB_LP_SESSION_API
    LabLpSession *session;
#endif
} Search;

static void outcome(Search *s, const char *status, const char *message) {
    s->result->status = status;
    snprintf(s->result->message, sizeof s->result->message, "%s", message);
}
static void release(Node *node) {
    while (node && --node->references == 0) {
        Node *parent = node->parent;
        free(node);
        node = parent;
    }
}
static Node *new_node(Search *s, Node *parent, size_t variable, double lo, double hi,
                      long double bound, int has_bound) {
    Node *node = lab_alloc(1, sizeof *node);
    node->parent = parent;
    if (parent)
        parent->references++;
    node->variable = variable;
    node->lower = lo;
    node->upper = hi;
    node->bound = bound;
    node->has_bound = has_bound;
    node->id = s->next_id++;
    node->depth = parent ? parent->depth + 1 : 0;
    node->references = 1;
    return node;
}
static void push(Search *s, Node *node) {
    node->next = s->pending;
    s->pending = node;
}
static void bounds(Search *s, const Node *node) {
    size_t n = s->model->cols;
    memcpy(s->lower, s->base_lower, n * sizeof(double));
    memcpy(s->upper, s->base_upper, n * sizeof(double));
    for (const Node *p = node; p; p = p->parent)
        if (p->variable != SIZE_MAX) {
            size_t j = p->variable;
            s->lower[j] = fmax(s->lower[j], p->lower);
            s->upper[j] = fmin(s->upper[j], p->upper);
        }
}
static double remaining(const Search *s) {
    return s->options->time_limit - (lab_now() - s->started);
}
/* Compare contributions relative to a fixed reference before summation.
 * Constant objective offsets and large fixed-variable contributions must not
 * hide a smaller but material change in the integer assignment. */
static long double score(const Search *s, const double *point) {
    long double value = 0;
    for (size_t j = 0; j < s->model->cols; j++)
        value += (long double)s->model->objective[j] * ((long double)point[j] - s->reference[j]);
    return s->sense * value;
}
static long double difference(const Search *s, const double *left, const double *right) {
    long double value = 0;
    for (size_t j = 0; j < s->model->cols; j++)
        value += (long double)s->model->objective[j] * ((long double)left[j] - right[j]);
    return s->sense * value;
}
static void trace(Search *s, const Node *node, const char *event, size_t variable, double value,
                  int finite) {
    if (!s->trace)
        return;
    fprintf(s->trace, "{\"node\":%zu,\"depth\":%zu,\"event\":", node->id, node->depth);
    lab_json_string(s->trace, event);
    fprintf(s->trace, ",\"variable\":");
    if (variable == SIZE_MAX)
        fputs("null", s->trace);
    else
        lab_json_string(s->trace, s->model->col_name[variable]);
    fprintf(s->trace, ",\"value\":");
    if (finite && isfinite(value))
        fprintf(s->trace, "%.17g", value);
    else
        fputs("null", s->trace);
    fprintf(s->trace, ",\"has_incumbent\":%s}\n", s->result->has_primal ? "true" : "false");
}
static void global_bound(Search *s, const Node *current) {
    long double bound = INFINITY;
    int any = 0, unknown = 0;
    if (current) {
        any = 1;
        if (current->has_bound)
            bound = current->bound;
        else
            unknown = 1;
    }
    for (Node *p = s->pending; p; p = p->next) {
        any = 1;
        if (!p->has_bound)
            unknown = 1;
        else
            bound = fminl(bound, p->bound);
    }
    if (any && !unknown && isfinite(bound)) {
        if (s->result->has_primal)
            bound = fminl(bound, s->incumbent);
        double published = (double)(s->constant + s->sense * bound);
        s->result->has_bound = isfinite(published);
        if (s->result->has_bound)
            s->result->best_bound = published;
    } else {
        s->result->has_bound = 0;
    }
}
static int call_lp(Search *s, const double *lower, const double *upper, LabLpResult *lp) {
    memset(lp, 0, sizeof *lp);
    if (remaining(s) <= 0) {
        outcome(s, "time_limit", "MILP time budget exhausted");
        return -1;
    }
    LabLpOptions opt = s->options->lp;
    opt.time_limit = fmin(opt.time_limit, remaining(s));
    opt.basis_in = opt.basis_out = NULL;
    s->result->lp_solves++;
#ifdef LAB_LP_SESSION_API
    int warm = s->session != NULL;
    if (warm) {
        for (size_t j = 0; j < s->model->cols; j++)
            if (lab_lp_session_set_bounds(s->session, j, lower[j], upper[j])) {
                warm = 0;
                break;
            }
    }
    if (warm) {
        lab_lp_session_solve(s->session, lp);
    } else
#endif
    {
        LabModel view = *s->model;
        view.col_lower = (double *)lower;
        view.col_upper = (double *)upper;
        lab_lp_solve(&view, &opt, lp);
    }
#ifdef LAB_LP_SESSION_API
    if (lp->warm_start)
        s->result->warm_starts++;
    else
#endif
        s->result->cold_starts++;
#ifdef LAB_LP_SESSION_API
    s->result->dual_iterations += lp->dual_iterations;
    s->result->cold_fallbacks += lp->cold_fallbacks;
    s->result->exact_attempts += lp->exact_attempts;
    s->result->exact_infeasibilities += lp->exact_infeasibilities;
    s->result->rejected_pivots += lp->rejected_pivots;
    s->result->harris_rescues += lp->harris_rescues;
    s->result->refinement_steps += lp->refinement_steps;
    if (lp->scaled_rows > s->result->scaled_rows)
        s->result->scaled_rows = lp->scaled_rows;
    if (lp->scaled_columns > s->result->scaled_columns)
        s->result->scaled_columns = lp->scaled_columns;
    s->result->initial_backward_error =
        fmax(s->result->initial_backward_error, lp->initial_backward_error);
    s->result->max_backward_error = fmax(s->result->max_backward_error, lp->max_backward_error);
    if (lp->factor_nnz_peak > s->result->factor_nnz_peak)
        s->result->factor_nnz_peak = lp->factor_nnz_peak;
    if (lp->standard_nnz > s->result->standard_nnz_peak)
        s->result->standard_nnz_peak = lp->standard_nnz;
#endif
    s->result->lp_iterations += lp->iterations;
    s->result->refactorizations += lp->refactorizations;
    if (s->trace && lp->farkas) {
        fprintf(s->trace, "{\"event\":\"lp_infeasibility_certificate\",\"lp_solve\":%zu,"
                         "\"scope\":\"LP with recorded column bounds; not a tree proof\","
                         "\"bounds\":[", s->result->lp_solves);
        mpq_t q;
        mpq_init(q);
        for (size_t j = 0; j < s->model->cols; j++) {
            if (j)
                fputc(',', s->trace);
            fputs("{\"column\":", s->trace);
            lab_json_string(s->trace, s->model->col_name[j]);
            for (int side = 0; side < 2; side++) {
                fputs(side ? ",\"upper\":" : ",\"lower\":", s->trace);
                double v = side ? upper[j] : lower[j];
                if (!isfinite(v))
                    fputs("null", s->trace);
                else {
                    const LabExactModel *e = s->model->exact;
                    if (e && v == (side ? e->col_upper_source[j] : e->col_lower_source[j]))
                        mpq_set(q, side ? e->col_upper[j] : e->col_lower[j]);
                    else
                        mpq_set_d(q, v);
                    char *text = mpq_get_str(NULL, 10, q);
                    lab_json_string(s->trace, text);
                    free(text);
                }
            }
            fputc('}', s->trace);
        }
        mpq_clear(q);
        fputs("],\"multipliers\":{", s->trace);
        for (size_t i = 0; i < lp->farkas_rows; i++) {
            if (i)
                fputc(',', s->trace);
            lab_json_string(s->trace, s->model->row_name[i]);
            fputc(':', s->trace);
            lab_json_string(s->trace, lp->farkas[i]);
        }
        fputs("},\"contradiction_margin\":", s->trace);
        lab_json_string(s->trace, lp->farkas_margin);
        fputs("}\n", s->trace);
    }
    if (remaining(s) <= 0) {
        outcome(s, "time_limit", "MILP time budget exhausted during an LP solve");
        return -1;
    }
    return 0;
}
static int unresolved(Search *s, const LabLpResult *lp) {
    if (!strcmp(lp->status, "optimal") || !strcmp(lp->status, "infeasible"))
        return 0;
    if (!strcmp(lp->status, "unbounded"))
        outcome(s, "lp_unbounded", "An LP relaxation is unbounded; MILP unboundedness is unproved");
    else if (!strcmp(lp->status, "time_limit"))
        outcome(s, "time_limit", "LP solve exhausted the remaining MILP time budget");
    else if (!strcmp(lp->status, "iteration_limit"))
        outcome(s, "iteration_limit", "A node LP exhausted its pivot budget");
    else if (!strcmp(lp->status, "unsupported"))
        outcome(s, "unsupported", lp->message);
    else
        outcome(s, "numerical_failure", lp->message[0] ? lp->message : "Unresolved node LP");
    return -1;
}
static int valid_candidate(Search *s, const double *value, int require_integer, double *objective) {
    const LabModel *m = s->model;
    long double result = m->offset;
    long double *activity = lab_alloc(m->rows, sizeof(long double));
    int valid = 1;
    for (size_t j = 0; j < m->cols; j++) {
        double x = value[j];
        if (!isfinite(x) || x < s->lower[j] - 1e-7 || x > s->upper[j] + 1e-7 ||
            (require_integer && m->integer[j] && x != round(x))) {
            valid = 0;
            break;
        }
        result += (long double)m->objective[j] * x;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; k++)
            activity[m->row_index[k]] += (long double)m->value[k] * x;
    }
    for (size_t i = 0; i < m->rows && valid; i++)
        if (!isfinite(activity[i]) || activity[i] < (long double)m->row_lower[i] - 1e-7L ||
            activity[i] > (long double)m->row_upper[i] + 1e-7L)
            valid = 0;
    free(activity);
    if (!isfinite((double)result))
        valid = 0;
    *objective = (double)result;
    return valid;
}
/* Fix and resolve integer values. Simply rounding the LP vector can violate
 * rows that couple its continuous variables to the integers. */
static int repair(Search *s, const LabLpResult *lp, const Node *node) {
    const LabModel *m = s->model;
    memcpy(s->repair_lower, s->lower, m->cols * sizeof(double));
    memcpy(s->repair_upper, s->upper, m->cols * sizeof(double));
    for (size_t j = 0; j < m->cols; j++)
        if (m->integer[j]) {
            double rounded = round(lp->primal[j]);
            if (rounded < s->lower[j] || rounded > s->upper[j])
                return 0;
            s->repair_lower[j] = s->repair_upper[j] = rounded;
        }
    LabLpResult fixed;
    if (call_lp(s, s->repair_lower, s->repair_upper, &fixed)) {
        lab_lp_result_free(&fixed);
        return -1;
    }
    if (unresolved(s, &fixed)) {
        lab_lp_result_free(&fixed);
        return -1;
    }
    if (!strcmp(fixed.status, "infeasible")) {
        lab_lp_result_free(&fixed);
        return 0;
    }
    if (!fixed.has_primal) {
        outcome(s, "numerical_failure", "Fixed-integer LP returned no primal witness");
        lab_lp_result_free(&fixed);
        return -1;
    }
    for (size_t j = 0; j < m->cols; j++)
        if (m->integer[j])
            fixed.primal[j] = s->repair_lower[j];
    double objective;
    if (!valid_candidate(s, fixed.primal, 1, &objective)) {
        outcome(s, "numerical_failure", "Repaired integer witness fails original-model checks");
        lab_lp_result_free(&fixed);
        return -1;
    }
    long double normalized = score(s, fixed.primal);
    long double delta = s->result->has_primal ? difference(s, fixed.primal, s->result->primal) : 0;
    long double repair_gap = difference(s, fixed.primal, lp->primal);
    if (!isfinite(normalized) || !isfinite(delta) || !isfinite(repair_gap)) {
        outcome(s, "numerical_failure", "Objective comparison overflow");
        lab_lp_result_free(&fixed);
        return -1;
    }
    if (!s->result->has_primal || delta < 0) {
        if (!s->result->primal)
            s->result->primal = lab_alloc(m->cols, sizeof(double));
        memcpy(s->result->primal, fixed.primal, m->cols * sizeof(double));
        s->result->objective = objective;
        s->result->has_primal = 1;
        s->incumbent = normalized;
        trace(s, node, "incumbent", SIZE_MAX, objective, 1);
    }
    /* A feasible rounded assignment need not close a fractional LP node.
     * Test the absolute cost-weighted point difference, never a relative
     * tolerance scaled by constant or fixed objective contributions. */
    int closes_node = fabsl(repair_gap) <= 1e-7L;
    lab_lp_result_free(&fixed);
    return closes_node;
}
static size_t branch_variable(Search *s, const double *primal, double tolerance) {
    for (size_t j = 0; j < s->model->cols; j++)
        if (s->model->integer[j] && fabs(primal[j] - round(primal[j])) > tolerance)
            return j;
    return SIZE_MAX;
}
static int branches(Search *s, Node *node, size_t j, double x) {
    double left = floor(x), right = ceil(x);
    int has_left = left >= s->lower[j] && left < s->upper[j];
    int has_right = right <= s->upper[j] && right > s->lower[j];
    if (!isfinite(x) || fabs(x) > 9007199254740991.0 || left == right ||
        (!has_left && !has_right)) {
        outcome(s, "numerical_failure", "Cannot make a progressing integer branch");
        return -1;
    }
    /* Push right first so the lower branch is explored first. Both children
     * inherit the solved parent's numerical lower bound. */
    if (has_right)
        push(s, new_node(s, node, j, right, s->upper[j], node->bound, node->has_bound));
    if (has_left)
        push(s, new_node(s, node, j, s->lower[j], left, node->bound, node->has_bound));
    trace(s, node, "branch", j, x, 1);
    return 0;
}
void lab_bb_result_free(LabBbResult *result) {
    free(result->primal);
    memset(result, 0, sizeof *result);
}
void lab_bb_solve(const LabModel *model, const LabBbOptions *options, LabBbResult *result) {
    memset(result, 0, sizeof *result);
    Search s = {0};
    s.model = model;
    s.options = options;
    s.result = result;
    s.sense = model->maximize ? -1 : 1;
    s.incumbent = INFINITY;
    s.started = lab_now();
    outcome(&s, "numerical_failure", "MILP search did not complete");
    if (!isfinite(options->time_limit) || !options->node_limit) {
        outcome(&s, "input_error", "Positive finite time and node limits are required");
        goto cleanup;
    }
    if (options->time_limit <= 0) {
        outcome(&s, "time_limit", "MILP time budget exhausted before search");
        goto cleanup;
    }
    if (options->trace_path) {
        s.trace = fopen(options->trace_path, "wx");
        if (!s.trace) {
            outcome(&s, "input_error", "Cannot create fresh search trace");
            goto cleanup;
        }
    }
    size_t n = model->cols;
    s.base_lower = lab_alloc(n, sizeof(double));
    s.base_upper = lab_alloc(n, sizeof(double));
    s.reference = lab_alloc(n, sizeof(double));
    s.lower = lab_alloc(n, sizeof(double));
    s.upper = lab_alloc(n, sizeof(double));
    s.repair_lower = lab_alloc(n, sizeof(double));
    s.repair_upper = lab_alloc(n, sizeof(double));
    s.constant = model->offset;
    for (size_t j = 0; j < n; j++) {
        double lo = model->col_lower[j], hi = model->col_upper[j];
        if (model->integer[j]) {
            if ((isfinite(lo) && fabs(lo) > 9007199254740991.0) ||
                (isfinite(hi) && fabs(hi) > 9007199254740991.0)) {
                outcome(&s, "unsupported", "Integer bounds exceed the exactly spaced double range");
                goto cleanup;
            }
            lo = ceil(lo);
            hi = floor(hi);
        }
        if (isnan(lo) || isnan(hi)) {
            outcome(&s, "input_error", "Nonfinite column-bound encoding");
            goto cleanup;
        }
        if (lo > hi) {
            outcome(&s, "infeasible", "A column has no permitted integer or continuous value");
            goto cleanup;
        }
        s.base_lower[j] = lo;
        s.base_upper[j] = hi;
        s.reference[j] = isfinite(lo) ? lo : isfinite(hi) ? hi : 0;
        s.constant += (long double)model->objective[j] * s.reference[j];
    }
#ifdef LAB_LP_SESSION_API
    if (!options->cold) {
        LabLpOptions opt = options->lp;
        opt.basis_in = opt.basis_out = NULL;
        opt.time_limit = fmin(opt.time_limit, remaining(&s));
        s.session = lab_lp_session_create(model, &opt);
    }
#endif
    push(&s, new_node(&s, NULL, SIZE_MAX, 0, 0, 0, 0));
    while (s.pending) {
        if (remaining(&s) <= 0) {
            global_bound(&s, NULL);
            outcome(&s, "time_limit", "MILP time budget exhausted");
            goto cleanup;
        }
        if (result->nodes >= options->node_limit) {
            global_bound(&s, NULL);
            outcome(&s, "node_limit", "MILP node budget exhausted");
            goto cleanup;
        }
        Node *node = s.pending;
        s.pending = node->next;
        node->next = NULL;
        bounds(&s, node);
        if (result->has_primal && node->has_bound && node->bound >= s.incumbent) {
            result->pruned_nodes++;
            trace(&s, node, "bound_prune", SIZE_MAX, (double)(s.constant + s.sense * node->bound),
                  1);
            release(node);
            continue;
        }
        LabLpResult lp;
        result->nodes++;
        if (call_lp(&s, s.lower, s.upper, &lp)) {
            global_bound(&s, node);
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
        trace(&s, node, lp.status, SIZE_MAX, lp.objective, lp.has_primal);
        if (unresolved(&s, &lp)) {
            global_bound(&s, node);
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
        if (!strcmp(lp.status, "infeasible")) {
            lab_lp_result_free(&lp);
            release(node);
            continue;
        }
        if (!lp.has_primal || !isfinite(lp.objective)) {
            outcome(&s, "numerical_failure", "Node LP returned no finite primal objective");
            global_bound(&s, node);
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
        double recomputed;
        if (!valid_candidate(&s, lp.primal, 0, &recomputed) ||
            fabs(recomputed - lp.objective) > 1e-7 * (1 + fabs(recomputed))) {
            outcome(&s, "numerical_failure", "Node LP witness fails original-model checks");
            global_bound(&s, node);
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
        lp.objective = recomputed;
        int huge_integer = 0;
        for (size_t j = 0; j < model->cols; j++)
            if (model->integer[j] && fabs(lp.primal[j]) > 9007199254740991.0)
                huge_integer = 1;
        if (huge_integer) {
            outcome(&s, "unsupported", "Integer LP values exceed the exactly spaced double range");
            global_bound(&s, node);
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
        long double normalized = score(&s, lp.primal);
        long double lower = normalized - 1e-7L * (1 + fabsl(normalized));
        if (!isfinite(normalized) || !isfinite(lower)) {
            outcome(&s, "numerical_failure", "Node objective comparison overflow");
            global_bound(&s, node);
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
        node->bound = node->has_bound ? fmaxl(node->bound, lower) : lower;
        node->has_bound = 1;
        if (result->has_primal && node->bound >= s.incumbent) {
            result->pruned_nodes++;
            trace(&s, node, "bound_prune", SIZE_MAX, (double)(s.constant + s.sense * node->bound),
                  1);
            lab_lp_result_free(&lp);
            release(node);
            continue;
        }
        size_t variable = branch_variable(&s, lp.primal, 1e-7);
        if (variable == SIZE_MAX) {
            int repaired = repair(&s, &lp, node);
            if (repaired < 0) {
                global_bound(&s, node);
                lab_lp_result_free(&lp);
                release(node);
                goto cleanup;
            }
            if (repaired > 0) {
                lab_lp_result_free(&lp);
                release(node);
                continue;
            }
            /* A failed rounded assignment does not establish this node's
             * infeasibility. Branch on any remaining exact fraction or stop. */
            variable = branch_variable(&s, lp.primal, 0);
            if (variable == SIZE_MAX) {
                outcome(&s, "numerical_failure", "Integral node failed fixed-integer repair");
                global_bound(&s, node);
                lab_lp_result_free(&lp);
                release(node);
                goto cleanup;
            }
        }
        if (branches(&s, node, variable, lp.primal[variable])) {
            global_bound(&s, node);
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
        lab_lp_result_free(&lp);
        release(node);
    }
    if (result->has_primal) {
        outcome(&s, "optimal",
                "Numerical branch-and-bound tree exhausted; no independent tree proof");
        result->has_bound = 1;
        result->best_bound = result->objective;
    } else {
        outcome(&s, "infeasible", "Numerical branch-and-bound tree exhausted without an incumbent");
    }
cleanup:
    while (s.pending) {
        Node *next = s.pending->next;
        release(s.pending);
        s.pending = next;
    }
#ifdef LAB_LP_SESSION_API
    lab_lp_session_free(s.session);
#endif
    if (s.trace) {
        int failed = ferror(s.trace);
        if (fclose(s.trace))
            failed = 1;
        if (failed)
            outcome(&s, "input_error", "Cannot write search trace");
    }
    free(s.base_lower);
    free(s.base_upper);
    free(s.reference);
    free(s.lower);
    free(s.upper);
    free(s.repair_lower);
    free(s.repair_upper);
    result->elapsed = lab_now() - s.started;
}

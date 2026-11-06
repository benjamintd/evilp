#include "bb.h"
#include "common.h"
#include "exact.h"
#include "presolve.h"
#include "search.h"
#include "augment.h"
#include "root_cut.h"
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
    long double parent_lp_score;
    double branch_distance;
    LabSearchDirection direction;
} Node;
typedef struct {
    const LabModel *model;
    LabModel lp_model;
    LabModel augmented;
    LabRootCutBatch cuts;
    double *lp_row_lower, *lp_row_upper;
    LabSearch branching;
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
static void current_bound(const Search *s, size_t j, double value, int upper, mpq_t q) {
    const LabExactModel *e = s->model->exact;
    if (e && value == (upper ? e->col_upper_source[j] : e->col_lower_source[j]))
        mpq_set(q, upper ? e->col_upper[j] : e->col_lower[j]);
    else
        mpq_set_d(q, value);
}
/* Every active cut is conditional on the exported root box. Verify this
 * premise in exact source units for node and fixed-integer repair LPs. */
static int cut_box_contains(const Search *s, const LabRootCutBatch *cuts,
                            const double *lower, const double *upper) {
    if (cuts->columns != s->model->cols || (cuts->columns && !cuts->box))
        return 0;
    mpq_t q;
    mpq_init(q);
    int valid = 1;
    for (size_t j = 0; valid && j < cuts->columns; j++) {
        const LabNativeBox *b = &cuts->box[j];
        if (b->has_lower) {
            if (!isfinite(lower[j]))
                valid = 0;
            else {
                current_bound(s, j, lower[j], 0, q);
                valid = mpq_cmp(q, b->lower) >= 0;
            }
        }
        if (valid && b->has_upper) {
            if (!isfinite(upper[j]))
                valid = 0;
            else {
                current_bound(s, j, upper[j], 1, q);
                valid = mpq_cmp(q, b->upper) <= 0;
            }
        }
    }
    mpq_clear(q);
    return valid;
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
    if (s->cuts.count && !cut_box_contains(s, &s->cuts, lower, upper)) {
        outcome(s, "numerical_failure", "LP bounds leave the certified root-cut box");
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
        LabModel view = s->lp_model;
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
    s->result->source_infeasibility_attempts += lp->source_infeasibility_attempts;
    s->result->source_infeasibilities += lp->source_infeasibilities;
    s->result->source_projected_rows += lp->source_projected_rows;
    s->result->source_removed_terms += lp->source_removed_terms;
    s->result->source_reconstructed_terms += lp->source_reconstructed_terms;
    s->result->rejected_pivots += lp->rejected_pivots;
    s->result->harris_rescues += lp->harris_rescues;
    s->result->primal_accepted_pivots += lp->primal_accepted_pivots;
    s->result->primal_zero_steps += lp->primal_zero_steps;
    s->result->primal_short_steps += lp->primal_short_steps;
    s->result->pricing_bland_fallbacks += lp->pricing_bland_fallbacks;
    s->result->primal_candidate_rejections += lp->primal_candidate_rejections;
    s->result->primal_repairs += lp->primal_repairs;
    s->result->primal_repair_iterations += lp->primal_repair_iterations;
    s->result->exact_primal_attempts += lp->exact_primal_attempts;
    s->result->exact_primal_recoveries += lp->exact_primal_recoveries;
    s->result->basis_updates += lp->basis_updates;
    s->result->basis_update_rejections += lp->basis_update_rejections;
    if (lp->factor_core_peak > s->result->factor_core_peak)
        s->result->factor_core_peak = lp->factor_core_peak;
    if (lp->singleton_eliminations_peak > s->result->singleton_eliminations_peak)
        s->result->singleton_eliminations_peak = lp->singleton_eliminations_peak;
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
        fprintf(s->trace, "{\"event\":\"%s\",\"lp_solve\":%zu,\"scope\":",
                s->cuts.count ? "mip_in_box_infeasibility_certificate"
                              : "lp_infeasibility_certificate", s->result->lp_solves);
        lab_json_string(s->trace, s->cuts.count
            ? "Original MIP with recorded bounds and GMI cuts; not original LP or a tree proof"
            : "LP with recorded column bounds; not a tree proof");
        if (s->cuts.count) {
            fputs(",\"cut_ids\":[", s->trace);
            for (size_t r = 0; r < s->cuts.count; r++) {
                if (r)
                    fputc(',', s->trace);
                lab_json_string(s->trace, s->lp_model.row_name[s->model->rows + r]);
            }
            fputc(']', s->trace);
        }
        fputs(",\"bounds\":[", s->trace);
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
            lab_json_string(s->trace, s->lp_model.row_name[i]);
            fputc(':', s->trace);
            lab_json_string(s->trace, lp->farkas[i]);
        }
        fputs(lp->farkas_margin_is_lower_bound
                  ? "},\"margin_kind\":\"source_enclosure\",\"contradiction_margin_lower_bound\":"
                  : "},\"contradiction_margin\":", s->trace);
        lab_json_string(s->trace, lp->farkas_margin);
        if (lp->farkas_margin_is_lower_bound)
            fprintf(s->trace, ",\"projected_rows\":%zu,\"removed_terms\":%zu,\"reconstructed_terms\":%zu",
                    lp->source_projected_rows, lp->source_removed_terms, lp->source_reconstructed_terms);
        fputs("}\n", s->trace);
    }
    if (remaining(s) <= 0) {
        outcome(s, "time_limit", "MILP time budget exhausted during an LP solve");
        return -1;
    }
    return 0;
}
/* A single bounded root round. New rows and their session are constructed
 * before changing the live search view; the original candidate model never
 * changes. Failure to derive/allocate usable cuts leaves the old LP intact. */
static int root_cuts(Search *s) {
    s->result->root_cut_attempts++;
    if (!s->session) {
        s->result->root_cut_status = "session_unavailable";
        return 0;
    }
    LabRootCutBatch candidate;
    lab_root_cut_batch_init(&candidate);
    double budget = fmin(1, .05 * remaining(s));
    LabRootCutLimits cut_limits = lab_root_cut_limits(lab_now() + budget);
    LabRootCutStatus status = lab_root_cut_derive(s->session, s->model, &cut_limits, &candidate);
    s->result->root_cut_status = lab_root_cut_status_name(status);
    if (remaining(s) <= 0) {
        lab_root_cut_batch_free(&candidate);
        outcome(s, "time_limit", "MILP time budget exhausted deriving root cuts");
        return -1;
    }
    if (status != LAB_ROOT_CUT_OK || !candidate.count) {
        lab_root_cut_batch_free(&candidate);
        return 0;
    }
    if (!cut_box_contains(s, &candidate, s->lower, s->upper)) {
        s->result->root_cut_status = "box_mismatch";
        lab_root_cut_batch_free(&candidate);
        return 0;
    }
    LabAddedRow *rows = lab_alloc(candidate.count, sizeof *rows);
    char **names = lab_alloc(candidate.count, sizeof *names);
    size_t serial = 0;
    for (size_t r = 0; r < candidate.count; r++) {
        char name[80];
        int collision;
        do {
            snprintf(name, sizeof name, "__milplab_gmi_%zu", serial++);
            collision = 0;
            for (size_t i = 0; i < s->lp_model.rows; i++)
                if (!strcmp(name, s->lp_model.row_name[i]))
                    collision = 1;
        } while (collision);
        names[r] = lab_strdup(name);
        rows[r] = (LabAddedRow){names[r], candidate.cut[r].coefficient,
                                candidate.cut[r].lower, INFINITY};
    }
    LabAugmentLimits limits = {.max_rows = 8000000, .max_nonzeros = 128000000,
                               .max_added_rows = 4};
    LabAugmentStatus augmented = lab_model_append_rows(&s->lp_model, candidate.count,
                                                       rows, &limits, &s->augmented);
    for (size_t r = 0; r < candidate.count; r++)
        free(names[r]);
    free(names);
    free(rows);
    if (augmented != LAB_AUGMENT_OK) {
        s->result->root_cut_status = augmented == LAB_AUGMENT_RESOURCE
            ? "augmentation_resource_limit" : "augmentation_invalid";
        lab_root_cut_batch_free(&candidate);
        return 0;
    }
    LabLpSession *next = NULL;
    if (!s->options->cold) {
        LabLpOptions options = s->options->lp;
        options.basis_in = options.basis_out = NULL;
        options.time_limit = fmin(options.time_limit, remaining(s));
        next = lab_lp_session_create(&s->augmented, &options);
        if (!next) {
            s->result->root_cut_status = "augmented_session_unavailable";
            lab_model_free(&s->augmented);
            lab_root_cut_batch_free(&candidate);
            return 0;
        }
    }
    lab_lp_session_free(s->session);
    s->session = next;
    s->lp_model = s->augmented;
    s->cuts = candidate;
    s->result->root_cuts = candidate.count;
    s->result->root_cut_status = "applied";
    if (s->trace)
        for (size_t r = 0; r < candidate.count; r++) {
            fputs("{\"event\":\"root_gmi_cut\",\"cut_id\":", s->trace);
            lab_json_string(s->trace, s->lp_model.row_name[s->model->rows + r]);
            fputs(",\"scope\":\"ORIGINAL_MIP_IN_BOX\",\"certificate\":", s->trace);
            lab_json_string(s->trace, candidate.cut[r].certificate);
            fputs("}\n", s->trace);
        }
    return 1;
}
static int unresolved(Search *s, const LabLpResult *lp) {
    if (!strcmp(lp->status, "optimal") || !strcmp(lp->status, "infeasible"))
        return 0;
    if (lp->primal)
        for (size_t j = 0; j < s->model->cols; j++)
            if (s->model->integer[j] && isfinite(lp->primal[j]) &&
                fabs(lp->primal[j]) > 9007199254740991.0) {
                outcome(s, "unsupported", "Integer values exceed the exactly spaced double range");
                return -1;
            }
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
    for (size_t j = 0; j < m->cols; j++) {
        double x = value[j];
        if (!isfinite(x) || x < s->lower[j] - 1e-7 || x > s->upper[j] + 1e-7 ||
            (require_integer && m->integer[j] && x != round(x)))
            return 0;
    }
    LabModel view = *m;
    view.col_lower = s->lower;
    view.col_upper = s->upper;
    double residual;
    return !lab_exact_evaluate(&view, value, objective, &residual) && residual <= 1e-7;
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
static LabSearchChoice branch_variable(Search *s, const double *primal, double tolerance,
                                       size_t *variable) {
    if (!s->options->first_fractional)
        return lab_search_choose(s->model, primal, s->lower, s->upper,
                                 &s->branching, tolerance, variable);
    *variable = SIZE_MAX;
    for (size_t j = 0; j < s->model->cols; j++)
        if (s->model->integer[j] && fabs(primal[j] - round(primal[j])) > tolerance) {
            *variable = j;
            return LAB_SEARCH_CHOSEN;
        }
    return LAB_SEARCH_INTEGRAL;
}
static int branches(Search *s, Node *node, size_t j, double x, long double lp_score) {
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
    if (has_right) {
        Node *child = new_node(s, node, j, right, s->upper[j], node->bound, node->has_bound);
        child->parent_lp_score = lp_score;
        child->branch_distance = right - x;
        child->direction = LAB_SEARCH_UP;
        push(s, child);
    }
    if (has_left) {
        Node *child = new_node(s, node, j, s->lower[j], left, node->bound, node->has_bound);
        child->parent_lp_score = lp_score;
        child->branch_distance = x - left;
        child->direction = LAB_SEARCH_DOWN;
        push(s, child);
    }
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
    s.lp_model = *model;
    s.options = options;
    s.result = result;
    s.sense = model->maximize ? -1 : 1;
    s.incumbent = INFINITY;
    s.started = lab_now();
    result->root_cut_status = options->no_root_cuts ? "disabled" : "not_attempted";
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
    }
    if (!options->no_presolve) {
        double available = remaining(&s);
        if (available <= 0) {
            outcome(&s, "time_limit", "MILP time budget exhausted before presolve");
            goto cleanup;
        }
        LabPresolveOptions opt = {.time_limit = fmin(2, .05 * available),
                                  .max_passes = 32, .relax = 0};
        LabPresolveResult presolved;
        lab_presolve(model, s.base_lower, s.base_upper, &opt, &presolved);
        result->presolve_passes = presolved.passes;
        result->presolve_tightenings = presolved.tightenings;
        result->presolve_fixed = presolved.fixed;
        result->presolve_rows = presolved.redundant_rows;
        if (remaining(&s) <= 0) {
            outcome(&s, "time_limit", "MILP time budget exhausted during root presolve");
            lab_presolve_result_free(&presolved);
            goto cleanup;
        }
        if (!strcmp(presolved.status, "input_error") || !strcmp(presolved.status, "infeasible")) {
            outcome(&s, presolved.status, presolved.message);
            lab_presolve_result_free(&presolved);
            goto cleanup;
        }
        memcpy(s.base_lower, presolved.lower, n * sizeof(double));
        memcpy(s.base_upper, presolved.upper, n * sizeof(double));
        if (presolved.redundant_rows) {
            s.lp_row_lower = lab_alloc(model->rows, sizeof(double));
            s.lp_row_upper = lab_alloc(model->rows, sizeof(double));
            memcpy(s.lp_row_lower, model->row_lower, model->rows * sizeof(double));
            memcpy(s.lp_row_upper, model->row_upper, model->rows * sizeof(double));
            for (size_t i = 0; i < model->rows; i++)
                if (presolved.redundant[i]) {
                    s.lp_row_lower[i] = -INFINITY;
                    s.lp_row_upper[i] = INFINITY;
                }
            s.lp_model.row_lower = s.lp_row_lower;
            s.lp_model.row_upper = s.lp_row_upper;
        }
        lab_presolve_result_free(&presolved);
    }
    s.lp_model.col_lower = s.base_lower;
    s.lp_model.col_upper = s.base_upper;
    s.constant = model->offset;
    for (size_t j = 0; j < n; j++) {
        double lo = s.base_lower[j], hi = s.base_upper[j];
        if (model->integer[j] &&
            ((isfinite(lo) && fabs(lo) > 9007199254740991.0) ||
             (isfinite(hi) && fabs(hi) > 9007199254740991.0))) {
            outcome(&s, "unsupported", "Integer bounds exceed the exactly spaced double range");
            goto cleanup;
        }
        s.reference[j] = isfinite(lo) ? lo : isfinite(hi) ? hi : 0;
        s.constant += (long double)model->objective[j] * s.reference[j];
    }
    if (!options->first_fractional)
        (void)lab_search_init(&s.branching, n); /* allocation failure keeps deterministic fallback */
#ifdef LAB_LP_SESSION_API
    if (!options->cold || !options->no_root_cuts) {
        LabLpOptions opt = options->lp;
        opt.basis_in = opt.basis_out = NULL;
        opt.time_limit = fmin(opt.time_limit, remaining(&s));
        s.session = lab_lp_session_create(&s.lp_model, &opt);
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
node_lp_ready:
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
        if (!node->parent && !options->no_root_cuts && !result->root_cut_attempts) {
            int fractional = 0;
            for (size_t j = 0; j < model->cols; j++)
                if (model->integer[j] && fabs(lp.primal[j] - round(lp.primal[j])) > 1e-7)
                    fractional = 1;
            int added = fractional ? root_cuts(&s) : 0;
            if (added < 0) {
                global_bound(&s, node);
                lab_lp_result_free(&lp);
                release(node);
                goto cleanup;
            }
            if (options->cold) {
                lab_lp_session_free(s.session);
                s.session = NULL;
            }
            if (added) {
                lab_lp_result_free(&lp);
                result->root_cut_lp_solves++;
                if (call_lp(&s, s.lower, s.upper, &lp)) {
                    global_bound(&s, node);
                    lab_lp_result_free(&lp);
                    release(node);
                    goto cleanup;
                }
                goto node_lp_ready;
            }
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
        if (!options->first_fractional && node->parent && node->variable != SIZE_MAX &&
            lab_search_observe(&s.branching, node->variable, node->direction,
                               node->branch_distance, normalized - node->parent_lp_score) == 1)
            result->branch_observations++;
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
        size_t variable = SIZE_MAX;
        LabSearchChoice selection = branch_variable(&s, lp.primal, 1e-7, &variable);
        if (selection < 0) {
            outcome(&s, selection == LAB_SEARCH_UNSUPPORTED ? "unsupported" : "numerical_failure",
                    "Cannot select a progressing integer branch");
            global_bound(&s, node);
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
        if (selection == LAB_SEARCH_INTEGRAL) {
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
            selection = branch_variable(&s, lp.primal, 0, &variable);
            if (selection != LAB_SEARCH_CHOSEN) {
                outcome(&s, selection == LAB_SEARCH_UNSUPPORTED ? "unsupported" : "numerical_failure",
                        "Integral node failed fixed-integer repair");
                global_bound(&s, node);
                lab_lp_result_free(&lp);
                release(node);
                goto cleanup;
            }
        }
        if (branches(&s, node, variable, lp.primal[variable], normalized)) {
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
    lab_model_free(&s.augmented);
    lab_root_cut_batch_free(&s.cuts);
    lab_search_free(&s.branching);
    if (s.trace) {
        int failed = ferror(s.trace);
        if (fclose(s.trace))
            failed = 1;
        if (failed)
            outcome(&s, "input_error", "Cannot write search trace");
    }
    free(s.base_lower);
    free(s.base_upper);
    free(s.lp_row_lower);
    free(s.lp_row_upper);
    free(s.reference);
    free(s.lower);
    free(s.upper);
    free(s.repair_lower);
    free(s.repair_upper);
    result->elapsed = lab_now() - s.started;
}

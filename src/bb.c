#include "bb.h"
#include "binary_search.h"
#include "assignment.h"
#include "pair_assignment.h"
#include "binary_sat.h"
#include "flow.h"
#include "lot_sizing.h"
#ifdef LAB_SCIP_INTEGER
#include "integer_mip.h"
#include "scip_mip.h"
#endif
#ifdef LAB_GLPK_BASIS
#include "activation_cover.h"
#include "pooled_assignment.h"
#endif
#include "packing.h"
#include "partition.h"
#include "common.h"
#include "exact.h"
#include "dependency.h"
#include "presolve.h"
#include "search.h"
#include "augment.h"
#include "root_cut.h"
#include "dive.h"
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
    LabRootCutBatch cuts;
    double *lp_row_lower, *lp_row_upper;
    LabSearchBound frontier;
    LabSearch branching;
    const LabBbOptions *options;
    LabBbResult *result;
    double started, sense;
    long double incumbent, constant;
    double objective_step;
    int checked_objective_step;
    double *base_lower, *base_upper, *reference, *lower, *upper, *repair_lower, *repair_upper;
    Node **heap;
    size_t heap_size, heap_cap;
    Node *plunge;
    size_t plunge_steps;
    size_t next_id;
    FILE *trace;
#ifdef LAB_LP_SESSION_API
    LabLpSession *session;
    int session_initialized;
#endif
} Search;

/* The search owns the primal storage. Only terminal observations may publish
 * a global bound, and all public bounds are in original objective units. */
static void publish(Search *s, LabSearchResult verdict, const char *message) {
    LabBbResult *r = s->result;
    verdict.point = r->has_primal ? r->primal : NULL;
    r->search = verdict;
    r->status = lab_mip_status_name(verdict.status);
    r->has_bound = verdict.bound.origin != LAB_BOUND_NONE;
    r->best_bound = verdict.bound.value;
    if (message)
        snprintf(r->message, sizeof r->message, "%s", message);
}

static void outcome(Search *s, LabMipStatus status, const char *message) {
    publish(s, (LabSearchResult){.route = "native", .status = status}, message);
}

static void close_at_incumbent(Search *s, const char *route, LabClosure closure,
                               LabBoundOrigin origin, const char *message) {
    publish(s, (LabSearchResult){.route = route, .status = LAB_MIP_OPTIMAL, .closure = closure,
                                 .bound = {origin, s->result->objective}}, message);
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
static int node_better(const Node *a, const Node *b) {
    if (a->has_bound != b->has_bound)
        return !a->has_bound;
    if (a->has_bound && a->bound != b->bound)
        return a->bound < b->bound;
    if (a->depth != b->depth)
        return a->depth > b->depth;
    return a->id < b->id;
}
static void push(Search *s, Node *node) {
    if (s->heap_size == s->heap_cap) {
        size_t next_cap = s->heap_cap ? s->heap_cap * 2 : 64;
        Node **next = lab_alloc(next_cap, sizeof *next);
        if (s->heap_size)
            memcpy(next, s->heap, s->heap_size * sizeof *next);
        free(s->heap);
        s->heap = next;
        s->heap_cap = next_cap;
    }
    size_t i = s->heap_size++;
    while (i > 0) {
        size_t p = (i - 1) / 2;
        if (!node_better(node, s->heap[p]))
            break;
        s->heap[i] = s->heap[p];
        i = p;
    }
    s->heap[i] = node;
}
static Node *pop_heap(Search *s) {
    if (!s->heap_size)
        return NULL;
    Node *top = s->heap[0];
    Node *last = s->heap[--s->heap_size];
    if (s->heap_size) {
        size_t i = 0;
        for (;;) {
            size_t l = 2 * i + 1, r = l + 1, best = i;
            if (l < s->heap_size && node_better(s->heap[l], last))
                best = l;
            if (r < s->heap_size &&
                node_better(s->heap[r], best == i ? last : s->heap[best]))
                best = r;
            if (best == i)
                break;
            s->heap[i] = s->heap[best];
            i = best;
        }
        s->heap[i] = last;
    }
    return top;
}
static Node *next_node(Search *s) {
    if (s->plunge) {
        if (s->options->first_fractional || s->plunge_steps < 16) {
            Node *n = s->plunge;
            s->plunge = NULL;
            s->plunge_steps++;
            return n;
        }
        push(s, s->plunge);
        s->plunge = NULL;
    }
    s->plunge_steps = 0;
    return pop_heap(s);
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
/* If only integer variables carry a cost, any two integer solutions differ
 * by a multiple of gcd(numerators) / lcm(denominators). Offsets cancel.
 * Work from the source rationals, including coefficients that round to zero. */
static double objective_step(const LabModel *m, double deadline) {
    mpq_t cost, step;
    mpq_inits(cost, step, NULL);
    double result = 0;
    for (size_t j = 0; j < m->cols; ++j) {
        if (((j & 255U) == 0 && lab_now() >= deadline) ||
            lab_source_exact_objective(m, j, cost))
            goto done;
        if (!mpq_sgn(cost))
            continue;
        if (!m->integer[j])
            goto done;
        mpz_gcd(mpq_numref(step), mpq_numref(step), mpq_numref(cost));
        mpz_lcm(mpq_denref(step), mpq_denref(step), mpq_denref(cost));
        if (mpz_sizeinbase(mpq_denref(step), 2) > 256)
            goto done;
    }
    /* GMP truncates toward zero, giving a conservative positive step. */
    result = mpq_get_d(step);
done:
    mpq_clears(cost, step, NULL);
    return isfinite(result) ? result : 0;
}
static int cannot_improve(Search *s, long double lower) {
    if (!s->result->has_primal)
        return 0;
    if (lower >= s->incumbent)
        return 1;
    if (!s->checked_objective_step) {
        s->checked_objective_step = 1;
        if (!s->options->lp.no_exact && !s->options->no_objective_lattice)
            s->objective_step = objective_step(s->model, lab_now() + fmin(.1, remaining(s)));
    }
    if (s->objective_step <= 0)
        return 0;
    /* The normalized incumbent is a source sum converted to double. Bound
     * it upward, then round the subtraction upward too. Equality leaves a
     * possible one-step improvement and must keep the node open. */
    long double upper = nextafter((double)s->incumbent, INFINITY);
    if (lower <= nextafterl(upper - s->objective_step, INFINITY))
        return 0;
    s->result->objective_lattice_prunes++;
    return 1;
}
static void current_bound(const Search *s,size_t j,double value,int upper,mpq_t q) {
    if (lab_source_exact_column_value(s->model,j,upper,value,q))
        mpq_set_d(q,value);
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
    double value;
    if (lab_source_objective_difference(s->model, point, s->reference, &value))
        return NAN;
    return s->sense * (long double)value;
}
static long double difference(const Search *s, const double *left, const double *right) {
    double value;
    if (lab_source_objective_difference(s->model, left, right, &value))
        return NAN;
    return s->sense * (long double)value;
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
    const Node *extra[2] = {current, s->plunge};
    for (size_t k = 0; k < 2; k++) {
        if (!extra[k])
            continue;
        any = 1;
        if (extra[k]->has_bound)
            bound = fminl(bound, extra[k]->bound);
        else
            unknown = 1;
    }
    if (s->heap_size) {
        any = 1;
        /* Min-heap orders !has_bound first, then smallest bound. */
        if (!s->heap[0]->has_bound)
            unknown = 1;
        else
            bound = fminl(bound, s->heap[0]->bound);
    }
    if (any && !unknown && isfinite(bound)) {
        if (s->result->has_primal)
            bound = fminl(bound, s->incumbent);
        double published = (double)(s->constant + s->sense * bound);
        s->frontier = isfinite(published) ? (LabSearchBound){LAB_BOUND_FRONTIER, published}
                                          : (LabSearchBound){0};
    } else {
        s->frontier = (LabSearchBound){0};
    }
}
static int call_lp(Search *s, const double *lower, const double *upper, LabLpResult *lp) {
    memset(lp, 0, sizeof *lp);
    if (remaining(s) <= 0) {
        outcome(s, LAB_MIP_TIME_LIMIT, "MILP time budget exhausted");
        return -1;
    }
    if (s->cuts.count && !cut_box_contains(s, &s->cuts, lower, upper)) {
        outcome(s, LAB_MIP_NUMERICAL_FAILURE, "LP bounds leave the certified root-cut box");
        return -1;
    }
    LabLpOptions opt = s->options->lp;
    opt.time_limit = fmin(opt.time_limit, remaining(s));
    opt.basis_in = opt.basis_out = NULL;
    s->result->lp_solves++;
#ifdef LAB_LP_SESSION_API
    if (!s->session_initialized && (!s->options->cold || !s->options->no_root_cuts)) {
        s->session_initialized = 1;
        s->session = lab_lp_session_create(&s->lp_model, &opt);
    }
    int warm = s->session && (!s->options->cold || !s->result->root_cut_attempts);
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
    s->result->glpk_basis_iterations += lp->glpk_basis_iterations;
    s->result->glpk_bases_accepted += lp->glpk_bases_accepted;
    s->result->primal_accepted_pivots += lp->primal_accepted_pivots;
    s->result->primal_zero_steps += lp->primal_zero_steps;
    s->result->primal_short_steps += lp->primal_short_steps;
    s->result->pricing_bland_fallbacks += lp->pricing_bland_fallbacks;
    s->result->primal_candidate_rejections += lp->primal_candidate_rejections;
    s->result->source_pricing_attempts += lp->source_pricing_attempts;
    s->result->source_pricing_pivots += lp->source_pricing_pivots;
    s->result->source_pricing_completions += lp->source_pricing_completions;
    s->result->source_dual_closures += lp->source_dual_closures;
    s->result->source_support_scans += lp->source_support_scans;
    s->result->source_support_binary_paths += lp->source_support_binary_paths;
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
                    current_bound(s,j,v,side,q);
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
        outcome(s, LAB_MIP_TIME_LIMIT, "MILP time budget exhausted during an LP solve");
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
        outcome(s, LAB_MIP_TIME_LIMIT, "MILP time budget exhausted deriving root cuts");
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
    LabAugmentStatus augmented = lab_lp_session_append_dense_rows(s->session, candidate.count, rows, &limits);
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
    s->lp_model = *lab_lp_session_domain(s->session);
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
                outcome(s, LAB_MIP_UNSUPPORTED, "Integer values exceed the exactly spaced double range");
                return -1;
            }
    if (!strcmp(lp->status, "unbounded"))
        outcome(s, LAB_MIP_LP_UNBOUNDED, "An LP relaxation is unbounded; MILP unboundedness is unproved");
    else if (!strcmp(lp->status, "time_limit"))
        outcome(s, LAB_MIP_TIME_LIMIT, "LP solve exhausted the remaining MILP time budget");
    else if (!strcmp(lp->status, "iteration_limit"))
        outcome(s, LAB_MIP_ITERATION_LIMIT, "A node LP exhausted its pivot budget");
    else if (!strcmp(lp->status, "unsupported"))
        outcome(s, LAB_MIP_UNSUPPORTED, lp->message);
    else
        outcome(s, LAB_MIP_NUMERICAL_FAILURE, lp->message[0] ? lp->message : "Unresolved node LP");
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
#if defined(LAB_GLPK_BASIS) || defined(LAB_SCIP_INTEGER)
static void retain_search_incumbent(Search *s, const Node *node) {
    double objective;
    if (!valid_candidate(s, s->repair_lower, 1, &objective) ||
        (s->result->has_primal && difference(s, s->repair_lower, s->result->primal) >= 0))
        return;
    long double value = score(s, s->repair_lower);
    if (!isfinite(value))
        return;
    size_t n = s->model->cols;
    if (!s->result->primal)
        s->result->primal = lab_alloc(n, sizeof(double));
    memcpy(s->result->primal, s->repair_lower, n * sizeof(double));
    s->result->has_primal = 1;
    s->result->objective = objective;
    s->incumbent = value;
    trace(s, node, "incumbent", SIZE_MAX, objective, 1);
}
#endif
static int store_terminal_point(Search *s) {
    double objective;
    if (!valid_candidate(s, s->repair_lower, 1, &objective) || remaining(s) <= 0) return 0;
    size_t n = s->model->cols;
    if (!s->result->primal)
        s->result->primal = lab_alloc(n, sizeof(double));
    memcpy(s->result->primal, s->repair_lower, n * sizeof(double));
    s->result->has_primal = 1;
    s->result->objective = objective;
    return 1;
}
static int complete_point_search(Search *s, const Node *node, const char *route, LabClosure closure,
                                 const char *message, const char *event) {
    if (!store_terminal_point(s)) return 0;
    close_at_incumbent(s, route, closure,
                       closure == LAB_CLOSE_EXHAUSTED ? LAB_BOUND_EXHAUSTION : LAB_BOUND_STRUCTURAL, message);
    trace(s, node, event, SIZE_MAX, s->result->objective, 1);
    return 1;
}
#define LAB_PHASE(sp, name) do { if (getenv("LAB_PHASE_TRACE")) \
    fprintf(stderr, "phase %s remaining=%.3f\n", (name), remaining(sp)); } while (0)
#ifdef LAB_SCIP_INTEGER
/* Returns 1 if search terminated (optimal/infeasible/time_limit/node_limit), 0 to continue. */
static int run_scip_mip(Search *s, const Node *node) {
    const LabModel *model = s->model;
    const LabBbOptions *options = s->options;
    LabBbResult *result = s->result;
    if (model->cols < 32 ||
        remaining(s) <= 0.05 || result->nodes >= options->node_limit)
        return 0;
    result->scip_route = "general";
    trace(s, node, "scip_mip_start", SIZE_MAX, 0, 0);
    LabScipMipStats sstats;
    LabSearchResult attempt = lab_scip_mip(model, s->base_lower, s->base_upper,
                              s->lp_row_lower, s->lp_row_upper,
                              result->has_primal ? result->primal : NULL,
                              lab_now() + remaining(s),
                              options->node_limit - result->nodes,
                              options->lp.iteration_limit, s->repair_lower, &sstats);
    result->scip_path_cliques += sstats.path_cliques;
    result->nodes += sstats.nodes;
    result->lp_solves += sstats.lp_solves;
    result->lp_iterations += sstats.lp_iterations;
    if (attempt.status == LAB_MIP_INFEASIBLE && remaining(s) > 0) {
        publish(s, attempt, "SCIP branch-and-cut reported infeasibility; no independent tree proof");
        return 1;
    }
    const char *message;
    if (attempt.status == LAB_MIP_OPTIMAL && store_terminal_point(s)) {
        message = "SCIP branch-and-cut reported optimality; original witness checked";
        trace(s, node, "scip_mip_optimum", SIZE_MAX, result->objective, 1);
    } else {
        if (attempt.point)
            retain_search_incumbent(s, node);
        attempt.closure = LAB_CLOSE_NONE;
        if (attempt.status == LAB_MIP_GAP_LIMIT) {
            message = "SCIP reached its configured MIP gap limit";
        } else if (attempt.status == LAB_MIP_TIME_LIMIT || remaining(s) <= 0.05) {
            attempt.status = LAB_MIP_TIME_LIMIT;
            message = "MILP time budget exhausted during branch-and-cut";
        } else if (attempt.status == LAB_MIP_NODE_LIMIT || (sstats.nodes && result->nodes >= options->node_limit)) {
            attempt.status = LAB_MIP_NODE_LIMIT;
            message = "MILP node budget exhausted during branch-and-cut";
        } else {
            result->scip_route = "general_fallback";
            result->scip_fallback_reason = sstats.fallback_reason ? sstats.fallback_reason :
                attempt.point ? "incumbent_only" : "no_verified_completion";
            trace(s, node, "scip_mip_fallback", SIZE_MAX, 0, 0);
            return 0;
        }
    }
    publish(s, attempt, message);
    return 1;
}
#endif
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
    int probe = s->session && !s->options->cold;
    if (probe && lab_lp_session_begin_probe(s->session)) {
        outcome(s, LAB_MIP_NUMERICAL_FAILURE, "Could not start fixed-integer LP probe");
        return -1;
    }
    int failed = call_lp(s, s->repair_lower, s->repair_upper, &fixed);
    if (probe)
        lab_lp_session_rollback_probe(s->session);
    if (failed) {
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
        outcome(s, LAB_MIP_NUMERICAL_FAILURE, "Fixed-integer LP returned no primal witness");
        lab_lp_result_free(&fixed);
        return -1;
    }
    for (size_t j = 0; j < m->cols; j++)
        if (m->integer[j])
            fixed.primal[j] = s->repair_lower[j];
    double objective;
    if (!valid_candidate(s, fixed.primal, 1, &objective)) {
        outcome(s, LAB_MIP_NUMERICAL_FAILURE, "Repaired integer witness fails original-model checks");
        lab_lp_result_free(&fixed);
        return -1;
    }
    long double normalized = score(s, fixed.primal);
    long double delta = s->result->has_primal ? difference(s, fixed.primal, s->result->primal) : 0;
    long double repair_gap = difference(s, fixed.primal, lp->primal);
    if (!isfinite(normalized) || !isfinite(delta) || !isfinite(repair_gap)) {
        outcome(s, LAB_MIP_NUMERICAL_FAILURE, "Objective comparison overflow");
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
/* The root LP prices only guide this exact search. The integer subsolver
 * checks the original model structure and restores all dual inequalities
 * in integer arithmetic before excluding a branch. */
static int root_partition(Search *s, const LabLpResult *lp, const Node *node) {
    const LabModel *m = s->model;
    /* Charge the exact search to the global node budget. Reserve one node
     * for the ordinary tree when a cover cannot be completed. */
    size_t left = s->options->node_limit - s->result->nodes;
    size_t limit = left > 1 ? left - 1 : 0;
    if (m->rows > 256 || m->cols < 64 || !limit || s->options->no_presolve ||
        s->options->lp.no_exact)
        return 0;
    for (size_t i = 0; i < m->rows; i++)
        if (m->row_lower[i] != 1 || m->row_upper[i] != 1)
            return 0;
    int fractional = 0;
    for (size_t j = 0; j < m->cols && !fractional; j++)
        fractional = m->integer[j] && fabs(lp->primal[j] - round(lp->primal[j])) > 1e-7;
    if (!fractional)
        return 0;
    LabPartitionResult cover;
    double *point = s->repair_lower;
    lab_partition_search(m, lp->dual, lab_now() + fmin(2, .5 * remaining(s)), limit, point,
                         &cover);
    s->result->nodes += cover.nodes;
    s->result->partition_nodes += cover.nodes;
    double objective;
    if (!cover.has_point || !valid_candidate(s, point, 1, &objective))
        return 0;
    long double normalized = score(s, point);
    if (!isfinite(normalized))
        return 0;
    if (s->result->has_primal && difference(s, point, s->result->primal) > 0) return 0;
    if (!s->result->primal) s->result->primal = lab_alloc(m->cols, sizeof(double));
    memcpy(s->result->primal, point, m->cols * sizeof(double));
    s->result->has_primal = 1;
    s->result->objective = objective;
    s->incumbent = normalized;
    if (!cover.optimal || remaining(s) <= 0) {
        trace(s, node, "incumbent", SIZE_MAX, objective, 1);
        return 0;
    }
    close_at_incumbent(s, "partition", LAB_CLOSE_EXHAUSTED, LAB_BOUND_EXHAUSTION,
                       "Set-partitioning search exhausted; original witness checked");
    trace(s, node, "partition_search_optimum", SIZE_MAX, objective, 1);
    return 1;
}
/* A failed heuristic establishes no node status or bound. Only a candidate
 * checked against the original model can update the ordinary incumbent. */
static void root_diving(Search *s, const LabLpResult *lp, const Node *node) {
    s->result->root_dive_attempts++;
    if (s->cuts.count && !cut_box_contains(s, &s->cuts, s->lower, s->upper)) {
        s->result->root_dive_status = "box_mismatch";
        return;
    }
    double seconds = fmin(1, .05 * remaining(s));
    if (seconds <= 0) {
        s->result->root_dive_status = "time_limit";
        return;
    }
    LabDiveOptions options = lab_dive_options(seconds);
    options.cold = s->options->cold;
    LabDiveResult dive = {0};
    lab_primal_dive(s->model, &s->lp_model, s->lower, s->upper, lp->primal,
                    &s->options->lp, &options, &dive);
    s->result->root_dive_status = lab_dive_status_name(dive.status);
    s->result->root_dive_lp_solves += dive.lp_solves;
    s->result->root_dive_iterations += dive.lp_iterations;
    s->result->root_dive_steps += dive.steps;
    s->result->root_dive_backtracks += dive.backtracks;
    s->result->root_dive_elapsed += dive.elapsed;
    if (dive.primal) {
        s->result->root_dive_candidates++;
        double objective;
        if (!valid_candidate(s, dive.primal, 1, &objective))
            s->result->root_dive_status = "candidate_rejected";
        else {
            long double normalized = score(s, dive.primal);
            long double delta = s->result->has_primal
                ? difference(s, dive.primal, s->result->primal) : 0;
            if (!isfinite(normalized) || !isfinite(delta))
                s->result->root_dive_status = "objective_comparison_unavailable";
            else if (!s->result->has_primal || delta < 0) {
                if (!s->result->primal)
                    s->result->primal = lab_alloc(s->model->cols, sizeof(double));
                memcpy(s->result->primal, dive.primal, s->model->cols * sizeof(double));
                s->result->objective = objective;
                s->result->has_primal = 1;
                s->incumbent = normalized;
                s->result->root_dive_incumbents++;
                trace(s, node, "incumbent", SIZE_MAX, objective, 1);
            }
        }
    }
    lab_dive_result_free(&dive);
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
        outcome(s, LAB_MIP_NUMERICAL_FAILURE, "Cannot make a progressing integer branch");
        return -1;
    }
    /* Hybrid best-bound + plunging: plunge down the closer child immediately
     * and place the sibling into the best-bound heap. */
    Node *left_child = NULL, *right_child = NULL;
    if (has_right) {
        right_child = new_node(s, node, j, right, s->upper[j], node->bound, node->has_bound);
        right_child->parent_lp_score = lp_score;
        right_child->branch_distance = right - x;
        right_child->direction = LAB_SEARCH_UP;
    }
    if (has_left) {
        left_child = new_node(s, node, j, s->lower[j], left, node->bound, node->has_bound);
        left_child->parent_lp_score = lp_score;
        left_child->branch_distance = x - left;
        left_child->direction = LAB_SEARCH_DOWN;
    }
    if (s->plunge) {
        push(s, s->plunge);
        s->plunge = NULL;
    }
    int prefer_right = !s->options->first_fractional && has_left && has_right &&
                       (right - x) < (x - left);
    if (prefer_right) {
        s->plunge = right_child;
        push(s, left_child);
    } else if (has_left) {
        s->plunge = left_child;
        if (right_child)
            push(s, right_child);
    } else if (has_right) {
        s->plunge = right_child;
    }
    trace(s, node, "branch", j, x, 1);
    return 0;
}
void lab_bb_result_free(LabBbResult *result) {
    free(result->primal);
    lab_assignment_free_decimals(result->primal_decimal, result->primal_columns);
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
#ifdef LAB_SCIP_INTEGER
    result->scip_route = options->no_scip ? "disabled" : "not_dispatched";
#else
    result->scip_route = "unavailable";
#endif
    result->root_cut_status = options->no_root_cuts ? "disabled" : "not_attempted";
    result->root_dive_status = options->no_root_diving ? "disabled" : "not_attempted";
    outcome(&s, LAB_MIP_NUMERICAL_FAILURE, "MILP search did not complete");
    if (!isfinite(options->time_limit) || !options->node_limit) {
        outcome(&s, LAB_MIP_INPUT_ERROR, "Positive finite time and node limits are required");
        goto cleanup;
    }
    if (options->time_limit <= 0) {
        outcome(&s, LAB_MIP_TIME_LIMIT, "MILP time budget exhausted before search");
        goto cleanup;
    }
    if (options->trace_path) {
        s.trace = fopen(options->trace_path, "wx");
        if (!s.trace) {
            outcome(&s, LAB_MIP_INPUT_ERROR, "Cannot create fresh search trace");
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
                outcome(&s, LAB_MIP_UNSUPPORTED, "Integer bounds exceed the exactly spaced double range");
                goto cleanup;
            }
            lo = ceil(lo);
            hi = floor(hi);
        }
        if (isnan(lo) || isnan(hi)) {
            outcome(&s, LAB_MIP_INPUT_ERROR, "Nonfinite column-bound encoding");
            goto cleanup;
        }
        if (lo > hi) {
            publish(&s, (LabSearchResult){.route = "domain", .status = LAB_MIP_INFEASIBLE, .closure = LAB_CLOSE_DOMAIN},
                    "A column has no permitted integer or continuous value");
            goto cleanup;
        }
        s.base_lower[j] = lo;
        s.base_upper[j] = hi;
    }
    if (!options->no_presolve) {
        double available = remaining(&s);
        if (available <= 0) {
            outcome(&s, LAB_MIP_TIME_LIMIT, "MILP time budget exhausted before presolve");
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
            outcome(&s, LAB_MIP_TIME_LIMIT, "MILP time budget exhausted during root presolve");
            lab_presolve_result_free(&presolved);
            goto cleanup;
        }
        if (!strcmp(presolved.status, "input_error") || !strcmp(presolved.status, "infeasible")) {
            int infeasible = !strcmp(presolved.status, "infeasible");
            publish(&s, (LabSearchResult){.route = "presolve",
                    .status = infeasible ? LAB_MIP_INFEASIBLE : LAB_MIP_INPUT_ERROR,
                    .closure = infeasible ? LAB_CLOSE_DOMAIN : LAB_CLOSE_NONE}, presolved.message);
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
    for (size_t j = 0; j < n; j++) {
        double lo = s.base_lower[j], hi = s.base_upper[j];
        if (model->integer[j] &&
            ((isfinite(lo) && fabs(lo) > 9007199254740991.0) ||
             (isfinite(hi) && fabs(hi) > 9007199254740991.0))) {
            outcome(&s, LAB_MIP_UNSUPPORTED, "Integer bounds exceed the exactly spaced double range");
            goto cleanup;
        }
        s.reference[j] = isfinite(lo) ? lo : isfinite(hi) ? hi : 0;
    }
    double source_constant;
    if (lab_source_objective_value(model, s.reference, &source_constant)) {
        outcome(&s, LAB_MIP_NUMERICAL_FAILURE, "Source objective reference unavailable");
        goto cleanup;
    }
    s.constant = source_constant;
    if (!options->first_fractional)
        (void)lab_search_init(&s.branching, n); /* allocation failure keeps deterministic fallback */
    push(&s, new_node(&s, NULL, SIZE_MAX, 0, 0, 0, 0));
    while (s.plunge || s.heap_size) {
        if (remaining(&s) <= 0) {
            global_bound(&s, NULL);
            outcome(&s, LAB_MIP_TIME_LIMIT, "MILP time budget exhausted");
            goto cleanup;
        }
        if (result->nodes >= options->node_limit) {
            global_bound(&s, NULL);
            outcome(&s, LAB_MIP_NODE_LIMIT, "MILP node budget exhausted");
            goto cleanup;
        }
        Node *node = next_node(&s);
        node->next = NULL;
        bounds(&s, node);
        if (node->has_bound && cannot_improve(&s, node->bound)) {
            result->pruned_nodes++;
            trace(&s, node, "bound_prune", SIZE_MAX, (double)(s.constant + s.sense * node->bound),
                  1);
            release(node);
            continue;
        }
        LabLpResult lp;
        result->nodes++;
        LAB_PHASE(&s, "specialized");
        if (!node->parent && !options->no_presolve && !options->lp.no_exact) {
            size_t pair_nodes = 0;
            int paired = lab_pair_assignment(model, lab_now() + fmin(2, .4 * remaining(&s)),
                                             (options->node_limit - result->nodes) / 2,
                                             s.repair_lower, &pair_nodes);
            result->nodes += pair_nodes;
            if (paired && complete_point_search(&s, node, "pair_assignment", LAB_CLOSE_EXHAUSTED,
                    "Pair assignment search exhausted; original witness checked",
                    "pair_assignment_optimum")) {
                release(node);
                goto cleanup;
            }
            size_t trials = 0;
            char **decimal = NULL;
            int matched = lab_assignment_solve(model, s.base_lower, s.base_upper,
                              lab_now() + fmin(1, .4 * remaining(&s)),
                              (options->node_limit - result->nodes) / 2, s.repair_lower,
                              &decimal, &trials);
            result->nodes += trials;
            if (matched && complete_point_search(&s, node, "assignment", LAB_CLOSE_MATCHING_BOUND,
                    "Assignment met an exact matching bound with feasible source potentials",
                    "assignment_optimum")) {
                result->primal_decimal = decimal;
                result->primal_columns = model->cols;
                release(node);
                goto cleanup;
            }
            lab_assignment_free_decimals(decimal, model->cols);
#ifdef LAB_SCIP_INTEGER
            if (trials && !options->no_scip && !options->first_fractional &&
                result->nodes < options->node_limit && remaining(&s) > 0.1) {
                LabScipMipStats compact;
                double cut_off = lab_now() + fmin(5.0, 0.35 * remaining(&s));
                int solved_compact = lab_assignment_quotient(model, s.base_lower, s.base_upper, cut_off,
                    (options->node_limit - result->nodes) / 2, options->lp.iteration_limit,
                    s.repair_lower, &decimal, &compact);
                result->nodes += compact.nodes;
                result->lp_solves += compact.lp_solves;
                result->lp_iterations += compact.lp_iterations;
                result->scip_path_cliques += compact.path_cliques;
                if (solved_compact && store_terminal_point(&s)) {
                    /* The quotient preserves objective coefficients, sense and offset. */
                    LabSearchResult quotient = {.route = "assignment_quotient", .status = LAB_MIP_OPTIMAL,
                                                .closure = LAB_CLOSE_BACKEND};
                    if (compact.has_bound)
                        quotient.bound = (LabSearchBound){LAB_BOUND_BACKEND, compact.best_bound};
                    publish(&s, quotient, "Assignment quotient reported optimality; lifted source witness checked");
                    trace(&s, node, "assignment_quotient_optimum", SIZE_MAX, result->objective, 1);
                    result->primal_decimal = decimal;
                    result->primal_columns = model->cols;
                    result->scip_route = "assignment_quotient";
                    release(node);
                    goto cleanup;
                }
                lab_assignment_free_decimals(decimal, model->cols);
            }
#endif
            int solved = lab_fixed_charge_flow(model, lab_now() + fmin(2, .5 * remaining(&s)),
                              (options->node_limit - result->nodes) / 2, s.repair_lower,
                              &result->flow_eliminations, &result->flow_terminal_states);
            result->nodes += result->flow_eliminations + result->flow_terminal_states;
            if (solved && complete_point_search(&s, node, "flow", LAB_CLOSE_EXHAUSTED,
                    "Fixed-charge network search exhausted; original witness checked", "flow_optimum")) {
                release(node);
                goto cleanup;
            }
            size_t lot_nodes = 0;
            solved = lab_lot_sizing(model, lab_now() + fmin(1, .35 * remaining(&s)),
                                    (options->node_limit - result->nodes) / 2,
                                    s.repair_lower, &lot_nodes);
            result->nodes += lot_nodes;
            if (solved && complete_point_search(&s, node, "lot_sizing", LAB_CLOSE_EXHAUSTED,
                    "Line production and shared opening search exhausted; original witness checked",
                        "lot_sizing_optimum")) {
                release(node);
                goto cleanup;
            }
            solved = lab_resource_packing(model, lab_now() + fmin(1, .35 * remaining(&s)),
                                          (options->node_limit - result->nodes) / 2, s.repair_lower,
                                          &result->packing_attempts);
            result->nodes += result->packing_attempts;
            if (solved && complete_point_search(&s, node, "packing", LAB_CLOSE_MATCHING_BOUND,
                    "Resource packing met the capacity lower bound; original witness checked", "packing_optimum")) {
                release(node);
                goto cleanup;
            }
#ifdef LAB_GLPK_BASIS
            LabPooledStats pool_stats;
            LabSearchResult pooled = lab_pooled_assignment(model, lab_now() + fmin(2, .4 * remaining(&s)),
                             (options->node_limit - result->nodes) / 2,
                             options->lp.iteration_limit, s.repair_lower, &pool_stats);
            result->nodes += pool_stats.nodes;
            result->lp_solves += pool_stats.lp_solves;
            result->lp_iterations += pool_stats.iterations;
            if (pooled.status == LAB_MIP_OPTIMAL && complete_point_search(&s, node, pooled.route, pooled.closure,
                    "Pooled assignment search closed with exact network bounds and an original witness",
                    "pooled_assignment_optimum")) {
                release(node);
                goto cleanup;
            }
            if (pooled.point)
                retain_search_incumbent(&s, node);
            size_t cover_iterations = 0;
            solved = result->nodes < options->node_limit &&
                     lab_activation_cover(model, lab_now() + fmin(2, .4 * remaining(&s)),
                                          options->lp.iteration_limit, s.repair_lower,
                                          &cover_iterations);
            result->lp_solves += solved || cover_iterations > 0;
            result->lp_iterations += cover_iterations;
            if (solved && complete_point_search(&s, node, "activation_cover", LAB_CLOSE_MATCHING_BOUND,
                    "Activated covering met a verified integer lower bound; original witness checked",
                    "activation_cover_optimum")) {
                release(node);
                goto cleanup;
            }
            size_t pattern_work = 0;
            solved = lab_packing_patterns(model, lab_now() + fmin(1, .35 * remaining(&s)),
                                          (options->node_limit - result->nodes) / 2,
                                          s.repair_lower, &pattern_work);
            result->packing_attempts += pattern_work;
            result->nodes += pattern_work;
            if (solved && complete_point_search(&s, node, "pattern", LAB_CLOSE_MATCHING_BOUND,
                    "Integer packing met its pattern lower bound; original witness checked", "pattern_optimum")) {
                release(node);
                goto cleanup;
            }
#endif
            if (lab_binary_search(model, lab_now() + fmin(2, .5 * remaining(&s)), options->node_limit,
                                s.repair_lower, &result->binary_search_queries) &&
                complete_point_search(&s, node, "binary_search", LAB_CLOSE_EXHAUSTED,
                    "Binary assignment search exhausted; original witness checked", "binary_search_optimum")) {
                release(node);
                goto cleanup;
            }
        }
#ifdef LAB_KISSAT
        LAB_PHASE(&s, "sat");
        if (!node->parent && !options->no_presolve && !options->lp.no_exact &&
            options->binary_sat >= 0 && (options->binary_sat || (n >= 64 && model->rows > 64))) {
            /* SAT either closes quickly or not at all on the benchmark corpus; cap its
             * share so a stalled encoding cannot starve the branch-and-cut engine. */
            double sat_budget = fmin(.85 * remaining(&s), fmax(25.0, .25 * remaining(&s)));
            LabSearchResult sat = lab_binary_sat(model, s.lower, s.upper, lab_now() + sat_budget,
                                     options->sat_decision_limit ? options->sat_decision_limit : options->node_limit,
                                     s.repair_lower, &result->sat_calls);
            double objective;
            if (sat.status == LAB_MIP_INFEASIBLE && remaining(&s) > 0) {
                publish(&s, sat, "Exact binary constraint encoding reported unsatisfiable");
                release(node);
                goto cleanup;
            }
            if (sat.point && valid_candidate(&s, s.repair_lower, 1, &objective) && remaining(&s) > 0) {
                result->primal = lab_alloc(n, sizeof(double));
                memcpy(result->primal, s.repair_lower, n * sizeof(double));
                result->has_primal = 1;
                result->objective = objective;
                s.incumbent = score(&s,result->primal);
                if (sat.status == LAB_MIP_OPTIMAL) {
                    publish(&s, sat, "Exact binary constraint search closed; original witness checked");
                    trace(&s, node, "binary_sat_optimum", SIZE_MAX, objective, 1);
                    release(node);
                    goto cleanup;
                }
                trace(&s,node,"incumbent",SIZE_MAX,objective,1);
            }
        }
#endif
#ifdef LAB_SCIP_INTEGER
        if (!node->parent && !options->no_presolve && !options->lp.no_exact && !options->cold &&
            !options->first_fractional && !options->no_scip) {
            LAB_PHASE(&s, "integer_mip");
            LabIntegerMipStats stats;
            LabSearchResult contracted = lab_integer_mip(model, s.base_lower, s.base_upper,
                                        lab_now() + .75 * remaining(&s),
                                        options->node_limit - result->nodes,
                                        options->lp.iteration_limit, s.repair_lower, &stats);
            result->nodes += stats.nodes;
            result->lp_solves += stats.lp_solves;
            result->lp_iterations += stats.lp_iterations;
            if (contracted.status == LAB_MIP_OPTIMAL && store_terminal_point(&s)) {
                close_at_incumbent(&s, contracted.route, contracted.closure, LAB_BOUND_INTEGER_LATTICE,
                                   "SCIP reported optimality after exact integer equality contraction");
                trace(&s, node, "integer_search_optimum", SIZE_MAX, result->objective, 1);
                result->scip_route = "contracted";
                release(node);
                goto cleanup;
            }
            if (contracted.point)
                retain_search_incumbent(&s, node);
            if (stats.nodes && result->nodes >= options->node_limit) {
                global_bound(&s, node);
                outcome(&s, LAB_MIP_NODE_LIMIT, "MILP node budget exhausted during integer search");
                release(node);
                goto cleanup;
            }
            int is_partition = model->rows <= 256 && model->cols >= 64;
            for (size_t i = 0; is_partition && i < model->rows; i++)
                if (model->row_lower[i] != 1 || model->row_upper[i] != 1)
                    is_partition = 0;
            LAB_PHASE(&s, "scip_mip");
            if (!is_partition && run_scip_mip(&s, node)) {
                release(node);
                goto cleanup;
            }
        }
#endif
        if (call_lp(&s, s.lower, s.upper, &lp)) {
            global_bound(&s, node);
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
node_lp_ready:
        trace(&s, node, lp.status, SIZE_MAX, lp.objective, lp.has_primal);
#ifdef LAB_SCIP_INTEGER
        if (!node->parent && !options->no_presolve && !options->lp.no_exact && !options->cold &&
            !options->first_fractional && !options->no_scip &&
            model->cols >= 32 && strcmp(result->scip_route, "general_fallback") &&
            strcmp(lp.status, "optimal") && strcmp(lp.status, "infeasible") &&
            strcmp(lp.status, "time_limit") && strcmp(lp.status, "unsupported") &&
            run_scip_mip(&s, node)) {
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
#endif
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
            outcome(&s, LAB_MIP_NUMERICAL_FAILURE, "Node LP returned no finite primal objective");
            global_bound(&s, node);
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
        double recomputed;
        if (!valid_candidate(&s, lp.primal, 0, &recomputed) ||
            fabs(recomputed - lp.objective) > 1e-7 * (1 + fabs(recomputed))) {
            outcome(&s, LAB_MIP_NUMERICAL_FAILURE, "Node LP witness fails original-model checks");
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
            outcome(&s, LAB_MIP_UNSUPPORTED, "Integer LP values exceed the exactly spaced double range");
            global_bound(&s, node);
            lab_lp_result_free(&lp);
            release(node);
            goto cleanup;
        }
        if (!node->parent && !result->root_cut_lp_solves && root_partition(&s, &lp, node)) {
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
        if (!node->parent && !options->no_root_diving && !result->root_dive_attempts) {
            int fractional = 0;
            for (size_t j = 0; j < model->cols; ++j)
                if (model->integer[j] && fabs(lp.primal[j] - round(lp.primal[j])) > 1e-7)
                    fractional = 1;
            if (fractional)
                root_diving(&s, &lp, node);
        }
        long double normalized = score(&s, lp.primal);
        long double lower = normalized - 1e-7L * (1 + fabsl(normalized));
        if (!isfinite(normalized) || !isfinite(lower)) {
            outcome(&s, LAB_MIP_NUMERICAL_FAILURE, "Node objective comparison overflow");
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
        if (cannot_improve(&s, node->bound)) {
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
            outcome(&s, selection == LAB_SEARCH_UNSUPPORTED ? LAB_MIP_UNSUPPORTED : LAB_MIP_NUMERICAL_FAILURE,
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
                outcome(&s, selection == LAB_SEARCH_UNSUPPORTED ? LAB_MIP_UNSUPPORTED : LAB_MIP_NUMERICAL_FAILURE,
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
        close_at_incumbent(&s, "native", LAB_CLOSE_EXHAUSTED, LAB_BOUND_EXHAUSTION,
                           "Numerical branch-and-bound tree exhausted; no independent tree proof");
    } else {
        publish(&s, (LabSearchResult){.route = "native", .status = LAB_MIP_INFEASIBLE, .closure = LAB_CLOSE_EXHAUSTED},
                "Numerical branch-and-bound tree exhausted without an incumbent");
    }
cleanup:
    if (s.frontier.origin && result->search.route && !strcmp(result->search.route, "native") &&
        result->search.status != LAB_MIP_OPTIMAL && result->search.status != LAB_MIP_INFEASIBLE) {
        LabSearchResult stopped = result->search;
        stopped.bound = s.frontier;
        /* The stop reason can be set before the frontier is evaluated. */
        publish(&s, stopped, NULL);
    }
    if (!options->lp.no_exact && !result->has_primal && s.repair_lower &&
        (!strcmp(result->status, "infeasible") || !strcmp(result->status, "numerical_failure"))) {
        double objective;
        int recovered = lab_recover_square_equalities(model, s.repair_lower, &objective,
                                                      s.started + options->time_limit);
        if (recovered > 0) {
            free(result->primal);
            result->primal = s.repair_lower;
            s.repair_lower = NULL;
            result->has_primal = 1;
            result->objective = objective;
            result->exact_primal_recoveries++;
            close_at_incumbent(&s, "equality_recovery", LAB_CLOSE_UNIQUE_POINT, LAB_BOUND_UNIQUE_POINT,
                               "Recovered a unique source equality solution and checked its witness");
            trace(&s, &(Node){0}, "exact_equality_recovery", SIZE_MAX, objective, 1);
        } else if (recovered < 0) {
            outcome(&s, LAB_MIP_NUMERICAL_FAILURE, "Exact equality solution is feasible but its witness cannot be exported");
        }
    }
    if (s.plunge) {
        release(s.plunge);
        s.plunge = NULL;
    }
    for (size_t i = 0; i < s.heap_size; i++)
        release(s.heap[i]);
    free(s.heap);
#ifdef LAB_LP_SESSION_API
    lab_lp_session_free(s.session);
#endif
    lab_root_cut_batch_free(&s.cuts);
    lab_search_free(&s.branching);
    if (s.trace) {
        int failed = ferror(s.trace);
        if (fclose(s.trace))
            failed = 1;
        if (failed)
            outcome(&s, LAB_MIP_INPUT_ERROR, "Cannot write search trace");
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

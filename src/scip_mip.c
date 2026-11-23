#include "scip_mip.h"
#include "scip_path_cliques.h"
#include "scip_feasjump.h"
#include "scip_presol_coefftight.h"
#include "scip_presol_parallel.h"
#include "common.h"
#include "exact.h"
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static double subscip_fraction(void) {
    const char *value = getenv("LAB_SUBSCIP_FRAC");
    if (!value)
        return 0.25;
    char *end;
    errno = 0;
    double fraction = strtod(value, &end);
    if (end != value && !*end && !errno && isfinite(fraction) &&
        fraction >= 0.0 && fraction <= 1.0)
        return fraction;
    fprintf(stderr, "Invalid LAB_SUBSCIP_FRAC, using the default of 0.25\n");
    return 0.25;
}

static int extract_verified(const LabModel *m, const double *lower, const double *upper,
                            SCIP *s, SCIP_VAR *const *vars, double *point, double *out_obj) {
    int nsols = SCIPgetNSols(s);
    if (nsols <= 0)
        return 0;
    SCIP_SOL **sols = SCIPgetSols(s);
    int limit = nsols < 8 ? nsols : 8;
    double *cand = lab_alloc(m->cols ? m->cols : 1, sizeof *cand);
    LabModel view = *m;
    view.col_lower = (double *)lower;
    view.col_upper = (double *)upper;
    int found = 0;
    for (int idx = 0; idx < limit; ++idx) {
        SCIP_SOL *sol = sols[idx];
        if (!sol)
            continue;
        int ok = 1;
        for (size_t j = 0; j < m->cols; ++j) {
            double v = SCIPgetSolVal(s, sol, vars[j]);
            if (!isfinite(v)) {
                ok = 0;
                break;
            }
            if (m->integer && m->integer[j]) {
                if (fabs(v) > 9007199254740991.0 || fabs(v - round(v)) > 1e-5) {
                    ok = 0;
                    break;
                }
                v = round(v);
            }
            if (v < lower[j] - 1e-7 || v > upper[j] + 1e-7) {
                ok = 0;
                break;
            }
            if (v < lower[j])
                v = lower[j];
            if (v > upper[j])
                v = upper[j];
            cand[j] = v;
        }
        if (!ok)
            continue;
        double obj = 0, viol = 0;
        if (!lab_exact_evaluate(&view, cand, &obj, &viol) && viol <= 1e-7) {
            memcpy(point, cand, m->cols * sizeof *point);
            *out_obj = obj;
            found = 1;
            break;
        }
    }
    free(cand);
    return found;
}

/* A signed sum of integer variables is integral. With a common source
 * coefficient magnitude, round the divided side in the inequality direction. */
static int scale_uniform_integer_row(const LabModel *m, size_t i, int upper, int len,
                                     SCIP_VAR **var, const int *entry, double *val,
                                     double *side, double deadline) {
    double magnitude = fabs(val[0]);
    for (int k = 0; k < len; ++k)
        if (fabs(val[k]) != magnitude || SCIPvarGetType(var[k]) == SCIP_VARTYPE_CONTINUOUS)
            return 0;

    mpq_t unit, q;
    mpq_inits(unit, q, NULL);
    int scaled = 0;
    for (int k = 0; k < len; ++k) {
        if (((k & 255) == 0 && lab_now() >= deadline) ||
            lab_source_exact_coefficient(m, (size_t)entry[k], q))
            goto done;
        if (val[k] < 0)
            mpq_neg(q, q);
        if (mpq_sgn(q) <= 0 || (k && !mpq_equal(q, unit)))
            goto done;
        if (!k)
            mpq_set(unit, q);
    }
    if (lab_source_exact_row_bound(m, i, upper, q))
        goto done;
    mpq_div(q, q, unit);
    mpz_t rounded;
    mpz_init(rounded);
    if (upper)
        mpz_fdiv_q(rounded, mpq_numref(q), mpq_denref(q));
    else
        mpz_cdiv_q(rounded, mpq_numref(q), mpq_denref(q));
    if (mpz_sizeinbase(rounded, 2) <= 53) {
        *side = mpz_get_d(rounded);
        for (int k = 0; k < len; ++k)
            val[k] = copysign(1.0, val[k]);
        scaled = 1;
    }
    mpz_clear(rounded);
done:
    mpq_clears(unit, q, NULL);
    return scaled;
}

/* For a positive binary covering inequality with uniform light coefficients,
 * either one heavy item satisfies the row or at least ceil(rhs/light) light
 * items are needed. Compute that ceiling from the source rationals. */
static int scale_binary_cover(const LabModel *m, size_t i, int upper, int len, SCIP_VAR **var,
                                     const int *entry, double *val, double *side, double deadline) {
    if (upper ? *side >= 0 : *side <= 0)
        return 0;
    int scaled = 0;
    double rhs = fabs(*side), light = 0.0;
    int heavy_count = 0;
    for (int k = 0; k < len; ++k) {
        double a = upper ? -val[k] : val[k];
        if (!SCIPvarIsBinary(var[k]) || a <= 0)
            return 0;
        if (a >= rhs)
            ++heavy_count;
        else if (light == 0.0)
            light = a;
        else if (a != light)
            return 0;
    }
    if (light == 0.0)
        return 0;
    mpq_t target, unit, coef, ratio;
    mpq_inits(target, unit, coef, ratio, NULL);
    if (lab_source_exact_row_bound(m, i, upper, target))
        goto done;
    if (upper)
        mpq_neg(target, target);
    if (mpq_sgn(target) <= 0)
        goto done;
    for (int k = 0; k < len; ++k) {
        if (((k & 255) == 0 && lab_now() >= deadline) ||
            lab_source_exact_coefficient(m, (size_t)entry[k], coef))
            goto done;
        if (upper)
            mpq_neg(coef, coef);
        if (mpq_sgn(coef) <= 0)
            goto done;
        int heavy = fabs(val[k]) >= rhs;
        if (heavy != (mpq_cmp(coef, target) >= 0))
            goto done;
        if (!heavy) {
            if (mpq_sgn(unit) && !mpq_equal(coef, unit))
                goto done;
            mpq_set(unit, coef);
        }
    }
    mpq_div(ratio, target, unit);
    mpz_t count;
    mpz_init(count);
    mpz_cdiv_q(count, mpq_numref(ratio), mpq_denref(ratio));
    /* If even all light items fall short, one heavy is required. */
    unsigned long max_needed = (unsigned long)(len - heavy_count) + 1;
    unsigned long needed = mpz_cmp_ui(count, max_needed) > 0 ?
                           max_needed : mpz_get_ui(count);
    mpz_clear(count);
    double sign = upper ? -1.0 : 1.0;
    for (int k = 0; k < len; ++k)
        val[k] = sign * (fabs(val[k]) >= rhs ? (double)needed : 1.0);
    *side = sign * (double)needed;
    scaled = 1;
done:
    mpq_clears(target, unit, coef, ratio, NULL);
    return scaled;
}

/* Integer variables make the activity a multiple of the rational gcd of all
 * coefficients. Divide by this step and round only the one-sided row bound.
 * This also handles a decimal bound that rounds onto the lattice in binary64. */
static void tighten_integer_lattice(const LabModel *m, size_t i, int upper, int len,
                                    SCIP_VAR **var, const int *entry, double *val,
                                    double *side, double deadline) {
    enum { MAX_TERMS = 64, MAX_BITS = 256 };
    if (len < 2 || len > MAX_TERMS || lab_now() >= deadline)
        return;
    for (int k = 0; k < len; ++k)
        if (SCIPvarGetType(var[k]) == SCIP_VARTYPE_CONTINUOUS)
            return;

    mpq_t step, q, target;
    mpq_inits(step, q, target, NULL);
    if (lab_source_exact_row_bound(m, i, upper, target) || !mpq_sgn(target) ||
        mpz_sizeinbase(mpq_numref(target), 2) > MAX_BITS ||
        mpz_sizeinbase(mpq_denref(target), 2) > MAX_BITS)
        goto done;
    double normalized[MAX_TERMS];
    for (int k = 0; k < len; ++k) {
        if (lab_source_exact_coefficient(m, (size_t)entry[k], q) || !mpq_sgn(q) ||
            mpz_sizeinbase(mpq_numref(q), 2) > MAX_BITS ||
            mpz_sizeinbase(mpq_denref(q), 2) > MAX_BITS)
            goto done;
        if (k == 0) {
            mpq_abs(step, q);
        } else {
            /* Inputs are reduced, so gcd(numerators)/lcm(denominators)
             * is already in lowest terms. */
            mpz_gcd(mpq_numref(step), mpq_numref(step), mpq_numref(q));
            mpz_lcm(mpq_denref(step), mpq_denref(step), mpq_denref(q));
            if (mpz_sizeinbase(mpq_denref(step), 2) > MAX_BITS)
                goto done;
        }
        /* Once the side belongs to this lattice it also belongs to the
         * finer lattice after including any remaining coefficients. */
        mpq_div(q, target, step);
        if (mpz_cmp_ui(mpq_denref(q), 1) == 0)
            goto done;
    }
    mpz_t rounded;
    mpz_init(rounded);
    if (upper)
        mpz_fdiv_q(rounded, mpq_numref(q), mpq_denref(q));
    else
        mpz_cdiv_q(rounded, mpq_numref(q), mpq_denref(q));
    int fits = mpz_sizeinbase(rounded, 2) <= 45;
    double bound = fits ? mpz_get_d(rounded) : 0;
    mpz_clear(rounded);
    if (!fits)
        goto done;
    for (int k = 0; k < len; ++k) {
        if (lab_source_exact_coefficient(m, (size_t)entry[k], q))
            goto done;
        mpq_div(q, q, step);
        if (mpz_cmp_ui(mpq_denref(q), 1) || mpz_sizeinbase(mpq_numref(q), 2) > 20)
            goto done;
        normalized[k] = mpz_get_d(mpq_numref(q));
    }
    memcpy(val, normalized, (size_t)len * sizeof *val);
    *side = bound;
done:
    mpq_clears(step, q, target, NULL);
}

static void scale_integer_inequality(const LabModel *m, size_t i, int upper, int len,
                                     SCIP_VAR **var, const int *entry, double *val,
                                     double *side, double deadline) {
    if (len && !scale_uniform_integer_row(m, i, upper, len, var, entry, val, side, deadline) &&
        !scale_binary_cover(m, i, upper, len, var, entry, val, side, deadline))
        tighten_integer_lattice(m, i, upper, len, var, entry, val, side, deadline);
}

enum { LOSSY_IMPORT = -3 };

static int compare_cost(const void *left, const void *right) {
    double a = *(const double *)left, b = *(const double *)right;
    return (a > b) - (a < b);
}

/* A per-unit price difference erased during dual presolve can change the
 * objective materially on a large activity. Check adjacent distinct prices
 * in both directions before handing the original model to that presolver. */
static int merges_distinct_costs(const LabModel *m, SCIP *s) {
    double *cost = lab_alloc(m->cols, sizeof *cost);
    memcpy(cost, m->objective, m->cols * sizeof *cost);
    qsort(cost, m->cols, sizeof *cost, compare_cost);
    int merged = 0;
    for (size_t j = 1; j < m->cols; ++j)
        if (cost[j] != cost[j - 1] && SCIPisEQ(s, cost[j], cost[j - 1])) {
            merged = 1;
            break;
        }
    free(cost);
    return merged;
}

static int solve_once(const LabModel *m, const double *lower, const double *upper,
                 const double *row_lower, const double *row_upper, const double *incumbent,
                 double deadline, size_t node_limit, size_t lp_iteration_limit, double *point,
                 LabScipMipStats *stats, int strict) {
    *stats = (LabScipMipStats){0};
    if (!m || !m->integer || !m->rows || !m->cols || !node_limit || !lp_iteration_limit ||
        node_limit > (size_t)LLONG_MAX || lp_iteration_limit > (size_t)LLONG_MAX ||
        m->rows > (size_t)INT_MAX || m->cols > (size_t)INT_MAX || m->nnz > (size_t)INT_MAX ||
        lab_now() >= deadline)
        return 0;

    int has_int = 0;
    if (!isfinite(m->offset) || fabs(m->offset) > 1e12)
        return 0;
    for (size_t j = 0; j < m->cols; ++j) {
        if (lower[j] > upper[j])
            return 0;
        if ((isfinite(lower[j]) && fabs(lower[j]) > 1e12) ||
            (isfinite(upper[j]) && fabs(upper[j]) > 1e12) ||
            fabs(m->objective[j]) > 1e12)
            return 0;
        if (m->integer[j])
            has_int = 1;
    }
    if (!has_int)
        return 0;

    SCIP *s = NULL;
    SCIP_VAR **vars = NULL;
    SCIP_CONS **cons = NULL;
    size_t *rstart = NULL;
    size_t *rcur = NULL;
    SCIP_VAR **rvars = NULL;
    SCIP_Real *rvals = NULL;
    int *source_entry = NULL;
    unsigned char *active = NULL;
    int result = 0;

    if (SCIPcreate(&s) != SCIP_OKAY)
        goto done;
    SCIPsetMessagehdlrQuiet(s, TRUE);
    if (SCIPincludeDefaultPlugins(s) != SCIP_OKAY ||
        lab_scip_include_feasjump(s) != SCIP_OKAY ||
        lab_scip_include_presol_parallel(s) != SCIP_OKAY ||
        lab_scip_include_presol_coefftight(s) != SCIP_OKAY ||
        SCIPaddRealParam(s, "heuristics/subsciptimefrac",
                         "maximal fraction of the remaining time given to a heuristic sub-SCIP",
                         NULL, FALSE, 0.25, 0.0, 1.0, NULL, NULL) != SCIP_OKAY ||
        SCIPsetRealParam(s, "heuristics/subsciptimefrac", subscip_fraction()) != SCIP_OKAY ||
        SCIPcreateProbBasic(s, "milp_lab") != SCIP_OKAY ||
        SCIPsetIntParam(s, "timing/clocktype", 2) != SCIP_OKAY ||
        SCIPsetIntParam(s, "parallel/maxnthreads", 1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "randomization/randomseedshift", 1) != SCIP_OKAY ||
        SCIPsetLongintParam(s, "limits/totalnodes", (SCIP_Longint)node_limit) != SCIP_OKAY ||
        SCIPsetLongintParam(s, "lp/iterlim", (SCIP_Longint)lp_iteration_limit) != SCIP_OKAY ||
        SCIPsetLongintParam(s, "lp/rootiterlim", (SCIP_Longint)lp_iteration_limit) != SCIP_OKAY ||
        SCIPsetRealParam(s, "limits/gap", 0.0) != SCIP_OKAY ||
        SCIPsetRealParam(s, "limits/absgap", 0.0) != SCIP_OKAY ||
        /* A small reduced cost on a large possible activity can still buy a
         * material objective improvement. Use the tighter LP pricing tolerance
         * for both initial and reimported general searches. */
        SCIPsetRealParam(s, "numerics/dualfeastol", 1e-9) != SCIP_OKAY ||
        /* Enable built-in C presolvers that default to 0 when PaPILO is expected,
         * and cap root aggregation-cut rounds to avoid dense LP blowup. */
        SCIPsetIntParam(s, "presolving/dualagg/maxrounds", -1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "presolving/dualinfer/maxrounds", -1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "presolving/stuffing/maxrounds", -1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "presolving/tworowbnd/maxrounds", -1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "presolving/redvub/maxrounds", -1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "presolving/sparsify/maxnonzeros", 64) != SCIP_OKAY ||
        SCIPsetRealParam(s, "presolving/restartfac", 0.02) != SCIP_OKAY ||
        SCIPsetRealParam(s, "presolving/immrestartfac", 0.08) != SCIP_OKAY ||
        SCIPsetRealParam(s, "presolving/restartminred", 0.05) != SCIP_OKAY ||
        SCIPsetIntParam(s, "propagating/probing/maxprerounds", 2) != SCIP_OKAY ||
        SCIPsetIntParam(s, "propagating/probing/maxuseless", 200) != SCIP_OKAY ||
        SCIPsetIntParam(s, "propagating/probing/maxtotaluseless", 400) != SCIP_OKAY ||
        SCIPsetIntParam(s, "separating/rapidlearning/freq", 0) != SCIP_OKAY ||
        SCIPsetIntParam(s, "separating/aggregation/maxroundsroot", 10) != SCIP_OKAY ||
        SCIPsetRealParam(s, "heuristics/rens/minfixingrate", 0.3) != SCIP_OKAY ||
        SCIPsetIntParam(s, "heuristics/rins/freq", 15) != SCIP_OKAY ||
        SCIPsetIntParam(s, "heuristics/farkasdiving/freq", -1) != SCIP_OKAY ||
        SCIPsetIntParam(s, "heuristics/conflictdiving/freq", 20) != SCIP_OKAY)
        goto done;
    if (SCIPgetParam(s, "heuristics/locks/maxlpiters") != NULL &&
        SCIPsetIntParam(s, "heuristics/locks/maxlpiters", 5000) != SCIP_OKAY)
        goto done;
    if ((m->rows >= 2000 || m->nnz >= 100000) &&
        SCIPsetBoolParam(s, "constraints/knapsack/disaggregation", FALSE) != SCIP_OKAY)
        goto done;
    /* Experiment hook: overlay a SCIP parameter file (unset in benchmark configurations). */
    if (getenv("LAB_SCIP_SET") && SCIPreadParams(s, getenv("LAB_SCIP_SET")) != SCIP_OKAY)
        goto done;
    if (strict && SCIPepsilon(s) > 1e-15 && SCIPsetRealParam(s, "numerics/epsilon", 1e-15) != SCIP_OKAY)
        goto done;

    if (merges_distinct_costs(m, s)) {
        result = LOSSY_IMPORT;
        goto done;
    }

    if (m->maximize && SCIPsetObjsense(s, SCIP_OBJSENSE_MAXIMIZE) != SCIP_OKAY)
        goto done;
    if (m->offset != 0.0 &&
        SCIPaddOrigObjoffset(s, m->offset) != SCIP_OKAY)
        goto done;

    vars = lab_alloc(m->cols, sizeof *vars);
    double inf = SCIPinfinity(s);
    for (size_t j = 0; j < m->cols; ++j) {
        if ((j & 1023U) == 0 && lab_now() >= deadline)
            goto done;
        double lo = isfinite(lower[j]) ? lower[j] : -inf;
        double hi = isfinite(upper[j]) ? upper[j] : inf;
        if ((isfinite(lower[j]) && fabs(lo) >= inf) || (isfinite(upper[j]) && fabs(hi) >= inf))
            goto done;
        if (m->objective[j] != 0.0 && SCIPisZero(s, m->objective[j])) {
            result = LOSSY_IMPORT;
            goto done;
        }
        SCIP_VARTYPE vt = SCIP_VARTYPE_CONTINUOUS;
        if (m->integer[j])
            vt = (lo >= 0.0 && hi <= 1.0) ? SCIP_VARTYPE_BINARY : SCIP_VARTYPE_INTEGER;
        char name[32];
        snprintf(name, sizeof name, "C%zu", j);
        if (SCIPcreateVarBasic(s, &vars[j], name, lo, hi, m->objective[j], vt) != SCIP_OKAY ||
            SCIPaddVar(s, vars[j]) != SCIP_OKAY)
            goto done;
    }

    active = lab_alloc(m->rows, sizeof *active);
    rstart = lab_alloc(m->rows + 1, sizeof *rstart);
    const double *rl = row_lower ? row_lower : m->row_lower;
    const double *ru = row_upper ? row_upper : m->row_upper;
    for (size_t i = 0; i < m->rows; ++i)
        if (isfinite(rl[i]) || isfinite(ru[i]))
            active[i] = 1;
    for (size_t j = 0; j < m->cols; ++j) {
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            if (active[i] && m->value[k] != 0.0)
                rstart[i + 1]++;
        }
    }
    for (size_t i = 0; i < m->rows; ++i)
        rstart[i + 1] += rstart[i];
    size_t total_nz = rstart[m->rows];
    rcur = lab_alloc(m->rows, sizeof *rcur);
    memcpy(rcur, rstart, m->rows * sizeof *rcur);
    rvars = lab_alloc(total_nz ? total_nz : 1, sizeof *rvars);
    source_entry = lab_alloc(total_nz ? total_nz : 1, sizeof *source_entry);
    rvals = lab_alloc(total_nz ? total_nz : 1, sizeof *rvals);
    for (size_t j = 0; j < m->cols; ++j) {
        if ((j & 1023U) == 0 && lab_now() >= deadline)
            goto done;
        for (size_t k = m->col_start[j]; k < m->col_start[j + 1]; ++k) {
            size_t i = m->row_index[k];
            double v = m->value[k];
            if (active[i] && v != 0.0) {
                size_t pos = rcur[i]++;
                rvars[pos] = vars[j];
                rvals[pos] = v;
                source_entry[pos] = (int)k;
            }
        }
    }
    cons = lab_alloc(m->rows, sizeof *cons);
    for (size_t i = 0; i < m->rows; ++i) {
        if ((i & 1023U) == 0 && lab_now() >= deadline)
            goto done;
        if (!active[i])
            continue;
        double l = isfinite(rl[i]) ? rl[i] : -inf;
        double u = isfinite(ru[i]) ? ru[i] : inf;
        if ((isfinite(rl[i]) && fabs(l) >= inf) || (isfinite(ru[i]) && fabs(u) >= inf))
            goto done;
        int len = (int)(rstart[i + 1] - rstart[i]);
        char name[32];
        snprintf(name, sizeof name, "R%zu", i);
        if (len && u == inf && l == m->row_lower[i] && !isfinite(m->row_upper[i]))
            scale_integer_inequality(m, i, 0, len, rvars + rstart[i], source_entry + rstart[i],
                               rvals + rstart[i], &l, deadline);
        else if (len && l == -inf && u == m->row_upper[i] && !isfinite(m->row_lower[i]))
            scale_integer_inequality(m, i, 1, len, rvars + rstart[i], source_entry + rstart[i],
                               rvals + rstart[i], &u, deadline);
        if (SCIPcreateConsBasicLinear(s, &cons[i], name, len, len ? &rvars[rstart[i]] : NULL,
                                      len ? &rvals[rstart[i]] : NULL, l, u) != SCIP_OKAY)
            goto done;
        if (SCIPgetNVarsLinear(s, cons[i]) != len ||
            SCIPgetLhsLinear(s, cons[i]) != l || SCIPgetRhsLinear(s, cons[i]) != u) {
            result = LOSSY_IMPORT;
            goto done;
        }
        SCIP_VAR **stored_vars = SCIPgetVarsLinear(s, cons[i]);
        SCIP_Real *stored_vals = SCIPgetValsLinear(s, cons[i]);
        for (int k = 0; k < len; ++k) {
            size_t pos = rstart[i] + (size_t)k;
            if (stored_vars[k] != rvars[pos] || stored_vals[k] != rvals[pos]) {
                result = LOSSY_IMPORT;
                goto done;
            }
        }
        if (SCIPaddCons(s, cons[i]) != SCIP_OKAY)
            goto done;
    }

    double now = lab_now();
    double cut_deadline = now + fmin(0.2, 0.02 * fmax(0.0, deadline - now));
    if (lab_scip_path_cliques(s, m, vars, rstart, source_entry, lower, upper, rl, ru,
                              cut_deadline, &stats->path_cliques) != SCIP_OKAY)
        goto done;

    if (incumbent) {
        SCIP_SOL *init_sol = NULL;
        if (SCIPcreateOrigSol(s, &init_sol, NULL) == SCIP_OKAY && init_sol) {
            int ok = 1;
            for (size_t j = 0; j < m->cols && ok; ++j)
                if (SCIPsetSolVal(s, init_sol, vars[j], incumbent[j]) != SCIP_OKAY)
                    ok = 0;
            SCIP_Bool stored = FALSE;
            if (ok)
                (void)SCIPaddSolFree(s, &init_sol, &stored);
            else
                (void)SCIPfreeSol(s, &init_sol);
        }
    }

    double remaining = deadline - lab_now();
    if (remaining <= 0.01) {
        stats->status = LAB_MIP_TIME_LIMIT;
        goto done;
    }
    if (SCIPsetRealParam(s, "limits/time", remaining * 0.96) != SCIP_OKAY)
        goto done;

    SCIP_RETCODE rc = SCIPsolve(s);
    if (rc != SCIP_OKAY)
        goto done;

    stats->nodes = (size_t)SCIPgetNTotalNodes(s);
    stats->lp_iterations = (size_t)SCIPgetNLPIterations(s);
    stats->lp_solves = (size_t)SCIPgetNLPs(s);

    SCIP_STATUS st = SCIPgetStatus(s);
    double dual = SCIPgetDualbound(s);
    if (isfinite(dual) && fabs(dual) < 0.5 * inf) {
        stats->has_bound = 1;
        stats->best_bound = dual;
    }

    if (st == SCIP_STATUS_INFEASIBLE) {
        stats->status = LAB_MIP_INFEASIBLE;
        result = -2;
        goto done;
    }

    double orig_obj = 0.0;
    int verified = extract_verified(m, lower, upper, s, vars, point, &orig_obj);
    if (st == SCIP_STATUS_OPTIMAL && verified) {
        double primal = SCIPgetPrimalbound(s);
        double scale = 1.0 + fabs(orig_obj);
        if (isfinite(primal) && stats->has_bound &&
            fabs(primal - orig_obj) <= 1e-4 * scale &&
            fabs(dual - orig_obj) <= 1e-4 * scale) {
            stats->status = LAB_MIP_OPTIMAL;
            result = 2;
            goto done;
        }
    }
    if (verified)
        result = 1;
    if (st == SCIP_STATUS_GAPLIMIT)
        stats->status = LAB_MIP_GAP_LIMIT;
    else if (st == SCIP_STATUS_TIMELIMIT)
        stats->status = LAB_MIP_TIME_LIMIT;
    else if (st == SCIP_STATUS_NODELIMIT || st == SCIP_STATUS_TOTALNODELIMIT)
        stats->status = LAB_MIP_NODE_LIMIT;

done:
    if (s) {
        if (cons)
            for (size_t i = 0; i < m->rows; ++i)
                if (cons[i])
                    (void)SCIPreleaseCons(s, &cons[i]);
        if (vars)
            for (size_t j = 0; j < m->cols; ++j)
                if (vars[j])
                    (void)SCIPreleaseVar(s, &vars[j]);
        if (SCIPfree(&s) != SCIP_OKAY) {
            result = 0;
            stats->has_bound = 0;
            stats->status = LAB_MIP_CONTINUE;
        }
    }
    free(cons);
    free(rvars);
    free(source_entry);
    free(rvals);
    free(rcur);
    free(rstart);
    free(active);
    free(vars);
    return result;
}

static int solve_with_retry(const LabModel *m, const double *lower, const double *upper,
                 const double *row_lower, const double *row_upper, const double *incumbent,
                 double deadline, size_t node_limit, size_t lp_iteration_limit, double *point,
                 LabScipMipStats *stats) {
    int found = solve_once(m, lower, upper, row_lower, row_upper, incumbent, deadline,
                           node_limit, lp_iteration_limit, point, stats, 0);
    if (found != -2 && found != LOSSY_IMPORT)
        return found;
    /* Reimport the untransformed problem with a tighter zero threshold.
     * A feasible point on the original model can contradict the first
     * search's cutoff. Discard that search's bound even on retry failure. */
    LabScipMipStats first = *stats;
    *stats = (LabScipMipStats){.nodes = first.nodes, .lp_solves = first.lp_solves,
                               .lp_iterations = first.lp_iterations};
    if (first.nodes >= node_limit || lab_now() >= deadline) {
        stats->status = first.nodes >= node_limit ? LAB_MIP_NODE_LIMIT : LAB_MIP_TIME_LIMIT;
        return 0;
    }
    found = solve_once(m, lower, upper, row_lower, row_upper, incumbent, deadline,
                       node_limit - first.nodes, lp_iteration_limit, point, stats, 1);
    stats->nodes += first.nodes;
    stats->lp_solves += first.lp_solves;
    stats->lp_iterations += first.lp_iterations;
    /* A second floating-point cutoff on a mixed model is not corroboration.
     * Let native search examine its relaxation against source decimals. */
    int crosscheck = incumbent != NULL;
    if (found == -2) {
        for (size_t j = 0; j < m->cols && !crosscheck; ++j)
            crosscheck = !m->integer[j] && m->col_start[j] < m->col_start[j + 1];
    }
    if (found == -2 && crosscheck) {
        stats->fallback_reason = incumbent ? "cutoff_conflicts_with_incumbent" : "mixed_cutoff_requires_source_check";
        stats->has_bound = 0;
        stats->status = LAB_MIP_CONTINUE;
        return 0;
    }
    if (found == LOSSY_IMPORT) {
        stats->fallback_reason = "lossy_model_import";
        return 0;
    }
    return found;
}

LabSearchResult lab_scip_mip(const LabModel *m, const double *lower, const double *upper,
                 const double *row_lower, const double *row_upper, const double *incumbent,
                 double deadline, size_t node_limit, size_t lp_iteration_limit, double *point,
                 LabScipMipStats *stats) {
    int found = solve_with_retry(m, lower, upper, row_lower, row_upper, incumbent, deadline,
                                 node_limit, lp_iteration_limit, point, stats);
    LabSearchResult r = {.route = "general_scip"};
    if (found == 1 || found == 2)
        r.point = point;
    if (found == 2 || found == -2) {
        r.status = found == 2 ? LAB_MIP_OPTIMAL : LAB_MIP_INFEASIBLE;
        r.closure = LAB_CLOSE_BACKEND;
    } else {
        r.status = stats->status;
    }
    if (stats->has_bound && r.status != LAB_MIP_CONTINUE)
        r.bound = (LabSearchBound){LAB_BOUND_BACKEND, stats->best_bound};
    return r;
}

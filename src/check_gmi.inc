/* Independent original-model GMI derivation checker. This translation-unit
 * component shares only check.c's exact reader and sparse algebra accumulator.
 * No optimizer parser, basis API, or cut-generation implementation is used. */
#include <float.h>
#include <math.h>

#if defined(__FAST_MATH__)
#error "Cut checking requires finite IEEE double semantics"
#endif
_Static_assert(FLT_RADIX == 2 && DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024,
               "Published cut checker targets binary64");

enum { GMI_MAX_TERMS = 4096, GMI_MAX_COLUMNS = 16384, GMI_MAX_ROWS = 16384, GMI_MAX_BITS = 65536 };
typedef struct {
    int lower_finite, upper_finite, seen;
    mpq_t lower, upper;
} GmiBox;
enum GmiTermKind { GMI_VARIABLE, GMI_SOURCE_ROW, GMI_BOUND_SLACK };
typedef struct {
    enum GmiTermKind kind;
    size_t index;
    int lower, integer;
    mpq_t coefficient, anchor, gamma;
} GmiTerm;
typedef struct {
    size_t basic, count, capacity;
    GmiBox *box;
    GmiTerm *term;
    mpq_t *equality, *fixed, *published, *exact;
    unsigned char *equality_seen, *fixed_seen, *published_seen;
    mpq_t beta, published_rhs, exact_rhs, error_lower, maximum_rhs, rounding_slack;
    int scope_seen, basic_seen, beta_seen, rhs_seen, identity_valid, premises_valid;
    int unavailable, cut_valid, mip_box_infeasible;
    double deadline;
} GmiProof;

static void gmi_proof_init(GmiProof *p) {
    memset(p, 0, sizeof *p);
    p->basic = SIZE_MAX;
    p->deadline = lab_now() + 20;
    mpq_inits(p->beta, p->published_rhs, p->exact_rhs, p->error_lower, p->maximum_rhs,
              p->rounding_slack, NULL);
}

static mpq_t *gmi_qvector(size_t n) {
    mpq_t *q = lab_alloc(n, sizeof *q);
    for (size_t j = 0; j < n; ++j)
        mpq_init(q[j]);
    return q;
}

static void gmi_proof_allocate(GmiProof *p, const Model *m) {
    p->box = lab_alloc(m->nc, sizeof *p->box);
    for (size_t j = 0; j < m->nc; ++j)
        mpq_inits(p->box[j].lower, p->box[j].upper, NULL);
    p->equality = gmi_qvector(m->nr);
    p->fixed = gmi_qvector(m->nc);
    p->published = gmi_qvector(m->nc);
    p->exact = gmi_qvector(m->nc);
    p->equality_seen = lab_alloc(m->nr, 1);
    p->fixed_seen = lab_alloc(m->nc, 1);
    p->published_seen = lab_alloc(m->nc, 1);
}

static void gmi_free_qvector(mpq_t *q, size_t n) {
    if (q)
        for (size_t j = 0; j < n; ++j)
            mpq_clear(q[j]);
    free(q);
}

static void gmi_proof_clear(GmiProof *p, const Model *m) {
    if (p->box)
        for (size_t j = 0; j < m->nc; ++j)
            mpq_clears(p->box[j].lower, p->box[j].upper, NULL);
    free(p->box);
    for (size_t k = 0; k < p->count; ++k)
        mpq_clears(p->term[k].coefficient, p->term[k].anchor, p->term[k].gamma, NULL);
    free(p->term);
    gmi_free_qvector(p->equality, m->nr);
    gmi_free_qvector(p->fixed, m->nc);
    gmi_free_qvector(p->published, m->nc);
    gmi_free_qvector(p->exact, m->nc);
    free(p->equality_seen);
    free(p->fixed_seen);
    free(p->published_seen);
    mpq_clears(p->beta, p->published_rhs, p->exact_rhs, p->error_lower, p->maximum_rhs,
               p->rounding_slack, NULL);
}

static int gmi_bounded_q(Model *m, GmiProof *p, const mpq_t q) {
    if (mpz_sizeinbase(mpq_numref(q), 2) > GMI_MAX_BITS ||
        mpz_sizeinbase(mpq_denref(q), 2) > GMI_MAX_BITS || lab_now() >= p->deadline) {
        p->unavailable = 1;
        fail(m, "GMI rational resource/deadline limit");
        return -1;
    }
    return 0;
}

static int gmi_parse_q(Model *m, GmiProof *p, mpq_t q, const char *s) {
    if (rational(q, s)) {
        fail(m, "invalid rational certificate field");
        return -1;
    }
    return gmi_bounded_q(m, p, q);
}

static int gmi_side(const char *s) {
    return !strcmp(s, "L") ? 1 : !strcmp(s, "U") ? 0 : -1;
}

static int gmi_finite_binary(const mpq_t q) {
    double d = mpq_get_d(q);
    if (!isfinite(d))
        return 0;
    mpq_t back;
    mpq_init(back);
    mpq_set_d(back, d);
    int equal = !mpq_cmp(back, q);
    mpq_clear(back);
    return equal;
}

static int gmi_tokens(char *line, char **fields, int capacity) {
    int count = 0;
    char *save;
    for (char *token = strtok_r(line, " \t\r\n", &save); token;
         token = strtok_r(NULL, " \t\r\n", &save)) {
        if (count >= capacity)
            return -1;
        fields[count++] = token;
    }
    return count;
}

static void gmi_read_certificate(Model *m, GmiProof *p, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fail(m, "cannot open GMI derivation certificate");
        return;
    }
    char line[32768], *v[12];
    int header = 0, ended = 0, ch = 0;
    size_t bytes = 0;
    while (!m->error[0]) {
        size_t length = 0;
        while ((ch = fgetc(f)) != EOF && ch != '\n') {
            if (!ch || length == sizeof line - 1) {
                fail(m, !ch ? "NUL byte in GMI certificate" : "GMI certificate line too long");
                break;
            }
            line[length++] = (char)ch;
        }
        if (m->error[0] || (!length && ch == EOF))
            break;
        bytes += length + 1;
        if (bytes > (1U << 24) || lab_now() >= p->deadline) {
            p->unavailable = 1;
            fail(m, "GMI certificate size/deadline limit");
            break;
        }
        line[length] = 0;
        char *start = line;
        while (isspace((unsigned char)*start))
            ++start;
        if (*start == '#')
            continue;
        int n = gmi_tokens(line, v, 12);
        if (!n)
            continue;
        if (ended || n < 0) {
            fail(m, "extra/malformed record after END or excess fields");
            break;
        }
        if (!header) {
            if (n != 1 || strcmp(v[0], "MILPLAB_GMI_1"))
                fail(m, "missing MILPLAB_GMI_1 header");
            header = 1;
            continue;
        }
        if (!strcmp(v[0], "SCOPE") && n == 2 && !p->scope_seen &&
            !strcmp(v[1], "ORIGINAL_MIP_IN_BOX")) {
            p->scope_seen = 1;
        } else if (!strcmp(v[0], "BASIC") && n == 2 && !p->basic_seen) {
            p->basic = lookup(&m->cols, v[1]);
            p->basic_seen = 1;
            if (p->basic == SIZE_MAX)
                fail(m, "unknown basic original column");
        } else if (!strcmp(v[0], "BETA") && n == 2 && !p->beta_seen) {
            p->beta_seen = 1;
            gmi_parse_q(m, p, p->beta, v[1]);
        } else if (!strcmp(v[0], "BOX") && n == 4) {
            size_t j = lookup(&m->cols, v[1]);
            if (j == SIZE_MAX || p->box[j].seen) {
                fail(m, "unknown/duplicate box column");
                continue;
            }
            GmiBox *box = &p->box[j];
            box->seen = 1;
            box->lower_finite = strcmp(v[2], "*") != 0;
            box->upper_finite = strcmp(v[3], "*") != 0;
            if (box->lower_finite)
                gmi_parse_q(m, p, box->lower, v[2]);
            if (box->upper_finite)
                gmi_parse_q(m, p, box->upper, v[3]);
            if (box->lower_finite && box->upper_finite && mpq_cmp(box->lower, box->upper) > 0)
                fail(m, "contradictory recorded box");
        } else if (!strcmp(v[0], "TERM") && n == 8) {
            char expected[32];
            snprintf(expected, sizeof expected, "%zu", p->count);
            int lower = gmi_side(v[6]);
            int kind = !strcmp(v[4], "VAR") ? GMI_VARIABLE : !strcmp(v[4], "ROW") ? GMI_SOURCE_ROW
                       : !strcmp(v[4], "BOUND") ? GMI_BOUND_SLACK : -1;
            int integer = !strcmp(v[3], "I") ? 1 : !strcmp(v[3], "C") ? 0 : -1;
            size_t j = lookup(kind == GMI_SOURCE_ROW ? &m->rows : &m->cols, v[5]);
            if (strcmp(expected, v[1]) || lower < 0 || kind < 0 || integer < 0 ||
                j == SIZE_MAX || (kind == GMI_SOURCE_ROW && m->r[j].type == 'N')) {
                fail(m, "invalid/unknown/nonconsecutive native term descriptor");
                continue;
            }
            if (p->count >= GMI_MAX_TERMS) {
                p->unavailable = 1;
                fail(m, "GMI term capacity limit");
                continue;
            }
            if (p->count == p->capacity) {
                p->capacity = p->capacity ? p->capacity * 2 : 16;
                GmiTerm *grown = realloc(p->term, p->capacity * sizeof *grown);
                if (!grown)
                    exit(70);
                p->term = grown;
            }
            GmiTerm *term = &p->term[p->count++];
            memset(term, 0, sizeof *term);
            term->kind = (enum GmiTermKind)kind;
            term->index = j; term->lower = lower; term->integer = integer;
            mpq_inits(term->coefficient, term->anchor, term->gamma, NULL);
            gmi_parse_q(m, p, term->coefficient, v[2]);
            gmi_parse_q(m, p, term->anchor, v[7]);
        } else if (!strcmp(v[0], "EQ") && n == 3) {
            size_t r = lookup(&m->rows, v[1]);
            if (r == SIZE_MAX || m->r[r].type == 'N' || p->equality_seen[r]) {
                fail(m, "unknown/duplicate/nonconstraint equality multiplier");
                continue;
            }
            p->equality_seen[r] = 1;
            gmi_parse_q(m, p, p->equality[r], v[2]);
        } else if (!strcmp(v[0], "FIX") && n == 3) {
            size_t j = lookup(&m->cols, v[1]);
            if (j == SIZE_MAX || p->fixed_seen[j]) {
                fail(m, "unknown/duplicate fixed-box multiplier");
                continue;
            }
            p->fixed_seen[j] = 1;
            gmi_parse_q(m, p, p->fixed[j], v[2]);
        } else if (!strcmp(v[0], "PUB_RHS") && n == 2 && !p->rhs_seen) {
            p->rhs_seen = 1;
            if (!gmi_parse_q(m, p, p->published_rhs, v[1]) && !gmi_finite_binary(p->published_rhs))
                fail(m, "published rhs is not an exact finite binary64 value");
        } else if (!strcmp(v[0], "PUB") && n == 3) {
            size_t j = lookup(&m->cols, v[1]);
            if (j == SIZE_MAX || p->published_seen[j]) {
                fail(m, "unknown/duplicate published column");
                continue;
            }
            p->published_seen[j] = 1;
            if (!gmi_parse_q(m, p, p->published[j], v[2]) && !gmi_finite_binary(p->published[j]))
                fail(m, "published coefficient is not an exact finite binary64 value");
        } else if (!strcmp(v[0], "END") && n == 1) {
            ended = 1;
        } else {
            fail(m, "unknown/duplicate/malformed GMI record");
        }
    }
    if (ferror(f))
        fail(m, "GMI certificate read error");
    fclose(f);
    if (!header || !ended || !p->scope_seen || !p->basic_seen || !p->beta_seen || !p->rhs_seen)
        fail(m, "incomplete GMI certificate");
    for (size_t j = 0; j < m->nc; ++j)
        if (!p->box[j].seen || !p->published_seen[j])
            fail(m, "incomplete box/published original-column vector");
}

static void gmi_check_premises(Model *m, GmiProof *p) {
    if (!m->c[p->basic].integer) {
        fail(m, "basic native variable is not an original integer column");
        return;
    }
    mpq_t lo, hi, scratch;
    mpq_inits(lo, hi, scratch, NULL);
    for (size_t k = 0; k < p->count && !m->error[0]; ++k) {
        GmiTerm *term = &p->term[k];
        if (term->kind == GMI_SOURCE_ROW) {
            int lf, hf;
            row_bounds(&m->r[term->index], lo, hi, &lf, &hf, scratch);
            if (!(term->lower ? lf : hf) || mpq_cmp(term->anchor, term->lower ? lo : hi))
                fail(m, "row slack anchor does not match an exact original finite row side");
            if (term->integer)
                fail(m, "row slacks must be continuous in this proof version");
        } else {
            Col *c = &m->c[term->index];
            GmiBox *box = &p->box[term->index];
            int source_side = (term->lower ? c->lf : c->hf) &&
                !mpq_cmp(term->anchor, term->lower ? c->lo : c->hi);
            int box_side = (term->lower ? box->lower_finite : box->upper_finite) &&
                !mpq_cmp(term->anchor, term->lower ? box->lower : box->upper);
            if (!source_side && !box_side)
                fail(m, "variable/bound anchor is not an exact source or recorded-box side");
            if (term->integer && (term->kind != GMI_VARIABLE || !c->integer ||
                                  mpz_cmp_ui(mpq_denref(term->anchor), 1)))
                fail(m, "native integer distance requires original integer column and integer anchor");
        }
    }
    for (size_t r = 0; r < m->nr && !m->error[0]; ++r) {
        if (!p->equality_seen[r])
            continue;
        int lf, hf;
        row_bounds(&m->r[r], lo, hi, &lf, &hf, scratch);
        if (!lf || !hf || mpq_cmp(lo, hi))
            fail(m, "EQ multiplier does not name an exact original equality");
    }
    for (size_t j = 0; j < m->nc && !m->error[0]; ++j)
        if (p->fixed_seen[j] && (!p->box[j].lower_finite || !p->box[j].upper_finite ||
                                 mpq_cmp(p->box[j].lower, p->box[j].upper)))
            fail(m, "FIX multiplier requires an exactly fixed recorded-box variable");
    p->premises_valid = !m->error[0];
    mpq_clears(lo, hi, scratch, NULL);
}

/* Aggregate row coefficients through the independent parser's A^T y pass.
 * No dense matrix or solver tableau/basis API is involved. */
static void gmi_reset_linear(Model *m) {
    for (size_t r = 0; r < m->nr; ++r)
        mpq_set_ui(m->r[r].dual, 0, 1);
    for (size_t j = 0; j < m->nc; ++j)
        mpq_set_ui(m->c[j].residual, 0, 1);
    m->farkas_enabled = 1;
    m->dual_enabled = 0;
}

static void gmi_add_affine(Model *m, GmiProof *p, const GmiTerm *term, const mpq_t weight,
                       mpq_t constant, mpq_t scratch) {
    if (term->kind == GMI_SOURCE_ROW) {
        if (term->lower)
            mpq_add(m->r[term->index].dual, m->r[term->index].dual, weight);
        else
            mpq_sub(m->r[term->index].dual, m->r[term->index].dual, weight);
        gmi_bounded_q(m, p, m->r[term->index].dual);
    } else {
        if (term->lower)
            mpq_add(m->c[term->index].residual, m->c[term->index].residual, weight);
        else
            mpq_sub(m->c[term->index].residual, m->c[term->index].residual, weight);
        gmi_bounded_q(m, p, m->c[term->index].residual);
    }
    mpq_mul(scratch, weight, term->anchor);
    if (term->lower)
        mpq_sub(constant, constant, scratch);
    else
        mpq_add(constant, constant, scratch);
    gmi_bounded_q(m, p, constant);
}

static int gmi_guard_linear(Model *m, GmiProof *p) {
    size_t bits = 0;
    for (size_t r = 0; r < m->nr && !m->error[0]; ++r) {
        if (gmi_bounded_q(m, p, m->r[r].dual))
            return -1;
        bits += mpz_sizeinbase(mpq_numref(m->r[r].dual), 2) +
                mpz_sizeinbase(mpq_denref(m->r[r].dual), 2);
    }
    for (size_t j = 0; j < m->nc && !m->error[0]; ++j) {
        if (gmi_bounded_q(m, p, m->c[j].residual))
            return -1;
        bits += mpz_sizeinbase(mpq_numref(m->c[j].residual), 2) +
                mpz_sizeinbase(mpq_denref(m->c[j].residual), 2);
    }
    if (bits > (1U << 20)) {
        p->unavailable = 1;
        fail(m, "GMI aggregate linear rational capacity limit");
    }
    return m->error[0] ? -1 : 0;
}

static void gmi_check_identity(Model *m, GmiProof *p, const char *model_path) {
    gmi_reset_linear(m);
    mpq_set_ui(m->c[p->basic].residual, 1, 1);
    mpq_t constant, scratch, lo, hi;
    mpq_inits(constant, scratch, lo, hi, NULL);
    mpq_neg(constant, p->beta);
    for (size_t k = 0; k < p->count && !m->error[0]; ++k)
        gmi_add_affine(m, p, &p->term[k], p->term[k].coefficient, constant, scratch);
    for (size_t r = 0; r < m->nr && !m->error[0]; ++r) {
        if (!p->equality_seen[r])
            continue;
        mpq_sub(m->r[r].dual, m->r[r].dual, p->equality[r]);
        int lf, hf;
        row_bounds(&m->r[r], lo, hi, &lf, &hf, scratch);
        mpq_mul(scratch, p->equality[r], lo);
        mpq_add(constant, constant, scratch);
        gmi_bounded_q(m, p, constant);
        gmi_bounded_q(m, p, m->r[r].dual);
    }
    for (size_t j = 0; j < m->nc && !m->error[0]; ++j) {
        if (!p->fixed_seen[j])
            continue;
        mpq_sub(m->c[j].residual, m->c[j].residual, p->fixed[j]);
        mpq_mul(scratch, p->fixed[j], p->box[j].lower);
        mpq_add(constant, constant, scratch);
        gmi_bounded_q(m, p, constant);
        gmi_bounded_q(m, p, m->c[j].residual);
    }
    if (m->error[0] || gmi_bounded_q(m, p, constant))
        goto done;
    if (mpq_sgn(constant)) {
        fail(m, "tableau identity has a nonzero exact constant residual");
        goto done;
    }
    if (!gmi_guard_linear(m, p))
        read_mps(m, model_path, 1);
    for (size_t j = 0; j < m->nc && !m->error[0]; ++j) {
        if (gmi_bounded_q(m, p, m->c[j].residual))
            break;
        if (mpq_sgn(m->c[j].residual))
            fail(m, "tableau identity has a nonzero exact original-column residual for %s",
                 m->c[j].name);
    }
    p->identity_valid = !m->error[0];
done:
    mpq_clears(constant, scratch, lo, hi, NULL);
}

/* Deliberately independent formula implementation; no gmi.c inclusion/link. */
static void gmi_derive_cut(Model *m, GmiProof *p, const char *model_path) {
    mpz_t whole;
    mpz_init(whole);
    mpq_t f, fractional_a, denominator, scratch, constant;
    mpq_inits(f, fractional_a, denominator, scratch, constant, NULL);
    mpz_fdiv_q(whole, mpq_numref(p->beta), mpq_denref(p->beta));
    mpq_set_z(scratch, whole);
    mpq_sub(f, p->beta, scratch);
    if (!mpq_sgn(f)) {
        p->unavailable = 1;
        fail(m, "integral native beta has no fractional GMI cut");
        goto done;
    }
    mpq_set_ui(denominator, 1, 1);
    mpq_sub(denominator, denominator, f);
    gmi_reset_linear(m);
    for (size_t k = 0; k < p->count && !m->error[0]; ++k) {
        GmiTerm *term = &p->term[k];
        if (term->integer) {
            mpz_fdiv_q(whole, mpq_numref(term->coefficient), mpq_denref(term->coefficient));
            mpq_set_z(scratch, whole);
            mpq_sub(fractional_a, term->coefficient, scratch);
            if (mpq_cmp(fractional_a, f) <= 0)
                mpq_set(term->gamma, fractional_a);
            else {
                mpq_set_ui(scratch, 1, 1);
                mpq_sub(scratch, scratch, fractional_a);
                mpq_div(term->gamma, scratch, denominator);
                mpq_mul(term->gamma, term->gamma, f);
            }
        } else if (mpq_sgn(term->coefficient) < 0) {
            mpq_div(term->gamma, term->coefficient, denominator);
            mpq_mul(term->gamma, term->gamma, f);
            mpq_neg(term->gamma, term->gamma);
        } else
            mpq_set(term->gamma, term->coefficient);
        if (gmi_bounded_q(m, p, term->gamma))
            break;
        gmi_add_affine(m, p, term, term->gamma, constant, scratch);
    }
    if (m->error[0])
        goto done;
    mpq_sub(p->exact_rhs, f, constant);
    if (!gmi_guard_linear(m, p))
        read_mps(m, model_path, 1);
    for (size_t j = 0; j < m->nc && !m->error[0]; ++j) {
        if (gmi_bounded_q(m, p, m->c[j].residual))
            break;
        mpq_set(p->exact[j], m->c[j].residual);
    }
    gmi_bounded_q(m, p, p->exact_rhs);
done:
    mpq_clears(f, fractional_a, denominator, scratch, constant, NULL);
    mpz_clear(whole);
}

static void gmi_check_binary_weakening(Model *m, GmiProof *p) {
    mpq_t difference, term, support;
    mpq_inits(difference, term, support, NULL);
    for (size_t j = 0; j < m->nc && !m->error[0]; ++j) {
        int original_sign = mpq_sgn(p->exact[j]), published_sign = mpq_sgn(p->published[j]);
        if (original_sign != published_sign) {
            fail(m, "published coefficient deleted, introduced, or reversed an exact coefficient");
            break;
        }
        mpq_sub(difference, p->published[j], p->exact[j]);
        int sign = mpq_sgn(difference);
        GmiBox *box = &p->box[j];
        if ((sign > 0 && !box->lower_finite) || (sign < 0 && !box->upper_finite)) {
            p->unavailable = 1;
            fail(m, "published coefficient error requires a missing recorded-box side");
            break;
        }
        if (sign) {
            mpq_mul(term, difference, sign > 0 ? box->lower : box->upper);
            mpq_add(p->error_lower, p->error_lower, term);
            gmi_bounded_q(m, p, p->error_lower);
        }
    }
    if (!m->error[0]) {
        mpq_add(p->maximum_rhs, p->exact_rhs, p->error_lower);
        mpq_sub(p->rounding_slack, p->maximum_rhs, p->published_rhs);
        if (gmi_bounded_q(m, p, p->rounding_slack))
            goto done;
        if (mpq_sgn(p->rounding_slack) < 0)
            fail(m, "published rhs exceeds exact box-compensated maximum rhs");
        else
            p->cut_valid = 1;
    }
    if (p->cut_valid) {
        int available = 1;
        for (size_t j = 0; j < m->nc; ++j) {
            int sign = mpq_sgn(p->exact[j]);
            GmiBox *box = &p->box[j];
            if ((sign > 0 && !box->upper_finite) || (sign < 0 && !box->lower_finite)) {
                available = 0;
                break;
            }
            if (sign) {
                mpq_mul(term, p->exact[j], sign > 0 ? box->upper : box->lower);
                mpq_add(support, support, term);
                if (gmi_bounded_q(m, p, support)) {
                    available = 0;
                    break;
                }
            }
        }
        p->mip_box_infeasible = available && mpq_cmp(support, p->exact_rhs) < 0;
        if (m->error[0])
            p->cut_valid = 0;
    }
done:
    mpq_clears(difference, term, support, NULL);
}

static int gmi_check_model(const char *model_path, const char *certificate_path) {
    Model m = {0};
    m.objective = SIZE_MAX;
    mpq_init(m.obj);
    GmiProof proof;
    gmi_proof_init(&proof);
    char model_hash[65] = "", certificate_hash[65] = "", observed_hash[65];
    if (lab_sha256(model_path, model_hash) || lab_sha256(certificate_path, certificate_hash))
        fail(&m, "cannot hash GMI model/certificate input");
    if (!m.error[0])
        read_mps(&m, model_path, 0);
    if (!m.error[0] && (m.nc > GMI_MAX_COLUMNS || m.nr > GMI_MAX_ROWS || m.nnz > (1U << 20))) {
        proof.unavailable = 1;
        fail(&m, "GMI original-model dimension limit");
    }
    if (!m.error[0]) {
        gmi_proof_allocate(&proof, &m);
        gmi_read_certificate(&m, &proof, certificate_path);
    }
    if (!m.error[0])
        gmi_check_premises(&m, &proof);
    if (!m.error[0])
        gmi_check_identity(&m, &proof, model_path);
    if (!m.error[0])
        gmi_derive_cut(&m, &proof, model_path);
    if (!m.error[0])
        gmi_check_binary_weakening(&m, &proof);
    if (!m.error[0] &&
        (lab_sha256(model_path, observed_hash) || strcmp(model_hash, observed_hash)))
        fail(&m, "model input changed during GMI checking");
    if (!m.error[0] &&
        (lab_sha256(certificate_path, observed_hash) || strcmp(certificate_hash, observed_hash)))
        fail(&m, "certificate input changed during GMI checking");
    int valid = proof.cut_valid && !m.error[0];
    printf("{\"schema\":1,\"status\":");
    lab_json_string(stdout, m.unsupported ? "unsupported" : proof.unavailable ? "unavailable"
                                  : valid ? "valid" : "invalid");
    printf(",\"certificate_scope\":\"original_mip_in_recorded_box\",\"certifies_cut_validity\":%s,"
           "\"native_integrality_premises\":%s,\"exact_tableau_identity\":%s,"
           "\"certifies_mip_in_box_infeasibility\":%s,"
           "\"certifies_original_lp_infeasibility\":false,"
           "\"certifies_original_mip_infeasibility\":false,\"certifies_full_tree\":false,"
           "\"recorded_box_proven_implied_by_original_model\":false,"
           "\"checked_primal_witness\":false,\"certifies_optimality\":false,\"message\":",
           valid ? "true" : "false", proof.premises_valid ? "true" : "false",
           proof.identity_valid ? "true" : "false",
           valid && proof.mip_box_infeasible ? "true" : "false");
    lab_json_string(stdout, m.error);
    printf(",\"model_sha256\":");
    if (*model_hash) lab_json_string(stdout, model_hash); else printf("null");
    printf(",\"certificate_sha256\":");
    if (*certificate_hash) lab_json_string(stdout, certificate_hash); else printf("null");
    printf(",\"exact_cut_rhs\":");
    if (valid) json_q(proof.exact_rhs); else printf("null");
    printf(",\"exact_cut\":{");
    if (valid)
        for (size_t j = 0; j < m.nc; ++j) {
            if (j) putchar(',');
            lab_json_string(stdout, m.c[j].name); putchar(':'); json_q(proof.exact[j]);
        }
    printf("},\"rounding_error_lower\":");
    if (valid) json_q(proof.error_lower); else printf("null");
    printf(",\"maximum_safe_published_rhs\":");
    if (valid) json_q(proof.maximum_rhs); else printf("null");
    printf(",\"rounding_slack\":");
    if (valid) json_q(proof.rounding_slack); else printf("null");
    printf("}\n");
    int status = m.unsupported ? 3 : proof.unavailable ? 4 : valid ? 0 : 1;
    gmi_proof_clear(&proof, &m);
    cleanup(&m);
    return status;
}

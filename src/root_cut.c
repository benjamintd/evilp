#include "root_cut.h"
#include "common.h"
#include "exact.h"
#include "gmi.h"

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const LabRootCutLimits *limits;
    size_t bits;
} CutGuard;

typedef struct {
    char *data;
    size_t length, capacity;
    CutGuard *guard;
} CutText;

static LabRootCutStatus cut_tick(const CutGuard *g) {
    return lab_now() >= g->limits->native.deadline ? LAB_ROOT_CUT_DEADLINE : LAB_ROOT_CUT_OK;
}

static LabRootCutStatus cut_check(CutGuard *g, const mpq_t q) {
    size_t a = mpz_sizeinbase(mpq_numref(q), 2);
    size_t b = mpz_sizeinbase(mpq_denref(q), 2);
    if (mpz_sgn(mpq_denref(q)) <= 0)
        return LAB_ROOT_CUT_INVALID;
    if (a > g->limits->native.max_bits || b > g->limits->native.max_bits ||
        a > SIZE_MAX - b || a + b > g->limits->native.max_total_bits - g->bits)
        return LAB_ROOT_CUT_RESOURCE;
    g->bits += a + b;
    return cut_tick(g);
}

#define CUT_TRY(expr) do { status = (expr); if (status != LAB_ROOT_CUT_OK) goto done; } while (0)

static mpq_t *cut_vector(size_t n) {
    if (n > SIZE_MAX / sizeof(mpq_t))
        return NULL;
    mpq_t *q = calloc(n ? n : 1, sizeof *q);
    if (q)
        for (size_t j = 0; j < n; ++j)
            mpq_init(q[j]);
    return q;
}

static void cut_vector_free(mpq_t *q, size_t n) {
    if (q)
        for (size_t j = 0; j < n; ++j)
            mpq_clear(q[j]);
    free(q);
}

void lab_root_cut_batch_init(LabRootCutBatch *b) {
    memset(b, 0, sizeof *b);
}

void lab_root_cut_batch_free(LabRootCutBatch *b) {
    if (!b)
        return;
    for (size_t i = 0; i < b->count; ++i) {
        free(b->cut[i].coefficient);
        free(b->cut[i].certificate);
    }
    free(b->cut);
    if (b->box)
        for (size_t j = 0; j < b->columns; ++j)
            mpq_clears(b->box[j].lower, b->box[j].upper, NULL);
    free(b->box);
    memset(b, 0, sizeof *b);
}

LabRootCutLimits lab_root_cut_limits(double deadline) {
    return (LabRootCutLimits){
        .native = {.max_rows = 4, .max_terms = 4096, .max_basis_dimension = 1024,
                   .max_basis_nnz = 100000, .max_canonical_nnz = 100000,
                   .max_bits = 8192, .max_total_bits = 16000000, .deadline = deadline},
        .max_certificate_bytes = 1U << 24};
}

const char *lab_root_cut_status_name(LabRootCutStatus s) {
    switch (s) {
    case LAB_ROOT_CUT_OK: return "ok";
    case LAB_ROOT_CUT_NONE: return "no_usable_cut";
    case LAB_ROOT_CUT_STALE: return "stale_or_unsolved_session";
    case LAB_ROOT_CUT_SOURCE_UNAVAILABLE: return "source_unavailable";
    case LAB_ROOT_CUT_RESOURCE: return "resource_limit";
    case LAB_ROOT_CUT_DEADLINE: return "deadline";
    case LAB_ROOT_CUT_INVALID: return "invalid_input";
    case LAB_ROOT_CUT_EXACT_FAILURE: return "exact_derivation_unavailable";
    }
    return "invalid_status";
}

static LabRootCutStatus native_status(LabNativeStatus s) {
    switch (s) {
    case LAB_NATIVE_OK: return LAB_ROOT_CUT_OK;
    case LAB_NATIVE_NO_ROWS: return LAB_ROOT_CUT_NONE;
    case LAB_NATIVE_STALE: return LAB_ROOT_CUT_STALE;
    case LAB_NATIVE_SOURCE_UNAVAILABLE: case LAB_NATIVE_FREE_SPLIT:
        return LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
    case LAB_NATIVE_RESOURCE: return LAB_ROOT_CUT_RESOURCE;
    case LAB_NATIVE_DEADLINE: return LAB_ROOT_CUT_DEADLINE;
    case LAB_NATIVE_INVALID: return LAB_ROOT_CUT_INVALID;
    case LAB_NATIVE_EXACT_FAILURE: case LAB_NATIVE_IDENTITY_FAILURE:
        return LAB_ROOT_CUT_EXACT_FAILURE;
    }
    return LAB_ROOT_CUT_INVALID;
}

static LabRootCutStatus gmi_status(LabGmiStatus s) {
    switch (s) {
    case LAB_GMI_OK: return LAB_ROOT_CUT_OK;
    case LAB_GMI_INTEGRAL: case LAB_GMI_UNBOUNDED_ROUNDING:
    case LAB_GMI_NONFINITE_ROUNDING: return LAB_ROOT_CUT_NONE;
    case LAB_GMI_RESOURCE: return LAB_ROOT_CUT_RESOURCE;
    case LAB_GMI_DEADLINE: return LAB_ROOT_CUT_DEADLINE;
    case LAB_GMI_INVALID: return LAB_ROOT_CUT_INVALID;
    }
    return LAB_ROOT_CUT_INVALID;
}

static double cut_now(void *context) {
    (void)context;
    return lab_now();
}

static int cut_name(const char *s) {
    if (!s || !*s || strlen(s) > 1024)
        return 0;
    for (; *s; ++s)
        if (isspace((unsigned char)*s) || (unsigned char)*s < 32)
            return 0;
    return 1;
}

static void source_row_side(const LabModel *m, size_t i, int upper, mpq_t q) {
    if (m->exact)
        mpq_set(q, upper ? m->exact->row_upper[i] : m->exact->row_lower[i]);
    else
        mpq_set_d(q, upper ? m->row_upper[i] : m->row_lower[i]);
}

/* The source identity is established by the exporter. Bind its exact matrix
 * and integrality to the caller's original source, and check every original
 * row premise used in that identity/cut. Relaxed unused rows may differ. */
static LabRootCutStatus bind_source(CutGuard *g, const LabModel *m,
                                   const LabNativeBatch *b) {
    if (m->rows != b->source_rows || m->cols != b->columns ||
        m->nnz != b->source_matrix.nnz || (m->parsed_decimal && !m->exact) ||
        !m->col_start || !m->integer || !m->col_name || !m->row_name ||
        !m->row_lower || !m->row_upper ||
        (m->nnz && (!m->row_index || !m->value)))
        return LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
    if (m->exact && (m->exact->rows != m->rows || m->exact->cols != m->cols ||
                     m->exact->nnz != m->nnz))
        return LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
    for (size_t j = 0; j <= m->cols; ++j)
        if (m->col_start[j] != b->source_matrix.start[j])
            return LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
    for (size_t j = 0; j < m->cols; ++j)
        if (!cut_name(m->col_name[j]) || m->integer[j] != b->integer[j])
            return LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
    for (size_t i = 0; i < m->rows; ++i)
        if (!cut_name(m->row_name[i]))
            return LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
    LabRootCutStatus status = LAB_ROOT_CUT_OK;
    mpq_t q;
    mpq_init(q);
    for (size_t k = 0; k < m->nnz; ++k) {
        if (m->row_index[k] != b->source_matrix.row[k] || !isfinite(m->value[k])) {
            status = LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
            goto done;
        }
        if (m->exact) mpq_set(q, m->exact->value[k]);
        else mpq_set_d(q, m->value[k]);
        CUT_TRY(cut_check(g, q));
        if (mpq_cmp(q, b->source_matrix.value[k])) {
            status = LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
            goto done;
        }
    }
    for (size_t r = 0; r < b->count; ++r) {
        const LabNativeTableau *t = &b->tableau[r];
        if (!t->basis_identity_verified || !t->native_identity_verified ||
            t->basic_column >= m->cols || !m->integer[t->basic_column]) {
            status = LAB_ROOT_CUT_EXACT_FAILURE;
            goto done;
        }
        for (size_t k = 0; k < t->terms; ++k) {
            const LabNativeTerm *term = &t->term[k];
            if (term->kind != LAB_NATIVE_ROW)
                continue;
            if (term->index >= m->rows ||
                !isfinite(term->upper ? m->row_upper[term->index] : m->row_lower[term->index])) {
                status = LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
                goto done;
            }
            source_row_side(m, term->index, term->upper, q);
            if (mpq_cmp(q, term->anchor)) {
                status = LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
                goto done;
            }
        }
        for (size_t k = 0; k < t->equalities; ++k) {
            size_t i = t->equality[k].index;
            if (i >= m->rows || !isfinite(m->row_lower[i]) || !isfinite(m->row_upper[i])) {
                status = LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
                goto done;
            }
            source_row_side(m, i, 0, q);
            if (mpq_cmp(q, b->source_row_box[i].lower)) {
                status = LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
                goto done;
            }
            source_row_side(m, i, 1, q);
            if (mpq_cmp(q, b->source_row_box[i].lower)) {
                status = LAB_ROOT_CUT_SOURCE_UNAVAILABLE;
                goto done;
            }
        }
    }
done:
    mpq_clear(q);
    return status;
}

static LabRootCutStatus text_add(CutText *s, const char *format, ...) {
    LabRootCutStatus status = cut_tick(s->guard);
    if (status != LAB_ROOT_CUT_OK)
        return status;
    va_list args, copy;
    va_start(args, format);
    va_copy(copy, args);
    int needed = vsnprintf(NULL, 0, format, copy);
    va_end(copy);
    if (needed < 0 || (size_t)needed > s->guard->limits->max_certificate_bytes - s->length) {
        va_end(args);
        return LAB_ROOT_CUT_RESOURCE;
    }
    size_t total = s->length + (size_t)needed + 1;
    if (total > s->capacity) {
        size_t capacity = s->capacity ? 2 * s->capacity : 1024;
        if (capacity < total)
            capacity = total;
        if (capacity > s->guard->limits->max_certificate_bytes + 1)
            capacity = s->guard->limits->max_certificate_bytes + 1;
        char *data = realloc(s->data, capacity);
        if (!data) {
            va_end(args);
            return LAB_ROOT_CUT_RESOURCE;
        }
        s->data = data;
        s->capacity = capacity;
    }
    int written = vsnprintf(s->data + s->length, s->capacity - s->length, format, args);
    va_end(args);
    if (written != needed)
        return LAB_ROOT_CUT_INVALID;
    s->length += (size_t)needed;
    return LAB_ROOT_CUT_OK;
}

static LabRootCutStatus text_q(CutText *s, const mpq_t q) {
    LabRootCutStatus status = cut_check(s->guard, q);
    if (status != LAB_ROOT_CUT_OK)
        return status;
    char *text = mpq_get_str(NULL, 10, q);
    status = text_add(s, "%s", text);
    free(text);
    return status;
}

static LabRootCutStatus certificate(CutGuard *g, const LabModel *m,
                                    const LabNativeBatch *b, const LabNativeTableau *t,
                                    const LabGmiBinaryCut *binary, char **out) {
    LabRootCutStatus status = LAB_ROOT_CUT_OK;
    CutText text = {.guard = g};
    mpq_t q;
    mpq_init(q);
    CUT_TRY(text_add(&text, "MILPLAB_GMI_1\nSCOPE ORIGINAL_MIP_IN_BOX\nBASIC %s\nBETA ",
                     m->col_name[t->basic_column]));
    CUT_TRY(text_q(&text, t->beta));
    CUT_TRY(text_add(&text, "\n"));
    for (size_t j = 0; j < b->columns; ++j) {
        CUT_TRY(text_add(&text, "BOX %s ", m->col_name[j]));
        if (b->box[j].has_lower) CUT_TRY(text_q(&text, b->box[j].lower));
        else CUT_TRY(text_add(&text, "*"));
        CUT_TRY(text_add(&text, " "));
        if (b->box[j].has_upper) CUT_TRY(text_q(&text, b->box[j].upper));
        else CUT_TRY(text_add(&text, "*"));
        CUT_TRY(text_add(&text, "\n"));
    }
    for (size_t k = 0; k < t->terms; ++k) {
        const LabNativeTerm *term = &t->term[k];
        CUT_TRY(text_add(&text, "TERM %zu ", k));
        CUT_TRY(text_q(&text, term->coefficient));
        CUT_TRY(text_add(&text, " %c %s %s %c ", term->is_integer ? 'I' : 'C',
                         term->kind == LAB_NATIVE_VAR ? "VAR" :
                         term->kind == LAB_NATIVE_ROW ? "ROW" : "BOUND",
                         term->kind == LAB_NATIVE_ROW ? m->row_name[term->index] :
                                                       m->col_name[term->index],
                         term->upper ? 'U' : 'L'));
        CUT_TRY(text_q(&text, term->anchor));
        CUT_TRY(text_add(&text, "\n"));
    }
    for (size_t k = 0; k < t->equalities; ++k) {
        CUT_TRY(text_add(&text, "EQ %s ", m->row_name[t->equality[k].index]));
        CUT_TRY(text_q(&text, t->equality[k].multiplier));
        CUT_TRY(text_add(&text, "\n"));
    }
    for (size_t k = 0; k < t->fixed; ++k) {
        CUT_TRY(text_add(&text, "FIX %s ", m->col_name[t->fix[k].index]));
        CUT_TRY(text_q(&text, t->fix[k].multiplier));
        CUT_TRY(text_add(&text, "\n"));
    }
    CUT_TRY(text_add(&text, "PUB_RHS "));
    mpq_set_d(q, binary->rhs);
    CUT_TRY(text_q(&text, q));
    CUT_TRY(text_add(&text, "\n"));
    for (size_t j = 0; j < b->columns; ++j) {
        CUT_TRY(text_add(&text, "PUB %s ", m->col_name[j]));
        mpq_set_d(q, binary->coefficient[j]);
        CUT_TRY(text_q(&text, q));
        CUT_TRY(text_add(&text, "\n"));
    }
    CUT_TRY(text_add(&text, "END\n"));
    *out = text.data;
    text.data = NULL;
done:
    mpq_clear(q);
    free(text.data);
    return status;
}

static LabRootCutStatus derive_one(CutGuard *g, const LabModel *m, const LabNativeBatch *b,
                                   const LabNativeTableau *t, LabRootCut *out) {
    LabRootCutStatus status = LAB_ROOT_CUT_RESOURCE;
    mpq_t *a = cut_vector(t->terms), *mapped = cut_vector(b->columns);
    unsigned char *integer = calloc(t->terms ? t->terms : 1, 1);
    LabGmiBox *box = calloc(b->columns ? b->columns : 1, sizeof *box);
    size_t initialized_boxes = 0;
    LabGmiCut exact;
    LabGmiBinaryCut binary;
    lab_gmi_cut_init(&exact);
    lab_gmi_binary_init(&binary);
    mpq_t rhs, factor, product;
    mpq_inits(rhs, factor, product, NULL);
    char *proof = NULL;
    if (!a || !mapped || !integer || !box)
        goto done;
    for (size_t j = 0; j < b->columns; ++j) {
        lab_gmi_box_init(&box[j]);
        initialized_boxes++;
        box[j].has_lower = b->box[j].has_lower;
        box[j].has_upper = b->box[j].has_upper;
        mpq_set(box[j].lower, b->box[j].lower);
        mpq_set(box[j].upper, b->box[j].upper);
    }
    for (size_t k = 0; k < t->terms; ++k) {
        mpq_set(a[k], t->term[k].coefficient);
        integer[k] = (unsigned char)t->term[k].is_integer;
    }
    size_t terms = g->limits->native.max_terms;
    if (terms < b->columns)
        terms = b->columns;
    LabGmiBudget budget = {.max_terms = terms, .max_bits = g->limits->native.max_bits,
                           .max_total_bits = g->limits->native.max_total_bits,
                           .deadline = g->limits->native.deadline, .now = cut_now};
    CUT_TRY(gmi_status(lab_gmi_build(&exact, t->terms, (const mpq_t *)a, integer,
                                    t->beta, &budget)));
    mpq_set(rhs, exact.rhs);
    for (size_t k = 0; k < t->terms; ++k) {
        const LabNativeTerm *term = &t->term[k];
        mpq_set(factor, exact.coefficient[k]);
        if (term->upper)
            mpq_neg(factor, factor);
        if (term->kind == LAB_NATIVE_ROW) {
            for (size_t j = 0; j < b->columns; ++j)
                for (size_t p = b->source_matrix.start[j]; p < b->source_matrix.start[j + 1]; ++p)
                    if (b->source_matrix.row[p] == term->index) {
                        mpq_mul(product, factor, b->source_matrix.value[p]);
                        CUT_TRY(cut_check(g, product));
                        mpq_add(mapped[j], mapped[j], product);
                        CUT_TRY(cut_check(g, mapped[j]));
                    }
        } else {
            mpq_add(mapped[term->index], mapped[term->index], factor);
            CUT_TRY(cut_check(g, mapped[term->index]));
        }
        mpq_mul(product, factor, term->anchor);
        CUT_TRY(cut_check(g, product));
        mpq_add(rhs, rhs, product);
        CUT_TRY(cut_check(g, rhs));
        CUT_TRY(cut_tick(g));
    }
    CUT_TRY(gmi_status(lab_gmi_round_box(&binary, b->columns, (const mpq_t *)mapped,
                                       rhs, box, &budget)));
    CUT_TRY(certificate(g, m, b, t, &binary, &proof));
    out->columns = b->columns;
    out->basic_column = t->basic_column;
    out->coefficient = binary.coefficient;
    binary.coefficient = NULL;
    out->lower = binary.rhs;
    out->certificate = proof;
    proof = NULL;
done:
    free(proof);
    cut_vector_free(a, t->terms);
    cut_vector_free(mapped, b->columns);
    free(integer);
    for (size_t j = 0; j < initialized_boxes; ++j)
        lab_gmi_box_clear(&box[j]);
    free(box);
    mpq_clears(rhs, factor, product, NULL);
    lab_gmi_cut_clear(&exact);
    lab_gmi_binary_clear(&binary);
    return status;
}

LabRootCutStatus lab_root_cut_derive(LabLpSession *session, const LabModel *source,
                                   const LabRootCutLimits *limits, LabRootCutBatch *out) {
    if (!session || !source || !limits || !out || !limits->max_certificate_bytes ||
        limits->max_certificate_bytes > (1U << 24))
        return LAB_ROOT_CUT_INVALID;
    LabNativeBatch native;
    lab_native_batch_init(&native);
    LabRootCutBatch result;
    lab_root_cut_batch_init(&result);
    LabRootCutStatus status = native_status(lab_lp_session_export_native(session, &limits->native,
                                                                       &native));
    if (status != LAB_ROOT_CUT_OK)
        goto done;
    CutGuard guard = {.limits = limits};
    CUT_TRY(bind_source(&guard, source, &native));
    result.columns = native.columns;
    result.cut = calloc(native.count ? native.count : 1, sizeof *result.cut);
    if (!result.cut) {
        status = LAB_ROOT_CUT_RESOURCE;
        goto done;
    }
    for (size_t r = 0; r < native.count; ++r) {
        status = derive_one(&guard, source, &native, &native.tableau[r], &result.cut[result.count]);
        if (status == LAB_ROOT_CUT_NONE)
            continue;
        if (status != LAB_ROOT_CUT_OK)
            goto done;
        result.count++;
    }
    if (!result.count) {
        status = LAB_ROOT_CUT_NONE;
        goto done;
    }
    CUT_TRY(cut_tick(&guard));
    result.box = native.box;
    native.box = NULL;
    {
        LabRootCutBatch old = *out;
        *out = result;
        result = old;
    }
done:
    lab_native_batch_free(&native);
    lab_root_cut_batch_free(&result);
    return status;
}

#undef CUT_TRY

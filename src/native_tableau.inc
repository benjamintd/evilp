/* Included here to reuse source-rational Work helpers.
 * This file exports algebra; it does not add constraints or prove cut validity. */
typedef struct {
    const LabNativeLimits *limit;
    double deadline;
    size_t stored_bits;
} NativeGuard;
#ifdef LAB_NATIVE_TABLEAU_TESTING
static size_t native_fail_after = SIZE_MAX;
void lab_native_test_fail_after(size_t allocations) { native_fail_after = allocations; }
#endif
static void *native_alloc(size_t count, size_t size) {
    if (size && count > SIZE_MAX / size)
        return NULL;
#ifdef LAB_NATIVE_TABLEAU_TESTING
    if (!native_fail_after)
        return NULL;
    if (native_fail_after != SIZE_MAX)
        native_fail_after--;
#endif
    return calloc(count ? count : 1, size ? size : 1);
}
static mpq_t *native_qvector(size_t n) {
    mpq_t *q = native_alloc(n, sizeof(mpq_t));
    if (q)
        for (size_t i = 0; i < n; i++)
            mpq_init(q[i]);
    return q;
}
static void native_qfree(mpq_t *q, size_t n) {
    if (q)
        for (size_t i = 0; i < n; i++)
            mpq_clear(q[i]);
    free(q);
}
static LabNativeStatus native_tick(const NativeGuard *g) {
    return lab_now() >= g->deadline ? LAB_NATIVE_DEADLINE : LAB_NATIVE_OK;
}
static LabNativeStatus native_check(NativeGuard *g, const mpq_t q, int stored) {
    size_t a = mpz_sizeinbase(mpq_numref(q), 2), b = mpz_sizeinbase(mpq_denref(q), 2);
    if (mpz_sgn(mpq_denref(q)) <= 0)
        return LAB_NATIVE_INVALID;
    if (a > g->limit->max_bits || b > g->limit->max_bits || a > SIZE_MAX - b)
        return LAB_NATIVE_RESOURCE;
    if (stored) {
        if (a + b > g->limit->max_total_bits - g->stored_bits)
            return LAB_NATIVE_RESOURCE;
        g->stored_bits += a + b;
    }
    return native_tick(g);
}
static void native_matrix_free(LabNativeMatrix *m) {
    free(m->start); free(m->row); native_qfree(m->value, m->nnz);
    memset(m, 0, sizeof *m);
}
static int native_matrix_alloc(LabNativeMatrix *m, size_t rows, size_t cols, size_t nnz) {
    m->rows = rows; m->cols = cols; m->nnz = nnz;
    m->start = native_alloc(cols + 1, sizeof(size_t));
    m->row = native_alloc(nnz, sizeof(size_t));
    m->value = native_qvector(nnz);
    return m->start && m->row && m->value ? 0 : -1;
}
static LabNativeBox *native_boxes(size_t n) {
    LabNativeBox *box = native_alloc(n, sizeof *box);
    if (box)
        for (size_t i = 0; i < n; i++)
            mpq_inits(box[i].lower, box[i].upper, NULL);
    return box;
}
static void native_boxes_free(LabNativeBox *box, size_t n) {
    if (box)
        for (size_t i = 0; i < n; i++)
            mpq_clears(box[i].lower, box[i].upper, NULL);
    free(box);
}
void lab_native_batch_init(LabNativeBatch *b) { memset(b, 0, sizeof *b); }
void lab_native_batch_free(LabNativeBatch *b) {
    if (!b)
        return;
    for (size_t i = 0; i < b->count; i++) {
        LabNativeTableau *t = &b->tableau[i];
        for (size_t j = 0; j < t->terms; j++)
            mpq_clears(t->term[j].coefficient, t->term[j].anchor,
                       t->term[j].canonical_scale, t->term[j].canonical_offset, NULL);
        for (size_t j = 0; j < t->equalities; j++)
            mpq_clear(t->equality[j].multiplier);
        for (size_t j = 0; j < t->fixed; j++)
            mpq_clear(t->fix[j].multiplier);
        free(t->term); free(t->equality); free(t->fix);
        native_qfree(t->weights, t->weight_count);
        mpq_clears(t->beta, t->signed_basic_scale, t->basic_shift, NULL);
    }
    free(b->tableau);
    native_boxes_free(b->box, b->columns);
    native_boxes_free(b->source_row_box, b->source_rows);
    free(b->integer); free(b->basis); free(b->row_id);
    native_qfree(b->rhs, b->dimension);
    native_qfree(b->basic_values, b->dimension);
    native_matrix_free(&b->source_matrix);
    native_matrix_free(&b->canonical_matrix);
    native_matrix_free(&b->basis_matrix);
    memset(b, 0, sizeof *b);
}
const char *lab_native_status_name(LabNativeStatus s) {
    switch (s) {
    case LAB_NATIVE_OK: return "ok";
    case LAB_NATIVE_NO_ROWS: return "no_fractional_integer_basic";
    case LAB_NATIVE_STALE: return "stale_or_unsolved_session";
    case LAB_NATIVE_SOURCE_UNAVAILABLE: return "exact_source_unavailable";
    case LAB_NATIVE_FREE_SPLIT: return "unsupported_free_split";
    case LAB_NATIVE_RESOURCE: return "resource_limit";
    case LAB_NATIVE_DEADLINE: return "deadline";
    case LAB_NATIVE_INVALID: return "invalid_input";
    case LAB_NATIVE_EXACT_FAILURE: return "exact_reconstruction_unavailable";
    case LAB_NATIVE_IDENTITY_FAILURE: return "exact_identity_rejected";
    }
    return "invalid_status";
}
#define NATIVE_TRY(expression) do { status = (expression); if (status != LAB_NATIVE_OK) goto done; } while (0)

static LabNativeStatus native_shift(const Work *w, size_t j, mpq_t q) {
    if (w->transform[j] == 'F') { mpq_set_ui(q,0,1); return LAB_NATIVE_OK; }
    return lab_source_exact_column_value(w->model,j,w->transform[j]=='U',w->shift[j],q)
        ? LAB_NATIVE_SOURCE_UNAVAILABLE : LAB_NATIVE_OK;
}
static LabNativeStatus native_basis_value(const Work *w, size_t j, size_t entry, mpq_t q) {
    const LabModel *m=w->model; const Column *c=&w->a[j]; size_t i=c->row[entry];
    if (j>=w->structural || w->source_row[i]==SIZE_MAX) { mpq_set_d(q,c->value[entry]); return LAB_NATIVE_OK; }
    size_t col=w->original[j],source=w->source_row[i],lo=m->col_start[col],hi=m->col_start[col+1];
    while(lo<hi){size_t mid=lo+(hi-lo)/2;if(m->row_index[mid]<source)lo=mid+1;else hi=mid;}
    if(lo==m->col_start[col+1]||m->row_index[lo]!=source) return LAB_NATIVE_EXACT_FAILURE;
    if(lab_source_exact_coefficient(m,lo,q)) return LAB_NATIVE_SOURCE_UNAVAILABLE;
    mpq_t scale;mpq_init(scale);mpq_set_d(scale,w->sign[j]);mpq_mul(q,q,scale);
    mpq_set_d(scale,w->row_scale[i]);mpq_mul(q,q,scale);mpq_clear(scale);return LAB_NATIVE_OK;
}
static LabNativeStatus native_rhs(const Work *w, NativeGuard *g, mpq_t *rhs) {
    const LabModel *m=w->model; LabNativeStatus status=LAB_NATIVE_OK;
    mpq_t *activity=native_qvector(m->rows);if(!activity)return LAB_NATIVE_RESOURCE;
    mpq_t shift,value,side,scale;mpq_inits(shift,value,side,scale,NULL);
    for(size_t j=0;j<m->cols;j++) {
        NATIVE_TRY(native_shift(w,j,shift));
        for(size_t k=m->col_start[j];k<m->col_start[j+1];k++) {
            if(lab_source_exact_coefficient(m,k,value)){status=LAB_NATIVE_SOURCE_UNAVAILABLE;goto done;}
            mpq_mul(value,value,shift);NATIVE_TRY(native_check(g,value,0));
            mpq_add(activity[m->row_index[k]],activity[m->row_index[k]],value);
            NATIVE_TRY(native_check(g,activity[m->row_index[k]],0));
        }
    }
    for(size_t i=0;i<w->m;i++) if(w->source_row[i]!=SIZE_MAX) {
        size_t source=w->source_row[i];
        if(lab_source_exact_row_bound(m,source,w->row_upper[i],side)){status=LAB_NATIVE_SOURCE_UNAVAILABLE;goto done;}
        mpq_sub(rhs[i],side,activity[source]);mpq_set_d(scale,w->row_scale[i]);mpq_mul(rhs[i],rhs[i],scale);
        NATIVE_TRY(native_check(g,rhs[i],0));
    }
    for(size_t j=0;j<m->cols;j++) {
        NATIVE_TRY(native_shift(w,j,shift));
        for(size_t i=0;i<w->m;i++) {
            int upper=w->row_id[i]==w->upper_row[j];
            if(!upper&&w->row_id[i]!=w->lower_row[j])continue;
            if(lab_source_exact_column_bound(m,j,upper,side)){status=LAB_NATIVE_SOURCE_UNAVAILABLE;goto done;}
            if(upper)mpq_sub(rhs[i],side,shift);else mpq_sub(rhs[i],shift,side);
            mpq_set_d(scale,w->equation_scale[i]);mpq_mul(rhs[i],rhs[i],scale);
            NATIVE_TRY(native_check(g,rhs[i],0));
        }
        NATIVE_TRY(native_tick(g));
    }
done:
    mpq_clears(shift,value,side,scale,NULL);native_qfree(activity,m->rows);return status;
}

/* A coalesced source coefficient may be nonzero even when its binary
 * sum is zero. Preserve these entries in the exact canonical matrix rather
 * than inheriting the floating matrix's omitted support. */
static LabNativeStatus native_missing_zero_count(const Work *w,NativeGuard *g,
    const LabNativeBatch *b,size_t j,size_t *count) {
    *count=0;if(j>=w->structural)return LAB_NATIVE_OK;
    size_t source=w->original[j];
    for(size_t k=w->model->col_start[source];k<w->model->col_start[source+1];k++) {
        if(w->model->value[k]!=0 || !mpq_sgn(b->source_matrix.value[k]))continue;
        for(size_t i=0;i<w->m;i++)if(w->source_row[i]==w->model->row_index[k])(*count)++;
        if(native_tick(g)!=LAB_NATIVE_OK)return LAB_NATIVE_DEADLINE;
    }
    return LAB_NATIVE_OK;
}
static LabNativeStatus native_missing_zeros(const Work *w,NativeGuard *g,
    LabNativeBatch *b,size_t j,size_t *at) {
    if(j>=w->structural)return LAB_NATIVE_OK;
    size_t source=w->original[j];LabNativeStatus status=LAB_NATIVE_OK;
    mpq_t scale;mpq_init(scale);
    for(size_t k=w->model->col_start[source];k<w->model->col_start[source+1];k++) {
        if(w->model->value[k]!=0 || !mpq_sgn(b->source_matrix.value[k]))continue;
        for(size_t i=0;i<w->m;i++)if(w->source_row[i]==w->model->row_index[k]) {
            b->canonical_matrix.row[*at]=i;
            mpq_set(b->canonical_matrix.value[*at],b->source_matrix.value[k]);
            mpq_set_d(scale,w->sign[j]);mpq_mul(b->canonical_matrix.value[*at],b->canonical_matrix.value[*at],scale);
            mpq_set_d(scale,w->row_scale[i]);mpq_mul(b->canonical_matrix.value[*at],b->canonical_matrix.value[*at],scale);
            NATIVE_TRY(native_check(g,b->canonical_matrix.value[(*at)++],1));
        }
        NATIVE_TRY(native_tick(g));
    }
done:mpq_clear(scale);return status;
}

static LabNativeStatus native_capture(const Work *w, NativeGuard *g, LabNativeBatch *b) {
    const LabModel *m = w->model;
    LabNativeStatus status = LAB_NATIVE_OK;
    b->columns = m->cols; b->source_rows = m->rows;
    b->dimension = w->m; b->total_rows = w->total_rows; b->canonical_columns = w->n;
    b->scope = "ORIGINAL_MIP_IN_BOX; algebra export only; no independent cut/tree proof";
    b->box = native_boxes(m->cols); b->source_row_box = native_boxes(m->rows);
    b->integer = native_alloc(m->cols, sizeof(unsigned char));
    b->basis = native_alloc(w->m, sizeof(size_t)); b->row_id = native_alloc(w->m, sizeof(size_t));
    b->rhs = native_qvector(w->m); b->basic_values = native_qvector(w->m);
    if (!b->box || !b->source_row_box || !b->integer || !b->basis || !b->row_id || !b->rhs || !b->basic_values)
        return LAB_NATIVE_RESOURCE;
    memcpy(b->integer, m->integer, m->cols * sizeof(unsigned char));
    memcpy(b->basis, w->basis, w->m * sizeof(size_t));
    memcpy(b->row_id, w->row_id, w->m * sizeof(size_t));
    for (size_t j = 0; j < m->cols; j++) {
        b->box[j].has_lower = isfinite(m->col_lower[j]);
        b->box[j].has_upper = isfinite(m->col_upper[j]);
        if (b->box[j].has_lower && lab_source_exact_column_bound(m,j,0,b->box[j].lower)) return LAB_NATIVE_SOURCE_UNAVAILABLE;
        if (b->box[j].has_upper && lab_source_exact_column_bound(m,j,1,b->box[j].upper)) return LAB_NATIVE_SOURCE_UNAVAILABLE;
        NATIVE_TRY(native_check(g, b->box[j].lower, 1));
        NATIVE_TRY(native_check(g, b->box[j].upper, 1));
        if (b->box[j].has_lower && b->box[j].has_upper && mpq_cmp(b->box[j].lower, b->box[j].upper) > 0) {
            status = LAB_NATIVE_STALE; goto done;
        }
    }
    for (size_t i = 0; i < m->rows; i++) {
        b->source_row_box[i].has_lower = isfinite(m->row_lower[i]);
        b->source_row_box[i].has_upper = isfinite(m->row_upper[i]);
        if (b->source_row_box[i].has_lower && lab_source_exact_row_bound(m,i,0,b->source_row_box[i].lower)) return LAB_NATIVE_SOURCE_UNAVAILABLE;
        if (b->source_row_box[i].has_upper && lab_source_exact_row_bound(m,i,1,b->source_row_box[i].upper)) return LAB_NATIVE_SOURCE_UNAVAILABLE;
        NATIVE_TRY(native_check(g, b->source_row_box[i].lower, 1));
        NATIVE_TRY(native_check(g, b->source_row_box[i].upper, 1));
    }
    if (native_matrix_alloc(&b->source_matrix, m->rows, m->cols, m->nnz)) { status = LAB_NATIVE_RESOURCE; goto done; }
    memcpy(b->source_matrix.start, m->col_start, (m->cols + 1) * sizeof(size_t));
    if (m->nnz) memcpy(b->source_matrix.row, m->row_index, m->nnz * sizeof(size_t));
    for (size_t k = 0; k < m->nnz; k++) {
        if(lab_source_exact_coefficient(m,k,b->source_matrix.value[k])) return LAB_NATIVE_SOURCE_UNAVAILABLE;
        NATIVE_TRY(native_check(g, b->source_matrix.value[k], 1));
    }
    size_t nnz = 0;
    for (size_t j = 0; j < w->n; j++)
        if (!w->artificial[j]) {
            size_t extra=0;NATIVE_TRY(native_missing_zero_count(w,g,b,j,&extra));
            if (w->a[j].count > g->limit->max_canonical_nnz - nnz) { status = LAB_NATIVE_RESOURCE; goto done; }
            nnz += w->a[j].count;
            if(extra>g->limit->max_canonical_nnz-nnz){status=LAB_NATIVE_RESOURCE;goto done;}
            nnz+=extra;
        }
    if (native_matrix_alloc(&b->canonical_matrix, w->m, w->n, nnz)) { status = LAB_NATIVE_RESOURCE; goto done; }
    size_t at = 0;
    for (size_t j = 0; j < w->n; j++) {
        b->canonical_matrix.start[j] = at;
        if (!w->artificial[j])
            for (size_t k = 0; k < w->a[j].count; k++) {
                b->canonical_matrix.row[at] = w->a[j].row[k];
                NATIVE_TRY(native_basis_value(w,j,k,b->canonical_matrix.value[at]));
                NATIVE_TRY(native_check(g, b->canonical_matrix.value[at++], 1));
            }
        if(!w->artificial[j])NATIVE_TRY(native_missing_zeros(w,g,b,j,&at));
    }
    b->canonical_matrix.start[w->n] = at;
    nnz = 0;
    for (size_t p = 0; p < w->m; p++) {
        if (w->basis[p] >= w->n || w->artificial[w->basis[p]]) { status = LAB_NATIVE_EXACT_FAILURE; goto done; }
        size_t count=b->canonical_matrix.start[w->basis[p]+1]-b->canonical_matrix.start[w->basis[p]];
        if (count > g->limit->max_basis_nnz - nnz) { status = LAB_NATIVE_RESOURCE; goto done; }
        nnz += count;
    }
    if (native_matrix_alloc(&b->basis_matrix, w->m, w->m, nnz)) { status = LAB_NATIVE_RESOURCE; goto done; }
    at = 0;
    for (size_t p = 0; p < w->m; p++) {
        size_t j = w->basis[p]; b->basis_matrix.start[p] = at;
        for (size_t k = b->canonical_matrix.start[j]; k < b->canonical_matrix.start[j + 1]; k++) {
            b->basis_matrix.row[at] = b->canonical_matrix.row[k];
            mpq_set(b->basis_matrix.value[at], b->canonical_matrix.value[k]);
            NATIVE_TRY(native_check(g, b->basis_matrix.value[at++], 1));
        }
    }
    b->basis_matrix.start[w->m] = at;
    NATIVE_TRY(native_rhs(w,g,b->rhs));
    for (size_t i = 0; i < w->m; i++) NATIVE_TRY(native_check(g, b->rhs[i], 1));
    if (lab_rational_sparse_solve(w->m, b->basis_matrix.start, b->basis_matrix.row,
                                 b->basis_matrix.value, b->rhs, b->basic_values, g->deadline)) {
        status = native_tick(g) == LAB_NATIVE_DEADLINE ? LAB_NATIVE_DEADLINE : LAB_NATIVE_EXACT_FAILURE;
        goto done;
    }
    for (size_t i = 0; i < w->m; i++) {
        NATIVE_TRY(native_check(g, b->basic_values[i], 1));
        if (mpq_sgn(b->basic_values[i]) < 0) { status = LAB_NATIVE_EXACT_FAILURE; goto done; }
    }
done:
    return status;
}
static LabNativeStatus native_transpose(const LabNativeMatrix *m, LabNativeMatrix *t) {
    if (native_matrix_alloc(t, m->cols, m->rows, m->nnz)) return LAB_NATIVE_RESOURCE;
    size_t *at = native_alloc(m->rows, sizeof(size_t));
    if (!at) return LAB_NATIVE_RESOURCE;
    for (size_t k = 0; k < m->nnz; k++) t->start[m->row[k] + 1]++;
    for (size_t i = 0; i < m->rows; i++) t->start[i + 1] += t->start[i];
    memcpy(at, t->start, m->rows * sizeof(size_t));
    for (size_t j = 0; j < m->cols; j++)
        for (size_t k = m->start[j]; k < m->start[j + 1]; k++) {
            size_t q = at[m->row[k]]++;
            t->row[q] = j; mpq_set(t->value[q], m->value[k]);
        }
    free(at); return LAB_NATIVE_OK;
}
static LabNativeStatus native_dot(NativeGuard *g, const LabNativeMatrix *a, size_t col,
                                  mpq_t *weight, mpq_t out) {
    LabNativeStatus status = LAB_NATIVE_OK;
    mpq_t product; mpq_init(product); mpq_set_ui(out, 0, 1);
    for (size_t k = a->start[col]; k < a->start[col + 1]; k++) {
        mpq_mul(product, a->value[k], weight[a->row[k]]);
        NATIVE_TRY(native_check(g, product, 0));
        mpq_add(out, out, product); NATIVE_TRY(native_check(g, out, 0));
    }
done:
    mpq_clear(product); return status;
}
static LabNativeStatus native_term(const Work *w, NativeGuard *g, LabNativeTableau *t,
                                   size_t j, const mpq_t alpha) {
    if (t->terms == g->limit->max_terms) return LAB_NATIVE_RESOURCE;
    LabNativeStatus status = LAB_NATIVE_OK;
    LabNativeTerm *term = &t->term[t->terms++];
    term->canonical_column = j;
    mpq_inits(term->coefficient, term->anchor, term->canonical_scale, term->canonical_offset, NULL);
    mpq_t scale, shift, delta, product;
    mpq_inits(scale, shift, delta, product, NULL);
    if (j < w->structural) {
        size_t original = w->original[j];
        if (w->transform[original] == 'F') { status = LAB_NATIVE_FREE_SPLIT; goto done; }
        term->kind = LAB_NATIVE_VAR; term->index = original;
        term->upper = w->transform[original] == 'U';
        if(lab_source_exact_column_bound(w->model,original,term->upper,term->anchor)){status=LAB_NATIVE_SOURCE_UNAVAILABLE;goto done;}
        NATIVE_TRY(native_shift(w,original,shift));
        mpq_set_d(scale, fabs(w->sign[j])); mpq_inv(term->canonical_scale, scale);
        if (term->upper) mpq_sub(delta, shift, term->anchor);
        else mpq_sub(delta, term->anchor, shift);
        mpq_mul(term->canonical_offset, delta, term->canonical_scale);
        term->is_integer = w->model->integer[original] && mpz_cmp_ui(mpq_denref(term->anchor), 1) == 0;
    } else {
        const Column *column = &w->a[j];
        if (column->count != 1) { status = LAB_NATIVE_EXACT_FAILURE; goto done; }
        size_t i = column->row[0];
        if (i >= w->m || j != w->structural + w->row_id[i] ||
            !(w->equation_scale[i] > 0) || !isfinite(w->equation_scale[i])) {
            status = LAB_NATIVE_EXACT_FAILURE; goto done;
        }
        mpq_set_d(term->canonical_scale, w->equation_scale[i]);
        if (w->source_row[i] != SIZE_MAX) {
            term->kind = LAB_NATIVE_ROW; term->index = w->source_row[i]; term->upper = w->row_upper[i];
            if(lab_source_exact_row_bound(w->model,term->index,term->upper,term->anchor)){status=LAB_NATIVE_SOURCE_UNAVAILABLE;goto done;}
        } else {
            term->kind = LAB_NATIVE_BOUND; term->index = SIZE_MAX;
            for (size_t c = 0; c < w->model->cols; c++) {
                if (w->row_id[i] == w->lower_row[c] || w->row_id[i] == w->upper_row[c]) {
                    term->index = c; term->upper = w->row_id[i] == w->upper_row[c]; break;
                }
            }
            if (term->index == SIZE_MAX) { status = LAB_NATIVE_EXACT_FAILURE; goto done; }
            if(lab_source_exact_column_bound(w->model,term->index,term->upper,term->anchor)){status=LAB_NATIVE_SOURCE_UNAVAILABLE;goto done;}
        }
    }
    mpq_mul(term->coefficient, t->signed_basic_scale, alpha);
    mpq_mul(term->coefficient, term->coefficient, term->canonical_scale);
    mpq_mul(product, t->signed_basic_scale, alpha);
    mpq_mul(product, product, term->canonical_offset);
    mpq_sub(t->beta, t->beta, product);
    NATIVE_TRY(native_check(g, term->coefficient, 1));
    NATIVE_TRY(native_check(g, term->anchor, 1));
    NATIVE_TRY(native_check(g, term->canonical_scale, 1));
    NATIVE_TRY(native_check(g, term->canonical_offset, 1));
    NATIVE_TRY(native_check(g, t->beta, 0));
done:
    mpq_clears(scale, shift, delta, product, NULL); return status;
}
static LabNativeStatus native_affine_row(NativeGuard *g, const LabNativeBatch *b, size_t i,
                                        const mpq_t factor, mpq_t *coefficient) {
    LabNativeStatus status = LAB_NATIVE_OK;
    mpq_t product; mpq_init(product);
    for (size_t j = 0; j < b->columns; j++)
        for (size_t k = b->source_matrix.start[j]; k < b->source_matrix.start[j + 1]; k++)
            if (b->source_matrix.row[k] == i) {
                mpq_mul(product, factor, b->source_matrix.value[k]);
                NATIVE_TRY(native_check(g, product, 0));
                mpq_add(coefficient[j], coefficient[j], product);
                NATIVE_TRY(native_check(g, coefficient[j], 0));
            }
done:
    mpq_clear(product); return status;
}
static LabNativeStatus native_identity(const Work *w, NativeGuard *g, const LabNativeBatch *b,
                                       LabNativeTableau *t) {
    LabNativeStatus status = LAB_NATIVE_OK;
    mpq_t *coef = native_qvector(b->columns);
    if (!coef) return LAB_NATIVE_RESOURCE;
    mpq_t constant, factor, product;
    mpq_inits(constant, factor, product, NULL);
    mpq_set_ui(coef[t->basic_column], 1, 1); mpq_neg(constant, t->beta);
    for (size_t k = 0; k < t->terms; k++) {
        const LabNativeTerm *term = &t->term[k];
        mpq_set(factor, term->coefficient);
        if (term->upper) mpq_neg(factor, factor);
        if (term->kind == LAB_NATIVE_ROW) {
            NATIVE_TRY(native_affine_row(g, b, term->index, factor, coef));
        } else {
            mpq_add(coef[term->index], coef[term->index], factor);
            NATIVE_TRY(native_check(g, coef[term->index], 0));
        }
        mpq_mul(product, factor, term->anchor); mpq_sub(constant, constant, product);
        NATIVE_TRY(native_check(g, constant, 0));
    }
    for (size_t i = 0; i < w->m; i++) {
        size_t source = w->source_row[i];
        if (source == SIZE_MAX || !b->source_row_box[source].has_lower || !b->source_row_box[source].has_upper ||
            mpq_cmp(b->source_row_box[source].lower, b->source_row_box[source].upper)) continue;
        mpq_mul(factor, t->signed_basic_scale, t->weights[i]);
        mpq_set_d(product, w->row_scale[i]); mpq_mul(factor, factor, product);
        if (!mpq_sgn(factor)) continue;
        LabNativeMultiplier *mu = &t->equality[t->equalities++];
        mu->index = source; mpq_init(mu->multiplier); mpq_set(mu->multiplier, factor);
        NATIVE_TRY(native_check(g, mu->multiplier, 1));
        mpq_neg(factor, factor); NATIVE_TRY(native_affine_row(g, b, source, factor, coef));
        mpq_mul(product, factor, b->source_row_box[source].lower); mpq_sub(constant, constant, product);
        NATIVE_TRY(native_check(g, constant, 0));
    }
    for (size_t j = 0; j < b->columns; j++) {
        if (!mpq_sgn(coef[j])) continue;
        if (!b->box[j].has_lower || !b->box[j].has_upper || mpq_cmp(b->box[j].lower, b->box[j].upper)) {
            status = LAB_NATIVE_IDENTITY_FAILURE; goto done;
        }
        LabNativeMultiplier *nu = &t->fix[t->fixed++];
        nu->index = j; mpq_init(nu->multiplier); mpq_set(nu->multiplier, coef[j]);
        NATIVE_TRY(native_check(g, nu->multiplier, 1));
        mpq_mul(product, coef[j], b->box[j].lower); mpq_add(constant, constant, product);
        mpq_set_ui(coef[j], 0, 1);
        NATIVE_TRY(native_check(g, constant, 0));
    }
    if (mpq_sgn(constant)) status = LAB_NATIVE_IDENTITY_FAILURE;
    else t->native_identity_verified = 1;
done:
    native_qfree(coef, b->columns); mpq_clears(constant, factor, product, NULL); return status;
}
static LabNativeStatus native_row(const Work *w, NativeGuard *g, LabNativeBatch *b,
                                  const LabNativeMatrix *transpose, size_t position) {
    LabNativeStatus status = LAB_NATIVE_OK;
    LabNativeTableau *t = &b->tableau[b->count++];
    mpq_inits(t->beta, t->signed_basic_scale, t->basic_shift, NULL);
    mpq_t alpha, beta, term, shift;
    mpq_inits(alpha, beta, term, shift, NULL);
    t->canonical_basic = w->basis[position]; t->basis_position = position;
    t->basic_column = w->original[t->canonical_basic];
    t->term = native_alloc(g->limit->max_terms, sizeof(LabNativeTerm));
    t->equality = native_alloc(b->source_rows, sizeof(LabNativeMultiplier));
    t->fix = native_alloc(b->columns, sizeof(LabNativeMultiplier));
    t->weight_count = w->m; t->weights = native_qvector(w->m);
    mpq_t *unit = native_qvector(w->m);
    if (!t->term || !t->equality || !t->fix || !t->weights || !unit) { status = LAB_NATIVE_RESOURCE; goto done; }
    mpq_set_ui(unit[position], 1, 1);
    if (lab_rational_sparse_solve(w->m, transpose->start, transpose->row, transpose->value,
                                 unit, t->weights, g->deadline)) {
        status = native_tick(g) == LAB_NATIVE_DEADLINE ? LAB_NATIVE_DEADLINE : LAB_NATIVE_EXACT_FAILURE;
        goto done;
    }
    for (size_t i = 0; i < w->m; i++) {
        NATIVE_TRY(native_check(g, t->weights[i], 1));
        NATIVE_TRY(native_dot(g, &b->basis_matrix, i, t->weights, alpha));
        if (mpq_cmp_ui(alpha, i == position ? 1UL : 0UL, 1)) { status = LAB_NATIVE_IDENTITY_FAILURE; goto done; }
        mpq_mul(term, t->weights[i], b->rhs[i]); mpq_add(beta, beta, term);
        NATIVE_TRY(native_check(g, beta, 0));
    }
    if (mpq_cmp(beta, b->basic_values[position])) { status = LAB_NATIVE_IDENTITY_FAILURE; goto done; }
    t->basis_identity_verified = 1;
    mpq_set_d(t->signed_basic_scale, w->sign[t->canonical_basic]);
    NATIVE_TRY(native_shift(w,t->basic_column,shift));
    mpq_set(t->basic_shift, shift);
    mpq_mul(t->beta, t->signed_basic_scale, beta); mpq_add(t->beta, t->beta, shift);
    for (size_t j = 0; j < w->n; j++) {
        if (w->artificial[j] || w->basic[j]) continue;
        NATIVE_TRY(native_dot(g, &b->canonical_matrix, j, t->weights, alpha));
        if (mpq_sgn(alpha)) NATIVE_TRY(native_term(w, g, t, j, alpha));
    }
    NATIVE_TRY(native_check(g, t->beta, 1));
    NATIVE_TRY(native_check(g, t->signed_basic_scale, 1));
    NATIVE_TRY(native_check(g, t->basic_shift, 1));
    NATIVE_TRY(native_identity(w, g, b, t));
done:
    mpq_clears(alpha, beta, term, shift, NULL);
    native_qfree(unit, w->m); return status;
}

LabNativeStatus lab_lp_session_export_native(LabLpSession *s, const LabNativeLimits *limit,
                                             LabNativeBatch *out) {
    if (!s || !limit || !out || !limit->max_rows || limit->max_rows > 4 || limit->max_terms > 4096 ||
        !limit->max_basis_dimension || limit->max_basis_dimension > 1024 || !limit->max_basis_nnz ||
        limit->max_basis_nnz > 100000 || !limit->max_canonical_nnz || limit->max_canonical_nnz > 100000 ||
        !limit->max_bits || limit->max_bits > 8192 || !limit->max_total_bits || !isfinite(limit->deadline))
        return LAB_NATIVE_INVALID;
    if (!s->export_ready || !s->ready || s->bound_generation != s->solved_generation)
        return LAB_NATIVE_STALE;
    const Work *w = &s->work;
    const LabModel *m = &s->model;
    if (s->options.no_exact || (m->parsed_decimal && !m->exact && !m->native_exact &&
        (!m->envelope || !m->envelope->complete))) return LAB_NATIVE_SOURCE_UNAVAILABLE;
    if (m->rows>16384 || m->cols>16384 || m->nnz>100000) return LAB_NATIVE_RESOURCE;
    if (w->m > limit->max_basis_dimension) return LAB_NATIVE_RESOURCE;
    if (!w->m) return LAB_NATIVE_NO_ROWS;
    for (size_t j = 0; j < m->cols; j++)
        if (w->transform[j] == 'F') return LAB_NATIVE_FREE_SPLIT;
    NativeGuard guard = {.limit = limit, .deadline = fmin(s->deadline, limit->deadline)};
    LabNativeStatus status = native_tick(&guard);
    LabNativeBatch result; lab_native_batch_init(&result);
    LabNativeMatrix transpose = {0};
    mpq_t value, shift, scale;
    mpq_inits(value, shift, scale, NULL);
    if (status != LAB_NATIVE_OK) goto done;
    NATIVE_TRY(native_capture(w, &guard, &result));
    NATIVE_TRY(native_transpose(&result.basis_matrix, &transpose));
    result.tableau = native_alloc(limit->max_rows, sizeof(LabNativeTableau));
    if (!result.tableau) { status = LAB_NATIVE_RESOURCE; goto done; }
    /* Deterministic original-column order, rather than canonical basis order. */
    for (size_t j = 0; j < m->cols; j++) {
        if (!m->integer[j]) continue;
        for (size_t p = 0; p < w->m; p++) {
            size_t c = w->basis[p];
            if (c >= w->structural || w->original[c] != j) continue;
            NATIVE_TRY(native_shift(w,j,shift)); mpq_set_d(scale,w->sign[c]);
            mpq_mul(value, scale, result.basic_values[p]); mpq_add(value, value, shift);
            NATIVE_TRY(native_check(&guard, value, 0));
            if (mpz_cmp_ui(mpq_denref(value), 1) == 0) continue;
            if (result.count == limit->max_rows) { result.row_limit_reached = 1; goto completed; }
            NATIVE_TRY(native_row(w, &guard, &result, &transpose, p));
        }
    }
completed:
    if (!result.count) { status = LAB_NATIVE_NO_ROWS; goto done; }
    NATIVE_TRY(native_tick(&guard));
    {
        LabNativeBatch old = *out; *out = result; result = old;
    }
done:
    mpq_clears(value, shift, scale, NULL);
    native_matrix_free(&transpose); lab_native_batch_free(&result);
    return status;
}
#undef NATIVE_TRY

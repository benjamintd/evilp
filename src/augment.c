#include "augment.h"
#include "common.h"
#include "exact.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int side_pair(double lower, double upper) {
    return !isnan(lower) && !isnan(upper) && lower != INFINITY && upper != -INFINITY &&
           lower <= upper;
}
static int valid_base(const LabModel *m) {
    if (!m || !m->col_start || (m->rows && (!m->row_name || !m->row_lower || !m->row_upper)) ||
        (m->cols && (!m->col_name || !m->col_lower || !m->col_upper || !m->objective ||
                      !m->integer)) || (m->nnz && (!m->row_index || !m->value)) ||
        !isfinite(m->offset) || m->col_start[0] || m->col_start[m->cols] != m->nnz)
        return 0;
    for (size_t i = 0; i < m->rows; i++)
        if (!m->row_name[i] || !side_pair(m->row_lower[i], m->row_upper[i]))
            return 0;
    for (size_t j = 0; j < m->cols; j++)
        if (!m->col_name[j] || !side_pair(m->col_lower[j], m->col_upper[j]) ||
            !isfinite(m->objective[j]) || m->integer[j] > 1 ||
            m->col_start[j] > m->col_start[j + 1])
            return 0;
    for (size_t k = 0; k < m->nnz; k++)
        if (m->row_index[k] >= m->rows || !isfinite(m->value[k]))
            return 0;
    const LabExactModel *e = m->exact;
    if (e && (e->rows != m->rows || e->cols != m->cols || e->nnz != m->nnz))
        return 0;
    const LabExactModel *native=m->native_exact;
    if(native && (native->rows!=m->rows || native->cols!=m->cols || native->nnz!=m->nnz))
        return 0;
    const LabSourceObjective *o=m->objective_exact;
    if(o && (o->cols!=m->cols || (m->cols&&(!o->value||!o->source))))return 0;
    const LabSourceEnvelope *s = m->envelope;
    if (s && (s->rows != m->rows || s->cols != m->cols || s->nnz != m->nnz))
        return 0;
    return 1;
}
static int storage_fits(size_t rows, size_t cols, size_t nnz) {
    /* These conservative maxima also cover double/name and interval arrays. */
    size_t unit = sizeof(mpq_t) > sizeof(LabInterval) ? sizeof(mpq_t) : sizeof(LabInterval);
    if (unit < sizeof(size_t))
        unit = sizeof(size_t);
    return rows <= SIZE_MAX / unit && cols < SIZE_MAX / unit && nnz <= SIZE_MAX / unit;
}
static void copy_column_metadata(const LabModel *base, LabModel *m) {
    for (size_t j = 0; j < base->cols; j++) {
        m->col_name[j] = lab_strdup(base->col_name[j]);
        m->col_lower[j] = base->col_lower[j];
        m->col_upper[j] = base->col_upper[j];
        m->objective[j] = base->objective[j];
        m->integer[j] = base->integer[j];
        if (m->exact) {
            mpq_set(m->exact->col_lower[j], base->exact->col_lower[j]);
            mpq_set(m->exact->col_upper[j], base->exact->col_upper[j]);
            mpq_set(m->exact->objective[j], base->exact->objective[j]);
            m->exact->col_lower_source[j] = base->exact->col_lower_source[j];
            m->exact->col_upper_source[j] = base->exact->col_upper_source[j];
        }
        if(m->native_exact) {
            mpq_set(m->native_exact->col_lower[j],base->native_exact->col_lower[j]);
            mpq_set(m->native_exact->col_upper[j],base->native_exact->col_upper[j]);
            mpq_set(m->native_exact->objective[j],base->native_exact->objective[j]);
            m->native_exact->col_lower_source[j]=base->native_exact->col_lower_source[j];
            m->native_exact->col_upper_source[j]=base->native_exact->col_upper_source[j];
        }
        if(m->objective_exact){
            mpq_set(m->objective_exact->value[j],base->objective_exact->value[j]);
            m->objective_exact->source[j]=base->objective_exact->source[j];
        }
        if (m->envelope) {
            m->envelope->objective[j]=base->envelope->objective[j];
            m->envelope->col_lower[j] = base->envelope->col_lower[j];
            m->envelope->col_upper[j] = base->envelope->col_upper[j];
            m->envelope->col_lower_source[j] = base->envelope->col_lower_source[j];
            m->envelope->col_upper_source[j] = base->envelope->col_upper_source[j];
        }
    }
}
static int compare_names(const void *left, const void *right) {
    return strcmp(*(const char *const *)left, *(const char *const *)right);
}
/* Sorting the added names avoids scanning all existing rows for each cut. */
static int distinct_row_names(const LabModel *base, size_t count, const LabSparseRow *rows) {
    if (!count)
        return 1;
    const char **names = lab_alloc(count, sizeof *names);
    for (size_t r = 0; r < count; ++r)
        names[r] = rows[r].name;
    qsort(names, count, sizeof *names, compare_names);
    int ok = 1;
    for (size_t r = 1; ok && r < count; ++r)
        ok = strcmp(names[r-1], names[r]) != 0;
    for (size_t i = 0; ok && i < base->rows; ++i)
        ok = !bsearch(&base->row_name[i], names, count, sizeof *names, compare_names);
    free(names);
    return ok;
}
LabAugmentStatus lab_model_append_sparse_rows(const LabModel *base, size_t count,
    const LabSparseRow *rows, const LabAugmentLimits *limits, LabModel *out) {
    if (!out || out == base || !limits || (count && !rows) || !valid_base(base))
        return LAB_AUGMENT_INVALID;
    if (count > limits->max_added_rows || base->rows > limits->max_rows ||
        count > limits->max_rows - base->rows || base->nnz > limits->max_nonzeros ||
        base->cols == SIZE_MAX || count > SIZE_MAX / sizeof *rows)
        return LAB_AUGMENT_RESOURCE;
    size_t nnz = base->nnz;
    for (size_t r = 0; r < count; r++) {
        if (!rows[r].name || !*rows[r].name || (rows[r].count && (!rows[r].column || !rows[r].value)) ||
            !side_pair(rows[r].lower, rows[r].upper))
            return LAB_AUGMENT_INVALID;
        if (rows[r].count > limits->max_nonzeros - nnz)
            return LAB_AUGMENT_RESOURCE;
        nnz += rows[r].count;
        for (size_t k = 0; k < rows[r].count; ++k)
            if (rows[r].column[k] >= base->cols || (k && rows[r].column[k-1] >= rows[r].column[k]) ||
                !isfinite(rows[r].value[k]) || rows[r].value[k] == 0)
                return LAB_AUGMENT_INVALID;
    }
    if (!storage_fits(base->rows + count, base->cols, nnz))
        return LAB_AUGMENT_RESOURCE;
    if (!distinct_row_names(base, count, rows))
        return LAB_AUGMENT_INVALID;
    LabModel m = {.rows = base->rows + count, .cols = base->cols, .nnz = nnz,
                  .offset = base->offset, .maximize = base->maximize,
                  .parsed_decimal = base->parsed_decimal};
    m.row_name = lab_alloc(m.rows, sizeof *m.row_name);
    m.col_name = lab_alloc(m.cols, sizeof *m.col_name);
    m.row_lower = lab_alloc(m.rows, sizeof *m.row_lower);
    m.row_upper = lab_alloc(m.rows, sizeof *m.row_upper);
    m.col_lower = lab_alloc(m.cols, sizeof *m.col_lower);
    m.col_upper = lab_alloc(m.cols, sizeof *m.col_upper);
    m.objective = lab_alloc(m.cols, sizeof *m.objective);
    m.integer = lab_alloc(m.cols, sizeof *m.integer);
    m.col_start = lab_alloc(m.cols + 1, sizeof *m.col_start);
    m.row_index = lab_alloc(m.nnz, sizeof *m.row_index);
    m.value = lab_alloc(m.nnz, sizeof *m.value);
    if (base->exact) {
        m.exact = lab_exact_model_new(m.rows, m.cols, m.nnz);
        mpq_set(m.exact->offset, base->exact->offset);
    }
    if(base->native_exact) {
        m.native_exact=lab_exact_model_new(m.rows,m.cols,m.nnz);
        mpq_set(m.native_exact->offset,base->native_exact->offset);
    }
    if(base->objective_exact){
        m.objective_exact=lab_source_objective_new(m.cols);
        mpq_set(m.objective_exact->offset,base->objective_exact->offset);
        m.objective_exact->offset_source=base->objective_exact->offset_source;
    }
    if (base->envelope) {
        m.envelope = lab_source_envelope_new(m.rows, m.cols, m.nnz);
        m.envelope->complete = base->envelope->complete;
        m.envelope->objective_complete=base->envelope->objective_complete;
        m.envelope->offset=base->envelope->offset;
        if (base->envelope->exception_count)
            m.envelope->exceptions = lab_alloc(base->envelope->exception_count,
                                               sizeof *m.envelope->exceptions);
    }
    copy_column_metadata(base, &m);
    for (size_t i = 0; i < m.rows; i++) {
        int original = i < base->rows;
        size_t r = original ? 0 : i - base->rows;
        m.row_name[i] = lab_strdup(original ? base->row_name[i] : rows[r].name);
        m.row_lower[i] = original ? base->row_lower[i] : rows[r].lower;
        m.row_upper[i] = original ? base->row_upper[i] : rows[r].upper;
        if (m.exact) {
            if (original) {
                mpq_set(m.exact->row_lower[i], base->exact->row_lower[i]);
                mpq_set(m.exact->row_upper[i], base->exact->row_upper[i]);
            } else {
                if (isfinite(m.row_lower[i]))
                    mpq_set_d(m.exact->row_lower[i], m.row_lower[i]);
                if (isfinite(m.row_upper[i]))
                    mpq_set_d(m.exact->row_upper[i], m.row_upper[i]);
            }
        }
        if(m.native_exact) {
            if(original) {
                mpq_set(m.native_exact->row_lower[i],base->native_exact->row_lower[i]);
                mpq_set(m.native_exact->row_upper[i],base->native_exact->row_upper[i]);
            } else {
                if(isfinite(m.row_lower[i]))mpq_set_d(m.native_exact->row_lower[i],m.row_lower[i]);
                if(isfinite(m.row_upper[i]))mpq_set_d(m.native_exact->row_upper[i],m.row_upper[i]);
            }
        }
        if (m.envelope) {
            m.envelope->row_lower[i] = original ? base->envelope->row_lower[i]
                : (LabInterval){m.row_lower[i], m.row_lower[i]};
            m.envelope->row_upper[i] = original ? base->envelope->row_upper[i]
                : (LabInterval){m.row_upper[i], m.row_upper[i]};
        }
    }
    /* Count the whole batch before filling CSC, including the decimal
     * exception remapping for existing coefficients. */
    size_t *cursor = lab_alloc(m.cols, sizeof *cursor);
    for (size_t r = 0; r < count; ++r)
        for (size_t k = 0; k < rows[r].count; ++k)
            m.col_start[rows[r].column[k] + 1]++;
    for (size_t j = 0; j < m.cols; ++j) {
        size_t next = m.col_start[j];
        m.col_start[j + 1] += next + (base->col_start[j + 1] - base->col_start[j]);
        for (size_t k = base->col_start[j]; k < base->col_start[j + 1]; k++, next++) {
            m.row_index[next] = base->row_index[k];
            m.value[next] = base->value[k];
            if (m.exact)
                mpq_set(m.exact->value[next], base->exact->value[k]);
            if(m.native_exact)
                mpq_set(m.native_exact->value[next],base->native_exact->value[k]);
            if (m.envelope) {
                m.envelope->value[next] = base->envelope->value[k];
                size_t t = m.envelope->exception_count;
                if (t < base->envelope->exception_count && base->envelope->exceptions[t].entry == k) {
                    const LabDecimalCoefficient *src = &base->envelope->exceptions[t];
                    LabDecimalCoefficient *dest = &m.envelope->exceptions[t];
                    dest->entry = next;
                    dest->rounded = src->rounded;
                    mpq_init(dest->value);
                    mpq_set(dest->value, src->value);
                    ++m.envelope->exception_count;
                }
            }
        }
        cursor[j] = next;
    }
    for (size_t r = 0; r < count; ++r) {
        for (size_t k = 0; k < rows[r].count; ++k) {
            size_t next = cursor[rows[r].column[k]]++;
            double value = rows[r].value[k];
            m.row_index[next] = base->rows + r;
            m.value[next] = value;
            if (m.exact)
                mpq_set_d(m.exact->value[next], value);
            if(m.native_exact)
                mpq_set_d(m.native_exact->value[next],value);
            if (m.envelope)
                m.envelope->value[next] = (LabInterval){value, value};
        }
    }
    free(cursor);
    lab_model_free(out);
    *out = m;
    return LAB_AUGMENT_OK;
}

/* Compatibility for callers which already have a dense cut vector. */
LabAugmentStatus lab_model_append_rows(const LabModel *base, size_t count,
    const LabAddedRow *rows, const LabAugmentLimits *limits, LabModel *out) {
    if (!limits || !out || base == out || (count && !rows) || !valid_base(base))
        return LAB_AUGMENT_INVALID;
    if (count > limits->max_added_rows || base->rows > limits->max_rows ||
        count > limits->max_rows - base->rows || base->nnz > limits->max_nonzeros ||
        count > SIZE_MAX / sizeof(LabSparseRow))
        return LAB_AUGMENT_RESOURCE;
    LabSparseRow *sparse = lab_alloc(count, sizeof *sparse);
    size_t total = 0;
    LabAugmentStatus status = LAB_AUGMENT_INVALID;
    for (size_t r = 0; r < count; ++r) {
        if (base->cols && !rows[r].coefficient)
            goto done;
        for (size_t j = 0; j < base->cols; ++j)
            sparse[r].count += rows[r].coefficient[j] != 0;
        if (sparse[r].count > limits->max_nonzeros - base->nnz - total) {
            status = LAB_AUGMENT_RESOURCE;
            goto done;
        }
        total += sparse[r].count;
    }
    if (!storage_fits(base->rows + count, base->cols, total)) {
        status = LAB_AUGMENT_RESOURCE;
        goto done;
    }
    size_t *columns = lab_alloc(total, sizeof *columns);
    double *values = lab_alloc(total, sizeof *values);
    size_t pos = 0;
    for (size_t r = 0; r < count; ++r) {
        sparse[r].name = rows[r].name;
        sparse[r].lower = rows[r].lower;
        sparse[r].upper = rows[r].upper;
        sparse[r].column = columns + pos;
        sparse[r].value = values + pos;
        for (size_t j = 0; j < base->cols; ++j)
            if (rows[r].coefficient[j] != 0) {
                columns[pos] = j;
                values[pos++] = rows[r].coefficient[j];
            }
    }
    status = lab_model_append_sparse_rows(base, count, sparse, limits, out);
    free(values);
    free(columns);
done:
    free(sparse);
    return status;
}

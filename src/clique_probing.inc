/* Certified binary cliques (sum a_k x_k <= 1 with a_k >= 1, x_k in [0,1]).
 * Supports both:
 *  1) Direct peer probing: setting x_j = 1 forces all clique peers to 0.
 *  2) Disjoint-clique activity strengthening: grouping positive (resp. negative)
 *     binary terms that share a certified clique replaces sum(a_k) by max(a_k)
 *     (resp. min(a_k)) in outward-rounded interval arithmetic. */

static int share_clique(const size_t *cstart, const size_t *cend, const size_t *crows,
                        size_t a, size_t b) {
    size_t pa = cstart[a], ea = cend[a];
    size_t pb = cstart[b], eb = cend[b];
    while (pa < ea && pb < eb) {
        if (crows[pa] == crows[pb])
            return 1;
        if (crows[pa] < crows[pb])
            ++pa;
        else
            ++pb;
    }
    return 0;
}

static void probe_cliques(const LabModel *m, const LabInterval *a, const LabInterval *rl,
                          const LabInterval *ru, const size_t *start, const size_t *entries,
                          const size_t *column, LabInterval *box, LabPresolveResult *r,
                          double deadline) {
    if (!m->integer || lab_now() >= deadline || m->nnz > 2000000)
        return;
    size_t *cstart = lab_alloc(m->cols + 1, sizeof *cstart);
    unsigned char *is_clique = lab_alloc(m->rows ? m->rows : 1, 1);
    size_t clique_entries = 0;
    for (size_t i = 0; i < m->rows; ++i) {
        if ((i & 63U) == 0 && lab_now() >= deadline)
            break;
        size_t len = start[i + 1] - start[i];
        if (!isfinite(ru[i].upper) || ru[i].upper > 1 || len < 2 || len > 512)
            continue;
        int eligible = 1;
        for (size_t p = start[i]; p < start[i + 1]; ++p) {
            size_t k = entries[p], j = column[k];
            if (!m->integer[j] || box[j].lower < 0 || box[j].upper > 1 || a[k].lower < 1 ||
                ((p & 255U) == 0 && lab_now() >= deadline)) {
                eligible = 0;
                break;
            }
        }
        if (!eligible)
            continue;
        if (clique_entries + len > 2000000)
            continue;
        is_clique[i] = 1;
        clique_entries += len;
        for (size_t p = start[i]; p < start[i + 1]; ++p)
            cstart[column[entries[p]] + 1]++;
    }
    if (!clique_entries) {
        free(is_clique);
        free(cstart);
        return;
    }
    for (size_t j = 0; j < m->cols; ++j)
        cstart[j + 1] += cstart[j];
    size_t *cursor = lab_alloc(m->cols ? m->cols : 1, sizeof *cursor);
    memcpy(cursor, cstart, m->cols * sizeof *cursor);
    size_t *crows = lab_alloc(clique_entries, sizeof *crows);
    for (size_t i = 0; i < m->rows; ++i) {
        if (!is_clique[i])
            continue;
        for (size_t p = start[i]; p < start[i + 1]; ++p)
            crows[cursor[column[entries[p]]]++] = i;
    }
    free(cursor);
    free(is_clique);

    /* Deduplicate per-column sorted clique row lists. */
    size_t *cend = lab_alloc(m->cols ? m->cols : 1, sizeof *cend);
    for (size_t j = 0; j < m->cols; ++j) {
        size_t w = cstart[j];
        for (size_t p = cstart[j]; p < cstart[j + 1]; ++p)
            if (w == cstart[j] || crows[p] != crows[w - 1])
                crows[w++] = crows[p];
        cend[j] = w;
    }
    /* Reuse cstart/cend via a compact helper: pack in place if needed. */
    for (size_t i = 0; i < m->rows; ++i) {
        if ((i & 63U) == 0 && lab_now() >= deadline)
            break;
        size_t begin = start[i], end = start[i + 1], len = end - begin;
        if (len < 2 || len > 32)
            continue;
        /* Detect duplicate columns in this row so clique grouping is exact. */
        int has_dup = 0;
        for (size_t u = begin; u < end && !has_dup; ++u)
            for (size_t v = u + 1; v < end; ++v)
                if (column[entries[u]] == column[entries[v]]) {
                    has_dup = 1;
                    break;
                }
        for (size_t p = begin; p < end; ++p) {
            size_t k = entries[p], j = column[k];
            if (!m->integer[j] || box[j].lower != 0 || box[j].upper != 1)
                continue;
            int j_has_clique = cstart[j] < cend[j];
            /* Test x_j = 1: start with a[k] (plus any duplicate of j). */
            LabInterval sum = a[k];
            int implied = 0;
            unsigned char used_hi[32] = {0}, used_lo[32] = {0};
            for (size_t u = begin; u < end; ++u) {
                if (u == p)
                    continue;
                size_t h = entries[u], c = column[h];
                if (c == j) {
                    sum = lab_interval_add(sum, a[h]);
                    continue;
                }
                if (j_has_clique && cstart[c] < cend[c] &&
                    share_clique(cstart, cend, crows, j, c)) {
                    implied = 1;
                    used_hi[u - begin] = used_lo[u - begin] = 1;
                    continue;
                }
            }
            /* Bound remaining terms with disjoint-clique grouping (no dup cols). */
            if (!has_dup) {
                for (size_t u = begin; u < end; ++u) {
                    if (u == p || used_hi[u - begin])
                        continue;
                    size_t h = entries[u], c = column[h];
                    int bin01 = m->integer[c] && box[c].lower >= 0 && box[c].upper <= 1 &&
                                cstart[c] < cend[c];
                    if (bin01 && a[h].lower >= 0) {
                        /* Find a clique of c shared with other unused positive terms. */
                        size_t best_clique = SIZE_MAX, best_cnt = 1;
                        for (size_t q = cstart[c]; q < cend[c]; ++q) {
                            size_t cr = crows[q], cnt = 1;
                            for (size_t v = u + 1; v < end; ++v) {
                                if (v == p || used_hi[v - begin])
                                    continue;
                                size_t hv = entries[v], cv = column[hv];
                                if (!m->integer[cv] || box[cv].lower < 0 || box[cv].upper > 1 ||
                                    a[hv].lower < 0 || cstart[cv] >= cend[cv])
                                    continue;
                                /* Binary search cr in crows[cstart[cv]..cend[cv]). */
                                size_t lo_i = cstart[cv], hi_i = cend[cv];
                                while (lo_i < hi_i) {
                                    size_t mid = lo_i + (hi_i - lo_i) / 2;
                                    if (crows[mid] < cr) lo_i = mid + 1;
                                    else hi_i = mid;
                                }
                                if (lo_i < cend[cv] && crows[lo_i] == cr)
                                    ++cnt;
                            }
                            if (cnt > best_cnt) {
                                best_cnt = cnt;
                                best_clique = cr;
                            }
                        }
                        if (best_clique != SIZE_MAX) {
                            double gmax = lab_interval_multiply(a[h], box[c]).upper;
                            used_hi[u - begin] = 1;
                            implied = 1;
                            for (size_t v = u + 1; v < end; ++v) {
                                if (v == p || used_hi[v - begin])
                                    continue;
                                size_t hv = entries[v], cv = column[hv];
                                if (!m->integer[cv] || box[cv].lower < 0 || box[cv].upper > 1 ||
                                    a[hv].lower < 0 || cstart[cv] >= cend[cv])
                                    continue;
                                size_t lo_i = cstart[cv], hi_i = cend[cv];
                                while (lo_i < hi_i) {
                                    size_t mid = lo_i + (hi_i - lo_i) / 2;
                                    if (crows[mid] < best_clique) lo_i = mid + 1;
                                    else hi_i = mid;
                                }
                                if (lo_i < cend[cv] && crows[lo_i] == best_clique) {
                                    double tv = lab_interval_multiply(a[hv], box[cv]).upper;
                                    if (tv > gmax)
                                        gmax = tv;
                                    used_hi[v - begin] = 1;
                                }
                            }
                            if (gmax < 0)
                                gmax = 0;
                            sum.upper = lab_interval_add(point(sum.upper), point(gmax)).upper;
                        }
                    }
                    if (!used_hi[u - begin]) {
                        LabInterval prod = lab_interval_multiply(a[h], box[c]);
                        sum.upper = lab_interval_add(point(sum.upper), point(prod.upper)).upper;
                        used_hi[u - begin] = 1;
                    }
                }
                for (size_t u = begin; u < end; ++u) {
                    if (u == p || used_lo[u - begin])
                        continue;
                    size_t h = entries[u], c = column[h];
                    int bin01 = m->integer[c] && box[c].lower >= 0 && box[c].upper <= 1 &&
                                cstart[c] < cend[c];
                    if (bin01 && a[h].upper <= 0) {
                        size_t best_clique = SIZE_MAX, best_cnt = 1;
                        for (size_t q = cstart[c]; q < cend[c]; ++q) {
                            size_t cr = crows[q], cnt = 1;
                            for (size_t v = u + 1; v < end; ++v) {
                                if (v == p || used_lo[v - begin])
                                    continue;
                                size_t hv = entries[v], cv = column[hv];
                                if (!m->integer[cv] || box[cv].lower < 0 || box[cv].upper > 1 ||
                                    a[hv].upper > 0 || cstart[cv] >= cend[cv])
                                    continue;
                                size_t lo_i = cstart[cv], hi_i = cend[cv];
                                while (lo_i < hi_i) {
                                    size_t mid = lo_i + (hi_i - lo_i) / 2;
                                    if (crows[mid] < cr) lo_i = mid + 1;
                                    else hi_i = mid;
                                }
                                if (lo_i < cend[cv] && crows[lo_i] == cr)
                                    ++cnt;
                            }
                            if (cnt > best_cnt) {
                                best_cnt = cnt;
                                best_clique = cr;
                            }
                        }
                        if (best_clique != SIZE_MAX) {
                            double gmin = lab_interval_multiply(a[h], box[c]).lower;
                            used_lo[u - begin] = 1;
                            implied = 1;
                            for (size_t v = u + 1; v < end; ++v) {
                                if (v == p || used_lo[v - begin])
                                    continue;
                                size_t hv = entries[v], cv = column[hv];
                                if (!m->integer[cv] || box[cv].lower < 0 || box[cv].upper > 1 ||
                                    a[hv].upper > 0 || cstart[cv] >= cend[cv])
                                    continue;
                                size_t lo_i = cstart[cv], hi_i = cend[cv];
                                while (lo_i < hi_i) {
                                    size_t mid = lo_i + (hi_i - lo_i) / 2;
                                    if (crows[mid] < best_clique) lo_i = mid + 1;
                                    else hi_i = mid;
                                }
                                if (lo_i < cend[cv] && crows[lo_i] == best_clique) {
                                    double tv = lab_interval_multiply(a[hv], box[cv]).lower;
                                    if (tv < gmin)
                                        gmin = tv;
                                    used_lo[v - begin] = 1;
                                }
                            }
                            if (gmin > 0)
                                gmin = 0;
                            sum.lower = lab_interval_add(point(sum.lower), point(gmin)).lower;
                        }
                    }
                    if (!used_lo[u - begin]) {
                        LabInterval prod = lab_interval_multiply(a[h], box[c]);
                        sum.lower = lab_interval_add(point(sum.lower), point(prod.lower)).lower;
                        used_lo[u - begin] = 1;
                    }
                }
            } else {
                for (size_t u = begin; u < end; ++u) {
                    if (u == p || used_hi[u - begin])
                        continue;
                    size_t h = entries[u], c = column[h];
                    if (c == j)
                        continue;
                    sum = lab_interval_add(sum, lab_interval_multiply(a[h], box[c]));
                }
            }
            if (implied && (sum.upper < rl[i].lower || sum.lower > ru[i].upper))
                tighten(&box[j], 0, 1, 1, r);
        }
    }
    free(cend);
    free(crows);
    free(cstart);
}

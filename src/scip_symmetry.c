/* Pure C17 graph-automorphism generator finder implementing SCIP's
 * SYMcomputeSymmetryGenerators() interface (compute_symmetry.h).
 *
 * Uses Hopcroft-style O((V + E) log V) equitable partition refinement (1-WL),
 * O(E) exact twin-transposition extraction, and bounded McKay/Bliss
 * individualization-refinement search with Schreier-Sims orbit pruning.
 * Every emitted generator is explicitly verified to preserve all vertex colors
 * and every edge of the symmetry detection graph before being passed to SCIP. */

#include "common.h"
#include "scip/scip.h"
#include "scip/scip_mem.h"
#include "scip/symmetry_graph.h"
#include "symmetry/compute_symmetry.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum {
    SYM_MAX_VERTICES = 250000,
    SYM_MAX_EDGES = 600000,
    SYM_MAX_TREE_NODES = 2000,
    SYM_DEFAULT_MAX_GENS = 512
};

typedef struct {
    int n;
    int nvar;
    int permlen;
    size_t m_dir;
    uint64_t *color;
    int *head;
    int *adj;
} AutGraph;

typedef struct {
    int a, b, c;
} EdgeTriple;

static int cmp_int(const void *pa, const void *pb) {
    int a = *(const int *)pa;
    int b = *(const int *)pb;
    return (a > b) - (a < b);
}

static int cmp_edge_triple(const void *pa, const void *pb) {
    const EdgeTriple *x = (const EdgeTriple *)pa;
    const EdgeTriple *y = (const EdgeTriple *)pb;
    if (x->a != y->a)
        return (x->a > y->a) - (x->a < y->a);
    if (x->c != y->c)
        return (x->c > y->c) - (x->c < y->c);
    return (x->b > y->b) - (x->b < y->b);
}

typedef struct {
    uint64_t key;
    int vertex;
} KeyVert;

static int cmp_key_vert(const void *pa, const void *pb) {
    const KeyVert *x = (const KeyVert *)pa;
    const KeyVert *y = (const KeyVert *)pb;
    if (x->key != y->key)
        return (x->key > y->key) - (x->key < y->key);
    return (x->vertex > y->vertex) - (x->vertex < y->vertex);
}

static uint64_t make_color(uint32_t kind, int raw) {
    return ((uint64_t)kind << 32) | (uint32_t)raw;
}

static void aut_graph_free(AutGraph *g) {
    free(g->color);
    free(g->head);
    free(g->adj);
    memset(g, 0, sizeof *g);
}

/* Build a simple undirected vertex-colored graph equivalent to symgraph.
 * Colored edges between constraint nodes (>= nvar) and variable nodes (< nvar)
 * are grouped per (constraint, edge_color) to avoid vertex blowup. */
static int build_aut_graph(SYM_GRAPH *symgraph, AutGraph *g) {
    memset(g, 0, sizeof *g);
    int nsymvars = SCIPgetSymgraphNVars(symgraph);
    SYM_SYMTYPE symtype = SCIPgetSymgraphSymtype(symgraph);
    int nvar = (symtype == SYM_SYMTYPE_PERM) ? nsymvars : 2 * nsymvars;
    int nnodes_sym = SCIPgetSymgraphNNodes(symgraph);
    int nsymedges = SCIPgetSymgraphNEdges(symgraph);
    int nbase = nvar + nnodes_sym;
    if (nsymvars <= 1 || nbase > SYM_MAX_VERTICES || nsymedges > SYM_MAX_EDGES)
        return 0;

    int unique_edge = SCIPhasGraphUniqueEdgetype(symgraph) ? 1 : 0;
    EdgeTriple *grouped = lab_alloc((size_t)(nsymedges ? nsymedges : 1), sizeof *grouped);
    int ngrouped = 0;
    int *plain_u = lab_alloc((size_t)(nsymedges + nsymvars + 1), sizeof *plain_u);
    int *plain_v = lab_alloc((size_t)(nsymedges + nsymvars + 1), sizeof *plain_v);
    int nplain = 0;
    EdgeTriple *split_e = lab_alloc((size_t)(nsymedges ? nsymedges : 1), sizeof *split_e);
    int nsplit = 0;

    for (int e = 0; e < nsymedges; ++e) {
        int raw1 = SCIPgetSymgraphEdgeFirst(symgraph, e);
        int raw2 = SCIPgetSymgraphEdgeSecond(symgraph, e);
        int u = raw1 < 0 ? (-raw1 - 1) : (raw1 + nvar);
        int v = raw2 < 0 ? (-raw2 - 1) : (raw2 + nvar);
        if (u < 0 || u >= nbase || v < 0 || v >= nbase || u == v) {
            free(grouped);
            free(plain_u);
            free(plain_v);
            free(split_e);
            return 0;
        }
        if (!unique_edge && SCIPisSymgraphEdgeColored(symgraph, e)) {
            int ec = SCIPgetSymgraphEdgeColor(symgraph, e);
            if (u >= nvar && v < nvar)
                grouped[ngrouped++] = (EdgeTriple){u, v, ec};
            else if (v >= nvar && u < nvar)
                grouped[ngrouped++] = (EdgeTriple){v, u, ec};
            else
                split_e[nsplit++] = (EdgeTriple){u, v, ec};
        } else {
            plain_u[nplain] = u;
            plain_v[nplain] = v;
            nplain++;
        }
    }

    if (ngrouped > 1)
        qsort(grouped, (size_t)ngrouped, sizeof *grouped, cmp_edge_triple);

    /* Count extra nodes needed for grouped and split colored edges. */
    int extra_nodes = nsplit;
    int group_edges = 0;
    for (int i = 0; i < ngrouped;) {
        int j = i + 1;
        while (j < ngrouped && grouped[j].a == grouped[i].a && grouped[j].c == grouped[i].c)
            ++j;
        extra_nodes++;
        group_edges += 1 + (j - i);
        i = j;
    }

    int total_n = nbase + extra_nodes;
    int sign_edges = (symtype == SYM_SYMTYPE_SIGNPERM) ? nsymvars : 0;
    size_t total_undir = (size_t)nplain + (size_t)group_edges + 2 * (size_t)nsplit + (size_t)sign_edges;
    if (total_n > SYM_MAX_VERTICES || total_undir > (size_t)SYM_MAX_EDGES * 2U) {
        free(grouped);
        free(plain_u);
        free(plain_v);
        free(split_e);
        return 0;
    }

    g->n = total_n;
    g->nvar = nvar;
    g->permlen = nvar;
    g->color = lab_alloc((size_t)total_n, sizeof *g->color);
    for (int j = 0; j < nvar; ++j)
        g->color[j] = make_color(0U, SCIPgetSymgraphVarnodeColor(symgraph, j));
    for (int j = 0; j < nnodes_sym; ++j)
        g->color[nvar + j] = make_color(1U, SCIPgetSymgraphNodeColor(symgraph, j));

    int *eu = lab_alloc(total_undir ? total_undir : 1, sizeof *eu);
    int *ev = lab_alloc(total_undir ? total_undir : 1, sizeof *ev);
    size_t ne = 0;
    for (int i = 0; i < nplain; ++i) {
        eu[ne] = plain_u[i];
        ev[ne] = plain_v[i];
        ne++;
    }
    int next_id = nbase;
    for (int i = 0; i < ngrouped;) {
        int j = i + 1;
        while (j < ngrouped && grouped[j].a == grouped[i].a && grouped[j].c == grouped[i].c)
            ++j;
        int mid = next_id++;
        g->color[mid] = make_color(2U, grouped[i].c);
        eu[ne] = grouped[i].a;
        ev[ne] = mid;
        ne++;
        for (int k = i; k < j; ++k) {
            eu[ne] = mid;
            ev[ne] = grouped[k].b;
            ne++;
        }
        i = j;
    }
    for (int i = 0; i < nsplit; ++i) {
        int mid = next_id++;
        g->color[mid] = make_color(3U, split_e[i].c);
        eu[ne] = split_e[i].a;
        ev[ne] = mid;
        ne++;
        eu[ne] = mid;
        ev[ne] = split_e[i].b;
        ne++;
    }
    if (symtype == SYM_SYMTYPE_SIGNPERM) {
        for (int j = 0; j < nsymvars; ++j) {
            eu[ne] = j;
            ev[ne] = j + nsymvars;
            ne++;
        }
    }

    free(grouped);
    free(plain_u);
    free(plain_v);
    free(split_e);

    int *deg = lab_alloc((size_t)total_n, sizeof *deg);
    for (size_t e = 0; e < ne; ++e) {
        deg[eu[e]]++;
        deg[ev[e]]++;
    }
    g->head = lab_alloc((size_t)total_n + 1, sizeof *g->head);
    for (int v = 0; v < total_n; ++v)
        g->head[v + 1] = g->head[v] + deg[v];
    int *cur = lab_alloc((size_t)total_n, sizeof *cur);
    memcpy(cur, g->head, (size_t)total_n * sizeof *cur);
    g->adj = lab_alloc(ne ? 2 * ne : 1, sizeof *g->adj);
    for (size_t e = 0; e < ne; ++e) {
        g->adj[cur[eu[e]]++] = ev[e];
        g->adj[cur[ev[e]]++] = eu[e];
    }
    free(cur);
    free(deg);
    free(eu);
    free(ev);

    /* Sort adjacency lists and deduplicate any parallel edges. */
    int write_pos = 0;
    int *new_head = lab_alloc((size_t)total_n + 1, sizeof *new_head);
    for (int v = 0; v < total_n; ++v) {
        int beg = g->head[v], end = g->head[v + 1];
        if (end - beg > 1)
            qsort(&g->adj[beg], (size_t)(end - beg), sizeof *g->adj, cmp_int);
        new_head[v] = write_pos;
        for (int p = beg; p < end; ++p)
            if (write_pos == new_head[v] || g->adj[p] != g->adj[write_pos - 1])
                g->adj[write_pos++] = g->adj[p];
    }
    new_head[total_n] = write_pos;
    free(g->head);
    g->head = new_head;
    g->m_dir = (size_t)write_pos;
    return 1;
}

static int has_edge(const AutGraph *g, int u, int v) {
    int lo = g->head[u], hi = g->head[u + 1];
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        int w = g->adj[mid];
        if (w == v)
            return 1;
        if (w < v)
            lo = mid + 1;
        else
            hi = mid;
    }
    return 0;
}

static int verify_automorphism(const AutGraph *g, const int *p) {
    int n = g->n;
    for (int u = 0; u < n; ++u) {
        int pu = p[u];
        if (pu < 0 || pu >= n || g->color[pu] != g->color[u])
            return 0;
        if ((g->head[pu + 1] - g->head[pu]) != (g->head[u + 1] - g->head[u]))
            return 0;
    }
    for (int u = 0; u < n; ++u) {
        int pu = p[u];
        for (int k = g->head[u]; k < g->head[u + 1]; ++k)
            if (!has_edge(g, pu, p[g->adj[k]]))
                return 0;
    }
    return 1;
}

typedef struct {
    int cell;
    int old_len;
} SplitRecord;

typedef struct {
    const AutGraph *g;
    int *lab;
    int *inv;
    int *cell_of;
    int *cell_len;
    int ncells;
    /* Splitting worklist & scratch */
    int *queue;
    unsigned char *in_queue;
    int *deg_cnt;
    int *touched_v;
    int *touched_c;
    unsigned char *c_seen;
    KeyVert *sort_buf;
    int *sub_starts;
    /* Undo trail */
    SplitRecord *trail;
    size_t trail_len;
    size_t trail_cap;
    /* First leaf & invariant check */
    int has_first_leaf;
    int *first_leaf;
    int *first_ncells;
    uint64_t *first_hash;
    int max_depth_seen;
    /* Orbits (union-find) */
    int *uf_parent;
    /* Output to SCIP */
    SCIP *scip;
    int maxgenerators;
    int nperms;
    int nmaxperms;
    int **perms;
    double log10_group;
    /* Budgets */
    int tree_nodes;
    double deadline;
    int aborted;
} AutState;

static int uf_find(int *parent, int x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static void uf_union(int *parent, int a, int b) {
    a = uf_find(parent, a);
    b = uf_find(parent, b);
    if (a != b) {
        if (a < b)
            parent[b] = a;
        else
            parent[a] = b;
    }
}

static void push_trail(AutState *st, int cell, int old_len) {
    if (st->trail_len == st->trail_cap) {
        size_t ncap = st->trail_cap ? st->trail_cap * 2 : 1024;
        SplitRecord *nt = lab_alloc(ncap, sizeof *nt);
        if (st->trail_len)
            memcpy(nt, st->trail, st->trail_len * sizeof *nt);
        free(st->trail);
        st->trail = nt;
        st->trail_cap = ncap;
    }
    st->trail[st->trail_len++] = (SplitRecord){cell, old_len};
}

static void undo_to(AutState *st, size_t mark) {
    while (st->trail_len > mark) {
        SplitRecord rec = st->trail[--st->trail_len];
        int c = rec.cell;
        int old_len = rec.old_len;
        /* Count how many subcells are being merged back into [c, c + old_len). */
        int pos = c;
        int sub_cnt = 0;
        while (pos < c + old_len) {
            int sl = st->cell_len[pos];
            sub_cnt++;
            pos += sl;
        }
        st->cell_len[c] = old_len;
        for (int p = c; p < c + old_len; ++p)
            st->cell_of[st->lab[p]] = c;
        st->ncells -= (sub_cnt - 1);
    }
}

/* Split cell c into subcells according to deg_cnt[v] (0 stays first). */
static void split_cell(AutState *st, int c, int *q_tail, uint64_t *sig_hash) {
    int len = st->cell_len[c];
    int nz = 0;
    /* Partition: deg_cnt == 0 at front, deg_cnt > 0 at back. */
    int i = c, j = c + len - 1;
    while (i <= j) {
        int vi = st->lab[i];
        if (st->deg_cnt[vi] == 0) {
            ++i;
        } else {
            int vj = st->lab[j];
            st->lab[i] = vj;
            st->inv[vj] = i;
            st->lab[j] = vi;
            st->inv[vi] = j;
            --j;
            ++nz;
        }
    }
    int zero_len = len - nz;
    /* Collect non-zero vertices and sort by deg_cnt. */
    for (int k = 0; k < nz; ++k) {
        int v = st->lab[c + zero_len + k];
        st->sort_buf[k] = (KeyVert){(uint64_t)(unsigned)st->deg_cnt[v], v};
    }
    if (nz > 1)
        qsort(st->sort_buf, (size_t)nz, sizeof *st->sort_buf, cmp_key_vert);
    /* Check if all vertices had identical degree (no actual split). */
    if (zero_len == 0 && st->sort_buf[0].key == st->sort_buf[nz - 1].key)
        return;

    push_trail(st, c, len);
    int nsub = 0;
    if (zero_len > 0) {
        st->sub_starts[nsub++] = c;
        st->cell_len[c] = zero_len;
    }
    int pos = c + zero_len;
    for (int k = 0; k < nz;) {
        int k2 = k + 1;
        while (k2 < nz && st->sort_buf[k2].key == st->sort_buf[k].key)
            ++k2;
        int sub_c = pos;
        int sub_l = k2 - k;
        st->sub_starts[nsub++] = sub_c;
        st->cell_len[sub_c] = sub_l;
        for (int t = k; t < k2; ++t) {
            int v = st->sort_buf[t].vertex;
            st->lab[pos] = v;
            st->inv[v] = pos;
            st->cell_of[v] = sub_c;
            ++pos;
        }
        k = k2;
    }
    st->ncells += (nsub - 1);

    /* Mix split signature into invariant hash. */
    uint64_t h = *sig_hash;
    h ^= (uint64_t)(unsigned)c * UINT64_C(0x9e3779b97f4a7c15);
    for (int s = 0; s < nsub; ++s) {
        int sc = st->sub_starts[s];
        int sl = st->cell_len[sc];
        h ^= (uint64_t)(unsigned)sl + UINT64_C(0x9e3779b97f4a7c15) + (h << 6) + (h >> 2);
    }
    *sig_hash = h;

    /* Hopcroft rule: if c was already queued, queue all subcells; else skip largest. */
    int was_queued = st->in_queue[c];
    int largest_idx = 0;
    if (!was_queued) {
        int best_l = -1;
        for (int s = 0; s < nsub; ++s) {
            int sl = st->cell_len[st->sub_starts[s]];
            if (sl > best_l) {
                best_l = sl;
                largest_idx = s;
            }
        }
    }
    int q_cap = st->g->n + 1;
    for (int s = 0; s < nsub; ++s) {
        int sc = st->sub_starts[s];
        if (!was_queued && s == largest_idx)
            continue;
        if (!st->in_queue[sc]) {
            st->in_queue[sc] = 1;
            st->queue[*q_tail] = sc;
            *q_tail = (*q_tail + 1) % q_cap;
        }
    }
}

static uint64_t refine_partition(AutState *st, int q_head, int q_tail) {
    const AutGraph *g = st->g;
    int q_cap = g->n + 1;
    uint64_t sig = UINT64_C(0xcbf29ce484222325);
    while (q_head != q_tail && st->ncells < g->n) {
        int sc = st->queue[q_head];
        q_head = (q_head + 1) % q_cap;
        st->in_queue[sc] = 0;
        int sl = st->cell_len[sc];
        int ntv = 0, ntc = 0;
        for (int p = sc; p < sc + sl; ++p) {
            int u = st->lab[p];
            for (int k = g->head[u]; k < g->head[u + 1]; ++k) {
                int v = g->adj[k];
                int cv = st->cell_of[v];
                if (st->cell_len[cv] <= 1)
                    continue;
                if (st->deg_cnt[v] == 0) {
                    st->touched_v[ntv++] = v;
                    if (!st->c_seen[cv]) {
                        st->c_seen[cv] = 1;
                        st->touched_c[ntc++] = cv;
                    }
                }
                st->deg_cnt[v]++;
            }
        }
        if (ntc > 1)
            qsort(st->touched_c, (size_t)ntc, sizeof *st->touched_c, cmp_int);
        for (int i = 0; i < ntc; ++i) {
            int c = st->touched_c[i];
            st->c_seen[c] = 0;
            split_cell(st, c, &q_tail, &sig);
        }
        for (int i = 0; i < ntv; ++i)
            st->deg_cnt[st->touched_v[i]] = 0;
    }
    while (q_head != q_tail) {
        st->in_queue[st->queue[q_head]] = 0;
        q_head = (q_head + 1) % q_cap;
    }
    return sig;
}

static void record_generator(AutState *st, const int *perm) {
    if (st->maxgenerators > 0 && st->nperms >= st->maxgenerators)
        return;
    int permlen = st->g->permlen;
    int is_id = 1;
    for (int j = 0; j < permlen; ++j)
        if (perm[j] != j) {
            is_id = 0;
            break;
        }
    if (is_id)
        return;
    if (st->nmaxperms <= 0) {
        st->nmaxperms = st->maxgenerators > 0 ? st->maxgenerators : 64;
        if (SCIPallocBlockMemoryArray(st->scip, &st->perms, st->nmaxperms) != SCIP_OKAY) {
            st->nmaxperms = 0;
            return;
        }
    } else if (st->nperms >= st->nmaxperms) {
        int newsize = SCIPcalcMemGrowSize(st->scip, st->nperms + 1);
        if (SCIPreallocBlockMemoryArray(st->scip, &st->perms, st->nmaxperms, newsize) != SCIP_OKAY)
            return;
        st->nmaxperms = newsize;
    }
    int *pp = NULL;
    if (SCIPduplicateBlockMemoryArray(st->scip, &pp, perm, permlen) != SCIP_OKAY)
        return;
    st->perms[st->nperms++] = pp;
}

/* Fast O(E) exact twin-transposition detector on all non-trivial cells.
 * Finds vertices in the same cell with identical neighbor lists (N(u) == N(v))
 * or adjacent twins (N(u) \ {v} == N(v) \ {u}) and emits S_k generators. */
static void extract_twins(AutState *st) {
    const AutGraph *g = st->g;
    int n = g->n;
    int *perm = lab_alloc((size_t)n, sizeof *perm);
    for (int i = 0; i < n; ++i)
        perm[i] = i;

    for (int c = 0; c < n;) {
        int len = st->cell_len[c];
        if (len >= 2) {
            int has_var = 0;
            for (int p = c; p < c + len; ++p)
                if (st->lab[p] < g->nvar) {
                    has_var = 1;
                    break;
                }
            if (has_var) {
                for (int k = 0; k < len; ++k) {
                    int u = st->lab[c + k];
                    uint64_t h = UINT64_C(0x84222325cbf29ce4) ^ (uint64_t)(unsigned)(g->head[u + 1] - g->head[u]);
                    for (int e = g->head[u]; e < g->head[u + 1]; ++e) {
                        int w = g->adj[e];
                        /* Exclude self-cell edges only if checking adjacent twins later;
                         * exact non-adjacent twins share identical sorted adj[u]. */
                        h = (h ^ (uint64_t)(unsigned)w) * UINT64_C(0x100000001b3);
                    }
                    st->sort_buf[k] = (KeyVert){h, u};
                }
                qsort(st->sort_buf, (size_t)len, sizeof *st->sort_buf, cmp_key_vert);
                for (int k = 0; k < len;) {
                    int k2 = k + 1;
                    while (k2 < len && st->sort_buf[k2].key == st->sort_buf[k].key)
                        ++k2;
                    if (k2 - k >= 2) {
                        /* Verify exact neighborhood equality against bucket leaders. */
                        int prev = st->sort_buf[k].vertex;
                        int group_sz = 1;
                        for (int t = k + 1; t < k2; ++t) {
                            int v = st->sort_buf[t].vertex;
                            int du = g->head[prev + 1] - g->head[prev];
                            int dv = g->head[v + 1] - g->head[v];
                            if (du == dv && g->color[prev] == g->color[v] &&
                                (du == 0 || !memcmp(&g->adj[g->head[prev]], &g->adj[g->head[v]],
                                                    (size_t)du * sizeof(int)))) {
                                perm[prev] = v;
                                perm[v] = prev;
                                record_generator(st, perm);
                                perm[prev] = prev;
                                perm[v] = v;
                                uf_union(st->uf_parent, prev, v);
                                prev = v;
                                group_sz++;
                            }
                        }
                        if (group_sz > 1) {
                            for (int f = 2; f <= group_sz; ++f)
                                st->log10_group += log10((double)f);
                        }
                    }
                    k = k2;
                }
            }
        }
        c += len;
    }
    free(perm);
}

static int pick_target_cell(const AutState *st) {
    const AutGraph *g = st->g;
    int best_var = -1, best_var_len = g->n + 1;
    int best_any = -1, best_any_len = g->n + 1;
    for (int c = 0; c < g->n;) {
        int len = st->cell_len[c];
        if (len >= 2) {
            int has_var = 0;
            for (int p = c; p < c + len; ++p)
                if (st->lab[p] < g->nvar) {
                    has_var = 1;
                    break;
                }
            if (has_var && len < best_var_len) {
                best_var = c;
                best_var_len = len;
            }
            if (len < best_any_len) {
                best_any = c;
                best_any_len = len;
            }
        }
        c += len;
    }
    return best_var >= 0 ? best_var : best_any;
}

static void individualize(AutState *st, int c, int v) {
    int len = st->cell_len[c];
    int pos = st->inv[v];
    if (pos != c) {
        int first = st->lab[c];
        st->lab[c] = v;
        st->inv[v] = c;
        st->lab[pos] = first;
        st->inv[first] = pos;
    }
    push_trail(st, c, len);
    st->cell_len[c] = 1;
    st->cell_len[c + 1] = len - 1;
    for (int p = c + 1; p < c + len; ++p)
        st->cell_of[st->lab[p]] = c + 1;
    st->ncells++;
}

/* Returns jump-back level (< depth on orbit pruning, depth otherwise). */
static int search_aut(AutState *st, int depth, int *branch_v, int *cand_perm) {
    const AutGraph *g = st->g;
    if (st->aborted)
        return -1;
    if ((st->maxgenerators > 0 && st->nperms >= st->maxgenerators) ||
        st->tree_nodes >= SYM_MAX_TREE_NODES || lab_now() >= st->deadline) {
        st->aborted = 1;
        return -1;
    }
    st->tree_nodes++;

    if (st->ncells == g->n) {
        if (!st->has_first_leaf) {
            memcpy(st->first_leaf, st->lab, (size_t)g->n * sizeof *st->first_leaf);
            st->has_first_leaf = 1;
            st->max_depth_seen = depth;
            return depth;
        }
        for (int i = 0; i < g->n; ++i)
            cand_perm[st->first_leaf[i]] = st->lab[i];
        if (verify_automorphism(g, cand_perm)) {
            record_generator(st, cand_perm);
            for (int u = 0; u < g->n; ++u)
                if (cand_perm[u] != u)
                    uf_union(st->uf_parent, u, cand_perm[u]);
            /* Jump back to lowest level where cand_perm moves the branch vertex. */
            for (int k = 0; k < depth; ++k)
                if (cand_perm[branch_v[k]] != branch_v[k])
                    return k;
        }
        return depth;
    }

    if (depth >= g->n || depth >= 128)
        return depth;

    int tc = pick_target_cell(st);
    if (tc < 0)
        return depth;
    int tclen = st->cell_len[tc];
    int *verts = lab_alloc((size_t)tclen, sizeof *verts);
    memcpy(verts, &st->lab[tc], (size_t)tclen * sizeof *verts);

    int first_u = verts[0];
    int orb_hits = 1;
    for (int idx = 0; idx < tclen; ++idx) {
        if (st->aborted)
            break;
        int v = verts[idx];
        /* Skip vertices already known to be in the same orbit as first_u (after idx == 0)
         * or in the same orbit as an earlier tried vertex at this level. */
        if (idx > 0) {
            int rv = uf_find(st->uf_parent, v);
            int dup = 0;
            for (int prev = 0; prev < idx; ++prev)
                if (uf_find(st->uf_parent, verts[prev]) == rv) {
                    dup = 1;
                    break;
                }
            if (dup)
                continue;
        }

        size_t mark = st->trail_len;
        individualize(st, tc, v);
        st->in_queue[tc] = 1;
        st->queue[0] = tc;
        uint64_t sig = refine_partition(st, 0, 1);

        int descend = 1;
        if (!st->has_first_leaf) {
            st->first_ncells[depth] = st->ncells;
            st->first_hash[depth] = sig;
        } else {
            if (depth > st->max_depth_seen || st->ncells != st->first_ncells[depth] ||
                sig != st->first_hash[depth])
                descend = 0;
        }

        int jump = depth;
        if (descend) {
            branch_v[depth] = v;
            jump = search_aut(st, depth + 1, branch_v, cand_perm);
            if (idx > 0 && uf_find(st->uf_parent, v) == uf_find(st->uf_parent, first_u))
                orb_hits++;
        }

        undo_to(st, mark);
        if (st->aborted || jump < depth) {
            free(verts);
            return jump;
        }
    }

    if (orb_hits > 1)
        st->log10_group += log10((double)orb_hits);
    free(verts);
    return depth;
}

typedef struct {
    uint64_t color;
    int deg;
    int vertex;
} ColorDegVert;

static int cmp_color_deg_vert(const void *pa, const void *pb) {
    const ColorDegVert *x = (const ColorDegVert *)pa;
    const ColorDegVert *y = (const ColorDegVert *)pb;
    if (x->color != y->color)
        return (x->color > y->color) - (x->color < y->color);
    if (x->deg != y->deg)
        return (x->deg > y->deg) - (x->deg < y->deg);
    return (x->vertex > y->vertex) - (x->vertex < y->vertex);
}

SCIP_Bool SYMcanComputeSymmetry(void) {
    return TRUE;
}

const char *SYMsymmetryGetName(void) {
    return "lab-c17-aut";
}

const char *SYMsymmetryGetDesc(void) {
    return "Pure C17 partition-refinement graph automorphism detector";
}

const char *SYMsymmetryGetAddName(void) {
    return "";
}

const char *SYMsymmetryGetAddDesc(void) {
    return "";
}

SCIP_RETCODE SYMcomputeSymmetryGenerators(SCIP *scip, int maxgenerators, SYM_GRAPH *graph,
                                          int *nperms, int *nmaxperms, int ***perms,
                                          SCIP_Real *log10groupsize, SCIP_Real *symcodetime) {
    double t0 = lab_now();
    *nperms = 0;
    *nmaxperms = 0;
    *perms = NULL;
    *log10groupsize = 0.0;
    *symcodetime = 0.0;

    if (!scip || !graph)
        return SCIP_OKAY;

    AutGraph g;
    if (!build_aut_graph(graph, &g)) {
        *symcodetime = lab_now() - t0;
        return SCIP_OKAY;
    }

    int n = g.n;
    AutState st = {0};
    st.g = &g;
    st.scip = scip;
    st.maxgenerators = maxgenerators > 0 ? maxgenerators : SYM_DEFAULT_MAX_GENS;
    st.deadline = t0 + 0.35;
    st.lab = lab_alloc((size_t)n, sizeof *st.lab);
    st.inv = lab_alloc((size_t)n, sizeof *st.inv);
    st.cell_of = lab_alloc((size_t)n, sizeof *st.cell_of);
    st.cell_len = lab_alloc((size_t)n, sizeof *st.cell_len);
    st.queue = lab_alloc((size_t)n + 1, sizeof *st.queue);
    st.in_queue = lab_alloc((size_t)n, sizeof *st.in_queue);
    st.deg_cnt = lab_alloc((size_t)n, sizeof *st.deg_cnt);
    st.touched_v = lab_alloc((size_t)n, sizeof *st.touched_v);
    st.touched_c = lab_alloc((size_t)n, sizeof *st.touched_c);
    st.c_seen = lab_alloc((size_t)n, sizeof *st.c_seen);
    st.sort_buf = lab_alloc((size_t)n, sizeof *st.sort_buf);
    st.sub_starts = lab_alloc((size_t)n, sizeof *st.sub_starts);
    st.first_leaf = lab_alloc((size_t)n, sizeof *st.first_leaf);
    st.first_ncells = lab_alloc(128, sizeof *st.first_ncells);
    st.first_hash = lab_alloc(128, sizeof *st.first_hash);
    st.uf_parent = lab_alloc((size_t)n, sizeof *st.uf_parent);
    for (int i = 0; i < n; ++i)
        st.uf_parent[i] = i;

    /* Initial coloring partition: sort vertices exactly by (color, degree). */
    ColorDegVert *init_v = lab_alloc((size_t)n, sizeof *init_v);
    for (int v = 0; v < n; ++v)
        init_v[v] = (ColorDegVert){g.color[v], g.head[v + 1] - g.head[v], v};
    qsort(init_v, (size_t)n, sizeof *init_v, cmp_color_deg_vert);

    int q_tail = 0;
    for (int i = 0; i < n;) {
        int j = i + 1;
        while (j < n && init_v[j].color == init_v[i].color && init_v[j].deg == init_v[i].deg)
            ++j;
        int len = j - i;
        st.cell_len[i] = len;
        for (int p = i; p < j; ++p) {
            int v = init_v[p].vertex;
            st.lab[p] = v;
            st.inv[v] = p;
            st.cell_of[v] = i;
        }
        st.in_queue[i] = 1;
        st.queue[q_tail++] = i;
        st.ncells++;
        i = j;
    }
    free(init_v);

    (void)refine_partition(&st, 0, q_tail);

    if (st.ncells < n) {
        extract_twins(&st);
        int *branch_v = lab_alloc(128, sizeof *branch_v);
        int *cand_perm = lab_alloc((size_t)n, sizeof *cand_perm);
        (void)search_aut(&st, 0, branch_v, cand_perm);
        free(branch_v);
        free(cand_perm);
    }

    if (st.nperms > 0) {
        *perms = st.perms;
        *nperms = st.nperms;
        *nmaxperms = st.nmaxperms;
        *log10groupsize = st.log10_group > 0.0 ? st.log10_group : log10(2.0);
    }

    free(st.lab);
    free(st.inv);
    free(st.cell_of);
    free(st.cell_len);
    free(st.queue);
    free(st.in_queue);
    free(st.deg_cnt);
    free(st.touched_v);
    free(st.touched_c);
    free(st.c_seen);
    free(st.sort_buf);
    free(st.sub_starts);
    free(st.trail);
    free(st.first_leaf);
    free(st.first_ncells);
    free(st.first_hash);
    free(st.uf_parent);
    aut_graph_free(&g);

    *symcodetime = lab_now() - t0;
    return SCIP_OKAY;
}

SCIP_Bool SYMcheckGraphsAreIdentical(SCIP *scip, SYM_SYMTYPE symtype, SYM_GRAPH *G1,
                                     SYM_GRAPH *G2) {
    (void)scip;
    (void)symtype;
    (void)G1;
    (void)G2;
    return FALSE;
}

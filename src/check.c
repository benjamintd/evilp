#define _POSIX_C_SOURCE 200809L
#include "common.h"
#include <ctype.h>
#include <errno.h>
#include <gmp.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

typedef struct {
    char *name;
    char type;
    mpq_t rhs, range, activity, dual;
    int has_range, has_rhs, dual_seen;
} Row;
typedef struct {
    char *name;
    mpq_t lo, hi, value, residual;
    int lf, hf, integer, bounded, explicit_lo, seen;
} Col;
typedef struct {
    char **keys;
    size_t *values, capacity, count;
} Map;
typedef struct {
    Row *r;
    Col *c;
    size_t nr, nc, rc, cc, nnz;
    Map rows, cols;
    size_t objective;
    int maximize;
    char *objname;
    mpq_t obj;
    char error[512];
    int unsupported, dual_enabled;
} Model;
static void fail(Model *m, const char *fmt, ...) {
    if (m->error[0])
        return;
    va_list a;
    va_start(a, fmt);
    vsnprintf(m->error, sizeof m->error, fmt, a);
    va_end(a);
}
static uint64_t hash(const char *s) {
    uint64_t h = 14695981039346656037ULL;
    for (; *s; s++) {
        h ^= (unsigned char)*s;
        h *= 1099511628211ULL;
    }
    return h;
}
static void map_grow(Map *m) {
    size_t cap = m->capacity ? m->capacity * 2 : 1024;
    Map n = {0};
    n.capacity = cap;
    n.keys = lab_alloc(cap, sizeof(char *));
    n.values = lab_alloc(cap, sizeof(size_t));
    for (size_t i = 0; i < m->capacity; i++)
        if (m->keys[i]) {
            size_t k = (size_t)hash(m->keys[i]) & (cap - 1);
            while (n.keys[k])
                k = (k + 1) & (cap - 1);
            n.keys[k] = m->keys[i];
            n.values[k] = m->values[i];
            n.count++;
        }
    free(m->keys);
    free(m->values);
    *m = n;
}
static size_t lookup(const Map *m, const char *s) {
    if (!m->capacity)
        return SIZE_MAX;
    size_t k = (size_t)hash(s) & (m->capacity - 1);
    while (m->keys[k]) {
        if (!strcmp(m->keys[k], s))
            return m->values[k];
        k = (k + 1) & (m->capacity - 1);
    }
    return SIZE_MAX;
}
static void insert(Map *m, char *s, size_t v) {
    if ((m->count + 1) * 2 >= m->capacity)
        map_grow(m);
    size_t k = (size_t)hash(s) & (m->capacity - 1);
    while (m->keys[k])
        k = (k + 1) & (m->capacity - 1);
    m->keys[k] = s;
    m->values[k] = v;
    m->count++;
}
/* Decimal/scientific input -> exact rational. No binary floating point conversion. */
static int number(mpq_t q, const char *s) {
    if (!s || strlen(s) > 4096)
        return -1;
    char digits[4098];
    size_t n = 0;
    int neg = 0, frac = 0, dot = 0;
    long exponent = 0;
    if (*s == '+' || *s == '-') {
        neg = *s == '-';
        s++;
    }
    while (*s && *s != 'e' && *s != 'E' && *s != 'd' && *s != 'D') {
        if (*s == '.' && !dot) {
            dot = 1;
            s++;
            continue;
        }
        if (!isdigit((unsigned char)*s))
            return -1;
        digits[n++] = *s++;
        if (dot)
            frac++;
    }
    if (!n)
        return -1;
    if (*s) {
        char *end;
        errno = 0;
        exponent = strtol(s + 1, &end, 10);
        if (errno || end == s + 1 || *end || exponent > 10000 || exponent < -10000)
            return -1;
    }
    digits[n] = 0;
    mpz_t power;
    mpz_init(power);
    if (mpz_set_str(mpq_numref(q), digits, 10)) {
        mpz_clear(power);
        return -1;
    }
    if (neg)
        mpz_neg(mpq_numref(q), mpq_numref(q));
    long shift = exponent - frac;
    mpz_ui_pow_ui(power, 10, (unsigned long)(shift < 0 ? -shift : shift));
    mpz_set_ui(mpq_denref(q), 1);
    if (shift < 0)
        mpz_set(mpq_denref(q), power);
    else
        mpz_mul(mpq_numref(q), mpq_numref(q), power);
    mpq_canonicalize(q);
    mpz_clear(power);
    return 0;
}
/* Certificate files may contain fractions. MPS and native SOL remain decimal-only. */
static int rational(mpq_t q, const char *s) {
    if (!s || strlen(s) > 4096)
        return -1;
    const char *slash = strchr(s, '/');
    if (!slash)
        return number(q, s);
    const char *p = s;
    if (*p == '+' || *p == '-')
        p++;
    if (p == slash)
        return -1;
    while (p < slash)
        if (!isdigit((unsigned char)*p++))
            return -1;
    p = slash + 1;
    if (!*p)
        return -1;
    while (*p)
        if (!isdigit((unsigned char)*p++))
            return -1;
    if (mpq_set_str(q, *s == '+' ? s + 1 : s, 10) || mpz_sgn(mpq_denref(q)) == 0)
        return -1;
    mpq_canonicalize(q);
    return 0;
}
static int tokens(char *s, char **v, int cap) {
    int n = 0;
    char *save;
    for (char *t = strtok_r(s, " \t\r\n", &save); t; t = strtok_r(NULL, " \t\r\n", &save)) {
        if (t[0] == '$')
            break;
        if (n >= cap)
            return -1;
        v[n++] = t;
    }
    return n;
}
static size_t column(Model *m, const char *s, int integer) {
    size_t i = lookup(&m->cols, s);
    if (i != SIZE_MAX) {
        if (integer >= 0 && m->c[i].integer != integer)
            fail(m, "conflicting integrality markers for column %s", s);
        return i;
    }
    if (m->nc == m->cc) {
        size_t cap = m->cc ? m->cc * 2 : 256;
        Col *p = realloc(m->c, cap * sizeof *p);
        if (!p)
            exit(70);
        m->c = p;
        m->cc = cap;
    }
    i = m->nc++;
    Col *c = &m->c[i];
    memset(c, 0, sizeof *c);
    c->name = lab_strdup(s);
    mpq_inits(c->lo, c->hi, c->value, c->residual, NULL);
    c->lf = 1;
    c->integer = integer > 0;
    if (c->integer) {
        c->hf = 1;
        mpq_set_ui(c->hi, 1, 1);
    }
    insert(&m->cols, c->name, i);
    return i;
}
static int vector_name(Model *m, char **selected, const char *s) {
    if (!*selected)
        *selected = lab_strdup(s);
    if (strcmp(*selected, s)) {
        m->unsupported = 1;
        fail(m, "multiple named RHS/range/bound vectors are unsupported");
        return -1;
    }
    return 0;
}
enum Section { NONE, ROWS, COLUMNS, RHS, RANGES, BOUNDS, SENSE, OBJNAME };
static int read_mps(Model *m, const char *path, int evaluate) {
    gzFile f = gzopen(path, "rb");
    if (!f) {
        fail(m, "cannot open model");
        return -1;
    }
    gzbuffer(f, 1 << 20);
    char line[16384], *v[12], lastcol[16384] = "";
    enum Section sec = NONE;
    int integer = 0, ended = 0, sawrows = 0, sawcols = 0;
    size_t ln = 0;
    char *rhsset = NULL, *rangeset = NULL, *boundset = NULL;
    mpq_t a, t;
    mpq_inits(a, t, NULL);
    while (gzgets(f, line, sizeof line) && !m->error[0]) {
        ln++;
        size_t len = strlen(line);
        if (len == sizeof(line) - 1 && line[len - 1] != '\n') {
            fail(m, "line too long at %zu", ln);
            break;
        }
        if (line[0] == '*')
            continue;
        int header = !isspace((unsigned char)line[0]);
        int n = tokens(line, v, 12);
        if (!n)
            continue;
        if (n < 0) {
            fail(m, "too many fields at line %zu", ln);
            break;
        }
        if (header) {
            if (!strcmp(v[0], "NAME"))
                continue;
            if (!strcmp(v[0], "ENDATA")) {
                ended = 1;
                break;
            }
            if (!strcmp(v[0], "ROWS")) {
                sec = ROWS;
                sawrows++;
                continue;
            }
            if (!strcmp(v[0], "COLUMNS")) {
                sec = COLUMNS;
                sawcols++;
                continue;
            }
            if (!strcmp(v[0], "RHS")) {
                sec = RHS;
                continue;
            }
            if (!strcmp(v[0], "RANGES")) {
                sec = RANGES;
                continue;
            }
            if (!strcmp(v[0], "BOUNDS")) {
                sec = BOUNDS;
                continue;
            }
            if (!strcmp(v[0], "OBJSENSE")) {
                sec = SENSE;
                if (n == 1)
                    continue;
                v[0] = v[1];
                n = 1;
            } else if (!strcmp(v[0], "OBJNAME")) {
                sec = OBJNAME;
                if (n == 1)
                    continue;
                v[0] = v[1];
                n = 1;
            } else if (sec == NONE || !strcmp(v[0], "SOS") || !strncmp(v[0], "Q", 1) ||
                       !strcmp(v[0], "INDICATORS") || !strcmp(v[0], "CSECTION")) {
                m->unsupported = 1;
                fail(m, "unsupported section %s at line %zu", v[0], ln);
                break;
            }
        }
        if (sec == SENSE) {
            if (n != 1 || (strcmp(v[0], "MIN") && strcmp(v[0], "MAX") && strcmp(v[0], "MINIMIZE") &&
                           strcmp(v[0], "MAXIMIZE"))) {
                fail(m, "invalid objective sense");
                break;
            }
            if (!evaluate)
                m->maximize = !strncmp(v[0], "MAX", 3);
            continue;
        }
        if (sec == OBJNAME) {
            if (n != 1) {
                fail(m, "invalid OBJNAME");
                break;
            }
            if (!evaluate) {
                free(m->objname);
                m->objname = lab_strdup(v[0]);
            }
            continue;
        }
        if (sec == ROWS) {
            if (n != 2 || strlen(v[0]) != 1 || !strchr("NLEG", v[0][0])) {
                fail(m, "unsupported/malformed row at line %zu", ln);
                break;
            }
            if (evaluate)
                continue;
            if (lookup(&m->rows, v[1]) != SIZE_MAX) {
                fail(m, "duplicate row %s", v[1]);
                break;
            }
            if (m->nr == m->rc) {
                size_t cap = m->rc ? m->rc * 2 : 256;
                Row *p = realloc(m->r, cap * sizeof *p);
                if (!p)
                    exit(70);
                m->r = p;
                m->rc = cap;
            }
            Row *r = &m->r[m->nr];
            memset(r, 0, sizeof *r);
            r->name = lab_strdup(v[1]);
            r->type = v[0][0];
            mpq_inits(r->rhs, r->range, r->activity, r->dual, NULL);
            if (r->type == 'N' && m->objective == SIZE_MAX)
                m->objective = m->nr;
            insert(&m->rows, r->name, m->nr++);
            continue;
        }
        if (sec == COLUMNS) {
            if (n >= 3 && !strcmp(v[1], "'MARKER'")) {
                if (!strcmp(v[2], "'INTORG'") && !integer)
                    integer = 1;
                else if (!strcmp(v[2], "'INTEND'") && integer)
                    integer = 0;
                else
                    fail(m, "invalid integer marker at line %zu", ln);
                continue;
            }
            const char *col;
            int start;
            if (n == 3 || n == 5) {
                col = v[0];
                snprintf(lastcol, sizeof lastcol, "%s", col);
                start = 1;
            } else if ((n == 2 || n == 4) && lastcol[0]) {
                col = lastcol;
                start = 0;
            } else {
                fail(m, "malformed columns at line %zu", ln);
                break;
            }
            size_t c = evaluate ? lookup(&m->cols, col) : column(m, col, integer);
            if (c == SIZE_MAX) {
                fail(m, "unknown column");
                break;
            }
            for (int j = start; j < n; j += 2) {
                size_t r = lookup(&m->rows, v[j]);
                if (r == SIZE_MAX || number(a, v[j + 1])) {
                    fail(m, "unknown row or invalid coefficient at line %zu", ln);
                    break;
                }
                if (evaluate) {
                    mpq_mul(t, a, m->c[c].value);
                    mpq_add(m->r[r].activity, m->r[r].activity, t);
                    if (m->dual_enabled) {
                        if (r == m->objective) {
                            if (m->maximize)
                                mpq_sub(m->c[c].residual, m->c[c].residual, a);
                            else
                                mpq_add(m->c[c].residual, m->c[c].residual, a);
                        } else if (m->r[r].type != 'N') {
                            mpq_mul(t, a, m->r[r].dual);
                            mpq_sub(m->c[c].residual, m->c[c].residual, t);
                        }
                    }
                } else
                    m->nnz++;
            }
            continue;
        }
        if (sec == RHS || sec == RANGES) {
            if (n != 3 && n != 5) {
                fail(m, "blank/malformed RHS or RANGES vector at line %zu", ln);
                break;
            }
            if (vector_name(m, sec == RHS ? &rhsset : &rangeset, v[0]))
                break;
            for (int j = 1; j < n; j += 2) {
                size_t r = lookup(&m->rows, v[j]);
                if (r == SIZE_MAX || number(a, v[j + 1])) {
                    fail(m, "unknown row or invalid RHS at line %zu", ln);
                    break;
                }
                if (!evaluate) {
                    Row *row = &m->r[r];
                    if (sec == RHS) {
                        if (row->has_rhs) {
                            fail(m, "duplicate RHS entry");
                            break;
                        }
                        row->has_rhs = 1;
                        mpq_set(row->rhs, a);
                    } else {
                        if (row->has_range) {
                            fail(m, "duplicate range entry");
                            break;
                        }
                        row->has_range = 1;
                        mpq_set(row->range, a);
                    }
                }
            }
            continue;
        }
        if (sec == BOUNDS) {
            if (n < 3 || n > 4 || vector_name(m, &boundset, v[1])) {
                fail(m, "malformed bound at line %zu", ln);
                break;
            }
            const char *type = v[0];
            int needs = (!strcmp(type, "LO") || !strcmp(type, "UP") || !strcmp(type, "FX") ||
                         !strcmp(type, "LI") || !strcmp(type, "UI"));
            if (strcmp(type, "LO") && strcmp(type, "UP") && strcmp(type, "FX") &&
                strcmp(type, "LI") && strcmp(type, "UI") && strcmp(type, "BV") &&
                strcmp(type, "FR") && strcmp(type, "MI") && strcmp(type, "PL")) {
                m->unsupported = 1;
                fail(m, "unsupported bound type %s", type);
                break;
            }
            if ((needs && n != 4) || (!needs && n != 3 && n != 4)) {
                fail(m, "bound arity at line %zu", ln);
                break;
            }
            /* Some standard MPS writers include an ignored value on BV/FR/MI/PL. */
            if (!needs && n == 4 && number(t, v[3])) {
                fail(m, "invalid ignored bound value");
                break;
            }
            if (needs && number(a, v[3])) {
                fail(m, "invalid bound number");
                break;
            }
            if (evaluate)
                continue;
            size_t ci = column(m, v[2], -1);
            Col *c = &m->c[ci];
            if (!c->bounded) {
                if (c->integer)
                    c->hf = 0;
                c->bounded = 1;
            }
            if (!strcmp(type, "LO") || !strcmp(type, "LI")) {
                c->lf = 1;
                c->explicit_lo = 1;
                mpq_set(c->lo, a);
                if (type[0] == 'L' && type[1] == 'I')
                    c->integer = 1;
            } else if (!strcmp(type, "UP") || !strcmp(type, "UI")) {
                c->hf = 1;
                mpq_set(c->hi, a);
                if (!c->explicit_lo && mpq_sgn(a) < 0)
                    c->lf = 0;
                if (type[1] == 'I')
                    c->integer = 1;
            } else if (!strcmp(type, "FX")) {
                c->lf = c->hf = c->explicit_lo = 1;
                mpq_set(c->lo, a);
                mpq_set(c->hi, a);
            } else if (!strcmp(type, "BV")) {
                c->integer = c->lf = c->hf = c->explicit_lo = 1;
                mpq_set_ui(c->lo, 0, 1);
                mpq_set_ui(c->hi, 1, 1);
            } else if (!strcmp(type, "FR")) {
                c->lf = c->hf = 0;
                c->explicit_lo = 1;
            } else if (!strcmp(type, "MI")) {
                c->lf = 0;
                c->explicit_lo = 1;
            } else if (!strcmp(type, "PL"))
                c->hf = 0;
            continue;
        }
        fail(m, "unexpected record at line %zu", ln);
    }
    int gzerr;
    gzerror(f, &gzerr);
    if (gzerr != Z_OK && gzerr != Z_STREAM_END)
        fail(m, "compressed model read error");
    if (!ended || sawrows != 1 || sawcols != 1 || integer)
        fail(m, "incomplete or repeated model sections");
    gzclose(f);
    free(rhsset);
    free(rangeset);
    free(boundset);
    mpq_clears(a, t, NULL);
    if (!evaluate && m->objname) {
        m->objective = lookup(&m->rows, m->objname);
        if (m->objective == SIZE_MAX || m->r[m->objective].type != 'N')
            fail(m, "invalid objective row");
    }
    if (m->objective == SIZE_MAX)
        fail(m, "missing objective row");
    return m->error[0] ? -1 : 0;
}
static void cleanup(Model *m) {
    for (size_t i = 0; i < m->nr; i++) {
        free(m->r[i].name);
        mpq_clears(m->r[i].rhs, m->r[i].range, m->r[i].activity, m->r[i].dual, NULL);
    }
    for (size_t i = 0; i < m->nc; i++) {
        free(m->c[i].name);
        mpq_clears(m->c[i].lo, m->c[i].hi, m->c[i].value, m->c[i].residual, NULL);
    }
    free(m->r);
    free(m->c);
    free(m->objname);
    free(m->rows.keys);
    free(m->rows.values);
    free(m->cols.keys);
    free(m->cols.values);
    mpq_clear(m->obj);
}
static int solution(Model *m, const char *path, const char *format, mpq_t claimed,
                    int *has_claimed) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fail(m, "cannot open solution");
        return -1;
    }
    char line[16384], copy[16384], *v[12];
    int active = strcmp(format, "highs") != 0, marked = 0, no_solution = 0;
    int highs_primal = 0, highs_finished = 0;
    size_t assigned = 0;
    while (fgets(line, sizeof line, f) && !m->error[0]) {
        if (strlen(line) == sizeof line - 1 && !strchr(line, '\n')) {
            fail(m, "solution line too long");
            break;
        }
        if (!strcmp(format, "highs")) {
            /* HiGHS repeats Columns/Rows headers for primal, dual and basis
             * data. Only the first primal block defines this witness. */
            if (!strncmp(line, "# Primal solution values", 24)) {
                if (highs_primal || highs_finished) {
                    fail(m, "duplicate/out-of-order HiGHS primal block");
                    break;
                }
                highs_primal = 1;
                active = 0;
                continue;
            }
            if (!strncmp(line, "# Dual", 6) || !strncmp(line, "# Basis", 7)) {
                highs_finished = 1;
                active = 0;
                continue;
            }
            if (highs_finished)
                continue;
            if (!strncmp(line, "# Columns", 9)) {
                active = highs_primal;
                continue;
            }
            if (!strncmp(line, "# Rows", 6)) {
                active = 0;
                continue;
            }
        }
        if (!strncmp(line, "=infeas=", 8) ||
            (!strcmp(format, "scip") &&
             ((!strncmp(line, "solution status:", 16) &&
               (strstr(line, "infeasible") || strstr(line, "no solution"))) ||
              !strncmp(line, "no solution available", 21))) ||
            (!strcmp(format, "cbc") &&
             (!strncmp(line, "Infeasible", 10) || !strncmp(line, "Integer infeasible", 18) ||
              strstr(line, "no integer solution"))) ||
            (!strcmp(format, "highs") &&
             (!strcmp(line, "Infeasible\n") || !strcmp(line, "None\n")) && !marked)) {
            no_solution = 1;
            break;
        }
        snprintf(copy, sizeof copy, "%s", line);
        int n = tokens(copy, v, 12);
        if (n < 0) {
            fail(m, "too many solution fields");
            break;
        }
        if (!n)
            continue;
        const char *name = NULL, *value = NULL, *obj = NULL;
        if (!strcmp(format, "highs")) {
            if (highs_primal && !strcmp(v[0], "Feasible")) {
                marked = 1;
                continue;
            }
            if (highs_primal && n == 2 && !strcmp(v[0], "Objective"))
                obj = v[1];
            else if (active && n == 2 && v[0][0] != '#') {
                name = v[0];
                value = v[1];
            }
        } else if (!strcmp(format, "scip")) {
            if (!strncmp(line, "solution status:", 16)) {
                marked = 1;
                continue;
            }
            if (n >= 3 && !strcmp(v[0], "objective") && !strcmp(v[1], "value:"))
                obj = v[2];
            else if (n >= 2 && v[0][0] != '#') {
                name = v[0];
                value = v[1];
            }
        } else if (!strcmp(format, "cbc")) {
            char *o = strstr(line, "objective value ");
            if (o) {
                char *e;
                obj = o + 16;
                for (e = (char *)obj; *e && !isspace((unsigned char)*e); e++) {
                }
                *e = 0;
                marked = 1;
            } else if (n >= 3) {
                int j = !strcmp(v[0], "**") ? 1 : 0;
                if (n < j + 3) {
                    fail(m, "malformed Cbc solution");
                    break;
                }
                name = v[j + 1];
                value = v[j + 2];
            }
        } else {
            if (n == 2 && !strcmp(v[0], "=obj="))
                obj = v[1];
            else if (n == 2 && v[0][0] != '#') {
                name = v[0];
                value = v[1];
            } else if (v[0][0] != '#')
                fail(m, "malformed simple solution record");
        }
        if (obj) {
            if (*has_claimed ||
                (!strcmp(format, "rational") ? rational(claimed, obj) : number(claimed, obj))) {
                fail(m, "invalid/duplicate claimed objective");
                break;
            }
            *has_claimed = 1;
            marked = 1;
        }
        if (name) {
            size_t c = lookup(&m->cols, name);
            if (c == SIZE_MAX) {
                fail(m, "unknown solution variable %s", name);
                break;
            }
            if (m->c[c].seen || (!strcmp(format, "rational") ? rational(m->c[c].value, value)
                                                             : number(m->c[c].value, value))) {
                fail(m, "duplicate/invalid value for %s", name);
                break;
            }
            m->c[c].seen = 1;
            assigned++;
        }
    }
    if (ferror(f))
        fail(m, "solution read error");
    fclose(f);
    if (no_solution)
        return 1;
    if (!marked && !assigned)
        fail(m, "empty solution");
    return m->error[0] ? -1 : 0;
}
static int read_dual(Model *m, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fail(m, "cannot open dual certificate");
        return -1;
    }
    char line[16384], *v[4];
    while (fgets(line, sizeof line, f) && !m->error[0]) {
        if (strlen(line) == sizeof line - 1 && !strchr(line, '\n')) {
            fail(m, "dual certificate line too long");
            break;
        }
        char *start = line;
        while (isspace((unsigned char)*start))
            start++;
        if (*start == '#')
            continue;
        int n = tokens(line, v, 4);
        if (!n || (n > 0 && v[0][0] == '#'))
            continue;
        if (n != 2) {
            fail(m, "malformed dual certificate record");
            break;
        }
        size_t r = lookup(&m->rows, v[0]);
        if (r == SIZE_MAX || m->r[r].type == 'N') {
            fail(m, "unknown/nonconstraint dual row %s", v[0]);
            break;
        }
        if (m->r[r].dual_seen || rational(m->r[r].dual, v[1])) {
            fail(m, "duplicate/invalid dual multiplier for %s", v[0]);
            break;
        }
        m->r[r].dual_seen = 1;
    }
    if (ferror(f))
        fail(m, "dual certificate read error");
    fclose(f);
    m->dual_enabled = !m->error[0];
    return m->error[0] ? -1 : 0;
}
static void row_bounds(const Row *r, mpq_t lo, mpq_t hi, int *lf, int *hf, mpq_t scratch) {
    *lf = r->type != 'L';
    *hf = r->type != 'G';
    mpq_set(lo, r->rhs);
    mpq_set(hi, r->rhs);
    if (r->has_range) {
        mpq_abs(scratch, r->range);
        *lf = *hf = 1;
        if (r->type == 'L' || (r->type == 'E' && mpq_sgn(r->range) < 0))
            mpq_sub(lo, lo, scratch);
        else
            mpq_add(hi, hi, scratch);
    }
}
/* Any rational y gives a safe support-function bound when all selected sides
 * are finite. Reduced costs are recomputed from original input, never trusted. */
static int certified_bound(Model *m, mpq_t bound, char message[512]) {
    mpq_t lo, hi, t;
    mpq_inits(lo, hi, t, NULL);
    mpq_set(bound, m->r[m->objective].rhs);
    if (!m->maximize)
        mpq_neg(bound, bound);
    int available = 1;
    for (size_t i = 0; i < m->nr && available; i++) {
        Row *r = &m->r[i];
        int sign = mpq_sgn(r->dual);
        if (r->type == 'N' || !sign)
            continue;
        int lf, hf;
        row_bounds(r, lo, hi, &lf, &hf, t);
        if ((sign > 0 && !lf) || (sign < 0 && !hf)) {
            snprintf(message, 512, "multiplier selects infinite side of row %s", r->name);
            available = 0;
            break;
        }
        mpq_mul(t, r->dual, sign > 0 ? lo : hi);
        mpq_add(bound, bound, t);
    }
    for (size_t i = 0; i < m->nc && available; i++) {
        Col *c = &m->c[i];
        int sign = mpq_sgn(c->residual);
        if (!sign)
            continue;
        if ((sign > 0 && !c->lf) || (sign < 0 && !c->hf)) {
            snprintf(message, 512, "reduced cost selects infinite side of column %s", c->name);
            available = 0;
            break;
        }
        mpq_mul(t, c->residual, sign > 0 ? c->lo : c->hi);
        mpq_add(bound, bound, t);
    }
    if (m->maximize)
        mpq_neg(bound, bound);
    mpq_clears(lo, hi, t, NULL);
    return available;
}
static void maximum(mpq_t max, const mpq_t a) {
    if (mpq_cmp(a, max) > 0)
        mpq_set(max, a);
}
static void json_q(const mpq_t q) {
    char *s = mpq_get_str(NULL, 10, q);
    lab_json_string(stdout, s);
    free(s);
}
int main(int argc, char **argv) {
    if (argc == 2 && !strcmp(argv[1], "--help")) {
        puts("milpcheck --model FILE [--inspect | --solution FILE | --dual FILE]\n"
             "          [--format simple|rational|highs|scip|cbc] [--relax] [--dual FILE]\n"
             "          [--lin-tol DECIMAL] [--int-tol DECIMAL] [--obj-tol DECIMAL] [--exact]\n"
             "Independent primal and LP-bound checking. Exact matching primal/dual witnesses\n"
             "can certify optimality; no infeasibility/unboundedness proof is implemented.");
        return 0;
    }
    const char *model = NULL, *sol = NULL, *dual = NULL, *format = "simple", *lt = "1e-4",
               *it = "1e-5", *ot = "1e-4";
    int inspect = 0, relax = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--inspect"))
            inspect = 1;
        else if (!strcmp(argv[i], "--relax"))
            relax = 1;
        else if (!strcmp(argv[i], "--exact"))
            lt = it = ot = "0";
        else if (i + 1 < argc && !strcmp(argv[i], "--model"))
            model = argv[++i];
        else if (i + 1 < argc && !strcmp(argv[i], "--solution"))
            sol = argv[++i];
        else if (i + 1 < argc && !strcmp(argv[i], "--dual"))
            dual = argv[++i];
        else if (i + 1 < argc && !strcmp(argv[i], "--format"))
            format = argv[++i];
        else if (i + 1 < argc && !strcmp(argv[i], "--lin-tol"))
            lt = argv[++i];
        else if (i + 1 < argc && !strcmp(argv[i], "--int-tol"))
            it = argv[++i];
        else if (i + 1 < argc && !strcmp(argv[i], "--obj-tol"))
            ot = argv[++i];
        else {
            fprintf(stderr, "unknown/missing argument: %s\n", argv[i]);
            return 64;
        }
    }
    if (!model || (!inspect && !sol && !dual) || (inspect && (sol || dual)) ||
        (strcmp(format, "simple") && strcmp(format, "highs") && strcmp(format, "scip") &&
         strcmp(format, "cbc") && strcmp(format, "rational"))) {
        fputs("usage: milpcheck --model FILE [--inspect | --solution FILE | --dual FILE] "
              "[--format simple|rational|highs|scip|cbc] [--relax] [--dual FILE] [--exact]\n",
              stderr);
        return 64;
    }
    Model m = {0};
    m.objective = SIZE_MAX;
    mpq_init(m.obj);
    mpq_t lin, integ, objtol, claimed, lo, hi, t, u, maxrow, maxbound, maxint, objerr, bound;
    mpq_inits(lin, integ, objtol, claimed, lo, hi, t, u, maxrow, maxbound, maxint, objerr, bound,
              NULL);
    if (number(lin, lt) || number(integ, it) || number(objtol, ot) || mpq_sgn(lin) < 0 ||
        mpq_sgn(integ) < 0 || mpq_sgn(objtol) < 0)
        fail(&m, "invalid tolerance");
    if (!m.error[0])
        read_mps(&m, model, 0);
    if (!m.error[0] && dual)
        read_dual(&m, dual);
    int exitcode = 0, has_claimed = 0, nosol = 0, bound_valid = 0;
    char bound_message[512] = "";
    size_t rb = 0, bb = 0, ib = 0;
    if (!m.error[0] && sol)
        nosol = solution(&m, sol, format, claimed, &has_claimed);
    if (!m.error[0] && !inspect && ((!nosol && sol) || dual))
        read_mps(&m, model, 1);
    if (!m.error[0] && dual)
        bound_valid = certified_bound(&m, bound, bound_message);
    if (!m.error[0] && sol && !nosol) {
        for (size_t i = 0; i < m.nc; i++) {
            Col *c = &m.c[i];
            mpq_set_ui(t, 0, 1);
            if (c->lf) {
                mpq_sub(u, c->lo, c->value);
                maximum(t, u);
            }
            if (c->hf) {
                mpq_sub(u, c->value, c->hi);
                maximum(t, u);
            }
            maximum(maxbound, t);
            if (mpq_cmp(t, lin) > 0)
                bb++;
            if (c->integer && !relax) {
                mpz_t z;
                mpz_init(z);
                mpz_fdiv_q(z, mpq_numref(c->value), mpq_denref(c->value));
                mpq_set_z(t, z);
                mpq_sub(t, c->value, t);
                mpq_set_ui(u, 1, 1);
                mpq_sub(u, u, t);
                if (mpq_cmp(u, t) < 0)
                    mpq_set(t, u);
                maximum(maxint, t);
                if (mpq_cmp(t, integ) > 0)
                    ib++;
                mpz_clear(z);
            }
        }
        for (size_t i = 0; i < m.nr; i++) {
            Row *r = &m.r[i];
            if (r->type == 'N')
                continue;
            int lf, hf;
            row_bounds(r, lo, hi, &lf, &hf, t);
            mpq_set_ui(t, 0, 1);
            if (lf) {
                mpq_sub(u, lo, r->activity);
                maximum(t, u);
            }
            if (hf) {
                mpq_sub(u, r->activity, hi);
                maximum(t, u);
            }
            maximum(maxrow, t);
            if (mpq_cmp(t, lin) > 0)
                rb++;
        }
        mpq_sub(m.obj, m.r[m.objective].activity, m.r[m.objective].rhs);
        if (has_claimed) {
            mpq_sub(objerr, m.obj, claimed);
            mpq_abs(objerr, objerr);
        }
        if (rb || bb || ib || (has_claimed && mpq_cmp(objerr, objtol) > 0))
            exitcode = 1;
    }
    if (m.error[0])
        exitcode = m.unsupported ? 3 : 2;
    else if (nosol)
        exitcode = 4;
    else if (dual && !sol && !bound_valid)
        exitcode = 4;
    int primal_valid = sol && !inspect && !exitcode;
    int optimal = primal_valid && bound_valid && !mpq_sgn(lin) && !mpq_sgn(integ) &&
                  !mpq_sgn(objtol) && !mpq_cmp(m.obj, bound);
    printf("{\"schema\":1,\"status\":");
    lab_json_string(stdout, m.error[0] ? (m.unsupported ? "unsupported" : "error")
                            : nosol    ? "no_solution"
                            : inspect  ? "supported"
                            : dual && !sol && !bound_valid ? "unavailable"
                            : exitcode                     ? "invalid"
                                                           : "valid");
    printf(",\"message\":");
    lab_json_string(stdout, m.error);
    printf(",\"rows\":%zu,\"columns\":%zu,\"nonzeros\":%zu,\"objective_sense\":\"%s\",\"checked_"
           "feasible\":%s,\"relaxation\":%s,\"certifies_optimality\":%s",
           m.nr, m.nc, m.nnz, m.maximize ? "max" : "min", primal_valid ? "true" : "false",
           relax ? "true" : "false", optimal ? "true" : "false");
    printf(",\"bound_status\":");
    lab_json_string(stdout, m.error[0] || !dual ? "absent" : bound_valid ? "valid" : "unavailable");
    printf(",\"bound_direction\":\"%s\",\"certified_bound\":", m.maximize ? "upper" : "lower");
    if (!m.error[0] && bound_valid)
        json_q(bound);
    else
        printf("null");
    printf(",\"bound_message\":");
    lab_json_string(stdout, bound_message);
    printf(",\"objective\":");
    if (sol && !m.error[0] && !nosol)
        json_q(m.obj);
    else
        printf("null");
    printf(",\"row_violations\":%zu,\"bound_violations\":%zu,\"integrality_violations\":%zu,\"max_"
           "row_violation\":",
           rb, bb, ib);
    json_q(maxrow);
    printf(",\"max_bound_violation\":");
    json_q(maxbound);
    printf(",\"max_integrality_violation\":");
    json_q(maxint);
    printf(",\"objective_error\":");
    json_q(objerr);
    printf(",\"linear_tolerance\":");
    json_q(lin);
    printf(",\"integer_tolerance\":");
    json_q(integ);
    printf(",\"objective_tolerance\":");
    json_q(objtol);
    printf("}\n");
    cleanup(&m);
    mpq_clears(lin, integ, objtol, claimed, lo, hi, t, u, maxrow, maxbound, maxint, objerr, bound,
               NULL);
    return exitcode;
}

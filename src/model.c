#define _POSIX_C_SOURCE 200809L
#include "model.h"
#include "common.h"
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

/* This reader belongs to the optimizer. It deliberately does not call or link
 * the independent checker's model reader or arithmetic implementation. */
enum { MODEL_ITEMS = 5000000, MODEL_ENTRIES = 50000000, LINE_SIZE = 16384 };
typedef struct {
    const char *name;
    size_t index;
} NameSlot;
typedef struct {
    NameSlot *slots;
    size_t capacity, used;
} Names;
typedef struct {
    char *name;
    double rhs, range;
    char kind;
    unsigned rhs_seen : 1, range_seen : 1;
} InputRow;
typedef struct {
    char *name;
    double lower, upper;
    unsigned integer : 1, marker : 1, marker_seen : 1;
    unsigned bound_seen : 1, explicit_lower : 1;
} InputCol;
typedef struct {
    size_t column, row;
    double value;
} Entry;
typedef struct {
    LabModel *output;
    InputRow *rows;
    InputCol *cols;
    Entry *entries;
    size_t nr, nc, ne, row_cap, col_cap, entry_cap;
    Names row_names, col_names;
    char *objective_name, *rhs_name, *range_name, *bound_name;
    size_t line;
} Reader;

static int problem(Reader *p, const char *format, ...) {
    if (!p->output->error[0]) {
        size_t used =
            (size_t)snprintf(p->output->error, sizeof p->output->error, "MPS line %zu: ", p->line);
        if (used < sizeof p->output->error) {
            va_list ap;
            va_start(ap, format);
            vsnprintf(p->output->error + used, sizeof p->output->error - used, format, ap);
            va_end(ap);
        }
    }
    return -1;
}
static int unsupported(Reader *p, const char *what) {
    p->output->unsupported = 1;
    return problem(p, "unsupported %s", what);
}
static size_t name_hash(const char *text) {
    size_t h = 5381;
    while (*text)
        h = (h * 33) ^ (unsigned char)*text++;
    h ^= h >> 17;
    return h;
}
static size_t find_name(const Names *names, const char *text) {
    if (!names->capacity)
        return SIZE_MAX;
    size_t k = name_hash(text) % names->capacity;
    while (names->slots[k].name) {
        if (!strcmp(names->slots[k].name, text))
            return names->slots[k].index;
        k = (k + 1) % names->capacity;
    }
    return SIZE_MAX;
}
static void save_name(Names *names, const char *text, size_t value) {
    if (!names->capacity || names->used >= names->capacity / 2) {
        size_t cap = names->capacity ? 2 * names->capacity : 128;
        NameSlot *next = lab_alloc(cap, sizeof *next);
        for (size_t i = 0; i < names->capacity; ++i) {
            if (!names->slots[i].name)
                continue;
            size_t k = name_hash(names->slots[i].name) % cap;
            while (next[k].name)
                k = (k + 1) % cap;
            next[k] = names->slots[i];
        }
        free(names->slots);
        names->slots = next;
        names->capacity = cap;
    }
    size_t k = name_hash(text) % names->capacity;
    while (names->slots[k].name)
        k = (k + 1) % names->capacity;
    names->slots[k].name = text;
    names->slots[k].index = value;
    ++names->used;
}
static int reserve(Reader *p, void **data, size_t *capacity, size_t count, size_t item,
                   size_t limit) {
    if (count <= *capacity)
        return 0;
    if (count > limit)
        return problem(p, "model exceeds reader resource limit");
    size_t cap = *capacity ? *capacity : 64;
    while (cap < count)
        cap = cap > limit / 2 ? limit : 2 * cap;
    if (cap > SIZE_MAX / item)
        return problem(p, "allocation size overflow");
    void *next = realloc(*data, cap * item);
    if (!next)
        return problem(p, "cannot allocate model storage");
    *data = next;
    *capacity = cap;
    return 0;
}
static int split(char *line, char **field, size_t *count) {
    *count = 0;
    char *save = NULL;
    for (char *s = strtok_r(line, " \t\r\n", &save); s; s = strtok_r(NULL, " \t\r\n", &save)) {
        if (*s == '$')
            break;
        if (*count == 8)
            return -1;
        field[(*count)++] = s;
    }
    return 0;
}
static int finite_number(Reader *p, const char *text, double *value) {
    size_t n = strlen(text);
    if (!n || n >= LINE_SIZE)
        return problem(p, "invalid numeric field");
    char copy[LINE_SIZE];
    for (size_t i = 0; i <= n; ++i)
        copy[i] = text[i] == 'D' || text[i] == 'd' ? 'e' : text[i];
    /* strtod accepts hexadecimal floats and NaNs, neither of which is MPS
     * decimal syntax. Validate the grammar before converting. */
    const char *s = copy;
    if (*s == '-' || *s == '+')
        ++s;
    int digits = 0;
    while (isdigit((unsigned char)*s)) {
        ++s;
        ++digits;
    }
    if (*s == '.') {
        ++s;
        while (isdigit((unsigned char)*s)) {
            ++s;
            ++digits;
        }
    }
    if (!digits)
        return problem(p, "invalid decimal %s", text);
    if (*s == 'e' || *s == 'E') {
        ++s;
        if (*s == '-' || *s == '+')
            ++s;
        const char *first = s;
        while (isdigit((unsigned char)*s))
            ++s;
        if (s == first)
            return problem(p, "invalid exponent %s", text);
    }
    if (*s)
        return problem(p, "invalid decimal %s", text);
    errno = 0;
    char *end = NULL;
    double result = strtod(copy, &end);
    if (errno || !end || *end || !isfinite(result))
        return problem(p, "numeric field is outside finite double range: %s", text);
    *value = result;
    return 0;
}
static size_t get_column(Reader *p, const char *name, int marker) {
    size_t j = find_name(&p->col_names, name);
    if (j != SIZE_MAX) {
        if (marker >= 0 && p->cols[j].marker_seen && p->cols[j].marker != (unsigned)marker) {
            problem(p, "conflicting integer markers for %s", name);
            return SIZE_MAX;
        }
        return j;
    }
    if (reserve(p, (void **)&p->cols, &p->col_cap, p->nc + 1, sizeof *p->cols, MODEL_ITEMS))
        return SIZE_MAX;
    j = p->nc++;
    InputCol *col = &p->cols[j];
    memset(col, 0, sizeof *col);
    col->name = lab_strdup(name);
    col->lower = 0;
    col->upper = marker > 0 ? 1 : INFINITY;
    col->integer = marker > 0;
    col->marker = marker > 0;
    col->marker_seen = marker >= 0;
    save_name(&p->col_names, col->name, j);
    return j;
}
static int select_vector(Reader *p, char **selected, const char *name) {
    if (!*selected)
        *selected = lab_strdup(name);
    else if (strcmp(*selected, name))
        return unsupported(p, "multiple named RHS/range/bound vectors");
    return 0;
}
static int add_row(Reader *p, size_t n, char **field) {
    if (n != 2 || strlen(field[0]) != 1 || !strchr("NLEG", field[0][0]))
        return problem(p, "malformed or unsupported row");
    if (find_name(&p->row_names, field[1]) != SIZE_MAX)
        return problem(p, "duplicate row %s", field[1]);
    if (reserve(p, (void **)&p->rows, &p->row_cap, p->nr + 1, sizeof *p->rows, MODEL_ITEMS))
        return -1;
    InputRow *r = &p->rows[p->nr];
    memset(r, 0, sizeof *r);
    r->name = lab_strdup(field[1]);
    r->kind = field[0][0];
    save_name(&p->row_names, r->name, p->nr++);
    return 0;
}
static int add_values(Reader *p, size_t column, size_t start, size_t n, char **field) {
    for (size_t i = start; i < n; i += 2) {
        size_t row = find_name(&p->row_names, field[i]);
        double number;
        if (row == SIZE_MAX)
            return problem(p, "unknown coefficient row %s", field[i]);
        if (finite_number(p, field[i + 1], &number))
            return -1;
        if (reserve(p, (void **)&p->entries, &p->entry_cap, p->ne + 1, sizeof *p->entries,
                    MODEL_ENTRIES))
            return -1;
        p->entries[p->ne++] = (Entry){column, row, number};
    }
    return 0;
}
static int row_side(Reader *p, size_t n, char **field, int ranges) {
    if (n != 3 && n != 5)
        return problem(p, "RHS/RANGES requires one named vector and row/value pairs");
    if (select_vector(p, ranges ? &p->range_name : &p->rhs_name, field[0]))
        return -1;
    for (size_t k = 1; k < n; k += 2) {
        size_t i = find_name(&p->row_names, field[k]);
        double value;
        if (i == SIZE_MAX)
            return problem(p, "unknown RHS/RANGES row %s", field[k]);
        if (finite_number(p, field[k + 1], &value))
            return -1;
        InputRow *r = &p->rows[i];
        if (ranges ? r->range_seen : r->rhs_seen)
            return problem(p, "duplicate RHS/RANGES entry for %s", r->name);
        if (ranges) {
            r->range = value;
            r->range_seen = 1;
        } else {
            r->rhs = value;
            r->rhs_seen = 1;
        }
    }
    return 0;
}
static int bound(Reader *p, size_t n, char **field) {
    if (n < 3 || n > 4)
        return problem(p, "malformed BOUNDS record");
    if (select_vector(p, &p->bound_name, field[1]))
        return -1;
    const char *kind = field[0];
    int numeric = !strcmp(kind, "LO") || !strcmp(kind, "UP") || !strcmp(kind, "FX") ||
                  !strcmp(kind, "LI") || !strcmp(kind, "UI");
    if (!numeric && strcmp(kind, "FR") && strcmp(kind, "MI") && strcmp(kind, "PL") &&
        strcmp(kind, "BV"))
        return unsupported(p, "bound type");
    if (numeric && n != 4)
        return problem(p, "missing bound value");
    double value = 0;
    if (n == 4 && finite_number(p, field[3], &value))
        return -1;
    size_t j = get_column(p, field[2], -1);
    if (j == SIZE_MAX)
        return -1;
    InputCol *c = &p->cols[j];
    if (!c->bound_seen) {
        /* Integer marker columns have an implicit upper bound of one only
         * when the column has no explicit BOUNDS records. */
        if (c->integer)
            c->upper = INFINITY;
        c->bound_seen = 1;
    }
    if (!strcmp(kind, "LO") || !strcmp(kind, "LI")) {
        c->lower = value;
        c->explicit_lower = 1;
        if (!strcmp(kind, "LI"))
            c->integer = 1;
    } else if (!strcmp(kind, "UP") || !strcmp(kind, "UI")) {
        c->upper = value;
        if (!c->explicit_lower && value < 0)
            c->lower = -INFINITY;
        if (!strcmp(kind, "UI"))
            c->integer = 1;
    } else if (!strcmp(kind, "FX")) {
        c->lower = c->upper = value;
        c->explicit_lower = 1;
    } else if (!strcmp(kind, "BV")) {
        c->integer = 1;
        c->lower = 0;
        c->upper = 1;
        c->explicit_lower = 1;
    } else if (!strcmp(kind, "FR")) {
        c->lower = -INFINITY;
        c->upper = INFINITY;
        c->explicit_lower = 1;
    } else if (!strcmp(kind, "MI")) {
        c->lower = -INFINITY;
        c->explicit_lower = 1;
    } else if (!strcmp(kind, "PL")) {
        c->upper = INFINITY;
    }
    return 0;
}
static int entry_order(const void *left, const void *right) {
    const Entry *a = left, *b = right;
    if (a->column != b->column)
        return a->column < b->column ? -1 : 1;
    if (a->row != b->row)
        return a->row < b->row ? -1 : 1;
    return 0;
}
static int finish(Reader *p) {
    LabModel *m = p->output;
    size_t objective = SIZE_MAX;
    if (p->objective_name)
        objective = find_name(&p->row_names, p->objective_name);
    else
        for (size_t i = 0; i < p->nr; ++i)
            if (p->rows[i].kind == 'N') {
                objective = i;
                break;
            }
    if (objective == SIZE_MAX || p->rows[objective].kind != 'N')
        return problem(p, "missing or invalid objective row");
    size_t *index = lab_alloc(p->nr, sizeof *index);
    for (size_t i = 0; i < p->nr; ++i)
        index[i] = p->rows[i].kind == 'N' ? SIZE_MAX : m->rows++;
    m->cols = p->nc;
    m->offset = -p->rows[objective].rhs;
    m->row_name = lab_alloc(m->rows, sizeof *m->row_name);
    m->row_lower = lab_alloc(m->rows, sizeof *m->row_lower);
    m->row_upper = lab_alloc(m->rows, sizeof *m->row_upper);
    m->col_name = lab_alloc(m->cols, sizeof *m->col_name);
    m->col_lower = lab_alloc(m->cols, sizeof *m->col_lower);
    m->col_upper = lab_alloc(m->cols, sizeof *m->col_upper);
    m->objective = lab_alloc(m->cols, sizeof *m->objective);
    m->integer = lab_alloc(m->cols, sizeof *m->integer);
    m->col_start = lab_alloc(m->cols + 1, sizeof *m->col_start);
    for (size_t i = 0; i < p->nr; ++i) {
        InputRow *r = &p->rows[i];
        if (index[i] == SIZE_MAX)
            continue;
        size_t k = index[i];
        m->row_name[k] = r->name;
        r->name = NULL;
        double lo = r->kind == 'L' ? -INFINITY : r->rhs;
        double hi = r->kind == 'G' ? INFINITY : r->rhs;
        if (r->range_seen) {
            if (r->kind == 'L' || (r->kind == 'E' && r->range < 0))
                lo = r->rhs - fabs(r->range);
            else
                hi = r->rhs + fabs(r->range);
            if (!isfinite(lo) || !isfinite(hi)) {
                free(index);
                return problem(p, "row range overflows finite double range");
            }
        }
        m->row_lower[k] = lo;
        m->row_upper[k] = hi;
    }
    for (size_t j = 0; j < p->nc; ++j) {
        m->col_name[j] = p->cols[j].name;
        p->cols[j].name = NULL;
        m->col_lower[j] = p->cols[j].lower;
        m->col_upper[j] = p->cols[j].upper;
        m->integer[j] = (unsigned char)p->cols[j].integer;
    }
    /* Sorting gives deterministic coalescing even if a writer interleaves
     * column records. Rows within every final CSC column are sorted. */
    if (p->ne)
        qsort(p->entries, p->ne, sizeof *p->entries, entry_order);
    size_t unique = 0;
    for (size_t e = 0; e < p->ne;) {
        Entry a = p->entries[e++];
        while (e < p->ne && a.column == p->entries[e].column && a.row == p->entries[e].row) {
            a.value += p->entries[e++].value;
            if (!isfinite(a.value)) {
                free(index);
                return problem(p, "duplicate coefficients overflow finite double range");
            }
        }
        if (a.row == objective)
            m->objective[a.column] = a.value;
        else if (index[a.row] != SIZE_MAX && a.value != 0) {
            a.row = index[a.row];
            p->entries[unique++] = a;
            ++m->col_start[a.column + 1];
        }
    }
    m->nnz = unique;
    m->row_index = lab_alloc(unique, sizeof *m->row_index);
    m->value = lab_alloc(unique, sizeof *m->value);
    for (size_t e = 0; e < unique; ++e) {
        m->row_index[e] = p->entries[e].row;
        m->value[e] = p->entries[e].value;
    }
    for (size_t j = 0; j < m->cols; ++j)
        m->col_start[j + 1] += m->col_start[j];
    free(index);
    return 0;
}
static void cleanup(Reader *p) {
    for (size_t i = 0; i < p->nr; ++i)
        free(p->rows[i].name);
    for (size_t j = 0; j < p->nc; ++j)
        free(p->cols[j].name);
    free(p->rows);
    free(p->cols);
    free(p->entries);
    free(p->row_names.slots);
    free(p->col_names.slots);
    free(p->objective_name);
    free(p->rhs_name);
    free(p->range_name);
    free(p->bound_name);
}
int lab_model_read(LabModel *model, const char *path) {
    Reader p = {.output = model};
    enum { NONE, ROWS, COLUMNS, RHS, RANGES, BOUNDS, SENSE, OBJECTIVE } section = NONE;
    int seen_rows = 0, seen_columns = 0, ended = 0, marker = 0;
    int seen_rhs = 0, seen_ranges = 0, seen_bounds = 0;
    gzFile input = gzopen(path, "rb");
    if (!input)
        return problem(&p, "cannot open model %s", path);
    gzbuffer(input, 1 << 20);
    char line[LINE_SIZE], previous[LINE_SIZE] = "";
    z_off_t position = gztell(input);
    while (!model->error[0] && gzgets(input, line, sizeof line)) {
        ++p.line;
        size_t length = strlen(line);
        z_off_t next_position = gztell(input);
        if (position < 0 || next_position < position ||
            (uintmax_t)(next_position - position) != (uintmax_t)length) {
            problem(&p, "binary or unreadable MPS record");
            break;
        }
        position = next_position;
        if (length == sizeof line - 1 && line[sizeof line - 2] != '\n') {
            problem(&p, "record exceeds line resource limit");
            break;
        }
        if (line[0] == '*')
            continue;
        int header = !isspace((unsigned char)line[0]);
        char *field[8];
        size_t n;
        if (split(line, field, &n)) {
            problem(&p, "too many fields");
            break;
        }
        if (!n)
            continue;
        if (header) {
            const char *word = field[0];
            if (!strcmp(word, "NAME")) {
                if (seen_rows || n > 2)
                    problem(&p, "misplaced NAME section");
                continue;
            }
            if (!strcmp(word, "ENDATA")) {
                if (n != 1)
                    problem(&p, "malformed ENDATA");
                ended = 1;
                break;
            }
            if (!strcmp(word, "OBJSENSE") || !strcmp(word, "OBJNAME")) {
                if (seen_columns || n > 2) {
                    problem(&p, "misplaced objective header");
                    break;
                }
                section = !strcmp(word, "OBJSENSE") ? SENSE : OBJECTIVE;
                if (n == 1)
                    continue;
                field[0] = field[1];
                n = 1;
            } else if (!strcmp(word, "ROWS")) {
                if (n != 1 || seen_rows || seen_columns)
                    problem(&p, "repeated or misplaced ROWS section");
                section = ROWS;
                seen_rows = 1;
                continue;
            } else if (!strcmp(word, "COLUMNS")) {
                if (n != 1 || !seen_rows || seen_columns)
                    problem(&p, "repeated or misplaced COLUMNS section");
                section = COLUMNS;
                seen_columns = 1;
                continue;
            } else if (!strcmp(word, "RHS") || !strcmp(word, "RANGES") || !strcmp(word, "BOUNDS")) {
                if (n != 1 || !seen_columns || marker) {
                    problem(&p, "misplaced RHS/RANGES/BOUNDS section");
                    break;
                }
                int *seen = !strcmp(word, "RHS")      ? &seen_rhs
                            : !strcmp(word, "RANGES") ? &seen_ranges
                                                      : &seen_bounds;
                if (*seen) {
                    problem(&p, "repeated RHS/RANGES/BOUNDS section");
                    break;
                }
                *seen = 1;
                section = !strcmp(word, "RHS") ? RHS : !strcmp(word, "RANGES") ? RANGES : BOUNDS;
                continue;
            } else if (!strcmp(word, "SOS") || !strcmp(word, "INDICATORS") ||
                       !strcmp(word, "CSECTION") || !strcmp(word, "QMATRIX") ||
                       !strcmp(word, "QSECTION") || !strcmp(word, "QUADOBJ") ||
                       !strcmp(word, "QCMATRIX")) {
                unsupported(&p, "model section");
                break;
            }
        }
        if (section == SENSE) {
            if (n != 1 || (strcmp(field[0], "MIN") && strcmp(field[0], "MINIMIZE") &&
                           strcmp(field[0], "MAX") && strcmp(field[0], "MAXIMIZE")))
                problem(&p, "invalid objective sense");
            else
                model->maximize = !strncmp(field[0], "MAX", 3);
        } else if (section == OBJECTIVE) {
            if (n != 1 || p.objective_name)
                problem(&p, "duplicate or malformed OBJNAME");
            else
                p.objective_name = lab_strdup(field[0]);
        } else if (section == ROWS) {
            add_row(&p, n, field);
        } else if (section == COLUMNS) {
            if (n == 3 && !strcmp(field[1], "'MARKER'")) {
                if (!strcmp(field[2], "'INTORG'") && !marker)
                    marker = 1;
                else if (!strcmp(field[2], "'INTEND'") && marker)
                    marker = 0;
                else
                    problem(&p, "unbalanced integer marker");
            } else {
                size_t start;
                const char *name;
                if (n == 3 || n == 5) {
                    name = field[0];
                    snprintf(previous, sizeof previous, "%s", name);
                    start = 1;
                } else if ((n == 2 || n == 4) && previous[0]) {
                    name = previous;
                    start = 0;
                } else {
                    problem(&p, "malformed COLUMNS record");
                    break;
                }
                size_t j = get_column(&p, name, marker);
                if (j == SIZE_MAX)
                    break;
                add_values(&p, j, start, n, field);
            }
        } else if (section == RHS || section == RANGES) {
            row_side(&p, n, field, section == RANGES);
        } else if (section == BOUNDS) {
            bound(&p, n, field);
        } else {
            problem(&p, "record outside supported section");
        }
    }
    int zstatus;
    (void)gzerror(input, &zstatus);
    if (zstatus != Z_OK && zstatus != Z_STREAM_END)
        problem(&p, "compressed model read failed");
    gzclose(input);
    if (!seen_rows || !seen_columns || !ended || marker)
        problem(&p, "truncated or incomplete model");
    if (!model->error[0])
        finish(&p);
    cleanup(&p);
    return model->error[0] ? -1 : 0;
}
void lab_model_free(LabModel *model) {
    for (size_t i = 0; i < model->rows; ++i)
        free(model->row_name ? model->row_name[i] : NULL);
    for (size_t j = 0; j < model->cols; ++j)
        free(model->col_name ? model->col_name[j] : NULL);
    free(model->row_name);
    free(model->col_name);
    free(model->row_lower);
    free(model->row_upper);
    free(model->col_lower);
    free(model->col_upper);
    free(model->objective);
    free(model->integer);
    free(model->col_start);
    free(model->row_index);
    free(model->value);
    memset(model, 0, sizeof *model);
}

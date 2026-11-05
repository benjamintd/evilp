#define _POSIX_C_SOURCE 200809L
#include "model.h"
#include "exact.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void require(int condition, const char *text) {
    if (!condition) {
        fprintf(stderr, "exact model contract failed: %s\n", text);
        exit(1);
    }
}
static void read_text(LabModel *m, const char *text, int success) {
    char path[] = "/tmp/milp-exact-model-XXXXXX";
    int fd = mkstemp(path);
    require(fd >= 0, "temporary model");
    FILE *f = fdopen(fd, "w");
    require(f != NULL, "temporary stream");
    require(fputs(text, f) >= 0 && !fclose(f), "write model");
    int failed = lab_model_read(m, path);
    unlink(path);
    if (success && failed)
        fprintf(stderr, "%s\n", m->error);
    require(success ? !failed : failed != 0, "model read outcome");
    require(m->parsed_decimal, "parsed model provenance origin");
}
static void decimal_parser(void) {
    mpq_t q;
    mpq_init(q);
    require(!lab_exact_decimal(q, "+.100D+1") && mpq_cmp_ui(q, 1, 1) == 0,
            "D exponent parsed exactly");
    require(!lab_exact_decimal(q, "-12.5e-1") && mpq_cmp_si(q, -5, 4) == 0,
            "signed decimal parsed exactly");
    require(!lab_exact_decimal(q, "0.000") && mpq_sgn(q) == 0, "canonical zero");
    const char *invalid[] = {"", "+", ".", "1.2.3", "1e", "1e+", "1/3", "nan",
                             "0x1p0", "1e4097", "0e99999999999999999999999999"};
    for (size_t i = 0; i < sizeof invalid / sizeof *invalid; i++)
        require(lab_exact_decimal(q, invalid[i]) != 0, "invalid/bounded decimal rejected");
    mpq_clear(q);
}
static void source_arithmetic(void) {
    LabModel m = {0};
    read_text(&m,
        "NAME DECIMAL\nROWS\n N O\n L R\nCOLUMNS\n x O 1 R .1\n"
        " x R .2\n y O 0 R .7\nRHS\n b R .3\nRANGES\n g R .2\n"
        "BOUNDS\n LO b x .1\n UP b x .3\n FR b y\nENDATA\n", 1);
    require(m.exact != NULL, "small parsed model retains exact provenance");
    require(m.exact->rows == 1 && m.exact->cols == 2 && m.exact->nnz == 2,
            "exact dimensions match CSC");
    require(mpq_cmp_ui(m.exact->value[0], 3, 10) == 0,
            "duplicate decimal .1+.2 is exactly .3");
    require(m.value[0] != .3, "test exposes binary duplicate rounding");
    require(mpq_cmp_ui(m.exact->row_lower[0], 1, 10) == 0 &&
                mpq_cmp_ui(m.exact->row_upper[0], 3, 10) == 0,
            "ranged bounds retain exact subtraction");
    require(mpq_cmp_ui(m.exact->col_lower[0], 1, 10) == 0 &&
                mpq_cmp_ui(m.exact->col_upper[0], 3, 10) == 0,
            "column decimal sides exact");
    require(m.exact->col_lower_source[0] == m.col_lower[0] &&
                m.exact->col_upper_source[0] == m.col_upper[0],
            "source float bounds identify later overlays");
    require(!isfinite(m.col_lower[1]) && !isfinite(m.col_upper[1]) &&
                mpq_sgn(m.exact->col_lower[1]) == 0 && mpq_sgn(m.exact->col_upper[1]) == 0,
            "infinite sides use zero placeholders");
    lab_model_free(&m);
}
static void topology_disagreement(void) {
    const char *models[] = {
        "NAME LOST\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 10000000000000001\n"
        " x R -10000000000000000\nRHS\n b R 1\nENDATA\n",
        "NAME EXTRA\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R .1\n"
        " x R .2\n x R -.3\nRHS\n b R 1\nENDATA\n",
        "NAME RANGE\nROWS\n N O\n L R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R .3\nRANGES\n g R 1e-20\nENDATA\n",
        "NAME BOUND\nROWS\n N O\nCOLUMNS\n x O 0\n"
        "BOUNDS\n LO b x .1\n UP b x .100000000000000001\nENDATA\n"
    };
    for (size_t i = 0; i < sizeof models / sizeof *models; i++) {
        LabModel m = {0};
        read_text(&m, models[i], 1);
        require(m.exact == NULL, "float/exact topology mismatch disables certificates");
        lab_model_free(&m);
    }
}
static void defaults_and_nonconstraints(void) {
    LabModel m = {0};
    read_text(&m,
        "NAME DEFAULTS\nOBJNAME O\nROWS\n N IGNORE\n N O\n E R\nCOLUMNS\n"
        " mark 'MARKER' 'INTORG'\n x O 0 R 1\n y O 0 R 1\n"
        " x IGNORE 10000000000000001\n x IGNORE -10000000000000000\n"
        " end 'MARKER' 'INTEND'\nRHS\n b R .3\nRANGES\n g R -.2\n"
        "BOUNDS\n LI b y .1\nENDATA\n", 1);
    require(m.exact != NULL, "ignored nonconstraint coefficients do not affect support");
    require(mpq_cmp_ui(m.exact->row_lower[0], 1, 10) == 0 &&
                mpq_cmp_ui(m.exact->row_upper[0], 3, 10) == 0,
            "negative equality range exact");
    require(mpq_sgn(m.exact->col_lower[0]) == 0 &&
                mpq_cmp_ui(m.exact->col_upper[0], 1, 1) == 0,
            "implicit marker binary upper bound retained");
    require(isinf(m.col_upper[1]) && mpq_cmp_ui(m.exact->col_lower[1], 1, 10) == 0,
            "explicit integer bound clears implicit upper bound");
    lab_model_free(&m);
    read_text(&m,
        "NAME NEGATIVE\nROWS\n N O\nCOLUMNS\n x O 0\n"
        "BOUNDS\n UP b x -.1\nENDATA\n", 1);
    require(m.exact != NULL && isinf(m.col_lower[0]) &&
                mpq_cmp_si(m.exact->col_upper[0], -1, 10) == 0,
            "negative upper-bound default lower side remains infinite");
    lab_model_free(&m);
}
static void bounded_storage_and_cleanup(void) {
    LabModel m = {0};
    read_text(&m,
        "NAME HUGEEXP\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 0e99999999\n"
        "RHS\n b R 0\nENDATA\n", 1);
    require(m.exact == NULL, "exact exponent resource limit leaves ordinary reader intact");
    lab_model_free(&m);
    read_text(&m,
        "NAME BAD\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R .1\n"
        "RHS\n b R .2\n b R .3\nENDATA\n", 0);
    lab_model_free(&m);
    size_t capacity = 1000000;
    char *text = calloc(capacity, 1);
    require(text != NULL, "storage test allocation");
    size_t used = (size_t)snprintf(text, capacity,
        "NAME LIMIT\nROWS\n N O\n G R\nCOLUMNS\n x O 0\n");
    for (size_t k = 0; k < 20001; k++)
        used += (size_t)snprintf(text + used, capacity - used, " x R .1\n");
    snprintf(text + used, capacity - used, "RHS\n b R 1\nENDATA\n");
    read_text(&m, text, 1);
    require(m.exact == NULL && m.nnz == 1, "entry budget disables only exact provenance");
    lab_model_free(&m);
    used = (size_t)snprintf(text, capacity, "NAME COLUMNS\nROWS\n N O\nCOLUMNS\n");
    for (size_t k = 0; k < 257; k++)
        used += (size_t)snprintf(text + used, capacity - used, " x%zu O 0\n", k);
    snprintf(text + used, capacity - used, "ENDATA\n");
    read_text(&m, text, 1);
    require(m.exact == NULL && m.cols == 257, "column budget leaves float model intact");
    lab_model_free(&m);
    used = (size_t)snprintf(text, capacity, "NAME ROWS\nROWS\n N O\n");
    for (size_t k = 0; k < 129; k++)
        used += (size_t)snprintf(text + used, capacity - used, " G R%zu\n", k);
    snprintf(text + used, capacity - used, "COLUMNS\n x O 0\nENDATA\n");
    read_text(&m, text, 1);
    require(m.exact == NULL && m.rows == 129, "row budget leaves float model intact");
    lab_model_free(&m);
    used = (size_t)snprintf(text, capacity,
        "NAME LONGBOUND\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R .1\n"
        "RHS\n b R 0\nBOUNDS\n LO b x .2\n UP b x 0.");
    memset(text + used, '0', 1100);
    used += 1100;
    snprintf(text + used, capacity - used, "\nENDATA\n");
    read_text(&m, text, 1);
    require(m.exact == NULL, "oversized bound lexeme permanently disables provenance");
    lab_model_free(&m);
    free(text);
}
int main(void) {
    decimal_parser();
    source_arithmetic();
    topology_disagreement();
    defaults_and_nonconstraints();
    bounded_storage_and_cleanup();
    puts("Original-decimal provenance, exact arithmetic, topology and resource contracts passed.");
    return 0;
}

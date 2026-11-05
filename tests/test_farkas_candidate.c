#define _POSIX_C_SOURCE 200809L
#include "farkas_candidate.h"
#include "exact.h"
#include "common.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

static void require(int ok, const char *message) {
    if (!ok) {
        fprintf(stderr, "source Farkas candidate: %s\n", message);
        exit(1);
    }
}
static void read_model(LabModel *m, const char *text, char path[64]) {
    strcpy(path, "/tmp/milp-source-farkas-XXXXXX");
    int fd = mkstemp(path);
    require(fd >= 0, "temporary model");
    FILE *f = fdopen(fd, "w");
    require(f && fputs(text, f) >= 0 && !fclose(f), "write model");
    if (lab_model_read(m, path)) {
        fprintf(stderr, "%s\n", m->error);
        require(0, "parse model");
    }
}
static void checker_replay(const char *checker, const char *path, const LabModel *m,
                             const LabFarkasCandidate *r, int expect_strictly_larger) {
    require(r->rows == m->rows && r->multiplier && r->margin_lower_bound &&
            !strcmp(r->margin_kind, "source_enclosure"), "candidate ownership and margin scope");
    mpq_t lower, actual;
    mpq_inits(lower, actual, NULL);
    require(!mpq_set_str(lower, r->margin_lower_bound, 10) && mpq_sgn(lower) > 0,
            "strict positive rational lower margin");
    if (!checker) {
        mpq_clears(lower, actual, NULL);
        return;
    }
    char certificate[] = "/tmp/milp-source-farkas-cert-XXXXXX";
    int fd = mkstemp(certificate);
    require(fd >= 0, "temporary candidate");
    FILE *f = fdopen(fd, "w");
    require(f != NULL, "candidate stream");
    fputs("MILPLAB_FARKAS_1\n", f);
    for (size_t i = 0; i < m->rows; i++)
        fprintf(f, "%s %s\n", m->row_name[i], r->multiplier[i]);
    require(!fclose(f), "write candidate");
    char output[] = "/tmp/milp-source-farkas-output-XXXXXX";
    int out = mkstemp(output);
    require(out >= 0, "checker output");
    pid_t pid = fork();
    require(pid >= 0, "start checker");
    if (!pid) {
        dup2(out, STDOUT_FILENO);
        close(out);
        execl(checker, checker, "--model", path, "--farkas", certificate, (char *)NULL);
        _exit(127);
    }
    close(out);
    int status;
    require(waitpid(pid, &status, 0) == pid && WIFEXITED(status) && WEXITSTATUS(status) == 0,
            "independent checker accepts the original-source row combination");
    f = fopen(output, "r");
    require(f != NULL, "read independent result");
    char text[16384];
    size_t count = fread(text, 1, sizeof text - 1, f);
    text[count] = 0;
    fclose(f);
    require(strstr(text, "\"certifies_infeasibility\":true") != NULL,
            "independent checker certifies infeasibility");
    char *margin = strstr(text, "\"contradiction_margin\":\"");
    require(margin != NULL, "independent exact margin present");
    margin += strlen("\"contradiction_margin\":\"");
    char *end = strchr(margin, '"');
    require(end != NULL, "independent margin terminator");
    *end = 0;
    require(!mpq_set_str(actual, margin, 10) && mpq_cmp(actual, lower) >= 0,
            "independent actual margin is at least the claimed lower bound");
    if (expect_strictly_larger)
        require(mpq_cmp(actual, lower) > 0, "conservative lower margin is not an exact-margin claim");
    unlink(certificate); unlink(output);
    mpq_clears(lower, actual, NULL);
}
static void basic_and_dyadic_mapping(const char *checker) {
    LabModel m = {0};
    char path[64];
    read_model(&m,
        "NAME BASIC\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R 2\nBOUNDS\n UP b x 1\nENDATA\n", path);
    size_t source[] = {0, 0, SIZE_MAX};
    double scale[] = {2, 4, 1}, pi[] = {.5, .25, 0};
    LabFarkasCandidate r;
    require(lab_farkas_candidate(&m, 3, source, scale, pi, lab_now() + 5, &r) == 1,
            "bounded contradiction can be validated without a basis");
    require(!strcmp(r.multiplier[0], "2") && !strcmp(r.margin_lower_bound, "2"),
            "duplicate original rows accumulate exact scale*pi products");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    pi[0] = NAN;
    require(lab_farkas_candidate(&m, 3, source, scale, pi, lab_now() + 5, &r) == -1,
            "NaN multiplier is rejected before GMP conversion");
    lab_farkas_candidate_free(&r);
    pi[0] = .5; scale[2] = INFINITY;
    require(lab_farkas_candidate(&m, 3, source, scale, pi, lab_now() + 5, &r) == -1,
            "even ignored-row scale encodings must be finite");
    lab_farkas_candidate_free(&r);
    scale[2] = 1;
    require(lab_farkas_candidate(&m, 3, source, scale, pi, lab_now() - 1, &r) == 0 &&
            !r.multiplier, "expired deadline cannot emit a candidate");
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
    read_model(&m,
        "NAME CONSTANT\nROWS\n N O\n G R\nCOLUMNS\n x O 0\n"
        "RHS\n b R 1\nBOUNDS\n FR b x\nENDATA\n", path);
    for (int exponent = -900; exponent <= 900; exponent += 1800) {
        scale[0] = pi[0] = scalbn(1, exponent);
        require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 1,
                "mapping preserves products that would underflow or overflow binary64");
        checker_replay(checker, path, &m, &r, 0);
        lab_farkas_candidate_free(&r);
    }
    unlink(path); lab_model_free(&m);
}
static char *decimal_large(int feasible) {
    char *text = calloc(16384, 1);
    require(text != NULL, "large source fixture");
    size_t at = (size_t)snprintf(text, 16384, "NAME SOURCE\nROWS\n N O\n G R\n");
    for (size_t i = 0; i < 128; i++)
        at += (size_t)snprintf(text + at, 16384 - at, " L D%zu\n", i);
    snprintf(text + at, 16384 - at,
             "COLUMNS\n x O 0 R .1\nRHS\n b R %s\nBOUNDS\n UP b x 1\nENDATA\n",
             feasible ? ".1" : ".2");
    return text;
}
static void envelope_lower_margin_and_feasible_adversary(const char *checker) {
    size_t source[] = {0};
    double scale[] = {1}, pi[] = {1};
    LabModel m = {0};
    LabFarkasCandidate r;
    char path[64], *text = decimal_large(0);
    read_model(&m, text, path); free(text);
    require(!m.exact && m.envelope && m.rows == 129,
            "source-envelope fixture exceeds small exact provenance size");
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 1,
            "large decimal-source model obtains a strict interval contradiction");
    checker_replay(checker, path, &m, &r, 1);
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
    text = decimal_large(1);
    read_model(&m, text, path); free(text);
    for (int sign = -1; sign <= 1; sign++) {
        pi[0] = sign;
        require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
                "feasible source decimal boundary cannot be certified infeasible");
        lab_farkas_candidate_free(&r);
    }
    unlink(path); lab_model_free(&m);
}
static void free_columns_and_source_box(const char *checker) {
    LabModel m = {0};
    char path[64];
    read_model(&m,
        "NAME CANCEL\nROWS\n N O\n E A\n E B\nCOLUMNS\n x O 0 A 3\n x B 9\n"
        "RHS\n b A 6 B 9\nBOUNDS\n FR b x\nENDATA\n", path);
    size_t source[] = {0, 1};
    double scale[] = {1, 1}, pi[] = {1.0 / 3, -1.0 / 9};
    LabFarkasCandidate r;
    require(lab_farkas_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 1,
            "exact dyadic cancellation gives zero coefficient on a free column");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    pi[1] = nextafter(pi[1], -INFINITY);
    require(lab_farkas_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 0,
            "a genuinely nonzero free-column coefficient is never rounded away");
    lab_farkas_candidate_free(&r);
    require(lab_farkas_projected_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 0,
            "row projection never clips a nonzero combined free-column coefficient");
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
    read_model(&m,
        "NAME SOURCEBOX\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R 1.00000000000000001\nBOUNDS\n"
        " FX b x 1.00000000000000002\nENDATA\n", path);
    pi[0] = 1;
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
            "rounded fixed-looking box cannot fabricate an exact contradiction");
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
    read_model(&m,
        "NAME LOST\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 10000000000000001\n"
        " x R -10000000000000000\nRHS\n b R 1\nBOUNDS\n FX b x 1\nENDATA\n", path);
    require(!m.exact && m.value[0] == 0, "uncertain-zero source fixture");
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
            "uncertain source coefficient is not treated as a literal zero");
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
}
static void overlays_and_unknown_sources(const char *checker) {
    LabModel m = {0};
    char path[64];
    read_model(&m,
        "NAME OVERLAY\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R 1\nBOUNDS\n UP b x 2\nENDATA\n", path);
    size_t source[] = {0};
    double scale[] = {1}, pi[] = {1};
    LabFarkasCandidate r;
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
            "original box is feasible");
    lab_farkas_candidate_free(&r);
    m.col_upper[0] = 0;
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 1,
            "changed finite node-box overlay can justify a local contradiction");
    LabModel box = {0};
    char box_path[64];
    read_model(&box,
        "NAME OVERLAY\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R 1\nBOUNDS\n UP b x 0\nENDATA\n", box_path);
    checker_replay(checker, box_path, &box, &r, 0);
    unlink(box_path); lab_model_free(&box); lab_farkas_candidate_free(&r);
    LabExactModel *exact = m.exact;
    LabSourceEnvelope *envelope = m.envelope;
    m.exact = NULL; m.envelope = NULL;
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
            "parsed model without trustworthy source data abstains");
    lab_farkas_candidate_free(&r);
    m.exact = exact; m.envelope = envelope;
    unlink(path); lab_model_free(&m);
    read_model(&m,
        "NAME WIDENUP\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R 2\nBOUNDS\n UP b x 1\nENDATA\n", path);
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 1,
            "original narrow upper side makes this source model infeasible");
    lab_farkas_candidate_free(&r);
    m.col_upper[0] = 3;
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
            "a widened finite upper side replaces the original and admits a feasible point");
    lab_farkas_candidate_free(&r);
    m.col_upper[0] = INFINITY;
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
            "a widened infinite upper side cannot silently retain the original bound");
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
    read_model(&m,
        "NAME WIDENLOW\nROWS\n N O\n L R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R -2\nBOUNDS\n LO b x -1\nENDATA\n", path);
    pi[0] = -1;
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 1,
            "original narrow lower side makes this source model infeasible");
    lab_farkas_candidate_free(&r);
    m.col_lower[0] = -3;
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
            "a widened finite lower side replaces the source lower bound");
    lab_farkas_candidate_free(&r);
    m.col_lower[0] = -INFINITY;
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
            "a widened infinite lower side replaces the source lower bound");
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
}
static void feasible_oracles_and_malformed_input(void) {
    size_t start[] = {0, 4, 8, 12};
    size_t rows[] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    double a[12], lower[] = {-2, -2, -2}, upper[] = {3, 3, 3};
    double rl[4], ru[4], point[] = {1, -1, 2};
    LabModel m = {.rows = 4, .cols = 3, .nnz = 12, .col_start = start,
                  .row_index = rows, .value = a, .row_lower = rl, .row_upper = ru,
                  .col_lower = lower, .col_upper = upper};
    size_t source[] = {0, 1, 2, 3, 0, 1, 2, SIZE_MAX};
    double scale[8], pi[8];
    uint32_t state = UINT32_C(7176723);
    for (size_t sample = 0; sample < 1000; sample++) {
        for (size_t k = 0; k < 12; k++) {
            state = UINT32_C(1664525) * state + UINT32_C(1013904223);
            a[k] = (double)((long)(state % 7) - 3);
        }
        for (size_t i = 0; i < 4; i++) {
            double activity = 0;
            for (size_t j = 0; j < 3; j++)
                activity += a[4 * j + i] * point[j];
            rl[i] = activity - (double)(sample % 3);
            ru[i] = activity + (double)(sample % 5);
        }
        for (size_t i = 0; i < 8; i++) {
            state = UINT32_C(1664525) * state + UINT32_C(1013904223);
            pi[i] = (double)((long)(state % 17) - 8) / 8;
            scale[i] = scalbn(1, (int)(i % 7) - 3);
        }
        LabFarkasCandidate r;
        require(lab_farkas_candidate(&m, 8, source, scale, pi, lab_now() + 5, &r) == 0,
                "1000 independent exact feasible points forbid a contradiction for arbitrary dyadic y");
        require(!r.multiplier && !r.margin_lower_bound, "unavailable result does not leak a witness");
        lab_farkas_candidate_free(&r);
        require(lab_farkas_projected_candidate(&m, 8, source, scale, pi, lab_now() + 5, &r) == 0,
                "1000 exact feasible points also forbid a projected-candidate contradiction");
        lab_farkas_candidate_free(&r);
    }
    source[0] = 4;
    LabFarkasCandidate r;
    require(lab_farkas_candidate(&m, 8, source, scale, pi, lab_now() + 5, &r) == -1,
            "out-of-range source row rejected");
    lab_farkas_candidate_free(&r);
    double dummy = 0;
    LabModel huge = {.rows = 4000001, .row_lower = &dummy, .row_upper = &dummy};
    require(lab_farkas_candidate(&huge, 0, NULL, NULL, NULL, lab_now() + 5, &r) == 0,
            "storage limit abstains before allocating or reading fictitious rows");
    lab_farkas_candidate_free(&r);
}
static void masks_and_one_sided_infinity(const char *checker) {
    LabModel m = {0};
    char path[64];
    read_model(&m,
        "NAME MASK\nROWS\n N O\n G R\n G SKIP\nCOLUMNS\n x O 0 R 1\n"
        " x SKIP 1\nRHS\n b R 2 SKIP -100\nBOUNDS\n UP b x 1\nENDATA\n", path);
    m.row_lower[1] = -INFINITY;
    m.row_upper[1] = INFINITY;
    size_t source[] = {0, 1};
    double scale[] = {1, 1}, pi[] = {1, 0};
    LabFarkasCandidate r;
    require(lab_farkas_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 1,
            "zero multiplier on a masked row preserves original-row indexing");
    require(!strcmp(r.multiplier[1], "0"), "masked row exports a literal zero");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    pi[1] = 1;
    require(lab_farkas_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 0,
            "nonzero multiplier selecting an infinite current row side is unavailable");
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
    read_model(&m,
        "NAME ONESIDED\nROWS\n N O\n L R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R -1\nENDATA\n", path);
    pi[0] = -1;
    require(lab_farkas_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 1,
            "negative q on a lower-bounded unbounded column has finite upper support");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
}
static void projected_wrong_sign_control(const char *checker) {
    LabModel m = {0};
    char path[64];
    read_model(&m,
        "NAME PROJECT\nROWS\n N O\n G A\n G B\nCOLUMNS\n x O 0 A 1\n"
        " x B 1\nRHS\n b A 2 B -100\nBOUNDS\n UP b x 1\nENDATA\n", path);
    size_t source[] = {0, 1};
    double scale[] = {1, 1}, pi[] = {1, -0x1p-60};
    LabFarkasCandidate r;
    require(lab_farkas_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 0,
            "strict constructor still rejects a wrong-sign one-sided multiplier");
    lab_farkas_candidate_free(&r);
    require(lab_farkas_projected_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 1,
            "projected constructor proposes a different fully validated row combination");
    require(r.projected_rows == 1 && !strcmp(r.multiplier[1], "0") &&
            !strcmp(r.margin_lower_bound, "1"), "projection is exact and recorded");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
}
int main(int argc, char **argv) {
    const char *checker = argc > 1 ? argv[1] : NULL;
    basic_and_dyadic_mapping(checker);
    envelope_lower_margin_and_feasible_adversary(checker);
    free_columns_and_source_box(checker);
    overlays_and_unknown_sources(checker);
    feasible_oracles_and_malformed_input();
    masks_and_one_sided_infinity(checker);
    projected_wrong_sign_control(checker);
    puts("Source-envelope Farkas candidates, dyadic mapping and independent exact-margin replay pass.");
    return 0;
}

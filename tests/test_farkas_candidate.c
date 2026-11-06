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

typedef int (*Constructor)(const LabModel *, size_t, const size_t *, const double *,
                           const double *, double, LabFarkasCandidate *);
static Constructor constructors[] = {
    lab_farkas_candidate, lab_farkas_projected_candidate, lab_farkas_simplified_candidate,
    lab_farkas_reconstructed_candidate
};

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
    require(lab_farkas_simplified_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 0 &&
            r.removed_terms == 0,
            "simplification never clips a surviving nonzero combined free-column coefficient");
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
    require(m.exact && m.value[0] == 1, "exact coalescing restores source fixture coefficient");
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
    for (size_t sample = 0; sample < 2000; sample++) {
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
            if ((sample + i) % 7 == 0)
                rl[i] = -INFINITY;
            if ((sample + i) % 11 == 0)
                ru[i] = INFINITY;
        }
        for (size_t j = 0; j < 3; j++) {
            lower[j] = (sample + j) % 13 == 0 ? -INFINITY : -2;
            upper[j] = (sample + j) % 17 == 0 ? INFINITY : 3;
        }
        for (size_t i = 0; i < 8; i++) {
            state = UINT32_C(1664525) * state + UINT32_C(1013904223);
            pi[i] = (double)((long)(state % 17) - 8) / (double)(3 + 2 * (state % 7));
            if ((sample + i) % 3 == 0)
                pi[i] = scalbn(pi[i], -100);
            scale[i] = scalbn(1, (int)(i % 7) - 3);
        }
        LabFarkasCandidate r;
        double pi_snapshot[8], scale_snapshot[8];
        memcpy(pi_snapshot, pi, sizeof pi);
        memcpy(scale_snapshot, scale, sizeof scale);
        for (size_t method = 0; method < sizeof constructors / sizeof *constructors; method++) {
            require(constructors[method](&m, 8, source, scale, pi, lab_now() + 5, &r) == 0,
                    "2000 exact feasible points forbid all four constructors from claiming infeasibility");
            require(!r.multiplier && !r.margin_lower_bound,
                    "unavailable result does not leak a witness");
            require(!memcmp(pi_snapshot, pi, sizeof pi) &&
                    !memcmp(scale_snapshot, scale, sizeof scale), "candidate inputs are immutable");
            lab_farkas_candidate_free(&r);
        }
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
    require(lab_farkas_simplified_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 1 &&
            r.removed_terms == 1 && r.projected_rows == 0,
            "simplification records omitted canonical terms before row projection");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    pi[1] = -.5;
    require(lab_farkas_simplified_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 1 &&
            r.removed_terms == 0 && r.projected_rows == 1,
            "simplified constructor still projects retained invalid-side source multipliers");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
}
static void simplified_cutoffs_and_scope(const char *checker) {
    LabModel m = {0};
    char path[64];
    read_model(&m,
        "NAME SIMPLIFY\nROWS\n N O\n G PROOF\n E NOISE\nCOLUMNS\n x O 0 NOISE 1\n"
        "RHS\n b PROOF 1 NOISE 0\nBOUNDS\n FR b x\nENDATA\n", path);
    size_t source[] = {0, 1, SIZE_MAX};
    double scale[] = {1, 1, 1}, pi[] = {1, 0x1p-60, 0};
    LabFarkasCandidate r;
    for (size_t method = 0; method < 2; method++) {
        require(constructors[method](&m, 2, source, scale, pi, lab_now() + 5, &r) == 0,
                "strict and projected APIs preserve the original noisy-candidate rejection");
        require(!r.removed_terms && !r.reconstructed_terms,
                "legacy constructors never simplify or reconstruct canonical multipliers");
        lab_farkas_candidate_free(&r);
    }
    require(lab_farkas_simplified_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 1 &&
            r.removed_terms == 1 && !r.projected_rows && !strcmp(r.multiplier[0], "1") &&
            !strcmp(r.multiplier[1], "0") && !strcmp(r.margin_lower_bound, "1"),
            "generic simplification proposes an independently verifiable different candidate");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    for (int exponent = -1022; exponent <= 1023; exponent += 2045) {
        pi[0] = scalbn(1, exponent);
        pi[1] = scalbn(1, exponent - 40);
        require(lab_farkas_simplified_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 1 &&
                r.removed_terms == 1, "exact inclusive relative cutoff works at extreme exponents");
        checker_replay(checker, path, &m, &r, 0);
        lab_farkas_candidate_free(&r);
        pi[1] = nextafter(pi[1], INFINITY);
        require(lab_farkas_simplified_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 0 &&
                !r.removed_terms, "a surviving one-ulp-above-cutoff free-column residual still rejects");
        lab_farkas_candidate_free(&r);
    }
    pi[0] = 1; pi[1] = 0x1p-60; pi[2] = 0x1p1023;
    require(lab_farkas_simplified_candidate(&m, 3, source, scale, pi, lab_now() + 5, &r) == 0 &&
            r.removed_terms == 2,
            "cutoff maximum and removal count include ignored canonical source rows consistently");
    lab_farkas_candidate_free(&r);
    pi[0] = pi[1] = pi[2] = 0;
    for (size_t method = 0; method < sizeof constructors / sizeof *constructors; method++) {
        require(constructors[method](&m, 3, source, scale, pi, lab_now() + 5, &r) == 0 &&
                !r.removed_terms, "zero maximum emits no certificate and removes no nonzero terms");
        lab_farkas_candidate_free(&r);
        require(constructors[method](&m, 0, NULL, NULL, NULL, lab_now() + 5, &r) == 0,
                "empty canonical multiplier vectors do not prove infeasibility");
        lab_farkas_candidate_free(&r);
        require(constructors[method](&m, 3, source, scale, pi, lab_now() - 1, &r) == 0 &&
                !r.multiplier && strstr(r.message, "deadline"), "all constructors obey expired budgets");
        lab_farkas_candidate_free(&r);
        pi[0] = NAN;
        require(constructors[method](&m, 3, source, scale, pi, lab_now() + 5, &r) == -1,
                "all constructors reject nonfinite pi before exact conversion or simplification");
        lab_farkas_candidate_free(&r);
        pi[0] = 0;
    }
    for (int exponent = -1074; exponent <= 1023; exponent += 2097) {
        scale[0] = pi[0] = scalbn(1, exponent);
        for (size_t method = 0; method < sizeof constructors / sizeof *constructors; method++) {
            int reconstructed_to_zero = method == 3 && exponent < 0;
            require(constructors[method](&m, 1, source, scale, pi, lab_now() + 5, &r) ==
                    (reconstructed_to_zero ? 0 : 1) && !r.removed_terms,
                    "legacy constructors retain tiny pi; bounded-denominator reconstruction can change it to zero");
            if (reconstructed_to_zero)
                require(r.reconstructed_terms == 1 && !r.multiplier,
                        "changing a subnormal canonical pi to zero is counted and cannot emit a proof");
            else
                checker_replay(checker, path, &m, &r, 0);
            lab_farkas_candidate_free(&r);
        }
    }
    unlink(path); lab_model_free(&m);
}
static void all_constructor_source_and_overlay_guards(const char *checker) {
    size_t source[] = {0};
    double scale[] = {1}, pi[] = {1};
    LabFarkasCandidate r;
    LabModel m = {0};
    char path[64], *text = decimal_large(0);
    read_model(&m, text, path); free(text);
    for (size_t method = 0; method < sizeof constructors / sizeof *constructors; method++) {
        require(constructors[method](&m, 1, source, scale, pi, lab_now() + 5, &r) == 1,
                "all constructors prove against original decimal coefficient enclosures");
        checker_replay(checker, path, &m, &r, 1);
        lab_farkas_candidate_free(&r);
    }
    unlink(path); lab_model_free(&m);
    text = decimal_large(1);
    read_model(&m, text, path); free(text);
    for (size_t method = 0; method < sizeof constructors / sizeof *constructors; method++) {
        require(constructors[method](&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
                "all constructors preserve exact feasibility at an uncertain decimal boundary");
        lab_farkas_candidate_free(&r);
    }
    unlink(path); lab_model_free(&m);
    read_model(&m,
        "NAME REPLACE\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 1\n"
        "RHS\n b R 2\nBOUNDS\n UP b x 1\nENDATA\n", path);
    for (size_t method = 0; method < sizeof constructors / sizeof *constructors; method++) {
        m.col_upper[0] = 3;
        require(constructors[method](&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
                "widened finite overlays replace the original bound for all constructors");
        lab_farkas_candidate_free(&r);
        m.col_upper[0] = INFINITY;
        require(constructors[method](&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
                "widened infinite overlays replace the original bound for all constructors");
        lab_farkas_candidate_free(&r);
        m.col_upper[0] = 1;
        LabExactModel *exact = m.exact;
        LabSourceEnvelope *envelope = m.envelope;
        m.exact = NULL; m.envelope = NULL;
        require(constructors[method](&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
                "unknown parsed provenance abstains in all constructors");
        lab_farkas_candidate_free(&r);
        m.exact = exact; m.envelope = envelope;
        m.col_upper[0] = 0;
        require(constructors[method](&m, 1, source, scale, pi, lab_now() + 5, &r) == 1,
                "changed narrowed overlays can prove only the current node box");
        LabModel box = {0};
        char box_path[64];
        read_model(&box,
            "NAME REPLACE\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 1\n"
            "RHS\n b R 2\nBOUNDS\n UP b x 0\nENDATA\n", box_path);
        checker_replay(checker, box_path, &box, &r, 0);
        unlink(box_path); lab_model_free(&box); lab_farkas_candidate_free(&r);
    }
    unlink(path); lab_model_free(&m);
}
static void active_candidate_budget(void) {
    double row_lower[] = {1}, row_upper[] = {INFINITY};
    LabModel m = {.rows = 1, .row_lower = row_lower, .row_upper = row_upper};
    const size_t count = 200000;
    size_t *source = calloc(count, sizeof *source);
    double *scale = calloc(count, sizeof *scale), *pi = calloc(count, sizeof *pi);
    require(source && scale && pi, "budget fixture arrays");
    for (size_t i = 0; i < count; i++) {
        scale[i] = 1;
        pi[i] = 1.0 / 3;
    }
    for (size_t method = 0; method < sizeof constructors / sizeof *constructors; method++) {
        LabFarkasCandidate r;
        require(constructors[method](&m, count, source, scale, pi, lab_now() + .001, &r) == 0 &&
                !r.multiplier && strstr(r.message, "deadline"),
                "large otherwise-positive candidates stop on an active LP proof budget");
        lab_farkas_candidate_free(&r);
    }
    free(source); free(scale); free(pi);
}
static void reconstruction_probe(double pi, const char *expected) {
    double row_lower[] = {pi < 0 ? -INFINITY : 1};
    double row_upper[] = {pi < 0 ? -1 : INFINITY};
    LabModel m = {.rows = 1, .row_lower = row_lower, .row_upper = row_upper};
    size_t source[] = {0};
    double scale[] = {1};
    mpq_t want, actual, original;
    mpq_inits(want, actual, original, NULL);
    require(!mpq_set_str(want, expected, 10), "independent bounded-denominator oracle syntax");
    mpq_set_d(original, pi);
    LabFarkasCandidate r;
    int result = lab_farkas_reconstructed_candidate(&m, 1, source, scale, &pi, lab_now() + 5, &r);
    require(r.reconstructed_terms == (size_t)(mpq_cmp(want, original) != 0) &&
            !r.removed_terms && !r.projected_rows, "reconstruction changes, including zero, are counted exactly");
    if (!mpq_sgn(want))
        require(result == 0 && !r.multiplier, "nearest zero cannot certify the contradictory row by itself");
    else {
        require(result == 1 && !mpq_set_str(actual, r.multiplier[0], 10) &&
                !mpq_cmp(actual, want), "C nearest rational matches the independently specified oracle");
        require(mpz_cmp_ui(mpq_denref(actual), 1048576UL) <= 0,
                "reconstructed denominator respects the fixed bound");
    }
    lab_farkas_candidate_free(&r);
    mpq_clears(want, actual, original, NULL);
}
static void reconstruction_curated_oracle(void) {
    /* Expected values generated independently with Fraction.limit_denominator(2^20). */
    static const struct { double pi; const char *expected; } cases[] = {
    { 0x0.0p+0, "0" },
    { -0x0.0p+0, "0" },
    { 0x1.0000000000000p+0, "1" },
    { -0x1.0000000000000p+0, "-1" },
    { 0x1.5555555555555p-2, "1/3" },
    { -0x1.5555555555555p-2, "-1/3" },
    { 0x1.999999999999ap-4, "1/10" },
    { -0x1.999999999999ap-4, "-1/10" },
    { 0x1.921fb54442d18p+1, "3126535/995207" },
    { -0x1.921fb54442d18p+1, "-3126535/995207" },
    { 0x1.6a09e667f3bcdp+0, "665857/470832" },
    { -0x1.6a09e667f3bcdp+0, "-665857/470832" },
    { 0x1.921fb78121fb8p+1, "355/113" },
    { 0x1.9000132308d74p+9, "532867589/666084" },
    { 0x1.fffffffffffffp+1023, "179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368" },
    { -0x1.fffffffffffffp+1023, "-179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368" },
    { 0x1.0000000000000p-1022, "0" },
    { -0x1.0000000000000p-1022, "0" },
    { 0x0.0000000000001p-1022, "0" },
    { -0x0.0000000000001p-1022, "0" },
    { -0x1.73dbcfda8e78dp-3, "-131072/721875" },
    { 0x1.54ded3ddad442p-7, "2048/196875" },
    { 0x1.fffffffffffffp-22, "0" },
    { -0x1.fffffffffffffp-22, "0" },
    { 0x1.0000000000000p-21, "0" },
    { -0x1.0000000000000p-21, "0" },
    { 0x1.0000000000001p-21, "1/1048576" },
    { -0x1.0000000000001p-21, "-1/1048576" },
    { 0x1.000007fffffffp+0, "1" },
    { -0x1.000007fffffffp+0, "-1" },
    { 0x1.0000080000000p+0, "1" },
    { -0x1.0000080000000p+0, "-1" },
    { 0x1.0000080000001p+0, "1048577/1048576" },
    { -0x1.0000080000001p+0, "-1048577/1048576" },
    { 0x1.fffffffffffffp-21, "1/1048576" },
    { -0x1.fffffffffffffp-21, "-1/1048576" },
    { 0x1.0000000000000p-20, "1/1048576" },
    { -0x1.0000000000000p-20, "-1/1048576" },
    { 0x1.0000000000001p-20, "1/1048576" },
    { -0x1.0000000000001p-20, "-1/1048576" },
    };
    for (size_t i = 0; i < sizeof cases / sizeof *cases; i++)
        reconstruction_probe(cases[i].pi, cases[i].expected);
}
static void reconstruction_external_oracle(const char *path) {
    if (!path)
        return;
    FILE *f = fopen(path, "r");
    require(f != NULL, "open independent reconstruction oracle");
    double input;
    char expected[4096];
    size_t count = 0;
    while (fscanf(f, "%la %4095s", &input, expected) == 2) {
        require(isfinite(input), "oracle pi must be finite");
        reconstruction_probe(input, expected);
        count++;
    }
    require(feof(f) && count > 0, "complete nonempty reconstruction oracle");
    fclose(f);
    printf("%zu independent bounded-denominator oracle cases pass.\n", count);
}
static void reconstruction_proof_and_residual_guards(const char *checker) {
    LabModel m = {0};
    char path[64];
    read_model(&m,
        "NAME RECONSTRUCT\nROWS\n N O\n E A\n E B\nCOLUMNS\n x O 0 A 3\n"
        " x B 10\nRHS\n b A 6 B 10\nBOUNDS\n FR b x\nENDATA\n", path);
    size_t source[] = {0, 1, SIZE_MAX};
    double scale[] = {1, 1, 1}, pi[] = {1.0 / 3, -.1, 0};
    LabFarkasCandidate r;
    for (size_t method = 0; method < 3; method++) {
        require(constructors[method](&m, 2, source, scale, pi, lab_now() + 5, &r) == 0 &&
                !r.reconstructed_terms, "all legacy constructors preserve the meaningful-residual rejection");
        lab_farkas_candidate_free(&r);
    }
    require(lab_farkas_reconstructed_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 1 &&
            r.reconstructed_terms == 2 && !strcmp(r.multiplier[0], "1/3") &&
            !strcmp(r.multiplier[1], "-1/10") && !strcmp(r.margin_lower_bound, "1"),
            "reconstruction proposes different meaningful multipliers and validates the whole proof");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    pi[0] = .25; pi[1] = -.125;
    require(lab_farkas_reconstructed_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 0 &&
            !r.reconstructed_terms, "a surviving exactly represented free-column residual is never clipped");
    lab_farkas_candidate_free(&r);
    pi[0] = 1.0 / 3; pi[1] = -1.0 / 11;
    require(lab_farkas_reconstructed_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 0 &&
            r.reconstructed_terms == 2,
            "a nonzero free-column residual after meaningful rational reconstruction also rejects");
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
    read_model(&m,
        "NAME POSITIVEQ\nROWS\n N O\n G PROOF\n E NOISE\nCOLUMNS\n x O 0 NOISE 1\n"
        "RHS\n b PROOF 1 NOISE 0\nENDATA\n", path);
    pi[0] = 1; pi[1] = 0x1p-20;
    require(lab_farkas_reconstructed_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 0 &&
            !r.reconstructed_terms, "positive surviving q on a one-sided unbounded box remains unavailable");
    lab_farkas_candidate_free(&r);
    pi[1] = 0x1p-60; pi[2] = 1.0 / 3;
    require(lab_farkas_reconstructed_candidate(&m, 3, source, scale, pi, lab_now() + 5, &r) == 1 &&
            r.reconstructed_terms == 2 && !strcmp(r.multiplier[1], "0"),
            "zero changes and ignored canonical terms participate in reconstruction accounting");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    scale[0] = 0x1p-1074; pi[0] = 0x1p-20;
    require(lab_farkas_reconstructed_candidate(&m, 1, source, scale, pi, lab_now() + 5, &r) == 1 &&
            !r.reconstructed_terms, "exact scale multiplication preserves reconstructed dyadic products below binary64");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    unlink(path); lab_model_free(&m);
    read_model(&m,
        "NAME PROJECTAFTER\nROWS\n N O\n G PROOF\n G NOISE\nCOLUMNS\n x O 0 PROOF 1\n"
        " x NOISE 1\nRHS\n b PROOF 2 NOISE -100\nBOUNDS\n UP b x 1\nENDATA\n", path);
    scale[0] = 1; pi[0] = 1; pi[1] = -1.0 / 3;
    require(lab_farkas_reconstructed_candidate(&m, 2, source, scale, pi, lab_now() + 5, &r) == 1 &&
            r.reconstructed_terms == 1 && r.projected_rows == 1,
            "invalid-side row projection follows exact reconstructed aggregation");
    checker_replay(checker, path, &m, &r, 0);
    lab_farkas_candidate_free(&r);
    mpz_set_ui(mpq_numref(m.exact->row_lower[0]), 1);
    mpz_mul_2exp(mpq_numref(m.exact->row_lower[0]), mpq_numref(m.exact->row_lower[0]), 8192);
    pi[1] = 0;
    for (size_t method = 0; method < sizeof constructors / sizeof *constructors; method++) {
        require(constructors[method](&m, 1, source, scale, pi, lab_now() + 5, &r) == 0,
                "all constructors retain the 8192-bit selected-source-side guard");
        lab_farkas_candidate_free(&r);
    }
    unlink(path); lab_model_free(&m);
    double dummy = 0;
    size_t dummy_index = 0;
    LabModel huge = {.rows = 4000001, .row_lower = &dummy, .row_upper = &dummy};
    LabModel wide = {.cols = 2000001, .col_start = &dummy_index,
                     .col_lower = &dummy, .col_upper = &dummy};
    LabModel dense = {.rows = 1, .cols = 1, .nnz = 32000001, .col_start = &dummy_index,
                      .row_index = &dummy_index, .value = &dummy, .row_lower = &dummy,
                      .row_upper = &dummy, .col_lower = &dummy, .col_upper = &dummy};
    LabModel empty = {0};
    for (size_t method = 0; method < sizeof constructors / sizeof *constructors; method++) {
        require(constructors[method](&huge, 0, NULL, NULL, NULL, lab_now() + 5, &r) == 0,
                "all constructors enforce storage limits before reading unavailable rows");
        lab_farkas_candidate_free(&r);
        require(constructors[method](&wide, 0, NULL, NULL, NULL, lab_now() + 5, &r) == 0,
                "all constructors enforce the column cap before reading fictitious columns");
        lab_farkas_candidate_free(&r);
        require(constructors[method](&dense, 0, NULL, NULL, NULL, lab_now() + 5, &r) == 0,
                "all constructors enforce the entry cap before reading fictitious entries");
        lab_farkas_candidate_free(&r);
        require(constructors[method](&empty, 8000001, source, scale, pi, lab_now() + 5, &r) == 0,
                "all constructors enforce the canonical term cap before reading fictitious terms");
        lab_farkas_candidate_free(&r);
    }
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
    simplified_cutoffs_and_scope(checker);
    all_constructor_source_and_overlay_guards(checker);
    active_candidate_budget();
    reconstruction_curated_oracle();
    reconstruction_proof_and_residual_guards(checker);
    reconstruction_external_oracle(argc > 2 ? argv[2] : NULL);
    puts("Source-envelope Farkas candidates, dyadic mapping and independent exact-margin replay pass.");
    return 0;
}

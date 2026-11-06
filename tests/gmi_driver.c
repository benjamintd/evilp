#define _POSIX_C_SOURCE 200809L
#include "gmi.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now(void *context) {
    (void)context;
    struct timespec t;
    if (clock_gettime(CLOCK_MONOTONIC, &t))
        return 1e100;
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}

static int read_q(mpq_t q) {
    char token[65536];
    if (scanf("%65535s", token) != 1 || mpq_set_str(q, token, 10) ||
        mpz_sgn(mpq_denref(q)) == 0)
        return -1;
    mpq_canonicalize(q);
    return 0;
}

static int read_bound(int *present, mpq_t q) {
    char token[65536];
    if (scanf("%65535s", token) != 1)
        return -1;
    *present = strcmp(token, "*") != 0;
    if (!*present)
        return 0;
    if (mpq_set_str(q, token, 10) || mpz_sgn(mpq_denref(q)) == 0)
        return -1;
    mpq_canonicalize(q);
    return 0;
}

int main(void) {
    char mode[16];
    size_t n;
    if (scanf("%15s %zu", mode, &n) != 2 || n > 4096)
        return 2;
    int build = !strcmp(mode, "gmi"), round_box = !strcmp(mode, "round");
    if (!build && !round_box)
        return 2;
    mpq_t rhs;
    mpq_init(rhs);
    mpq_t *a = n ? calloc(n, sizeof *a) : NULL;
    unsigned char *integer = n ? calloc(n, sizeof *integer) : NULL;
    LabGmiBox *box = n ? calloc(n, sizeof *box) : NULL;
    if (n && (!a || !integer || !box))
        abort();
    for (size_t j = 0; j < n; ++j) {
        mpq_init(a[j]);
        lab_gmi_box_init(&box[j]);
    }
    int parse_error = read_q(rhs);
    for (size_t j = 0; j < n && !parse_error; ++j) {
        parse_error = read_q(a[j]);
        if (build) {
            unsigned int flag;
            if (scanf("%u", &flag) != 1 || flag > 1)
                parse_error = -1;
            else
                integer[j] = (unsigned char)flag;
        } else {
            parse_error |= read_bound(&box[j].has_lower, box[j].lower);
            parse_error |= read_bound(&box[j].has_upper, box[j].upper);
        }
    }
    LabGmiBudget budget = {4096, 16384, 1U << 24, now(NULL) + 10, now, NULL};
    LabGmiCut exact;
    LabGmiBinaryCut binary;
    lab_gmi_cut_init(&exact);
    lab_gmi_binary_init(&binary);
    LabGmiStatus status = LAB_GMI_INVALID;
    if (!parse_error)
        status = build ? lab_gmi_build(&exact, n, (const mpq_t *)a, integer, rhs, &budget)
                       : lab_gmi_round_box(&binary, n, (const mpq_t *)a, rhs, box, &budget);
    printf("{\"status\":\"%s\"", lab_gmi_status_name(status));
    if (status == LAB_GMI_OK) {
        if (build) {
            gmp_printf(",\"rhs\":\"%Qd\",\"coefficient\":[", exact.rhs);
            for (size_t j = 0; j < n; ++j)
                gmp_printf("%s\"%Qd\"", j ? "," : "", exact.coefficient[j]);
            printf("]");
        } else {
            mpq_t value;
            mpq_init(value);
            mpq_set_d(value, binary.rhs);
            gmp_printf(",\"rhs_binary\":\"%Qd\",\"error_lower\":\"%Qd\","
                       "\"adjusted_rhs\":\"%Qd\",\"coefficient_binary\":[",
                       value, binary.error_lower, binary.adjusted_rhs);
            for (size_t j = 0; j < n; ++j) {
                mpq_set_d(value, binary.coefficient[j]);
                gmp_printf("%s\"%Qd\"", j ? "," : "", value);
            }
            printf("]");
            mpq_clear(value);
        }
    }
    puts("}");
    lab_gmi_cut_clear(&exact);
    lab_gmi_binary_clear(&binary);
    for (size_t j = 0; j < n; ++j) {
        mpq_clear(a[j]);
        lab_gmi_box_clear(&box[j]);
    }
    free(a);
    free(integer);
    free(box);
    mpq_clear(rhs);
    return parse_error ? 2 : 0;
}

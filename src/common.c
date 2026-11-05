#define _POSIX_C_SOURCE 200809L
#include "common.h"
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
void *lab_alloc(size_t n, size_t size) {
    if (size && n > SIZE_MAX / size) {
        fputs("allocation overflow\n", stderr);
        exit(70);
    }
    void *p = calloc(n ? n : 1, size ? size : 1);
    if (!p) {
        fputs("out of memory\n", stderr);
        exit(70);
    }
    return p;
}
char *lab_strdup(const char *s) {
    size_t n = strlen(s) + 1;
    char *p = lab_alloc(n, 1);
    memcpy(p, s, n);
    return p;
}
void lab_json_string(FILE *f, const char *s) {
    fputc('"', f);
    for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
        if (*p == '"' || *p == '\\') {
            fputc('\\', f);
            fputc(*p, f);
        } else if (*p < 32)
            fprintf(f, "\\u%04x", (unsigned)*p);
        else
            fputc(*p, f);
    }
    fputc('"', f);
}
double lab_now(void) {
    struct timespec t;
    if (clock_gettime(CLOCK_MONOTONIC, &t))
        abort();
    return (double)t.tv_sec + (double)t.tv_nsec / 1e9;
}
/* SHA-256, FIPS 180-4. Unsigned 32-bit arithmetic is intentionally modular. */
static uint32_t ror(uint32_t a, unsigned n) {
    return (a >> n) | (a << (32 - n));
}
static void block(uint32_t h[8], const unsigned char p[64]) {
    static const uint32_t k[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2};
    uint32_t w[64];
    for (unsigned i = 0; i < 16; i++)
        w[i] = ((uint32_t)p[4 * i] << 24) | ((uint32_t)p[4 * i + 1] << 16) |
               ((uint32_t)p[4 * i + 2] << 8) | p[4 * i + 3];
    for (unsigned i = 16; i < 64; i++) {
        uint32_t x = w[i - 15], y = w[i - 2];
        w[i] = w[i - 16] + (ror(x, 7) ^ ror(x, 18) ^ (x >> 3)) + w[i - 7] +
               (ror(y, 17) ^ ror(y, 19) ^ (y >> 10));
    }
    uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], z = h[7];
    for (unsigned i = 0; i < 64; i++) {
        uint32_t t = z + (ror(e, 6) ^ ror(e, 11) ^ ror(e, 25)) + ((e & f) ^ (~e & g)) + k[i] + w[i];
        uint32_t u = (ror(a, 2) ^ ror(a, 13) ^ ror(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
        z = g;
        g = f;
        f = e;
        e = d + t;
        d = c;
        c = b;
        b = a;
        a = t + u;
    }
    h[0] += a;
    h[1] += b;
    h[2] += c;
    h[3] += d;
    h[4] += e;
    h[5] += f;
    h[6] += g;
    h[7] += z;
}
int lab_sha256(const char *path, char hex[65]) {
    FILE *f = fopen(path, "rb");
    if (!f)
        return -1;
    uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                     0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    unsigned char p[64];
    uint64_t bytes = 0;
    size_t n;
    while ((n = fread(p, 1, 64, f)) == 64) {
        block(h, p);
        bytes += 64;
    }
    if (ferror(f)) {
        fclose(f);
        return -1;
    }
    fclose(f);
    bytes += (uint64_t)n;
    p[n++] = 0x80;
    if (n > 56) {
        memset(p + n, 0, 64 - n);
        block(h, p);
        n = 0;
    }
    memset(p + n, 0, 56 - n);
    uint64_t bits = bytes * 8;
    for (unsigned i = 0; i < 8; i++)
        p[63 - i] = (unsigned char)(bits >> (8 * i));
    block(h, p);
    for (unsigned i = 0; i < 8; i++)
        snprintf(hex + 8 * i, 9, "%08x", h[i]);
    hex[64] = 0;
    return 0;
}

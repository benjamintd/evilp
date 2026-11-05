#ifndef LAB_COMMON_H
#define LAB_COMMON_H
#include <stddef.h>
#include <stdio.h>
void *lab_alloc(size_t n, size_t size);
char *lab_strdup(const char *s);
void lab_json_string(FILE *f, const char *s);
double lab_now(void);
int lab_sha256(const char *path, char hex[65]);
#endif

#include "dependency.h"
#include "common.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
int main(int argc, char **argv) {
    if (argc < 3)
        return 2;
    LabModel m = {0};
    if (lab_model_read(&m, argv[1]) || (size_t)(argc - 3) != m.rows)
        return 2;
    double *y = lab_alloc(m.rows, sizeof(double));
    for (size_t i = 0; i < m.rows; i++)
        y[i] = strtod(argv[i + 3], NULL);
    LabModel view = m;
    if (!strcmp(argv[2], "unknown")) {
        view.exact = view.native_exact = NULL;
        view.envelope = NULL;
    }
    double deadline = lab_now() + (!strcmp(argv[2], "expired") ? -1 : 5);
    printf("%d\n", lab_source_row_dependency(&view, y, deadline));
    free(y);
    lab_model_free(&m);
    return 0;
}

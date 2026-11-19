#include "common.h"
#include "packing.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    LabModel m = {0};
    if (argc != 3 || lab_model_read(&m, argv[1]))
        return 2;
    size_t steps = 0, limit = (size_t)strtoul(argv[2], NULL, 10);
    double *point = lab_alloc(m.cols, sizeof *point);
    int solved = lab_resource_packing(&m, lab_now()+.1, limit, point, &steps);
    printf("{\"solved\":%d,\"attempts\":%zu,\"columns\":[", solved, steps);
    int printed = 0;
    if (solved)
        for (size_t j = 0; j < m.cols; ++j)
            if (point[j]) {
                printf("%s%zu", printed ? "," : "", j);
                printed = 1;
            }
    puts("]}");
    free(point);
    lab_model_free(&m);
    return 0;
}

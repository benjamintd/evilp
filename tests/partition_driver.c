#include "partition.h"
#include "common.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc != 4)
        return 2;
    LabModel model = {0};
    if (lab_model_read(&model, argv[1]))
        return 2;
    size_t limit = (size_t)strtoull(argv[2], NULL, 10);
    int mode = atoi(argv[3]);
    double *price = lab_alloc(model.rows, sizeof *price);
    double *point = lab_alloc(model.cols, sizeof *point);
    for (size_t i = 0; i < model.rows; i++) {
        if (mode == 1)
            price[i] = (double)((int)(i % 7) - 3) * 7.125;
        else if (mode == 2)
            price[i] = i % 2 ? -0x1p30 : 0x1p30;
        else if (mode == 3)
            price[i] = i % 2 ? NAN : INFINITY;
    }
    LabPartitionResult result;
    lab_partition_search(&model, price, lab_now() + (mode == 4 ? -1 : 5), limit, point, &result);
    printf("{\"nodes\":%zu,\"point\":%d,\"optimal\":%d,\"columns\":[",
           result.nodes, result.has_point, result.optimal);
    int printed = 0;
    if (result.has_point)
        for (size_t j = 0; j < model.cols; j++)
            if (point[j]) {
                printf("%s%zu", printed ? "," : "", j);
                printed = 1;
            }
    puts("]}");
    free(point);
    free(price);
    lab_model_free(&model);
    return 0;
}

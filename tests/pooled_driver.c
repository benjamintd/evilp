#include "common.h"
#include "pooled_assignment.h"
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
    if (argc != 5)
        return 2;
    LabModel m = {0};
    if (lab_model_read_mip(&m, argv[1]))
        return 2;
    double *x = lab_alloc(m.cols, sizeof *x);
    LabPooledStats stats;
    int r = lab_pooled_assignment(&m, lab_now() + 3, strtoul(argv[2], NULL, 10),
                                  strtoul(argv[3], NULL, 10), x, &stats);
    if (r) {
        FILE *f = fopen(argv[4], "wx");
        if (!f)
            return 2;
        for (size_t j = 0; j < m.cols; j++)
            if (x[j] != 0)
                fprintf(f, "%s %.17g\n", m.col_name[j], x[j]);
        if (fclose(f))
            return 2;
    }
    printf("{\"result\":%d,\"nodes\":%zu,\"lp_solves\":%zu,\"iterations\":%zu}\n", r, stats.nodes,
           stats.lp_solves, stats.iterations);
    free(x);
    lab_model_free(&m);
    return 0;
}

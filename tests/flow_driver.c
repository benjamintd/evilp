#include "common.h"
#include "flow.h"
#include "solution.h"
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc != 4)
        return 2;
    LabModel m = {0};
    if (lab_model_read(&m, argv[1]))
        return 2;
    size_t count, states, limit = (size_t)strtoul(argv[3], NULL, 10);
    double *point = lab_alloc(m.cols, sizeof *point);
    int ok = lab_fixed_charge_flow(&m, lab_now() + 5, limit, point, &count, &states);
    if (ok) {
        FILE *f = fopen(argv[2], "wx");
        if (!f || lab_solution_write_point(f, &m, point) || fclose(f))
            return 2;
    }
    printf("{\"solved\":%d,\"eliminations\":%zu,\"states\":%zu}\n", ok, count, states);
    free(point);
    lab_model_free(&m);
    return 0;
}

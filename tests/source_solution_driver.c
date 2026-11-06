#include "exact.h"
#include "solution.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv) {
    if (argc != 4 && (argc != 5 || strcmp(argv[4], "without-objective-source")))
        return 2;
    LabModel model = {0};
    if (lab_model_read(&model, argv[1]) || model.cols != 1) {
        lab_model_free(&model);
        return 2;
    }
    char *end;
    double point = strtod(argv[3], &end);
    if (*end || !isfinite(point)) {
        lab_model_free(&model);
        return 2;
    }
    if (argc == 5) {
        lab_exact_model_free(model.exact);
        model.exact = NULL;
        lab_source_objective_free(model.objective_exact);
        model.objective_exact = NULL;
        if (model.envelope)
            model.envelope->objective_complete = 0;
    }
    FILE *f = fopen(argv[2], "wx");
    if (!f) {
        lab_model_free(&model);
        return 2;
    }
    int bad = lab_solution_write_point(f, &model, &point);
    bad |= fclose(f) != 0;
    lab_model_free(&model);
    return bad ? 1 : 0;
}

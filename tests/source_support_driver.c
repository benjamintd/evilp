#include "common.h"
#include "exact.h"
#include "model.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct {
    const LabModel *model;
} Work;
#include "../src/source_support.inc"
static int count;
#define CHECK(c)                                                                                   \
    do {                                                                                           \
        count++;                                                                                   \
        if (!(c)) {                                                                                \
            fprintf(stderr, "interval support check%d line%d: %s\n", count, __LINE__, #c);         \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)
static int verdict(Work *w, double x, double y, double deadline) {
    int av = 77;
    int bad = source_support_interval_gap(w, &x, &y, deadline, &av);
    return av == 1 ? (bad ? 1 : 0) : 2;
}
/* Cross the old row, column and nonzero limits at once. The last column
 * controls the gap, so a truncated scan would accept the second point. */
static int large_sparse_support(void) {
    const size_t rows = 5000, cols = 17000, per_column = 8, nnz = cols * per_column;
    LabModel m = {.rows = rows, .cols = cols, .nnz = nnz, .parsed_decimal = 1};
    m.envelope = lab_source_envelope_new(rows, cols, nnz);
    LabSourceEnvelope *e = m.envelope;
    e->objective_complete = 1;
    m.col_start = lab_alloc(cols + 1, sizeof *m.col_start);
    m.row_index = lab_alloc(nnz, sizeof *m.row_index);
    m.row_lower = lab_alloc(rows, sizeof(double));
    m.row_upper = lab_alloc(rows, sizeof(double));
    m.col_lower = lab_alloc(cols, sizeof(double));
    m.col_upper = lab_alloc(cols, sizeof(double));
    m.objective = lab_alloc(cols, sizeof(double));
    double *x = lab_alloc(cols, sizeof(double)), *y = lab_alloc(rows, sizeof(double));
    for (size_t i = 0; i < rows; i++) {
        m.row_lower[i] = -INFINITY;
        m.row_upper[i] = 2;
        e->row_lower[i] = (LabInterval){-INFINITY, -INFINITY};
        e->row_upper[i] = (LabInterval){2, 2};
    }
    for (size_t j = 0; j < cols; j++) {
        m.col_start[j] = j * per_column;
        m.col_upper[j] = e->col_upper_source[j] = 1;
        e->col_upper[j] = (LabInterval){1, 1};
        for (size_t k = j * per_column; k < (j + 1) * per_column; k++) {
            m.row_index[k] = k % rows;
            e->value[k] = (LabInterval){1, 1};
        }
    }
    m.col_start[cols] = nnz;
    m.objective[cols - 1] = 1;
    e->objective[cols - 1] = (LabInterval){1, 1};
    Work w = {&m};
    int available = 0;
    CHECK(!source_support_interval_gap(&w, x, y, lab_now() + 1, &available) && available == 1);
    x[cols - 1] = 1;
    CHECK(source_support_interval_gap(&w, x, y, lab_now() + 1, &available) && available == 1);
    x[cols - 1] = 0;
    e->value[nnz - 1] = (LabInterval){2, 1};
    CHECK(!source_support_interval_gap(&w, x, y, lab_now() + 1, &available) && !available);
    free(y);
    free(x);
    lab_model_free(&m);
    return 0;
}
int main(int argc, char **argv) {
    double rl = -INFINITY, ru = 1, cl = 0, cu = 2, cost = -1;
    size_t start[] = {0, 1}, row[] = {0};
    double value = 1;
    LabModel m = {.rows = 1,
                  .cols = 1,
                  .nnz = 1,
                  .col_start = start,
                  .row_index = row,
                  .value = &value,
                  .row_lower = &rl,
                  .row_upper = &ru,
                  .col_lower = &cl,
                  .col_upper = &cu,
                  .objective = &cost,
                  .parsed_decimal = 1};
    LabSourceEnvelope *e = lab_source_envelope_new(1, 1, 1);
    m.envelope = e;
    e->objective_complete = 1;
    e->value[0] = (LabInterval){1, 1};
    e->row_lower[0] = (LabInterval){-INFINITY, -INFINITY};
    e->row_upper[0] = (LabInterval){1, 1};
    e->col_lower[0] = (LabInterval){0, 0};
    e->col_upper[0] = (LabInterval){2, 2};
    e->objective[0] = (LabInterval){-1, -1};
    e->col_lower_source[0] = 0;
    e->col_upper_source[0] = 2;
    Work w = {&m};
    if (argc == 2 && !strcmp(argv[1], "--stream")) {
        double q[20];
        while (scanf("%lf", q) == 1) {
            for (int i = 1; i < 20; i++)
                if (scanf("%lf", q + i) != 1)
                    return 2;
            e->value[0] = (LabInterval){q[0], q[1]};
            e->objective[0] = (LabInterval){q[2], q[3]};
            value = q[0] / 2 + q[1] / 2;
            cost = q[2] / 2 + q[3] / 2;
            e->row_lower[0] = (LabInterval){q[6], q[7]};
            e->row_upper[0] = (LabInterval){q[8], q[9]};
            e->col_lower[0] = (LabInterval){q[10], q[11]};
            e->col_upper[0] = (LabInterval){q[12], q[13]};
            cl = q[14];
            cu = q[15];
            e->col_lower_source[0] = q[16];
            e->col_upper_source[0] = q[17];
            rl = q[18];
            ru = q[19];
            int av = 0;
            int bad = source_support_interval_gap(&w, q + 4, q + 5, lab_now() + 1, &av);
            printf("%d %d\n", av, bad);
        }
        lab_source_envelope_free(e);
        return 0;
    }
    CHECK(!large_sparse_support());
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 0);
    e->value[0] = (LabInterval){nextafter(1, -INFINITY), nextafter(1, INFINITY)};
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 0);
    e->row_upper[0] = (LabInterval){nextafter(1, -INFINITY), nextafter(1, INFINITY)};
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 0);
    cu = INFINITY;
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    cu = 2;
    cu = 1e15;
    e->col_upper_source[0] = cu;
    e->col_upper[0] = (LabInterval){cu, cu};
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 1);
    cu = 2;
    e->col_upper_source[0] = cu;
    e->col_upper[0] = (LabInterval){cu, cu};
    CHECK(verdict(&w, 1, -1, lab_now() - 1) == 2);
    e->complete = 0;
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    e->complete = 1;
    e->rows = 2;
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    e->rows = 1;
    e->objective_complete = 0;
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    e->objective_complete = 1;
    e->value[0] = (LabInterval){2, 1};
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    e->value[0] = (LabInterval){1, 1};
    e->value[0] = (LabInterval){NAN, 1};
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    e->value[0] = (LabInterval){1, 1};
    CHECK(verdict(&w, NAN, -1, lab_now() + 1) == 2);
    CHECK(verdict(&w, 1, INFINITY, lab_now() + 1) == 2);
    row[0] = 1;
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    row[0] = 0;
    start[1] = 0;
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    start[1] = 1;
    LabInterval *saved = e->objective;
    e->objective = NULL;
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    e->objective = saved;
    m.rows = 4097;
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    m.rows = 1;
    m.cols = 16385;
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    m.cols = 1;
    m.nnz = 100001;
    CHECK(verdict(&w, 1, -1, lab_now() + 1) == 2);
    m.nnz = 1;
    cost = 1;
    e->objective[0] = (LabInterval){1, 1};
    e->row_upper[0] = (LabInterval){1, 1};
    cl = .1;
    e->col_lower_source[0] = .1;
    e->col_lower[0] = (LabInterval){nextafter(.1, -INFINITY), nextafter(.1, INFINITY)};
    CHECK(verdict(&w, .1, 0, lab_now() + 1) == 0);
    cost = 1e15;
    e->objective[0] = (LabInterval){cost, cost};
    CHECK(verdict(&w, .1, 0, lab_now() + 1) == 1);
    cl = .25;
    CHECK(verdict(&w, .25, 0, lab_now() + 1) == 0);
    cost = 0;
    e->objective[0] = (LabInterval){0, 0};
    cl = -INFINITY;
    cu = INFINITY;
    rl = -INFINITY;
    ru = INFINITY;
    CHECK(verdict(&w, 0, 0, lab_now() + 1) == 0);
    CHECK(verdict(&w, 0, -1, lab_now() + 1) == 2);
    CHECK(!fesetround(FE_UPWARD));
    CHECK(verdict(&w, 0, 0, lab_now() + 1) == 2);
    CHECK(!fesetround(FE_TONEAREST));
    lab_source_envelope_free(e);
    printf("{\"passed\":true,\"checks\":%d}\n", count);
    return 0;
}

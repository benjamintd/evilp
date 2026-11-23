#ifndef LAB_POOLED_ASSIGNMENT_H
#define LAB_POOLED_ASSIGNMENT_H
#include "model.h"
#include "search_result.h"
typedef struct {
    size_t nodes, lp_solves, iterations;
} LabPooledStats;
LabSearchResult lab_pooled_assignment(const LabModel *model, double deadline, size_t nodes, size_t iterations,
                          double *point, LabPooledStats *stats);
#endif

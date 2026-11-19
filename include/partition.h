#ifndef LAB_PARTITION_H
#define LAB_PARTITION_H
#include "model.h"
typedef struct {
    size_t nodes;
    int has_point, optimal;
} LabPartitionResult;
/* Choose minimum normalized cost among exact covers of source-unit equality
 * rows by binary columns. Row duals only rank choices. A completed search
 * with a witness can report optimal; exhaustion without a witness abstains. */
void lab_partition_search(const LabModel *model, const double *row_dual, double deadline,
                          size_t node_limit, double *point, LabPartitionResult *result);
#endif

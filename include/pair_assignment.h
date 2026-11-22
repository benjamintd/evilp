#ifndef LAB_PAIR_ASSIGNMENT_H
#define LAB_PAIR_ASSIGNMENT_H
#include "model.h"

/* Search exact-cover assignments whose auxiliary columns are forced products.
 * Returns success only after the relaxation closes with an original witness. */
int lab_pair_assignment(const LabModel *model, double deadline, size_t limit, double *point,
                        size_t *nodes);
#endif

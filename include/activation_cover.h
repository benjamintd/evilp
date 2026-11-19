#ifndef LAB_ACTIVATION_COVER_H
#define LAB_ACTIVATION_COVER_H
#include "model.h"

/* Return a checked optimal assignment if an integral covering relaxation
 * on the activation variables meets a verified integer objective bound. */
int lab_activation_cover(const LabModel *model, double deadline, size_t iteration_limit,
                         double *point, size_t *iterations);
#endif

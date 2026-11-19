#ifndef LAB_LOT_SIZING_H
#define LAB_LOT_SIZING_H
#include "model.h"

/* Solve uncapacitated production on paths with exclusive opening decisions.
 * Return one only after exhaustive conflict search and witness reconstruction. */
int lab_lot_sizing(const LabModel *model, double deadline, size_t limit, double *point,
                   size_t *nodes);
#endif

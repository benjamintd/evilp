#ifndef LAB_FARKAS_CANDIDATE_H
#define LAB_FARKAS_CANDIDATE_H
#include "model.h"
typedef struct {
    size_t rows;
    size_t projected_rows;
    char **multiplier;
    char *margin_lower_bound;
    const char *margin_kind; /* source_enclosure; never an exact-margin claim */
    char message[256];
    double elapsed;
} LabFarkasCandidate;
/* Map each active canonical multiplier to its original row as scale*pi, in
 * exact dyadic arithmetic. SIZE_MAX source rows are ignored. Validate a strict
 * original-source Farkas contradiction using conservative source enclosures and
 * the current column box. A changed finite or infinite side REPLACES its source
 * side, matching LP-view semantics; unchanged identities use the source box.
 * Returns1 candidate,0 unavailable/limited/nonpositive margin,-1 invalid input.
 * The owned result contains rational original-row multipliers and a POSITIVE
 * LOWER BOUND on contradiction margin. An independent checker must recompute
 * its actual margin; a replay should require actual>=lower>0, not equality.
 * The original model and input arrays remain unchanged. */
int lab_farkas_candidate(const LabModel *model, size_t terms, const size_t *source_row,
                          const double *scale, const double *pi, double deadline,
                          LabFarkasCandidate *candidate);
/* Alternate candidate constructor: after exact original-row aggregation,
 * zero only multipliers whose sign selects an infinite CURRENT row side.
 * The entire unchanged source-enclosure/support proof must then pass again.
 * No combined-column residual is projected or rounded. */
int lab_farkas_projected_candidate(const LabModel *model, size_t terms,
                                    const size_t *source_row, const double *scale,
                                    const double *pi, double deadline,
                                    LabFarkasCandidate *candidate);
void lab_farkas_candidate_free(LabFarkasCandidate *candidate);
#endif

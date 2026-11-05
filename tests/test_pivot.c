#include "pivot.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
static void require(int c, const char *message) {
    if (!c) {
        fprintf(stderr, "pivot contract: %s\n", message);
        exit(1);
    }
}
int main(void) {
    size_t basis[] = {0, 1}, row;
    double xb[] = {0, 1}, d[] = {1e-8, 2};
    require(!lab_primal_ratio(2, xb, d, basis, 1, &row) && row == 0,
            "old degenerate rule selects the weak pivot");
    require(!lab_primal_ratio(2, xb, d, basis, 0, &row) && row == 1,
            "Harris selects a strong pivot within existing feasibility budget");
    /* A true value -1e-12 can be accepted and clamped to zero. Dividing by
     * the weak 1e-8 pivot amplifies it to -1e-4. The strong pivot leaves the
     * other basic value above the unchanged -1e-8 acceptance threshold. */
    require(-1e-12 / d[0] < -LAB_PRIMAL_TOLERANCE, "weak pivot amplifies tiny error");
    require(-1e-12 - xb[1] / d[1] * d[0] >= -LAB_PRIMAL_TOLERANCE,
            "strong pivot preserves the original guard");
    xb[0] = 0; xb[1] = 1e8; d[0] = 1e-13; d[1] = 1;
    require(lab_primal_ratio(2, xb, d, basis, 0, &row) == -1,
            "tiny positive direction still limits the step");
    d[0] = d[1] = 0;
    require(lab_primal_ratio(2, xb, d, basis, 0, &row) == 1, "genuine missing bound");
    d[0] = 1e-13;
    require(lab_primal_ratio(2, xb, d, basis, 0, &row) == -1, "tiny is not unbounded");
    d[0] = d[1] = 1; xb[0] = xb[1] = 1; basis[0] = 9; basis[1] = 2;
    require(!lab_primal_ratio(2, xb, d, basis, 0, &row) && row == 1,
            "equal pivot ties retain deterministic Bland order");
    d[0] = NAN;
    require(lab_primal_ratio(2, xb, d, basis, 0, &row) == -1, "invalid direction");
    /* Taking every relaxed strong pivot can create a tolerated negative that
     * the following weak pivot amplifies. Keep Bland as the normal path and
     * validate a Harris rescue only after its ordinary pivot is rejected. */
    xb[0] = 1; xb[1] = ldexp(4, -28);
    d[0] = 1; d[1] = ldexp(5, -28);
    require(!lab_primal_ratio(2, xb, d, basis, 0, &row) && row == 0,
            "unconditional Harris can cross a tiny-direction bound");
    require(!lab_primal_ratio(2, xb, d, basis, 1, &row) && row == 1,
            "normal Bland path avoids the constructed Harris trap");
    puts("Stable primal ratio contracts passed.");
    return 0;
}

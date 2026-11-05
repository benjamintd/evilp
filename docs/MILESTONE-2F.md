# Original-unit primal reconstruction and restoration

Version 0.6 addresses the remaining `mas74` failure and source-arithmetic
regressions found while investigating it. The independent checker is unchanged.
The original `1e-7` primal acceptance policy and existing basis residual guards
remain in force. Integer search and its global bounds are still numerical
claims; this milestone does not add a complete tree proof.

## Diagnosis

The previous search stopped at its 345th processed node, node ID 445 at depth
107. An exact analysis of its 101-row basis found a genuine negative slack
near `-6.99e-9`. That passed the canonical `1e-8` guard, but the row's change
of units was `2^-39`, making the original-row violation about 3840.92.
Original-decimal and binary-basis solutions agreed closely: residual refinement
and same-basis reconstruction could not turn this infeasible basis feasible.

The basis was exactly dual feasible. An independent rational diagnostic needed
two dual pivots to restore feasibility. The separate checker certified that
isolated LP's rational primal/dual optimum, approximately 19227.188241847,
also agreeing with HiGHS. That diagnostic certificate concerns the isolated
node LP, not the original integer model or the production search tree.

The regression investigation also found arithmetic failures. For example,
fixing `z=2^53` in `x + .1 z >= 900719925474099.25` requires `x>=.05` in the
original decimal model. The old optimizer returned `x=0`, objective zero and
residual zero because rounded fixed contributions lost the source gap.
Other generated cases exposed inaccurate feasible objectives, shift/recovery
cancellation and duplicate objective terms that rounded a nonzero cost to zero.

## Changes

For retained small-model provenance, shifted RHS values now use original
decimal coefficients, bounds and fixed contributions. The optimizer also
retains exact objective coefficients and offset, and initializes floating costs
from the exact coalesced objective. Nonzero objective underflow is explicit.

Before publishing a small-model primal, the optimizer reconstructs its basis
over original coefficients and exact cached shifts. It removes singleton rows
and columns, solves a remaining rational core of at most 128 rows, then
back-substitutes and checks every original equation. This avoids allocating a
full dense basis for the many unit slack and bound columns. The bounded sparse
reconstruction caps the full dimension at 1024, stored entries at 100,000 and
rational sizes at 8192 bits, within one second and the existing LP deadline.

Every reconstructed basic value must be nonnegative. The rational original
point must satisfy all original rows, including deleted rows, and the current
column box. Conversion to binary64 happens after shift cancellation, followed
by the unchanged original-unit check on the actual returned double point.
An unrepresentable, limited or infeasible reconstruction cannot manufacture a
primal witness. `--no-exact` disables this reconstruction as well as Phase I
witness reconstruction; original-model evaluation still applies.

Original-model evaluation uses GMP source arithmetic where available. Other
models use compensated FMA accumulation with conservative error bounds. The
existing support-gap rejection check also uses source coefficients and bounds
for retained models, so rounded row activity cannot hide material objective
error. These optimizer checks remain separate from `milpcheck` and its parser.

When cold recovery finds an infeasible basis, dual simplex can restore it
using basic feasibility thresholds mapped to original units. For a slack,
the threshold is no larger than `1e-7` times its equation scale; for a structural
column it accounts for the variable scale. No threshold exceeds the existing
canonical tolerance. Recovery and validation run again after repair. A warm
session that fails recovery keeps the existing cold-rebuild path, avoiding
expensive repairs in a needlessly large cached basis.

`--no-primal-repair` disables the restoration retry for comparison. JSON records
repair calls, repair pivots, exact primal attempts and successful reconstructions.
These diagnostics are not independent certificates. Source provenance remains
bounded as described in [milestone 2E](MILESTONE-2E.md); larger models still use
the numerical path and the separate checker.

## Validation scope

The new Python suite covers 58 LP/MILP runs across nine methods. It compares
points and objectives using independent fractions, checks submitted and
equivalent normalized rows, and calls the independent checker at `1e-7`.
The C session contract checks 176 warm/cold outcomes across both backends,
scaling controls and row units. Owned regressions use generated coefficients,
including a retained-basis negative-slack example.

The frozen previous solver fails 39 Python subcases. A passing new contract
means the assertions pass: some hard cases permit explicit numerical failure.
It does not mean 39 new proven optima. Required successful cases are separately
identified in the retained validation report.

JSON objectives remain binary64 values. Beside large fixed objective constants,
the regressions compare the exact cost of the exported point and check that
point without its rounded objective claim. They do not assert exact objective
serialization at magnitudes where binary64 cannot represent the small change.

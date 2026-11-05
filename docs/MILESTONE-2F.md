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

## Recorded results

The evaluated engine is frozen at `6d5298c`. The selected campaign's source
snapshots are clean and identical before and after all fourteen runs; binary
hashes, commands, raw witnesses and checker records are retained. All campaigns
completed successfully as validation procedures. This does not mean all of
their MILP instances were solved.

- [Release and sanitizer contracts](../results/milestone-2f-contracts-v1/summary.json):
  all eighteen suites pass in both builds, with zero failed or skipped tests.
- [Paired primal regressions](../results/primal-regressions-v1/summary.json):
  the previous frozen engine has 39 failed assertion subcases; the new contract
  has none. Of its 58 Python runs, eighteen require successful optimal outcomes
  and forty permit explicit numerical failure. All 176 C session outcomes pass.
- [Selected comparison](../results/primal-recovery-v1/report.html): all fourteen
  runs are retained and all predeclared assertions pass, including the old false optimum and its checker
  rejection as evidence. These are development cases, not a holdout.
- [Native MILP crosscheck](../results/mip-native-v5/summary.json): all forty
  expected outcomes match across eight models; all 35 feasible candidates pass
  independent checking.

On the full `mas74` model, the default engine reaches the 1,000-node limit and
retains a checked incumbent of 19563.4046944. The previous engine and the new
`--no-primal-repair` control stop numerically at node 345. Default search uses
two restoration calls and four dual repair pivots. The isolated node LP needs
two pivots, passes the original-model checker, and agrees with HiGHS within
`1e-7`. Its [separate dual check](../results/primal-recovery-v1/isolated_node/default/dual-checker.json)
has a valid rational LP lower bound, but does not certify exact optimality.
This isolated LP evidence does not prove the original integer model or tree.

Default `gen-ip054` also reaches 1,000 nodes with a checked incumbent,
8852.207808267. Its 217 recorded LP infeasibility certificates replay correctly.
The old, default and no-repair profiles together retain 651 valid individual
LP certificates, with matching exact contradiction margins; these are not
complete tree proofs. Disabling exact reconstruction still stops at node 35.
Neither full MILP reaches its pinned reference optimum in this campaign.

The serial comparison requests one thread and gives each MILP fifteen seconds
and 1,000 nodes. In this run, `mas74` default solver time is 9.640494 seconds,
versus 7.030875 with exact reconstruction disabled. The old engine takes
3.367941 seconds before failing at 345 nodes, so it performs different work.
Default `gen-ip054` takes 3.213296 seconds versus the old engine's 2.511301 for
the same node limit. Per-LP GMP reconstruction has a measurable cost; these
single development runs do not establish a performance ranking.

The [two-second MIPLIB smoke](../results/miplib-smoke-2s-v4/runs/summary.html)
retains four checked incumbents: `gen-ip054`, `mas76`, `neos5` and `pk1`.
Six cases exceed the configured size cap and six hit the time limit, with no
invalid candidates, checker errors or reached reference optima. `mas74` has
no incumbent at this budget. The short `gen-ip054` incumbent is 15019.62685667,
worse than the longer controlled run, so there is no uniform short-budget gain.

The [autoresearch gate](AUTORESEARCH.md) remains closed. Time, node and iteration
limits, unsupported cases and numerical failures are unresolved problem
outcomes, even when the harness correctly records them and its assertions pass.

# Numerical scaling and basis refinement

Version 0.4 adds changes of units for the LP and a compensated iterative
refinement layer for basis solves. The independent checker and its tolerance
policy are unchanged. The existing [integer-search contract](MILESTONE-2BC.md)
still applies; numerical completion claims are not independent tree proofs.

## Why this change

The earlier MIPLIB smoke stopped on `mas74` and `mas76` during root LP basis
solves. Debugger captures and independent rational analysis found genuine LU
accuracy loss: componentwise backward errors around 1.6e-9 and 1.9e-9 violated
the existing 1e-10 guard. One correction using a high-accuracy residual reduced
those captured-basis errors to roughly 2.2e-17 and 1.8e-17. No tolerance change
or coefficient perturbation was required.

The `gen-ip054` search failure is different. The captured auxiliary LP has a
positive Phase I optimum and exact zero reduced costs; floating error turns
some zero signs slightly negative. Scaling and refinement do not constitute an
infeasibility certificate. The conservative uncertainty guard remains intact.
Exact reconstruction or a checked infeasibility witness is future work.

## Changes of units

The optimizer chooses power-of-two structural-variable scales. Finite column
widths guide bounded-variable units; normalized column magnitudes guide other
variables. Every structural coefficient, including those in bound equations,
and its objective coefficient gets the variable scale. Recovered original
primal values use the same scale.

Rows are then normalized by structural coefficient magnitude, with unit slack
and artificial columns representing slack values in the new row units. Original
row dual multipliers include the inverse mapping back to original row units.
All selected changes of units are clipped toward unity when needed to preserve
finite, reversible stored coefficients and RHS values. Scaling never drops a
nonzero coefficient. The original model stays immutable.

Warm bound updates retain the cached equation scale. If a bound difference
would underflow or become irreversible in those units, the session rebuilds
cold. A regression covers a `DBL_MIN` tightening that previously vanished after
scaling. Row deletion preserves the corresponding scale/identity metadata.

Phase I and artificial-variable cleanup operate in normalized units. Their
thresholds do not replace original-model feasibility checks: recovery still
requires the existing absolute 1e-7 primal guard, and the separate rational
checker still reads the original MPS input under its recorded policy.

## Basis solves and refinement

`src/linear.c` retains an immutable CSC copy of each basis. It equilibrates the
basis with power-of-two row and column factors and uses the existing sparse or
dense LU backend. Normal and transpose solves map back to their input-basis
units before residual checking.

Residuals use fused multiply-add product errors and compensated summation.
The implementation does not assume `long double` has extra precision; on the
current ARM host it does not. Up to eight correction solves reuse the same LU,
retaining the best finite iterate and stopping on convergence or stagnation.
The target is 1e-14 componentwise backward error; acceptance still requires
the original 1e-10 guard. Sparse mode retains storage proportional to nonzeros
and fill, without a hidden dense fallback.

Backward error describes the input CSC basis, which has already undergone the
LP's changes of units. It is distinct from original-MPS primal residual and
from forward solution accuracy on an ill-conditioned system. Both original
primal checking layers remain necessary.

The error-free arithmetic requires binary64, round-to-nearest operations.
Unsafe fast-math, finite-only or excess-precision compilation modes are rejected;
unsupported runtime rounding modes fail explicitly. These restrictions apply
even when corrections are disabled because residual checking remains active.

## Controls and diagnostics

```sh
build/milpsolve MODEL --mip
build/milpsolve MODEL --relax --no-scaling
build/milpsolve MODEL --relax --no-refinement
build/milpsolve MODEL --relax --no-scaling --no-refinement
```

`--no-scaling` disables both LP changes of units and basis equilibration.
`--no-refinement` disables correction solves while retaining compensated
residual measurement and the unchanged acceptance guard. Neither flag bypasses
validation. The combined ablation is not bit-for-bit equivalent to the old
binary because the new residual implementation still applies.

JSON includes enable flags, correction counts, counts of canonically scaled
rows/structural columns and initial/returned basis backward-error diagnostics.
MILP counts aggregate corrections and retain the largest row/column count
over LP calls. Error diagnostics include attempted solves and cold fallbacks;
they are not independent certificates.

## Validation

All eleven CTest suites pass in Release and AddressSanitizer/UndefinedBehaviorSanitizer
builds. The existing 1,000-LP and 300-MILP exact oracles remain green, as do
the 200 bound-overlay sequences and independent certificate tests.

New tests add 160 scaled LPs and 60 scaled MILPs with exact oracles, original
and transformed row checks, free/fixed variables, column-unit spans, large
constants, dense/sparse backends and warm/cold search. The frozen old solver
has 42 explicit numerical failures on the shared numerical cases; the new
engine passes their unchanged assertions. Separate C tests exercise genuinely
improving refinement, normal/transpose/in-place solves, extreme units, rounding
mode rejection, malformed inputs and sparse storage.

The optional `tests/test_live_numerics.py` retains five old/new profiles and
HiGHS root references for `mas74`, `mas76` and `gen-ip054`. These are declared
development models selected from observed failures, not a holdout. All ablation
failures stay in the results. Objective agreement with HiGHS is a crosscheck;
only the independent rational checker establishes the scope of a certificate.

Earlier raw campaigns and evaluator tags remain unchanged. Future performance
work should confirm any gains on additional model families and longer budgets.

## Recorded results

The evaluator is frozen at `28aae81` (`milestone-2d-evaluator`). Both source
snapshots around the controlled root campaign are clean and identical. Binary
hashes, model hashes, commands, raw solutions and checker records accompany the
results; this documentation checkpoint does not change the evaluated engine.

- [Release and sanitizer evidence](../results/milestone-2d-contracts-v1/summary.json):
  all eleven suites pass in both builds, with zero skipped or failed tests.
- [Paired numerical oracle evidence](../results/numerics-oracles-v1/summary.json):
  the old frozen binary has 42 explicit numerical failures across the four
  shared numerical test methods. The new build passes all eight methods,
  including those same assertions. The original checker tolerances are unchanged.
- [Controlled root report](../results/numerics-roots-v1/report.html):
  eighteen runs cover three models, the old binary, four new configurations
  and HiGHS 1.15.1. All six intended numerical failures are retained.
- [Native MILP crosscheck](../results/mip-native-v3/summary.json):
  all forty reference outcomes match across HiGHS, SCIP, Cbc, warm lab search
  and cold lab search. All 35 feasible candidates pass independent checking.

For `mas74` and `mas76`, both the old binary and the new build without
refinement fail the basis-residual guard. Refinement alone and the default
combination finish the root LP, pass original-model primal checking and obtain
valid independent rational LP bounds. Scaling alone still fails. The default
root objectives are approximately 10482.7952803312 and 38893.9036405226,
respectively, agreeing with independently checked HiGHS solutions within 2e-11.
The checker does not certify exact optimality: its conservative lower bounds
leave gaps of approximately 1.5e-5 and 5.4e-5 to the checked candidate objectives.

The default `gen-ip054` root also passes primal checking and agrees with HiGHS
within 2e-11, but its independent dual-bound result is unavailable. Its later
Phase I sign uncertainty remains unresolved. A successful root solve does not
establish the integer model's optimum or the validity of every subsequent node.

The [twelve-instance MIPLIB smoke](../results/miplib-smoke-2s-v2/runs/summary.html)
uses the same two-second budget, seed and configured size cap as the earlier
campaign. It retains three independently checked incumbents (`mas74`, `mas76`,
`pk1`), with no invalid candidates or reference contradictions. Six instances
exceed the size cap, four hit their time limit and two terminate numerically.
`mas74` retains an incumbent before a later loss of basis primal feasibility;
`gen-ip054` still encounters the conservative Phase I sign guard. No pinned
reference optimum is reached. The earlier `neos5` incumbent is absent at this
budget, so these results do not establish a uniform search improvement or a
speed advantage over the native baselines.

The next numerical work is stable ratio and pivot selection for the later
primal-feasibility loss, followed by checked Phase I infeasibility witnesses or
exact reconstruction for uncertain signs. Larger models still require removing
the standard-form capacity bottleneck. These remain separate from search
heuristics, cutting planes and independently checkable tree proofs.

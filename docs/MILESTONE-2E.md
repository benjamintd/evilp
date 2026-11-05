# Pivot safeguards and exact infeasibility witnesses

Version 0.5 extends the numerical solver with bounded exact Phase I
reconstruction and an independent original-model Farkas checker. Numerical
MILP statuses and bounds remain claims: individual LP certificates do not
certify a complete search tree.

## Exact contradiction

An exported certificate contains `MILPLAB_FARKAS_1` followed by original
constraint row names and rational multipliers. Omitted rows have multiplier
zero. A positive multiplier selects a row's finite lower side; a negative
multiplier selects its finite upper side. The independent checker recomputes
`q = A^T y` from the original MPS decimals and proves

```
sum_i y_i * selected_row_side_i - sup_{column box} q^T x > 0.
```

All arithmetic in this test is exact GMP rational arithmetic. A required
infinite side makes the certificate unavailable; a nonpositive margin is
invalid. No feasibility tolerance can make the contradiction positive.
An explicitly empty original column box is also an exact certificate. Integer
bounds are not rounded in this mode: the certificate proves LP infeasibility,
which also excludes integer solutions.

```sh
build/milpsolve MODEL --relax --farkas work/infeasible.farkas
build/milpcheck --model MODEL --farkas work/infeasible.farkas
```

The solver creates a fresh file only when it has a candidate. Its JSON still
reports `certifies_infeasibility: false`: the separate checker establishes
the exported witness's scope. This capability does not cover every numerical
infeasibility outcome, unbounded rays or a whole MILP tree.

## Reconstruction and source data

The optimizer reader optionally preserves original decimal coefficients,
duplicate sums, row ranges and bounds independently of its floating data.
The checker's parser and proof arithmetic remain separate code. Provenance is
bounded to 128 constraint rows, 256 columns, 20,000 input entries and bounded
numeric spellings. Larger inputs retain ordinary floating optimization.
Disagreement between exact and floating zero/equality topology disables
reconstruction rather than filling gaps with rounded coefficients.

When a small Phase I solve stops with positive artificial objective, the
solver reconstructs its basis dual over these original coefficients. It maps
the multipliers out of the LP's row units, then checks the exact original-row
and current-column-box contradiction. It does not trust the rounded auxiliary
RHS or infer proof from tolerance-small reduced costs.

Exact elimination is explicitly capped at 128 basis rows, 8192 bits per stored
numerator/denominator and one second within the existing LP deadline. This is
a small dense rational reconstruction, separate from the normal sparse
floating factorization. Singular, oversized, expired or unproved attempts
cannot bypass the uncertain-sign guard. In-memory models without decimal
provenance use exact binary input values; parsed models that lost provenance
cannot silently take that route. Bound overlays preserve original decimal
sides when unchanged and use exact binary values for changed in-memory sides.

`--no-exact` disables reconstruction for controlled comparisons. JSON records
attempts, successful exact contradictions and provenance availability.

## Checking recorded node LPs

With `--trace`, MILP search records each reconstructed contradiction together
with its full current column box, including fixed-integer repair calls.
These records carry explicit LP scope and do not establish that every branch
or prune in the search is valid.

The trace replay tool preserves original matrix, objective, RHS and range
records, substitutes the recorded column box, and calls the independent
checker on each certificate. It retains every replay model, certificate,
checker record and margin. Checking these LPs is separate from validating an
incumbent on the original integer model.

```sh
build/milpsolve MODEL --mip --trace work/search.jsonl
python3 scripts/check_trace_certificates.py --model MODEL \
  --trace work/search.jsonl --checker build/milpcheck --out results/fresh-proof-replay
```

## Pivot selection

The earlier `mas74` failure amplified a tolerated basic value of about
`-1.15e-12` through a pivot near `1.23e-8`, producing an actual new basic value
near `-9.39e-5`. Iterative residual refinement cannot repair that infeasible
basis. A small dyadic regression reproduces the amplification without copying
benchmark coefficients.

The normal path preserves Bland's entering and leaving order. It factors the
proposed new basis and applies the existing residual and primal guards before
accepting that pivot. A rejected pivot restores the original basis entry and
can retry one stronger Harris alternative. The retry must pass those same
actual basis checks. The accepted factor is cached for the next iteration, so
the ordinary path still factors once per pivot. Diagnostics count rejected
pivots and successful rescues; attempts consume the existing pivot budget.

The alternative uses two Harris passes. Every positive direction, including
one too small to serve as a stable pivot, limits the step. The second pass
prefers a larger pivot within 75% of the existing `1e-8` canonical primal
allowance, with deterministic index ties. A fused multiply-add check rejects
an unsafe predicted step; prediction alone never accepts the new basis. The
actual basis guard and original-model primal tolerance remain unchanged.

An unconditional Harris policy was rejected during development: it crossed a
tiny-direction bound during `mas74` Phase I, leaving a tolerated negative
value that the following weak pivot amplified. A separate two-pivot dyadic
regression preserves this failure mechanism. Keeping the successful Bland
path avoids this new root regression.

`--bland-ratio` preserves the earlier policy without a Harris retry for
controlled comparisons. Bland's full anti-cycling theorem is not claimed for
paths that use the rescue. Pivot and time budgets remain explicit limits.

## Recorded validation and remaining failures

The evaluated engine is frozen at `eaf03d4` (`milestone-2e-evaluator`).
The [contract evidence](../results/milestone-2e-contracts-v1/summary.json)
records all sixteen suites passing in both Release and ASan/UBSan builds,
with zero skipped or failed tests. Existing exact LP/MILP oracles and scaled
regressions remain in the suite. New coverage includes 180 independent
Farkas-oracle cases, 36 generated decimal reconstruction cases, parser/resource
limits, malformed witnesses, replay adversaries and the two pivot mechanisms.

The [controlled comparison](../results/numerical-stability-v1/report.html)
retains eight MILP runs: two declared development models, the frozen 0.4
binary, the new default and two ablations. Each has a ten-second budget and
1,000-node cap. Source snapshots before and after are identical and clean;
binary/model hashes, exact evaluator source, commands and raw records are
retained. This is a correctness comparison on observed failures, not a holdout
or a speed ranking.

`gen-ip054` previously stopped at node 35 with uncertain Phase I signs. The
new default reaches the 1,000-node limit and retains a separately checked
integer incumbent with objective approximately 8852.207808267. All 217 recorded
LP certificates pass independent exact replay with matching contradiction
margins. Disabling reconstruction reproduces the failure at node 35. Disabling
Harris retries still reaches 1,000 nodes with the same 217 checked certificates,
isolating the effect of reconstruction. The pinned reference optimum remains
unreached; these individual contradictions do not certify the entire tree.

`mas74` remains numerically unresolved at node 345. The guarded policy makes
three successful Harris rescues, then stops at the unchanged original-model
`1e-7` primal check. The earlier policy stops at the basis primal guard. Both
retain the independently checked incumbent near 19563.4046944. Moving the stop
to a later check is not a completed numerical repair. Original-unit feasibility
restoration and accurate primal reconstruction are the next work here.

The [native crosscheck](../results/mip-native-v4/summary.json) again matches all
forty known reference outcomes across HiGHS, SCIP, Cbc and warm/cold lab search;
all 35 feasible candidates pass the independent checker.

The [twelve-instance two-second MIPLIB smoke](../results/miplib-smoke-2s-v3/runs/summary.html)
records six size limits and six time limits, with valid incumbents on
`gen-ip054`, `mas76`, `neos5` and `pk1`. There are no invalid candidates,
checker errors, reference contradictions or reached reference optima. The
short budget stops `mas74` before its known later failure and yields no
incumbent there. Absence of numerical failures in that short run is therefore
not evidence that the remaining numerical issue is fixed, nor a uniform
performance improvement over the prior version.

Constant-row and conflicting-bound presolve outcomes can still be numerical
claims without an exported solver witness. The checker can validate explicit
certificates for those cases, but solver export does not yet cover every
presolve path. Larger exact reconstructions, complete tree proofs and the
standard-form capacity bottleneck remain separate work.

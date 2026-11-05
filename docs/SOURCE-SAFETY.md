# Source-sign repair and checked infeasibility candidates

These checkpoints address numerical stops exposed by the complete corpus run.
The independent checker, original-model acceptance policies and proof guards
remain separate from optimizer statuses. The frozen full-corpus v2 outcome is
still zero correct completions, 235 time limits and five numerical stops.

## Repair a genuinely negative basic value

The `mad` basis contains a source-exact basic value near `-6.804e-17` that the
floating path has clamped to zero. Exact reconstruction correctly rejects it,
but the old tolerance-based dual leaving rule selected no repair row.
An independent Fraction reconstruction of all 460 basis values and shifted
right-hand sides confirms the sign and verifies the original rows, including
the 11 masked rows. The row-mask metadata is correct.

A completed exact solve now retains one deterministically selected negative
basic-position/column hint. Original-unit dual repair may use that proven sign
to select an ordinary guarded pivot. It does not perturb the right-hand side
or round the negative source value away. A fresh exact reconstruction must
establish the next state. Basis, row, loaded-basis and RHS changes invalidate
the hint. Existing dual, pivot, residual and original-model checks still apply.

Restoration has at most 32 attempts, exact repeated-basis detection, and the
existing time and pivot budgets. Failed or unavailable repairs remain explicit
failures; they never establish infeasibility. Generated examples cover source
coefficients `1 + 10^-18`, `1 + 10^-20` and `1 + 10^-24` that round to binary
one, including loaded and natural bases, both backends, multiple repairs,
maximization, upper-bound coordinates and failure controls.

The captured `mad` root needs one repair pivot. A ten-second development run
reaches 568 nodes and retains an independently checked incumbent near 1.6162;
the locked optimum is 0.0268. The model remains unresolved. Evidence is retained
in `work/mad-audit/source-sign-validation.json`.

## Validate a larger Phase I candidate

Small dense exact reconstruction remains bounded. Larger models can instead
propose original-row multipliers from the numerical Phase I dual. Mapping
`scale*pi` uses exact dyadic GMP multiplication, including products outside the
binary64 range. The validator uses exact source data or conservative outward
source enclosures, and the current column box.

For signed row combination `y`, it computes a conservative lower bound on
the selected row-side sum and an upper bound on the box support of `A^T y`.
Only a strictly positive difference accepts the candidate. Infinite sides,
unbounded support, missing provenance, arithmetic/resource limits and deadlines
cannot manufacture a certificate. Changed column bounds replace their source
sides as exact binary overlays; unchanged parsed bounds retain source data.

The first constructor preserves the proposed row combination. If unavailable,
a second may zero multipliers whose signs select infinite row sides, after
exact original-row aggregation. This proposes a different combination and
reruns the entire unchanged proof check. It never clips a nonzero combined
column coefficient. Both attempts share one deadline of at most one second
within the LP budget. `--no-exact` disables this retry path.

Exported multipliers are rational strings. A source enclosure supplies a
positive margin lower bound, not an exact-margin claim. Trace events therefore
use `margin_kind: source_enclosure` and `contradiction_margin_lower_bound`.
Independent replay requires `actual_margin >= lower_bound > 0`. Existing exact
events retain equality checks. Mixed, unknown, zero and excessive margin claims
are rejected. These are individual LP-box proofs, not complete search-tree
certificates.

## Validation and selected diagnostics

- `ade4924`: all 24 Release and sanitizer suites pass. The
  [source-safety summary](../results/source-safety-contracts-v1/summary.json)
  records hashes and logs. Its source lock includes every C contract source.
- `e02a45e`: all 24 Release suites and focused certificate sanitizer checks
  pass. The raw and projected constructors each reject 1,000 independently
  generated feasible adversaries. Positive candidates are independently replayed.
- A retained 129-row decimal contradiction exceeds the old reconstruction
  limits. Its exported candidate and typed trace pass the independent checker;
  corrupted margin controls fail replay.
- The selected ten-second projected-candidate rerun gets past the original
  stops in `fhnw-binpack4-4` and `neos859080`; both then hit time limits.
  All 31 emitted local certificates pass independent replay, including every
  positive margin lower bound. None of the four runs submits a solution.
  Evidence is in `work/source-stop-probes-v2-audit/independent-audit.json`.

These diagnostics ran during engineering and sanitizer work. They establish
specific correctness behavior and local proof validity; they are unsuitable
for a speed ranking and do not change full-corpus solved counts.

`glass4` and `graphdraw-domain` still stop numerically in those selected runs.
For the captured `glass4` node, an exact positive combined coefficient near
`6.731e-21` multiplies a variable with no finite upper bound. The independent
checker confirms unbounded support. That nonzero residual must not be clipped
to zero. Alternative candidate construction is being investigated separately.

The [autoresearch gate](AUTORESEARCH.md) stays closed. No research paper has
been selected. Full benchmark completion and a later matched HiGHS comparison
over the unchanged 240-instance denominator remain required.

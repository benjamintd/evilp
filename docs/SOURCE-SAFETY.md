# Source-sign repair and checked infeasibility candidates

These checkpoints address numerical stops exposed by the complete corpus run.
The independent checker, original-model acceptance policies and proof guards
remain separate from optimizer statuses. The frozen full-corpus v3 outcome is
zero correct completions, 238 time limits and two numerical stops. The earlier
v2 inventory retains its 235 time limits and five numerical stops.

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

## Later validated checkpoints

`af61788` adds phase-local primal entering-column ranking while retaining the
full reduced-cost scan, existing eligibility predicates and every pivot guard.
The default chooses the most negative eligible canonical reduced cost. Exact
ties use column order. A validated zero step, or eight consecutive positive
short steps, switches permanently to the original policy for that phase.
`--bland-pricing` requests the original policy directly. Requested mode and
accepted-pivot/stall/fallback counters are recorded separately. Deferred pivots
are counted only after their basis is validated; a rejected trial is not an
accepted step. Both builds pass all 26 suites, including 140 independent
solver/checker pairs per profile; the
[pricing record](../results/pricing-contracts-v1/summary.json) binds the source
to the checkpoint. A generated family needs one productive pivot versus 32
under the original rule. That establishes controlled work reduction only.

`390fcba` bypasses identity `scalbn` operations for finite zero-exponent values.
All nonzero-exponent and rejection paths retain their guards. The work-only
equivalence check compares status, output bits and floating exception flags on
101,404 scalar cases and 1,256 public solves per profile, including signed
zeros, subnormals, extreme magnitudes, both backends and in-place operations.
Generated operation counts show the avoided scaling calls; concurrent timings
are not used as a speed ranking.

`332afbd` adds a third candidate constructor. It omits canonical multipliers
with exact magnitude at most `max(abs(pi))*2^-40`, aggregates the remaining
multipliers exactly, applies invalid-side projection and repeats the complete
source/support proof. Source coefficients and combined-column residuals are
unchanged. All three attempts share the existing deadline. The
[combined record](../results/simplified-source-contracts-v1/summary.json) passes
all 26 suites in both builds. Each constructor rejects 2,000 independently
generated feasible adversaries. Captured `glass4` and `graphdraw-domain` node
proofs validate after this candidate construction, but a new default-pricing
`glass4` stop appears later, at node 78.

The five ten-second default-pricing development probes then retain four time
limits and that one numerical stop. All 40 emitted local proofs pass replay;
`mad` submits the only solution, accepted independently at `1e-7`. Evidence is
in `work/source-stop-probes-v3/independent-audit.json`. These selected outcomes
do not replace the frozen full-corpus counts.

The node-78 residual involves meaningful multipliers rather than tiny terms.
Offline nearest-rational reconstruction with denominator limit `2^20` gives a
new vector whose complete source-enclosure proof has margin `13/420`; the
independent checker confirms it. The fourth C constructor is committed at
`4eead21`; all 26 Release and sanitizer suites pass at that checkpoint. A fresh
`glass4` diagnostic advances to node 459, where another numerical stop remains.
All 148 emitted local proofs pass independent replay in
`work/reconstructed-glass-proof-replay-v1/summary.json`. Its acceptance still
requires the complete unchanged proof. No nonzero combined coefficient is
clipped.

## Source-unit primal ratios and bounded retries

A captured `mas76` basis exposed a separate primal issue: scaled row slacks
with equation scale `2^-39` could be clamped by a global canonical `1e-7`
allowance despite material violations in original row units. Basic feasibility
and both leaving-row choices now use the existing original-unit allowance.
Every positive direction limits the step, including directions too weak for a
stable pivot. Tiny-step ties are relative, and the chosen actual binary step
must stay inside each row's own allowance. Exact dyadic unit-change controls
preserve the true limiting row in all 50 cases; the old tie rule crosses a
material allowance in 25 comparison cases.

Rejected trial bases restore and revalidate the old basis and factor before
continuing. Up to 32 unsafe entering candidates can be excluded for that basis
only, with the existing budget guards. A validated basis change clears the
exclusions. If every eligible candidate has been excluded, the result is an
unresolved numerical failure, never an optimum or unboundedness claim. A weak
positive-direction LP explicitly tests this behavior with updates enabled and
disabled and under the Bland-ratio control. Retained dual-feasible bases use
the guarded dual repair when source-unit basic feasibility requires it.

The captured `mas76` node now solves with zero independently checked row and
bound violations. Selected full `mas76` and Nidda runs pass their prior
numerical stops and reach time limits with checked incumbents. A captured
Nidda node reprices once and produces an independently checked exact local
Farkas contradiction. The patch's 26-suite Release and sanitizer evidence,
generic cases and selected positive/negative controls are retained in
`work/pivot-recovery-next/handoff.json`. Later combined tests additionally
exercise these paths with active root cuts. None of these selected outcomes
establishes full MIPLIB completion or matched solver speed.

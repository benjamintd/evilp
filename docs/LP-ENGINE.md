# Native C LP engine

The current implementation is version 0.7. The frozen milestone-2A contract
and its original validation records are preserved below. Later numerical and
witness policies are described in [milestone 2D](MILESTONE-2D.md),
[milestone 2E](MILESTONE-2E.md) and [milestone 2F](MILESTONE-2F.md).

## Current 0.7 implementation

The optimizer and independent checker remain separate C17 programs with
separate model readers. The numerical engine uses sparse CSC storage and
bounded GMP reconstruction where original-decimal provenance is retained.
`milpcheck` evaluates original, unpresolved models; optimizer presolve,
branch scores and numerical stopping statuses do not replace that check.

The sparse factor backend first removes structural singleton rows and columns.
A singleton row fixes a variable; a singleton column defers an equation for
reverse substitution. These steps introduce no fill. Only the remaining
coupled core is factored, while original coefficients are retained for forward
and transpose substitution. Up to two guarded compensated core corrections
precede recovery of deferred variables. A large sparse basis with many bound
and slack identities can therefore have a small core. A large remaining
coupled core can still require expensive elimination and fill; structural
peeling does not make every large MIPLIB model cheap to solve.

Guarded product-form eta updates reuse a base LU factorization across column
replacements. Normal and transpose solves apply the update chain in the
appropriate order. The current original-unit basis CSC remains available for
compensated residuals and iterative refinement. The unchanged componentwise
backward-error guard is `1e-10`, with a refinement target of `1e-14`; published
primal points retain the original `1e-7` feasibility policy.

An update with a weak pivot, nonfinite arithmetic or an inaccurate direction
is rejected without changing the old linear basis. The LP driver then factors
the complete new basis. Updated solves that fail their residual checks also
retry through full refactorization. Rejected primal pivot trials refactor the
alternative basis before continuing. At most 32 eta columns and 4,000,000 eta
entries are retained; reaching either guard requests refactorization. Storage
is proportional to sparse entries, elimination fill and retained eta entries,
with no dense fallback in the sparse path. `--no-updates` selects the
full-refactor ablation. A maintained membership bitmap makes each pricing
membership lookup constant-time; pricing still examines candidate columns.

### Capacity and controls

The sparse standard-form guard defaults to 16,000,000 and the sparse basis
dimension ceiling is 8,000,000. For transformed structural count `z` and
expanded equation count `r`, the standard form reserves `z + 2r` columns.
`--max-dimension` controls that finite standard-form guard. Entry and fill
budgets can also stop a solve explicitly. These are resource limits, not
claims that all problems below them fit memory or complete within a budget.

`--backend dense` remains an explicitly selected reference backend. It defaults
to a 4096 standard-form limit and refuses larger requested limits before dense
allocation. Increasing sparse capacity does not increase the dense ceiling.
Scaling and refinement remain separately controlled by `--no-scaling` and
`--no-refinement`; none of these flags relaxes original-model validation.

MILP search applies source-aware root bound propagation before opening LP
sessions. It uses retained exact decimals or conservative outward source
envelopes, integer bound rounding, and redundant-row masks for the returned
box. Reductions require trustworthy source information. The original model
stays immutable, and retained candidates are checked against its unpresolved
rows. Presolve receives at most 32 passes and a budget of
`min(2 seconds, 5% of the remaining MILP budget)`, within the cumulative solve
budget. `--no-presolve` disables it. Presolve statuses remain optimizer claims;
this module does not emit an independent infeasibility certificate.

Branch selection uses deterministic reliability-based pseudocosts learned from
solved children. Two valid observations per direction qualify a learned
estimate; incomplete history uses fallback statistics or most fractionality.
Scores rank eligible fractional integer variables and never establish a bound,
feasibility or a pruning justification. `--first-fractional` restores the
original selection control. Both search controls are MILP-only. `--cold`
continues to rebuild LPs instead of reusing saved sessions.

JSON exposes update acceptance/rejection counts, refactorizations, peak core
dimension and singleton eliminations, alongside presolve and branch-observation
diagnostics. These are work counters and numerical diagnostics, not speed
results or independent certificates.

### Current validation and completion gate

The current [Release contract run](../results/milestone-2g-contracts-v1/release-ctest.log)
passes all 22 suites. The previous 18 remain, with additive presolve, branch
selection, benchmark-gate and scalability contracts. The generated scalability
tests check 2,500-variable LPs with 5,003 basis rows and a three-variable core,
known original-model points and objectives, explicit dense/sparse limits, and
refactor-count reduction against `--no-updates`; they assert no timing ranking.
The [full sanitizer run](../results/milestone-2g-contracts-v1/sanitizer-ctest.log)
also passes all 22 suites at engine checkpoint `82c1355`. The
[contract summary](../results/milestone-2g-contracts-v1/summary.json) records the
matching source/test hashes. The frozen
[18-suite milestone-2F evidence](../results/primal-regressions-v1/final-validation.json)
continues to describe its own Release and sanitizer builds.

The [full 240-instance v1 inventory](../results/miplib-candidate-v1/inventory.json)
records 209 unsupported cases, 31 time limits, seven independently checked
incumbents and zero correctly completed problems. It used the earlier candidate
at a two-second budget and is explicitly not a fair timing campaign. The
[completion and paired-speed gate](../results/miplib-candidate-v1/gate.json)
remains closed. Capacity guards and a large coupled core are separate blockers;
new engineering requires a fresh versioned campaign, not edits to that inventory.
No performance advantage over HiGHS or other solvers is established.

## Frozen milestone-2A contract

The remainder describes the original milestone-2A checkpoint, including its
historical defaults, limitations and next gate. Subsequent sparse storage,
reoptimization and MILP work is preserved in [milestone 2B/C](MILESTONE-2BC.md).

`milpsolve MODEL --relax` solves the continuous relaxation using project-owned
C17 code. It does not call HiGHS, SCIP, Cbc or another optimizer. GMP belongs
to the separate checker; the numerical solve uses floating-point arithmetic.
This is a correctness-focused reference engine, not a competitive MIP solver.

### Model and numerical architecture

`src/model.c` implements the optimizer's own plain/gzip MPS reader. It preserves
original row/column names, ranged row bounds, column bounds, integrality,
objective sense and offset in an immutable sparse CSC model. It does not reuse
the independent checker's parser. All 240 imported MIPLIB models were parsed
successfully; that does not imply they fit or solve in this first backend.
Duplicate coefficients are coalesced. Nonfinite, overflowing, underflowing,
malformed and unsupported inputs are rejected explicitly.

`src/lp.c` converts the model into nonnegative standard form. Lower-bounded
variables are shifted, upper-only variables reflected, free variables split,
and fixed variables substituted. Finite upper bounds become equations with
slacks. Signed rows, ranged rows, objective offsets and original-variable
reconstruction are preserved. A Phase I feasibility solve and artificial-column
cleanup precede revised primal simplex. Index-based entering and ratio-tie rules
make choices deterministic for a given binary and floating-point environment.

`src/factor.c` provides factorization, solve and transpose-solve operations.
The first backend uses dense partial-pivoted LU and refactors at every pivot.
It has no factor updates, sparse LU or dual reoptimization yet. The model is
sparse; the transformed LP and basis computations in this backend are dense.

Let `z` be the transformed structural variable count and `r` the expanded
equation count, including finite column-upper-bound equations. The reference
backend reserves `z + 2r` standard-form columns. It rejects problems exceeding
`--max-dimension` before dense allocation: default 512, hard ceiling 4096.
This explicit cap prevents a large MIPLIB model from triggering an accidental
dense allocation. Increasing it does not make the algorithm suitable for large
models; factorization cost is cubic in the equation count.

Original-model primal residuals are recomputed with long-double accumulation;
residuals exceeding 1e-7 cause numerical failure. Pivot selection accounts for
objective scale and solve uncertainty. A finite floating support-bound gap
check rejects materially inconsistent stopping points. These are numerical
safeguards, not formal certificates. Singular/unstable bases, overflow and
ambiguous cleanup have explicit failure paths. There is no advanced scaling,
perturbation, iterative refinement, presolve, cutting plane or integer search.

### Command and artifact contract

```sh
build/milpsolve MODEL --relax \
  --solution NEW.sol --dual NEW.dual --basis NEW.basis \
  --time-limit 60 --iteration-limit 100000 --max-dimension 512
```

Stdout is one JSON object. It records original and standard-form dimensions,
objective, iterations, refactorizations, redundant rows, primal/reduced-cost
residual diagnostics, model hash and elapsed times. Solver statuses are
`optimal`, `infeasible`, `unbounded`, `time_limit`, `iteration_limit`,
`numerical_failure`, `unsupported` or `input_error`. They remain numerical
solver claims. `certifies_optimality` is always false in optimizer output.
Exit codes are 0 for completed numerical statuses, 2 for input/output errors,
3 for unsupported dimensions/formats, and 4 for limits/numerical failure.

The budget includes hashing and reading, but is checked after reading and
between simplex operations; an in-progress read or LU operation is not
preempted. The benchmark runner supplies an external watchdog. Limits do not
currently export a partially solved LP state. Solution/dual/basis output uses
exclusive creation and never replaces an existing file.

The simple solution file contains original variable names and a claimed
objective. Output uses the shortest significant-digit decimal that roundtrips
to each computed double. The dual file contains original constraint-row
multipliers for normalized minimization (`-c` for original maximization).
These are proposed witnesses; their textual precision does not establish
correctness.

`--read-basis FILE` restores a saved Phase II basis for the same model hash,
dimensions and canonical row identities. Phase I and redundant-row discovery
are replayed before restoration, so a basis file cannot instruct the engine
to drop constraints. This first persistence format does not support bound
changes, cross-model warm starts or cross-version compatibility. Saved bases
are hints and are checked numerically before use.

The harness has an explicit `lab` adapter. It requires `milpbench --relax`,
uses the simple solution format, retains the proposed dual file and invokes
the checker with relaxation mode explicitly enabled. LP mode currently accepts
only that adapter. Baseline solvers are run separately on continuous MPS input.
Suite references for LP runs must be LP objectives, never MIPLIB's integer
reference objectives. Reports mark relaxation mode and separate checked bounds,
reference agreement and exact certificates.

### Independent rational bounds and exact witnesses

```sh
build/milpcheck --model MODEL --relax --solution NEW.sol --dual NEW.dual
build/milpcheck --model MODEL --relax --solution FRACTION.sol \
  --format rational --dual FRACTION.dual --exact
```

`--relax` explicitly suppresses integrality checks and records that scope.
The original integer-model check remains the default. `--format rational`
allows numerator/denominator values in primal witnesses; MPS and native SOL
formats retain their decimal grammar. Dual files accept decimals or fractions,
reject unknown/nonconstraint/duplicate names, and default missing multipliers
to zero.

For normalized minimization, the checker independently recomputes
`d = c - A^T y` in exact rational arithmetic. It sums the objective offset,
row-side supports selected by the signs of `y`, and column-bound supports
selected by the signs of `d`. Any rational multiplier vector produces a safe
bound when all selected sides are finite. Selected infinite sides make the
bound unavailable; they do not become zero or a fabricated certificate.
Maximization reverses the normalized bound into an upper bound.

JSON distinguishes `bound_status` (`valid`, `unavailable`, `absent`),
`certified_bound`, `bound_direction`, tolerance-based primal acceptance and
`certifies_optimality`. A bound-only check never claims primal feasibility.
A valid primal plus an unavailable bound can still pass primal checking; both
fields must be inspected. Exact optimality is certified only when an exactly
feasible rational primal matches the independently derived bound and all
tolerances are zero. A relaxation certificate proves the LP optimum, not the
integer optimum. An exactly integer-feasible witness matching a valid LP bound
can also establish an integer optimum.

Floating decimal output does not always admit exact certification. In particular,
nonterminating rational optima and small multiplier residuals on free columns
can require exact basis reconstruction, which this engine does not implement.
Unavailable certificates stay visible. Farkas infeasibility certificates,
unbounded-ray verification and complete MILP proof logging are not implemented.

### Validation and next gate

Completed validation at clean source revision `f94761e`:

| Check | Result |
| --- | --- |
| Release CTest | All six suites pass |
| AddressSanitizer / UndefinedBehaviorSanitizer CTest | All six suites pass |
| Generated exact LP oracle | 1,000 cases agree |
| Native continuous LP comparisons | 32/32 candidates pass and match exact references within the recorded comparison tolerance |
| Independently checked lab bounds | 8/8 valid and on the safe side of the exact optimum |
| Exact exported primal/dual example | LP optimum 853/100 independently certified |

The native comparison consists of eight fixtures × HiGHS, SCIP, Cbc and the
lab solver. Maximum objective error is about 5.24e-9, from Cbc's printed
solution precision. There are 31 clean process exits and one documented HiGHS
row-name output warning with exit code 1. Its raw `claimed_status` remains
`process_error`, while `reported_solver_status` is `optimal`; the exact-reference
and primal checks pass. The integration script accepts only that specific
known warning pattern and retains it separately. No timing ranking is claimed.

Artifacts:

- [Final native crosscheck summary](../results/lp-native-v2/summary.json),
  [native baseline report](../results/lp-native-v2/native/summary.html) and
  [lab report](../results/lp-native-v2/lab/summary.html).
- [Contract evidence](../results/lp-contracts-v1/summary.json) and
  [full sanitizer log](../results/lp-contracts-v1/sanitizer-ctest.log).
- [Exact LP certificate example](../results/lp-exact-demo-v1/checker.json).
- [All-240 reader audit](../results/lp-parser-audit-v1/summary.json).

The earlier [failed diagnostic campaign](../results/lp-native-v1/summary.json)
is preserved. It exposed the HiGHS primal/dual-section parsing bug; final
results use a fresh directory, fixed checker, frozen binaries and source lock.
Milestone-1 MIPLIB runs and their evaluator tags remain untouched.

CTest covers 1,000 generated bounded rational LPs using independent Fraction
vertex enumeration, along with free/fixed/upper-only variables, ranges,
min/max offsets, degeneracy, infeasibility, unboundedness, malformed input,
limits, basis tampering, original-model reconstruction and numerical regressions.
Separate tests exercise 120 randomized exact support bounds and the native
solver-to-checker boundary, including exact and unavailable certificates.
The optional `tests/test_live_lp.py` compares eight known continuous LPs across
all four engines with retained raw records and exact reference values.

The next gate is sparse factorization and reliable basis reoptimization, then
minimal branch-and-bound on exhaustively checked tiny MILPs. Do not use this
prototype's claimed LP bounds to make certified pruning decisions. Keep
evaluation changes separate from future optimizer trials, and preserve the
milestone-1 source tags and raw baseline records.

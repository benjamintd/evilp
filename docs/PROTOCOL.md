# Milestone 1 evaluation contract

The C runner launches native command-line solvers without a shell and validates
every solution file with a separate C process. The checker links GMP and zlib,
and links no optimization solver.

## Dataset

MIPLIB 2017 benchmark v2, all 240 instances; reference file v37.
config/miplib.lock.json pins downloads and individual compressed model hashes.
The importer refuses changes to an existing lock. MPS files stay compressed.
Dataset files are ignored by Git and fetched from their original publisher.
Per-instance redistribution rights remain unresolved in the manifest; do not
ship the fetched dataset under the project's MIT license. The MIPLIB publication
says CC BY-SA 4.0 is the submission default but allows exceptions.

The public dataset is a development/reference corpus. It is not a genuinely
unseen evaluation set for future agents. Hidden generators and prospective
real-world model families belong in a later milestone.

## Checker semantics

Exact GMP rational arithmetic evaluates original decimal/scientific MPS input
and decimal solution values. No floating-point presolve, solver model conversion,
or solver feasibility API is reused. Two streaming passes keep storage proportional
to rows and variables, rather than the number of matrix entries.

Supported: MIN/MAX, OBJNAME, N/L/E/G rows, objective RHS offset, COLUMNS
integer markers and continuation records, named RHS, signed RANGES, LO/UP/FX/
FR/MI/PL/BV/LI/UI bounds, gzip/plain MPS, and HiGHS/SCIP/Cbc/simple solution
files. Extra numeric values on BV/FR/MI/PL records are parsed and ignored.
Unknown rows/variables, duplicate solution assignments, malformed/nonfinite
numbers, empty or missing solutions, incomplete models, and unsupported extensions
fail closed. Only one RHS/range/bound vector is supported. Blank RHS vector
names, embedded spaces within fixed-width names, semicontinuous variables,
SOS, indicators, and quadratic/nonlinear extensions are unsupported.
Conflicting integer-marker declarations for the same column are rejected.

Default absolute tolerances: linear rows/bounds 1e-4; integrality 1e-5;
claimed-objective consistency 1e-4. All comparisons use rationals.
--exact sets all three tolerances to zero; this checks a supplied finite-decimal
witness exactly. By itself it does not prove optimality or infeasibility.
Missing variable assignments mean zero (sparse SOL convention).

Milestone 2A adds explicit --relax checking, --format rational for fractional
witnesses, and --dual for independently recomputed LP support bounds. Exact
feasible primal and matching dual witnesses can certify optimality; tolerant
acceptance and reference agreement cannot. Bound-only checks never claim primal
feasibility, and infinite selected sides make a bound unavailable. See the
[LP contract](LP-ENGINE.md) for proof scope and limitations. Historical baseline
evaluator tags and records retain their original primal-only semantics.

Exit codes: 0 valid/supported; 1 invalid candidate; 2 malformed/input error;
3 unsupported model; 4 no candidate. JSON includes every violation count,
largest residual, recomputed rational objective and tolerances.
The rows count includes N rows; nonzeros counts column entries, including
objective entries and explicitly stored zeros.

## Run policy

One solver process at a time, rotating solver order by instance index.
Positive seed 1, one requested thread, zero requested absolute/relative MIP gap.
Solver defaults otherwise apply. External numerical libraries receive one-thread
environment settings; requested thread limits do not establish bitwise determinism.
Internal time limits follow each solver's own clock and exclude different setup
work. Always compare recorded end-to-end wall time as well. A separate watchdog
signals the process group with SIGTERM after limit + grace, escalating to
SIGKILL one second later if needed. Checker time is excluded.
Do not treat a short-budget pilot as a fair speed ranking.

Every run records raw argv/settings/log/solution, exit/signal/watchdog outcome,
process wall/CPU/peak RSS, checker JSON, reported status/bound, references,
and SHA-256 identities of model, solver binary, runner and checker.
The environment record pins suite and solver-config hashes.
Release-archive digests and source-build provenance complement binary hashes;
hashing only an executable does not freeze its system libraries.

Reported dual bounds and optimal/infeasible statuses remain solver claims.
The lab's separately checked certified_bound field is distinct from a solver's
reported_dual_bound. The raw reported_solver_status is also distinct from the
effective claimed_status, which can be process_error after a nonzero exit.
Raw exit codes and signals are preserved even when a retained witness validates.
Cbc's explicit terminal solution header can refine an unknown/ambiguous log
status only after a clean process exit. Both reports and their source remain
recorded. Conflicting explicit claims become unknown; timeouts, signals and
nonzero exits keep priority. Reference objectives never influence status parsing.
Reference matching uses 1e-4 absolute + 1e-9 relative error on =opt= objectives.
It is not an optimality certificate. Contradictions remain visible in reports.
Reference metadata can be revised upstream; upgrades require a new lock
and a new experiment directory.

The original timed campaign uses the evaluator-v1 Git tag. Subsequent checker
updates can be applied to retained solutions with scripts/revalidate.py; this
writes new validation records and revalidated.jsonl, without changing run.json,
solver logs, solutions, or timings. Report these with report.py --runs
revalidated.jsonl --prefix validated. Freeze the original evaluator from its
Git tag and record the revalidation evaluator separately.

Each new run records the checker argument vector and SHA-256 digests for
the primal and dual files it retains. Replay verifies those digests and the
model hash, substitutes the selected checker, and resolves evidence paths in
the current run directory. For older records it reconstructs the original
adapter format and relaxation mode from the run metadata. A dual-only run can
reproduce a certified LP bound. A checker invocation failure or a missing or
changed recorded artifact produces an error row and a nonzero revalidation
exit. Runs that produced no evidence appear separately in the summary.

The Cbc adapter explicitly passes the MPS objective sense, because Cbc's reader
ignores OBJSENSE. Maximization regression tests guard this behavior. The pinned
MIPLIB benchmark models use minimization; this adapter correction does not
reinterpret their original timed runs.

The lab adapter requires explicit --relax with LP references, or --mip with
integer-model references. These modes are exclusive. LP mode excludes native
MILP adapters; native LP crosschecks provide continuous MPS models in a separate
run directory. MIP mode permits the lab alongside HiGHS, SCIP and Cbc and checks
integrality normally. Lab MILP best_bound values remain uncertified numerical
claims; no complete branch-and-bound proof is logged.

On the general SCIP route, a terminal optimal, gap-limit, time-limit or
node-limit result publishes SCIP's finite dual bound in original objective
units, or null when unavailable. A checked incumbent can differ from that
bound even when the solver reports optimality. The runner retains the number
in `reported_dual_bound`. Replay and reporting compare it with checked
witnesses under the evaluator's existing policy. A general SCIP attempt that
hands control back to native search contributes no terminal bound of its own.

The general importer can multiply a row and both finite sides by a power of
two to preserve small coefficients. It checks the stored sides, coefficient
values and variable identities. If a finite value would reach SCIP's infinity
sentinel, it rejects the import. When a rescued term or finite side falls
within feasibility tolerance in unit-normalized coordinates, a backend time
or node limit may return a checked point with a null bound. Such points also
undergo a residual check relative to the source row activities. A purported
closure on these rows returns control to native search with fallback reason
`ill_conditioned_row_requires_source_check`.

The MIP JSON, runner ledger, replay and report carry four provenance fields.

| Field | Interpretation |
| --- | --- |
| `search_route` | Search that supplied the final verdict. An abandoned SCIP attempt can precede a native verdict. |
| `closure_reason` | How that route closed, such as a backend status, exhaustion or a matching structural bound. A limit has no closure. |
| `bound_origin` | How the route obtained its bound, such as the backend dual, an open frontier or exhaustion. |
| `bound_scope` | `original_mip` for a published original-objective MIP bound and `none` when no bound is available. |

Missing provenance in older records stays null. The evaluator reports a relative
gap only when a checked objective and a compatible bound have usable signs, a
clean process result and consistent evidence. Equal values give zero. For two
nonzero values of the same sign, the gap is the directed difference divided by
the smaller absolute value. The gap is descriptive and carries no proof credit.

## Verifying a candidate and replaying fleet controls

The default optional SCIP build and feature CI use the lock named in
[`candidate-scip-lock.txt`](../config/candidate-scip-lock.txt). Configure a Release
build with that dependency, GLPK and Kissat, build all targets, then run:

```sh
python3 scripts/verify_candidate.py --build "$BUILD" --out "$OUT/release.json"
python3 scripts/check_fleet_release.py --build "$BUILD" \
  --verification "$OUT/release.json" --out "$OUT/fleet"
```

Use fresh output paths. Verification binds the tested executables and dependency
archives to the CTest results. The fleet check rechecks the historical controls
on the original models, runs the six locked snapshots at 60 seconds each, and
scores the results against the archived witnesses. Its `release-check.json`
lists the outcome of every attempt. A checked schedule from a timeout can
disprove another solver's claim while the timed-out search remains unresolved.
The archives are available only on hosts holding the private fleet inputs.

## Generating matched runs and development panels

Use a suite TSV with locked model hashes and references, and a solver TSV
with the paths and digests of the Lab and HiGHS executables. Use a fresh
directory for each invocation. The runner alternates solver order by instance
and records the campaign, budget, host, engine configurations and evaluator
hashes before launching the first pair. It keeps both solvers' witnesses.
The timed interval covers the solver process. Checking and file hashing happen
outside that interval.

```sh
milpbench --mip --paired --suite "$SUITE" --solvers "$SOLVERS" \
  --checker "$CHECKER" --limit 60 --grace 5 --out "$OUT"
python3 scripts/benchmark_gate.py --runs "$OUT/runs.jsonl" \
  --baseline-runs "$OUT/runs.jsonl" --campaign "$OUT/campaign.json" \
  --panel-suite "$SUITE" --out "$OUT/panel-score.json"
```

This command records a development panel. The score lists checked feasibility,
reference agreement, completion and timeouts for each engine. A complete panel
ledger means every selected model has one attempt per solver. For a formal run,
use the full locked suite, add `--quiescent` only on a reserved idle host, and
omit `--panel-suite` when scoring. The operator attests that the host is idle.
The gate also requires matching run identities and usable suspension readings.
A standalone candidate completion inventory passes
`--candidate-config "$OUT/lab-config.json"` to the gate.

Lab configuration records capture `LAB_SUBSCIP_FRAC`, phase tracing and the
contents of `LAB_SCIP_SET` when provided. The runner executes against the
saved parameter file. Paired runs use the built-in SCIP parameter profile
and require phase tracing to be off. The configuration digest binds the
runner's adapter defaults to the solver binary and the chosen budget.

For solver diagnostics, `milpsolve --mip --no-scip` opts out of both SCIP
routes. The MIP JSON field `scip_route` records dispatch and `scip_fallback_reason` identifies why general search returned
control to the native solver. Enabling `--trace`
records the general SCIP boundary and native events. The trace does not
contain an independent SCIP tree certificate. The `--no-root-cuts` and
`--no-root-diving` switches control the native heuristics.

## Reproducibility boundaries

This milestone records an evaluator freeze; it does not enforce a security
boundary against an agent able to edit the evaluator, binary or operating system.
Before autoresearch, keep evaluation in a separately controlled checkout/service.
Warm caches, laptop power management, thermal state, system load, and one seed
limit performance interpretation. Multi-seed, longer-budget runs are confirmation
work, not evidence supplied by this pilot.

The binary SAT route reports the argument that closes its objective bound.
`objective_box` means a feasible point reaches the best permitted column endpoints.
`source_equalities` means an exact combination of original equalities fixes the
objective. `row_penalty` identifies a lower bound on the penalty from satisfying
the rows. `sat_cutoff_unsat` means the SAT search excluded the remaining better
objective values. The same origin applies when an earlier unsatisfiable cutoff
raised the bound and a later witness attained it. A limited probe can supply a
witness. A subsequent route reports its own bound if it establishes one.
Independent proof status comes from the checker.

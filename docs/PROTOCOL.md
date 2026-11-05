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

The Cbc adapter explicitly passes the MPS objective sense, because Cbc's reader
ignores OBJSENSE. Maximization regression tests guard this behavior. The pinned
MIPLIB benchmark models use minimization; this adapter correction does not
reinterpret their original timed runs.

The milestone-2A lab adapter requires --relax and an LP-reference suite; this
mode currently excludes native MILP adapters. Native LP crosschecks instead
provide continuous MPS models to those adapters in a separate run directory.

## Reproducibility boundaries

This milestone records an evaluator freeze; it does not enforce a security
boundary against an agent able to edit the evaluator, binary or operating system.
Before autoresearch, keep evaluation in a separately controlled checkout/service.
Warm caches, laptop power management, thermal state, system load, and one seed
limit performance interpretation. Multi-seed, longer-budget runs are confirmation
work, not evidence supplied by this pilot.

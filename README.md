# MILP Lab

A C17 laboratory for building and evaluating a new MILP solver.
Milestone 1 implements the dataset import, independent primal checker,
native benchmark runner and HiGHS / SCIP / Cbc baseline adapters.
The current prototype adds sparse LP storage and LU factorization, saved-basis
dual reoptimization, and a minimal deterministic branch-and-bound MILP solver.
It also applies reversible power-of-two changes of units and compensated
iterative refinement, with the original feasibility policies unchanged.
Small uncertain Phase I solves can reconstruct exact original-decimal
infeasibility witnesses; guarded pivots can retry a stronger leaving choice.
Small-model primal recovery reconstructs original decimal arithmetic and can
restore infeasible bases using feasibility thresholds in original units.
Version 0.7 peels structural singleton rows and columns before factoring the
remaining coupled core. Guarded basis updates reuse that factorization across
pivots, and a membership bitmap avoids repeatedly scanning the basis during
pricing. Integer search adds conservative source-aware root propagation and
deterministic reliability-based pseudocost selection.
The independent checker validates original-model witnesses, rational LP bounds
and exact infeasibility certificates.

See [sparse basis updates and full-corpus coverage](docs/MILESTONE-2G.md),
[original-unit primal recovery](docs/MILESTONE-2F.md),
[pivot safeguards and infeasibility witnesses](docs/MILESTONE-2E.md),
[numerical scaling and refinement](docs/MILESTONE-2D.md), and
[the integer-search contract](docs/MILESTONE-2BC.md).
The [primal recovery comparison](results/primal-recovery-v1/report.html),
[40-case MILP crosscheck](results/mip-native-v5/summary.json), and
[MIPLIB smoke report](results/miplib-smoke-2s-v4/runs/summary.html) preserve
successful witnesses, numerical failures and capacity limits. The earlier
[pivot and proof comparison](results/numerical-stability-v1/report.html),
[numerical ablations](results/numerics-roots-v1/report.html),
[MILP report](results/mip-native-v2/report.html) and
[MIPLIB smoke](results/miplib-smoke-2s-v1/runs/summary.html) remain available.

The [autoresearch protocol](docs/AUTORESEARCH.md) starts only after every MILP
problem in the existing benchmarks passes. That gate is not yet complete.
The [full 240-instance v2 inventory](results/miplib-candidate-v2/inventory.json)
records zero unsupported cases, 235 time limits, five numerical stops, nine
independently checked incumbents and zero correctly completed instances.
All 240 attempts and all failures are retained. Concurrent engineering makes
this development coverage evidence unsuitable for a speed ranking.
The [gate result](results/miplib-candidate-v2/gate.json) remains closed. The
subsequent performance target is at least 120 strict paired wins against HiGHS
over the fixed 240-instance denominator under matched conditions.

The evaluated [Release contract log](results/milestone-2g-contracts-v1/release-ctest.log)
passes all 22 suites. These retain the previous 18 suites and add presolve,
branch selection, benchmark-gate and scalability contracts. All 22 suites also
pass in the [sanitizer run](results/milestone-2g-contracts-v1/sanitizer-ctest.log).
The tested engine is committed at `82c1355`; the earlier
[18-suite Release and sanitizer evidence](results/primal-regressions-v1/final-validation.json)
remains frozen and unchanged.
The [fresh eight-model native crosscheck](results/mip-native-v6/summary.json)
matches all 40 expected outcomes; all 35 feasible candidates pass independent
checking. Later column-norm caching (`303d9bf`) and transactional column
overlays (`4c02d7e`) also pass all 22 Release suites and their focused sanitizer
checks. They were added after the frozen v2 campaign; its results do not
measure those changes.

Milestone 2A is preserved in [the original LP engine contract](docs/LP-ENGINE.md),
[the lab's LP report](results/lp-native-v2/lab/summary.html), and
[the four-engine crosscheck](results/lp-native-v2/summary.json).

Milestone 1 is implemented and exercised on all 240 MIPLIB benchmark instances.
See [the baseline summary](docs/BASELINES.md),
[the complete 720-run report](results/miplib-2s-seed1/validated.html), and
[the 36-run smoke report](results/smoke-5s-seed1/validated-final.html).
Reports, downloaded models and native tools are present in this local project;
the source archive excludes these fetched/generated artifacts.

## Build and test

Requirements: C17 compiler, CMake 3.20+, GMP, zlib, Python 3.12+ for setup
(test/report scripts also run on 3.9), curl, and Git. Cbc's source build also
requires a C++ compiler and make.

~~~sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure
cmake -S . -B build-sanitize -DLAB_SANITIZE=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build-sanitize -j 4
ctest --test-dir build-sanitize --output-on-failure
~~~

## Setup

~~~sh
python3 scripts/setup.py solvers
python3 scripts/setup.py miplib
python3 scripts/audit_corpus.py
~~~

Pinned native assets target macOS arm64: HiGHS 1.15.1,
SCIP 10.1.0 / SoPlex 8.1.0, and Cbc 2.10.13. The Cbc published asset
requires a GCC runtime absent on this host; setup invokes build_cbc.py to select
a local build from the locked dependency tree instead. Linux users can supply native
executables and their SHA-256 digests in config/solvers.tsv.
Paths in suite/config TSVs resolve relative to the current directory:
run commands from this project root. No system-wide solver installation is needed.

## Check a solution

~~~sh
build/milpcheck --model tests/fixtures/tiny.mps \
  --solution tests/fixtures/valid.sol --format simple --exact
~~~

Use --format highs, scip, or cbc for native solution files.
See [the protocol](docs/PROTOCOL.md) for accepted MPS features, exact arithmetic,
absolute tolerances, and the difference between feasible witnesses and proof.

## Solve an LP relaxation

~~~sh
mkdir -p work
build/milpsolve tests/fixtures/tiny.mps --relax \
  --solution work/tiny-lp.sol --dual work/tiny-lp.dual --basis work/tiny-lp.basis
build/milpcheck --model tests/fixtures/tiny.mps --relax \
  --solution work/tiny-lp.sol --dual work/tiny-lp.dual --exact
~~~

The example's LP optimum is 8.53; the original integer optimum is 8.03.
The separate checker certifies this example's exact LP optimum. In general,
floating-point solutions are tolerance-accepted and their dual bounds can be
unavailable; inspect both statuses. Output files must be new. `--relax` is
required so fractional LP solutions cannot silently be presented as MILP results.
Sparse LU is the default. The standard-form capacity guard defaults to
16,000,000, with a separate sparse basis ceiling of 8,000,000. Entry and fill
budgets remain finite. `--max-dimension` selects a smaller standard-form guard;
`--backend dense` selects the reference backend and defaults to its hard 4096
limit. Large remaining coupled cores can still be expensive after singleton
peeling; capacity limits do not promise successful completion or sufficient
memory for every model.
Scaling and refinement are enabled by default. `--no-scaling` and
`--no-refinement` select controlled ablations; the residual guard still applies.
Basis updates are also enabled by default. At most 32 eta columns and 4,000,000
eta entries are retained before full refactorization. Weak, inaccurate or
nonfinite updates also fall back to refactorization. `--no-updates` provides
the full-refactor control. See [the current engine contract](docs/LP-ENGINE.md)
for storage and validation details.

~~~sh
build/milpsolve tests/fixtures/tiny.mps --relax --read-basis work/tiny-lp.basis
python3 tests/test_live_lp.py --out results/my-live-lp
~~~

The optional live LP test requires the three installed baseline solvers. It runs
eight continuous models through all four engines and checks known objectives
and the lab's independent bounds. Its generated lab configuration also shows
how to use `milpbench --relax` with an LP-reference suite. Native baseline
adapters operate on continuous MPS models in a separate run directory.

## Solve a MILP

~~~sh
build/milpsolve tests/fixtures/tiny.mps --mip \
  --solution work/tiny-mip.sol --trace work/tiny-mip-trace.jsonl
build/milpcheck --model tests/fixtures/tiny.mps --solution work/tiny-mip.sol --exact
build/milpsolve tests/fixtures/tiny.mps --mip --cold
python3 tests/test_live_mip.py --out results/my-live-mip
~~~

The tiny integer optimum is 8.03. Search defaults to warm LP sessions; `--cold`
rebuilds each LP for comparison. Use `--node-limit` and `--time-limit` to bound
search. Root presolve receives at most 32 propagation passes and a budget of
`min(2 seconds, 5% of the remaining MILP budget)`. It uses source-aware bounds
and redundant-row masks without changing the model checked independently.
`--no-presolve` disables it. Learned branch scores use solved-child observations;
incomplete history falls back deterministically, and `--first-fractional`
restores the original branch-selection control. These two flags are MILP-only.
Limits retain any incumbent found so far. Search statuses and global
bounds remain numerical claims: complete tree proof logging is not implemented.
An unbounded LP relaxation returns `lp_unbounded`, without asserting MILP
unboundedness. `--mip` and `--relax` are explicit, exclusive modes.

The `lab` benchmark adapter accepts `milpbench --mip` alongside the native
solvers. LP runs use `--relax` and LP reference objectives in separate suites.

## Run baselines

~~~sh
build/milpbench --suite config/smoke.tsv --out results/my-smoke --limit 5 --seed 1
python3 scripts/report.py results/my-smoke
build/milpbench --suite data/miplib2017/suite.tsv \
  --out results/my-miplib-pilot --limit 2 --seed 1
python3 scripts/report.py results/my-miplib-pilot
python3 scripts/freeze.py results/my-miplib-pilot/evaluator-lock.json
~~~

Retained solutions can be checked again without rerunning the solvers:

~~~sh
python3 scripts/revalidate.py results/my-miplib-pilot --checker build/milpcheck
python3 scripts/report.py results/my-miplib-pilot --runs revalidated.jsonl --prefix validated
~~~

For a subsequent validation version, use --name final-validation and report
--runs final-validation.jsonl. Previous records remain intact.

Each output directory must be new. Raw run artifacts and JSONL remain available
even when a solver times out or supplies no solution. Reports link to every log
and run record. Seed zero is rejected because some solvers interpret it as a
wall-clock seed. Benchmark runs are serial to avoid oversubscribing timed jobs.

The short-budget campaign is a pipeline baseline, not a solver ranking.
Longer campaigns can use --limit 300 and separate directories for seeds 1, 2,
and 3. Inspect all validation failures before accepting performance changes.

## Layout

- src/check.c: independent streaming MPS/solution checker using GMP rationals.
- src/bench.c: POSIX C process control, adapters, watchdog, provenance records.
- src/common.c: allocation, JSON escaping, monotonic timing and SHA-256.
- src/model.c: independent optimizer MPS reader and sparse CSC model.
- src/lp.c, src/factor.c, src/solve.c: LP engine, LU backend and command line.
- src/linear.c: guarded eta updates, equilibrated solves and compensated refinement.
- src/presolve.c: source-aware root bound propagation and redundant-row masks.
- src/search.c: deterministic reliability-based pseudocost branch ranking.
- src/exact.c, src/rational.c: original decimals and bounded exact reconstruction.
- src/pivot.c: deterministic ratio choices for validated pivot retries.
- src/bb.c: branch-and-bound, integer repair, limits and deterministic traces.
- config/: dataset locks, native assets, solver identities and smoke suite.
- tests/: exact arithmetic oracle, adversarial cases and runner contracts.
- data/, tools/, work/: locally fetched data/tools and build scratch; ignored.
- results/: raw experiments and generated HTML/JSON reports.

The new project code is MIT-licensed. Imported datasets, solver binaries,
GMP and upstream build dependencies retain their own licenses.

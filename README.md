# MILP Lab

A C17 laboratory for building and evaluating a new MILP solver.
Milestone 1 implements the dataset import, independent primal checker,
native benchmark runner and HiGHS / SCIP / Cbc baseline adapters.
It does not contain a new optimization solver yet.

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
- config/: dataset locks, native assets, solver identities and smoke suite.
- tests/: exact arithmetic oracle, adversarial cases and runner contracts.
- data/, tools/, work/: locally fetched data/tools and build scratch; ignored.
- results/: raw experiments and generated HTML/JSON reports.

The new project code is MIT-licensed. Imported datasets, solver binaries,
GMP and upstream build dependencies retain their own licenses.

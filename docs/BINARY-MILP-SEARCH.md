# Binary MILP search

The accelerated solver can translate an eligible MIP into Boolean clauses and solve it with Kissat 4.0.4. The normal MILP entry point tries this before the root relaxation. Eligibility follows the model coefficients and domains. The solver receives the original model and resource limits. The benchmark runner alone has access to reference objectives. [C17 implementation](../src/binary_sat.c) · [Kissat version, source hash and license](../config/kissat-source-v1.lock.json)

## Benchmark results

The C17 runner gave both solvers 60 seconds per problem on the same arm64 Mac, with one requested thread and seed 1. The seven-model panel checks the four structurally eligible models and the three previous completions. The independent checker used the original unpresolved MPS. Wall times include reading and process startup. [Run records](../results/binary-sat-v1/paired-60s/runs.jsonl)

| Model | Candidate | Seconds | HiGHS 1.15.1 | Seconds | Checked candidate objective |
|---|---|---:|---|---:|---:|
| bnatt400 | optimal | 0.873 | time_limit | 60.138 | 1 |
| bnatt500 | infeasible | 1.091 | time_limit | 60.058 | — |
| chromaticindex1024-7 | optimal | 10.826 | time_limit | 60.492 | 4 |
| chromaticindex512-7 | optimal | 3.633 | watchdog_timeout | 75.041 | 4 |
| enlight_hard | optimal | 0.030 | optimal | 9.950 | 37 |
| markshare_4_0 | optimal | 0.081 | time_limit | 60.029 | 1 |
| nw04 | optimal | 0.810 | optimal | 45.759 | 16862 |

The panel contains 6 checked, reference-matched candidate completions and 2 strictly faster completions where both solvers finished correctly. The candidate also reported bnatt500 infeasible. That case remains unresolved under the required independent-proof contract. The target remains 120 strict wins across the fixed 240 models. [Classifications](../results/binary-sat-v1/summary.json)

## Full inventory

At two seconds per model, the same release completed 4 of 240 correctly. The solver reported 1 infeasible, 18 numerical_failure, 4 optimal, 217 time_limit. The shorter inventory overlapped development checks and measures coverage only. [Gate record](../results/binary-sat-v1/gate.json) · [Full runs](../results/binary-sat-v1/full-2s/runs.jsonl)

## Exact translation and stopping

Every column must have an original integer domain of zero or one. The search supports up to sixteen nonzero objective coefficients. Source rationals determine the objective, including improvement cutoffs smaller than one floating-point unit. Positive penalties on signed literals handle either optimization direction. The objective offset cancels in comparisons.

For short inequalities with small dyadic coefficients, the encoder emits clauses that forbid overweight subsets. Long single clauses and at-most-one constraints have linear-size encodings. For short decimal rows the solver evaluates binary assignments with outward source intervals and consults exact source values when an interval overlaps a boundary. It removes a literal from a clause only if every newly covered assignment is also forbidden. An unsupported bound, uncertain assignment or exhausted encoding budget returns control to the ordinary search.

After finding a complete feasible assignment, the solver tries a strict objective cutoff. An unsatisfiable subset of the constraints suffices to rule out a cheaper assignment. It tries neighborhoods with uncommon patterns of row and variable counts first. Only the full constraint set can supply a feasible point. Incomplete searches can pass an incumbent to branch-and-bound after original-model validation. The CLI reports the number of SAT calls, including local cutoff checks, in `sat_calls`.

Kissat runs in a single child while the parent waits. The parent enforces a monotonic deadline during search and preprocessing and discards incomplete responses. Time spent reading the model counts toward the CLI budget. An explicit node limit also caps the SAT decision allowance. The default decision allowance across all calls is two million.

## Building and reproducing

The tested POSIX release uses the locked C sources for GLPK 5.0 and Kissat 4.0.4. GLPK supplies an optional simplex starting basis. Kissat retains its upstream MIT license, and GLPK retains GPL-3.0-or-later. The independent checker links GMP and zlib. The runner links zlib. The frozen release stores both third-party notices beside its executables. [GLPK source lock](../config/glpk-source-v1.lock.json) · [Build and source hashes](../results/binary-sat-v1/manifest.json)

From the repository root, with the two static libraries built, configure CMake with the source include paths and archive paths. These commands use the dependency locations in the measured workspace.

~~~sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DLAB_GLPK_BASIS=ON \
  -DGLPK_INCLUDE_DIR="$PWD/work/glpk-basis-v1/glpk-5.0/src" \
  -DGLPK_LIBRARY="$PWD/work/glpk-basis-v1/glpk-5.0/src/.libs/libglpk.a" \
  -DLAB_KISSAT=ON \
  -DKISSAT_INCLUDE_DIR="$PWD/work/binary-sat-v1/release-dep/kissat-rel-4.0.4/src" \
  -DKISSAT_LIBRARY="$PWD/work/binary-sat-v1/release-dep/kissat-rel-4.0.4/build/libkissat.a"
cmake --build build -j4
ctest --test-dir build --output-on-failure
~~~

On a fresh checkout, obtain the archives from the version locks, verify the SHA-256 checksums and extract them into an untracked build directory. Configure GLPK with `--disable-shared --enable-static --disable-dependency-tracking CFLAGS="-std=c17 -O3"` and run `make`. Configure Kissat with `--safe --quiet` and run `make CC="cc -std=c17 -D_DARWIN_C_SOURCE -O3 -DNDEBUG -DQUIET -DSAFE"` on macOS. The CMake cache variables above can also point to system-installed headers and libraries.

Use `--no-binary-sat` for an ablation on an accelerated build. Use `--binary-sat` to try the translator on smaller models. The default automatic size threshold is 64 columns and more than 64 rows. The committed regression suite includes exhaustive assignment oracles, signed costs and inequalities, non-dyadic source data, long cardinality rows, local relaxations, explicit budgets and checked original-model witnesses. [Contract tests](../tests/test_binary_sat.py) · [Release tests](../results/binary-sat-v1/release-final.log) · [ASan and UBSan tests](../results/binary-sat-v1/sanitizer-final.log) · [Build without external solvers](../results/binary-sat-v1/native-release.log)

Later releases apply a three-line pointer-arithmetic fix to the pinned Kissat source. The original sanitizer logs above used UBSan recovery. A strict run exposed the issue, and the patched library passed the same failing input. Use the v2 lock when building the current release. [Patch and source lock](../config/kissat-source-v2.lock.json) · [Current build and validation](NETWORKS-AND-PACKING.md)

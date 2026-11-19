# Integer search and sparse LP reconstruction

The combined C17 release recorded 9 strict wins in nine selected serial comparisons. The target is 120 on the fixed 240-model benchmark. Both solvers used one thread, seed 1 and a 90-second budget. Times include model input and solution output. A win required correct completion by both solvers, independently checked original-model witnesses and agreement with the locked references. [Build and source hashes](../results/discrete-search-v1/final/manifest.json) · [Measured pairs](../results/discrete-search-v1/final/paired-90s/runs.jsonl)

| Model | Candidate (s) | HiGHS (s) | Candidate status | HiGHS status |
|---|---:|---:|---|---|
| decomp2 | 0.094 | 4.145 | optimal | optimal |
| drayage-100-23 | 0.234 | 0.855 | optimal | optimal |
| enlight_hard | 0.023 | 10.521 | optimal | optimal |
| exp-1-500-5-5 | 0.028 | 4.128 | optimal | optimal |
| neos-3083819-nubu | 13.305 | 48.715 | optimal | optimal |
| neos-3381206-awhea | 0.022 | 1.476 | optimal | optimal |
| nw04 | 0.841 | 47.061 | optimal | optimal |
| p200x1188c | 0.030 | 0.729 | optimal | optimal |
| sp150x300d | 0.022 | 0.146 | optimal | optimal |

The candidate completed 13 models correctly in the full two-second inventory. 13 finished within two seconds measured from process startup to exit. The independent checker accepted 24 feasible witnesses across the 240 attempts, including 11 on unfinished models. The sweep used two concurrent processes. The recorded statuses were 209 time_limit, 13 optimal, 1 infeasible, 17 numerical_failure. The two selected longer panels produced 16 correct completions. [Inventory](../results/discrete-search-v1/final/gate.json) · [Score](../results/discrete-search-v1/final/summary.json)

## Bounded integer search

The solver recognizes zero-valued equalities of the form a·x − a·y = 0 and represents each connected group with one integer variable. It intersects member bounds, adds original costs and combines row coefficients. Scalar conversions, merges and row normalizations use exact integers. The method requires bounded integer variables, an integral objective and supported coefficients and row sides representable in thousandths. The branch-and-cut engine is SCIP 9.2.4 with the QSopt_ex double simplex, built from C17 sources. [Implementation](../src/integer_mip.c) · [Dependency lock](../config/scip-integer-v1.lock.json)

On neos-3083819-nubu the method removed 4,215 equality rows, reducing 8,644 variables to 4,429 and 4,725 rows to 510. A rational replay checked every coefficient, objective term, bound and component against the source. After search, the solver expands its integer vector, checks every original row and compares the integer cost with the numerical primal and dual bounds before accepting SCIP’s optimality claim. The independent checker evaluated the original-model witness at zero tolerance and obtained 6,307,996. The evidence establishes primal feasibility and reference agreement. No independent MILP tree proof was submitted. [Equivalence replay](../results/equality-projection-v1/equivalence.json) · [Original witness check](../results/discrete-search-v1/final/exact-nubu.json)

In development, SCIP solved the original formulation in 8.574 seconds. Two in-process runs on the reduced formulation took 15.067 and 13.650 seconds. Both produced the same exact original-model objective. With native search alone, the reduced model stopped on a numerical gap check after 11 nodes. A file-based experiment exercised an unsupported free-row operation in QSopt_ex during LP graph conflict analysis. The integrated configuration uses dual conflict analysis and disables strong-branching conflict analysis. [Original SCIP and HiGHS](../results/equality-projection-v1/original-60s/runs.jsonl) · [File-based failure](../results/equality-projection-v1/projected-1/solve.log) · [Follow-up runs](../results/equality-projection-v1/continued/plan.json)

## Sparse LP results

The exact sparse basis solver indexes rows and columns that have a single remaining coefficient. Each index enters its queue at most once. It now accepts up to 8,192 basis rows while keeping the existing 128-row limit on the dense remainder and the original bit, nonzero and time limits. The code checks every original equation before publishing the solution. [Rational solver](../src/rational.c)

In paired five-second LP runs, binkar10_1 and cost266-UUE changed from numerical failures to independently checked points and valid rational dual bounds. At ten seconds, cost266-UUE advanced from one failed root to 15 node LPs and exhausted its time budget. Both builds visited 247 nodes and found a checked incumbent on binkar10_1. Neither completed the MILP. [Root comparisons](../results/sparse-elimination-v1/root-comparison.jsonl) · [Sparse-only MILP timings](../results/sparse-elimination-v1/paired-60s/runs.jsonl)

## Build and verification

Build the optional C17 SCIP and QSopt_ex libraries from the two pinned archives using Python 3.12 or later. Supply a fresh output directory. The builder verifies hashes, applies the compatibility patches, checks the SCIP library’s compilation units, and retains the upstream licenses. The measured combination also enables GLPK 5.0 and Kissat 4.0.4. Supply their include and library paths as described in the preceding release. [Archive URLs and hashes](../config/scip-integer-v1.lock.json) · [GLPK and Kissat build notes](NETWORKS-AND-PACKING.md#reproducing-the-measured-build)

~~~sh
python3 scripts/build_scip_integer.py --scip-archive /path/to/scip.tar.gz --qsopt-archive /path/to/qsopt-ex.tar.gz --out work/scip-integer-deps
cmake -C work/scip-integer-deps/lab-cache.cmake -S . -B build-integer -DCMAKE_BUILD_TYPE=Release -DLAB_GLPK_BASIS=ON -DLAB_KISSAT=ON
cmake --build build-integer -j 4
ctest --test-dir build-integer --output-on-failure
~~~

All 60 accelerated release contracts, 60 accelerated ASan/UBSan contracts and 54 native-only contracts passed. Sanitizers covered the project’s integration code. The imported SCIP and QSopt_ex libraries used the same Release archives in both accelerated configurations. Tests cover permuted equality groups, signed objectives and offsets, sub-double coefficients, node budgets, and an 8,191-row nonsingular chain and its singular control. A former failure on 1,025 basis rows now has an exact point and bound. The regression also checks failure beyond the expanded limit. [Release tests](../results/discrete-search-v1/validation/release-contracts.log) · [Sanitizers](../results/discrete-search-v1/validation/san-contracts.log) · [Native-only tests](../results/discrete-search-v1/validation/native-contracts.log)

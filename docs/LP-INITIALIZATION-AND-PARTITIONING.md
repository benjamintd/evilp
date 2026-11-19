# LP initialization and set-partitioning search

The installed C17 solver completed **3 of 240** locked MIPLIB instances at two seconds per model. The preceding release completed 2. The 60-second comparisons cover those three MILPs. The full target remains 120 strict wins. [Benchmark gate](../results/partition-cover-v1/gate.json)

## Full-inventory coverage

| Two-second run | Correct completions | Time limits | Numerical failures | Watchdog stops | Checked feasible | Invalid witnesses | Timeouts at one node and one LP call |
|---|---:|---:|---:|---:|---:|---:|---:|
| Previous release | 2 | 227 | 11 | 0 | 13 | 0 | 166 |
| Current release | 3 | 219 | 18 | 0 | 15 | 0 | 113 |

The full inventory overlapped validation and development runs. Its timings describe coverage only. The native solver accepted at least one proposed basis on 111 models. The largest numerical failure groups were 7 unavailable original-source objective pricing decisions and 5 unconfirmed redundant rows. The corpus gate remains closed. [All runs](../results/partition-cover-v1/full-2s/runs.jsonl)

## Serial MILP comparisons

The C17 runner executed each pair on the same arm64 Mac with one thread and seed 1. Wall times include process startup. Both solvers received 60 seconds. The checker evaluated original, unpresolved MPS files.

| Model | Candidate status | Candidate wall seconds | HiGHS status | HiGHS wall seconds | Candidate checked objective |
|---|---|---:|---|---:|---:|
| enlight_hard | optimal | 0.028 | optimal | 9.857 | 37 |
| markshare_4_0 | optimal | 0.055 | time_limit | 60.025 | 1 |
| nw04 | optimal | 0.813 | optimal | 44.495 | 16862 |

Strict faster completions in these pairs were enlight_hard, nw04. The benchmark requires correct completion by both solvers for a strict win. An optimal status is the solver’s claim. Feasibility is checked independently and reference agreement is recorded separately. [60-second pairs](../results/partition-cover-v1/paired-mip-60s/runs.jsonl)

A separate ten-second comparison used the four LP-development models in MILP mode.

| Model | Candidate | HiGHS |
|---|---|---|
| air05 | time_limit, 10.026 s | time_limit, 10.093 s |
| nw04 | optimal, 0.805 s | time_limit, 10.519 s |
| ran14x18-disj-8 | time_limit, 10.026 s | time_limit, 10.036 s |
| sp150x300d | time_limit, 10.008 s | optimal, 0.149 s |

[Ten-second MILP runs](../results/partition-cover-v1/paired-mip-10s/runs.jsonl)

## Serial LP comparisons

Three runs per solver and model, ten seconds per run, with rotating order. Values are median solver-process CPU seconds. A timeout is labeled even when the measured CPU lies below the wall limit.

| Original LP relaxation | Previous release | Current release | HiGHS 1.15.1 | Current checked points and lower bounds |
|---|---:|---:|---:|---:|
| air05 | 9.971 (timeout) | 0.737 | 0.315 | 3/3 |
| nw04 | 1.321 | 0.556 | 1.356 | 3/3 |
| ran14x18-disj-8 | 2.201 | 0.075 | 0.054 | 3/3 |
| sp150x300d | 0.061 | 0.016 | 0.008 | 3/3 |

The checker evaluated primal points and dual bounds independently. All solver outcomes and checker records are retained. [Paired LP runs](../results/partition-cover-v1/paired-lp/records.json)

## Implementation and validation

The optional GNU GLPK 5.0 C library proposes a starting simplex basis on the row-bounded, variable-bounded relaxation. The native engine lifts the indices into its own equality form, refactors the basis, checks primal feasibility, prices columns and validates the original model. A rejected proposal falls back to ordinary simplex within the remaining budget. The checker and runner have no GLPK linkage. GLPK retains its GPL-3.0-or-later license. [Dependency lock](../config/glpk-source-v1.lock.json) · [Linkage check](../results/partition-cover-v1/linkage.json)

A native C17 set-partitioning search handles up to 64 unit equality rows, binary columns and source-exact integer costs. It removes dominated identical covers, repairs the proposed row prices to exact dual feasibility in integer arithmetic, and searches disjoint covers by increasing cost. The node budget includes the cover search. A partial search may supply a checked incumbent. Only exhaustion of every cheaper cover can close the search with a witness. No independent MILP proof file is emitted. [Implementation](../src/partition.c)

The interval support check now scans up to four million nonzeros with bounded memory and an independent deadline. It keeps the absolute gap threshold and checks every source enclosure, including entries whose contribution is zero. The LPs of bnatt400, drayage-25-23 and fastxgemm-n2r6s0t2 passed independent point and bound checks after this change. The source check on app1-1 still abstained. [Diagnostic LPs](../results/glpk-basis-v1/support-panel/runs)

All 50 suites passed in the accelerated release and the ASan/UBSan build. All 49 applicable suites passed in a release without GLPK. The tests include an independent subset oracle for the cover search, original-decimal counterexamples, budget exhaustion and a 136,000-entry support check with late malformed data. The sparse-update test checks the default optimum and every original variable, then compares update and refactor work from a common native starting basis. Its first run exposed that an imported basis could already contain the 64 active bounds. [Release tests](../results/partition-cover-v1/release-verified.log) · [Sanitizers](../results/partition-cover-v1/sanitizer.log) · [Native-only release](../results/partition-cover-v1/native-release.log) · [Source and binary hashes](../results/partition-cover-v1/verification-manifest.json)

## Running the accelerated build

The workspace executable is `build/milpsolve`. Configure another checkout with an installed GLPK development library using `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLAB_GLPK_BASIS=ON`, then run `cmake --build build`. Standard CMake builds keep GLPK optional. In an enabled build the CLI chooses it automatically for models with at least 100 rows and 100 columns, or at least one row and 1,024 columns. Explicit pricing policies and basis files retain their existing startup path. Use `--glpk-basis` or `--no-glpk-basis` to select the initializer explicitly. A supplied input basis takes precedence.

The fixed inventory still includes supportcase22. MIPLIB lists its incumbent objective of 110 and its optimization status as open. The locked best-known objective therefore remains unresolved for the completion gate. [MIPLIB instance record](https://miplib.zib.de/instance_details_supportcase22.html)

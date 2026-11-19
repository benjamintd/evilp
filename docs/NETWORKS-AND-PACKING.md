# MILP performance on network, routing and packing models

The combined C17 release adds three ways for the normal MILP command to finish a problem before an expensive branch-and-bound search. All use the original model, including its objective sense, row sides, bounds and decimal coefficients. Each path has a time and work budget. The solver selects these methods from coefficients and bounds. The runner compares their output with the locked reference objective. [Frozen release manifest](../results/compact-lp-v2/manifest.json) · [Experiment plan](../results/compact-lp-v2/plan.txt)

## Matched MILP runs

The C17 runner timed each solver process from startup to exit on the same arm64 Mac. Both received one thread, seed 1 and 60 seconds, with a 15-second watchdog grace. The checker evaluated the original unpresolved MPS after timing. Both solvers finished correctly and the candidate was faster on 4 of these 6 selected pairs. The corpus contains 240 problems and the target is 120 strict wins. [Raw pairs](../results/compact-lp-v2/paired-60s/runs.jsonl) · [Classifications](../results/compact-lp-v2/summary.json)

| Model | Candidate seconds | HiGHS seconds | HiGHS outcome | Candidate checked objective |
|---|---:|---:|---|---:|
| eil33-2 | 0.167 | 60.063 | time_limit | 233501979/250000 |
| enlight_hard | 0.025 | 12.218 | optimal | 37 |
| nw04 | 0.807 | 51.230 | optimal | 16862 |
| p200x1188c | 0.047 | 0.782 | optimal | 15078 |
| sp150x300d | 0.025 | 0.149 | optimal | 69 |
| supportcase18 | 0.044 | 60.115 | time_limit | 48 |

The same binary ran four additional cases with a 60-second budget. Across the two panels, 10 candidate MILPs finished correctly with an independently checked feasible solution at the locked reference objective. Only the six freshly paired cases enter the table above. Each MILP optimal status is the solver’s claim. [Additional runs](../results/compact-lp-v2/boolean-60s/runs.jsonl)

The two-second candidate-only coverage inventory completed 7 of 240 correctly. It used four concurrent lanes during validation and is unsuitable for comparing wall times. Reported statuses were 215 time_limit, 1 infeasible, 16 numerical_failure, 7 optimal, 1 watchdog_timeout. Checker outcomes were 221 unavailable, 19 valid. The completion and speed gates remain closed. [Inventory](../results/compact-lp-v2/full-2s/runs.jsonl) · [Full gate](../results/compact-lp-v2/gate.json)

## Standalone LP measurements

Medians of three solver-process wall times in seconds. Solver order rotated between repetitions. Each original LP had a ten-second limit. The previous build and this build used the same GLPK and Kissat versions. [Per-run times, point checks and dual checks](../results/compact-lp-v2/paired-lp/records.json)

| LP relaxation | Previous C17 | Current C17 | HiGHS | Current point and bound valid |
|---|---:|---:|---:|---:|
| air05 | 0.781 | 0.293 | 0.347 | 3/3 points, 3/3 bounds |
| nw04 | 0.594 | 0.632 | 1.471 | 3/3 points, 3/3 bounds |
| ran14x18-disj-8 | 0.098 | 0.073 | 0.071 | 3/3 points, 3/3 bounds |
| sp150x300d | 0.024 | 0.023 | 0.025 | 3/3 points, 3/3 bounds |

On sp150x300d the median solver CPU times were 0.016 seconds for the candidate and 0.009 for HiGHS. Process launch and scheduling are appreciable at these wall times.

All 55 accelerated contract suites passed in Release and in an AddressSanitizer build with immediate UndefinedBehaviorSanitizer termination. The scalability and source-pricing contracts exercise the direct check and the previous basis path on the same inputs. Their first runs flagged work counters specific to the basis path. Those original assertions remain in the expanded tests. The same sources passed all 52 suites available without either optional library. [Release](../results/compact-lp-v2/release-verified.log) · [Strict sanitizers](../results/compact-lp-v2/sanitizer.log) · [Native-only tests](../results/compact-lp-v2/native-release.log)

## Fixed-charge networks

The solver recognizes continuous directed flow with integer node balances and binary arc-opening decisions. Each arc must have its own gate, a nonnegative opening cost after normalization, a finite integer capacity, and zero per-unit flow cost. Fixing a gate assignment gives a bounded network with integer capacities and demands, so a feasible continuous flow has an integral representative.

For a vertex of degree two, the solver eliminates the balance equation and combines the costs on its incident arcs. It represents the minimum opening cost for each possible transfer by constant-cost intervals. Parallel arcs combine by minimizing the sum of the two costs over compatible transfers. Interval endpoints determine the representation size even for large capacities. Replaying the elimination choices reconstructs an original flow and gate assignment.

If eliminations stop and the original network has one supply, at most ten demands and sufficient capacity on each usable arc to carry the whole supply, the solver can instead optimize a directed tree connecting the terminals. For each terminal subset it tries splitting the subset at a vertex, then propagates distances along reversed arcs. A symmetric reversal handles a single demand and several supplies. It recovers the original flow and requires its opening cost to equal the computed minimum. [Network implementation](../src/flow.c) · [Enumeration and independent max-flow tests](../tests/test_flow.py)

## Decimal set-partitioning costs

The set-partitioning search now obtains a common denominator from original rational objective coefficients. If conversion fits its integer arithmetic bounds, it searches with scaled integer costs and computes its pruning bounds in the same units. The search supports signed costs and either optimization direction. Limits on the denominator and coefficient magnitude cause the search to defer to the general solver. On eil33-2 the recovered objective was exactly 233501979/250000. [Search](../src/partition.c) · [Fractional subset oracles](../tests/test_decimal_partition.py) · [Original-model check](../results/decimal-cover-v1/exact-witness.json)

## Capacity lower bounds for resource packing

This path recognizes unit item-assignment equalities and resource-capacity constraints with a common charge for opening a resource. It computes a candidate weight for each item from the cheapest singleton assignment in units of resource consumption. It verifies, for every original assignment column, that its total resource consumption covers the sum of the weights for all items it serves. This also checks bundled options that serve several items at once.

The ceiling of total item weight divided by maximum resource capacity is then a lower bound on the number of paid resources. A deterministic search builds an assignment using that many resources, moving or exchanging items to remove overloads. Only a feasible assignment meeting the bound can close the problem. On supportcase18 the verified item weights total 7,078 and each bin holds 150. The checker accepted an assignment using 48 bins at zero tolerance. [Implementation](../src/packing.c) · [Small assignment oracles and numerical counterexamples](../tests/test_packing.py) · [Original-model check](../results/resource-packing-v1/witness.json)

## Standalone LP solves

The optional GLPK initializer already computes a candidate primal point and dual prices. The default standalone LP path maps them through all row and variable transformations and checks the original feasible point and an outward-rounded absolute support gap of at most 1e-7. When both checks pass, the caller can return the point and dual immediately. If the check is uncertain, the proposed basis remains available to native simplex. Session solvers still retain their basis for subsequent bound changes. Explicit basis requests keep the basis path.

Disconnected columns with an exactly zero objective and a source-feasible zero lower bound can be fixed to zero for the standalone solve. This preserves the objective and every original dual support bound. It also lets the source-pricing fallback resolve a generated large LP with a sub-double cost difference and a material 0.01 objective gap. The original bound-only solve remains an unresolved negative control on that input. [Proposal verification](../src/glpk_basis.inc) · [Exact objective and signed-transform tests](../tests/test_standalone_proposal.py)

## Transfer to the MILP tree

The standalone LP shortcut returns a checked point and dual, while search nodes can also use a simplex basis to warm-start their descendants and derive cuts. An isolated ten-model test passed the checked proposal straight to the MILP root without retaining that basis. At three seconds per model the control processed 1,046 nodes and returned one checked incumbent. The experimental path processed 896 nodes and returned two checked incumbents, with one numerical stop. Neither build completed a model on this panel. The measured release keeps the standalone shortcut out of the MILP tree. [Paired development records](../results/mip-root-transfer-v1/records.jsonl)

An earlier experiment retained an original-row simplex problem across changes to branch bounds. At five seconds on another ten-model panel, it processed 11,829 nodes versus 3,161 for its control and returned six checked incumbents versus four. Both versions recorded nine time limits and one numerical failure. The node counts show how much work changed while a finished-MILP count remained at zero. [Node, incumbent and status records](../results/mip-reoptimization-v1/mip-panel/records.jsonl)

## Reproducing the measured build

Use the GLPK 5.0 and Kissat 4.0.4 archives and SHA-256 values in the dependency locks. Before building Kissat, apply the patch listed in its v2 lock with `patch -p1 -i /absolute/path/to/kissat-4.0.4-empty-original.patch` from the extracted source directory. The patch avoids adding a zero offset to an unallocated original-clause buffer. The v1 sanitizer suites had returned success with UBSan recovery enabled. A later run with immediate UBSan termination exposed this operation. The subsequent patched release passed the original failing input with immediate termination enabled. [Kissat v2 lock and patch digest](../config/kissat-source-v2.lock.json) · [GLPK lock](../config/glpk-source-v1.lock.json) · [Patched sanitizer run](../results/flow-network-v2/sanitized-bnatt.json)

Build the static libraries as described in the [binary search build notes](BINARY-MILP-SEARCH.md). Point the `KISSAT_INCLUDE_DIR` and `KISSAT_LIBRARY` cache entries at the patched build. The measured executable is `work/compact-lp-v2/validated/milpsolve`. The GMP-based independent checker and the C17 runner remain separate executables. Imported code retains its upstream licenses. [Build and source hashes](../results/compact-lp-v2/manifest.json)

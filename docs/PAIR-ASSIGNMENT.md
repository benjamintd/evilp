# Pair assignment search and MIPLIB comparison

## Serial comparison at 60 seconds

| Model | Lab outcome | Lab wall seconds | HiGHS outcome | HiGHS wall seconds |
|---|---|---:|---|---:|
| qap10 | Optimal, 340 | 0.139 | Time limit, incumbent 348 | 60.276 |
| drayage-100-23 | Optimal | 0.227 | Optimal | 0.793 |
| p200x1188c | Optimal | 0.025 | Optimal | 0.731 |
| neos17 | Optimal | 11.503 | Optimal | 6.010 |
| neos859080 | Infeasible | 2.298 | Infeasible | 0.969 |
| drayage-25-23 | Time limit, reference-matching incumbent | 57.767 | Optimal | 1.899 |

Each solver completed five of the six selected cases. Four completed in both engines and each engine was faster on two. The reference-status contract accepts the infeasibility report on `neos859080`. All optimal outcomes had checker-valid original witnesses matching their locked objectives. This selected panel gives no corpus-wide speed ranking. [Serial records](../results/pair-assignment-v1/serial-60s/runs.jsonl) · [Case classifications](../results/pair-assignment-v1/evaluator-v2/panel.json)

On `qap10`, the 4,150 original columns encode 100 primary choices and 4,050 auxiliary products. A separate run of the pre-change v10 executable timed out at the same 60-second limit with a checked incumbent of 362. The new search reconstructed a witness at exactly 340. The independent checker accepted it at zero tolerance for constraints, integrality, bounds and objective consistency. Optimality is the solver's search outcome. [Pre-change run](../results/pair-assignment-v1/control/qap10/lab/seed-1/run.json) · [Zero-tolerance check](../results/pair-assignment-v1/qap10-zero-tolerance.json)

## Two-second coverage across all 240 MIPLIB instances

| Outcome | Lab | HiGHS 1.15.1 |
|---|---:|---:|
| Correctly completed with a reference-matching checked witness | 18 | 3 |
| Independently checked feasible witnesses, including limits | 104 | 98 |
| Solver time limits | 222 | 208 |
| Input time limits | 0 | 26 |
| Watchdog timeouts | 0 | 3 |

The three HiGHS completions also completed in Lab. Of the 222 unresolved Lab runs, 86 produced a checked incumbent and 136 finished without one. Every submitted Lab witness passed the checker. The completion and speed gates remain closed. [All classifications](../results/pair-assignment-v1/coverage-2s/gate.json) · [Counts and solved instances](../results/pair-assignment-v1/coverage-2s/counts.json)

## Solver changes

The root search recognizes exact-cover equations with auxiliary product variables. Each primary integer column belongs to two cover equations with right-hand side one. The remaining equations must express that a sum of nonnegative auxiliary variables equals one primary variable. Each auxiliary variable has two distinct primary owners. All partners in a linking equation must share a cover equation and each pair must have a unique auxiliary column.

An auxiliary variable can be positive only when both owners are selected. When both are selected, their linking equation has exactly one eligible auxiliary variable, which must equal one. This reduces the search to picking a set of primary columns that covers every equation exactly once. The implementation verifies the coefficients, bounds and integral costs against their source rationals. Signed costs and both optimization senses are supported. Every complete cover selects the same number of primary columns and pairs, so shifting their costs to be nonnegative gives a lower bound at every partial assignment.

A missing link can admit extra covers in the relaxation. The solver checks every original equality on the reconstructed point before accepting the optimum. Exceeding the node or time budget returns control to the existing search. The algorithm has no dependency on a model name, a reference objective or a solver library. [C17 implementation](../src/pair_assignment.c) · [Permutation oracle and adversarial cases](../tests/test_pair_assignment.py)

The MPS reader reuses the current column index, recognizes literal unit coefficients and sorts rows within grouped columns only when necessary. Interleaved columns and duplicate coefficients retain the full sorting path and its summation order. The interval arithmetic used in presolve exploits monotonicity when multiplying by a singleton. Multiplication and division by an exact signed unit preserve endpoints, including subnormals, after the usual rounding-mode and validity checks. [MPS reader](../src/model.c) · [Interval arithmetic](../src/exact.c) · [Exact regression cases](../tests/test_exact_model.c)

## Reader and presolve timings

Median wall seconds over three serial repetitions on the same machine, alternating executable order. The reader measurement includes parsing and building the source arithmetic representation. The presolve measurement starts after model reading and runs one complete propagation pass on the LP relaxation. All three presolve pairs made the same number of tightenings within each model.

| Operation | Model | Before | After | Time reduction |
|---|---|---:|---:|---:|
| Read MPS | square41 | 5.539 | 3.527 | 36% |
| Read MPS | neos-5114902-kasavu | 3.175 | 2.266 | 29% |
| Read MPS | supportcase19 | 3.899 | 3.122 | 20% |
| Presolve pass | highschool1-aigio | 0.306 | 0.208 | 32% |
| Presolve pass | neos-827175 | 0.0198 | 0.0126 | 36% |
| Presolve pass | sorrell3 | 0.0582 | 0.0459 | 21% |

These measurements isolate two costs incurred before branch-and-bound. The earlier 240-instance profile identified repeated node LP solves as the largest cost over the full search. [Reader repetitions](../results/pair-assignment-v1/reader-pairs.jsonl) · [Presolve repetitions](../results/pair-assignment-v1/presolve-pairs.jsonl) · [SCIP component profile](SCIP-V9-TOLERANCES-AND-BUDGETS.md)

## Benchmark protocol

The accelerated release uses the pre-existing v10 SCIP/QSopt_ex configuration, GLPK and Kissat. HiGHS is the 1.15.1 binary in the solver lock. The C17 process runner requests one thread, seed 1 and zero MIP gaps. Its wall time includes process startup, model reading and solution writing. The independent checker reads the original unpresolved MPS after timing.

The coverage experiment uses four concurrent lanes, alternating solver order by model within each lane. Each solver receives two seconds and a 15-second watchdog grace. Contract tests overlapped some runs. The serial comparison uses 60 seconds per solver, the same grace and a separate output directory. Limits remain unresolved and strict paired timing wins require correct completion by both engines. [Dataset lock](../config/miplib.lock.json) · [Gate evaluator](../scripts/benchmark_gate.py)

For infeasible outcomes, the evaluator now recognizes the checker's explicit `no_solution` response (exit code 4) to a header-only solution file under the frozen reference-status contract. The original evaluator treated this as a checker process error for HiGHS on `neos859080`, while accepting the same status from a solver that emitted no file. Invalid witnesses, parser errors, missing required proofs and time limits continue to fail. Reclassifying the existing records with the corrected evaluator leaves both two-second completion counts unchanged. [Evaluator regressions](../tests/test_benchmark_gate.py) · [Reclassified coverage gate](../results/pair-assignment-v1/evaluator-v2/coverage-gate.json)

## Validation

The accelerated Release and ASan/UBSan builds passed all 63 suites. The native-only Release passed its 55 applicable suites. The sanitizer build instrumented the project's C17 sources and linked the existing release archives for optional engines. The product search tests cover randomized row and column orders, signed costs, both objective senses, continuous auxiliaries, the full mask width, node exhaustion and malformed product encodings. The updated evaluator passed all 15 scoring regressions. [Release contracts](../results/pair-assignment-v1/release-all.log) · [Sanitizers](../results/pair-assignment-v1/contracts-san-final.log) · [Additional sanitizer cases](../results/pair-assignment-v1/san-regressions.log) · [Native-only contracts](../results/pair-assignment-v1/contracts-native.log) · [Experiment identities](../results/pair-assignment-v1/manifest.json)

## Reference discrepancy

MIPLIB lists 6741.38 for the submitted solution of `binkar10_1`, with an integrality violation of 8e-7 and an “Exact” value of 6742.2. It also lists a solution at 6742.2 with zero recorded violations. The earlier Lab and HiGHS optimality claims at 6742.200024 therefore remain disagreements with the locked v37 objective. The published values do not independently establish MILP optimality. [MIPLIB solution table](https://miplib.zib.de/instance_details_binkar10_1.html) · [Historical screen](SCIP-V9-TOLERANCES-AND-BUDGETS.md)

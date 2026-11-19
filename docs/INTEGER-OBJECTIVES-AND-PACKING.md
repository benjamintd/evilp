# Integer objectives and repeated-item packing

One thread and seed 1, 60 seconds per solver process. The runner measures process startup through exit and the independent checker evaluates the original model afterward. A strict win requires both solvers to complete correctly and the candidate to use less wall time. All pairs here used the same frozen binary. [Experiment manifest](../results/packing-pattern-v1/manifest.json) · [Raw pairs](../results/packing-pattern-v1/paired-60s/runs.jsonl)

| Locked model | Candidate seconds | HiGHS seconds | Checked objective |
|---|---:|---:|---:|
| enlight_hard | 0.030 | 10.286 | 37 |
| neos-3381206-awhea | 0.025 | 1.456 | 453 |
| nw04 | 0.812 | 46.746 | 16862 |
| p200x1188c | 0.057 | 0.723 | 15078 |
| sp150x300d | 0.021 | 0.161 | 69 |

The release recorded 5 strict wins toward the 120-instance target. Across these pairs and six additional candidate runs, 11 selected models completed with a checked feasible witness and reference-objective agreement. The MILP completion status remains a solver claim. [Classifications](../results/packing-pattern-v1/summary.json) · [Additional candidate runs](../results/packing-pattern-v1/boolean-60s/runs.jsonl)

The two-second full inventory completed 8 of 240. The four coverage lanes overlapped validation. These runs measure coverage only. The fixed inventory and all failures are in the [gate record](../results/packing-pattern-v1/gate.json).

## Repeated items in interchangeable bins

The solver checks source coefficients, item demands, upper bounds, and the opening cost of every bin before it groups them. It enumerates the feasible nonempty contents of one bin. GLPK proposes prices for the item types. The solver converts those prices into an integer inequality and checks it against every enumerated pattern. It then asks for integer pattern counts at the resulting lower bound, expands them into the original variables, and checks the witness. A time limit, incomplete enumeration, or a mismatch returns control to ordinary search. The GLPK dependency uses the existing version lock.

On neos-3381206-awhea the full model has 475 interchangeable bins, four item types and 2,375 integer variables. The solver found a packing using 453 bins and matched its pattern bound. The original-model checker accepted the solution under exact arithmetic. [Exact checker output](../results/packing-pattern-v1/dev-check.json)

## Pruning between attainable objectives

When every nonzero objective coefficient belongs to an integer variable, the solver computes a guaranteed objective increment using exact source rationals. Once it has an incumbent, a node whose lower bound excludes a whole increment of improvement can close. Unknown cost provenance and nonzero continuous costs disable the rule. The CLI switch `--no-objective-lattice` reproduces search without this pruning.

For a generated 24-variable case, the LP objective was 2.875 and every integer objective was a multiple of 0.25. The previous binary timed out after 30,625 nodes with objective 3. The new rule closed the same model at the root. This generated case contributes zero benchmark wins. [Previous run](../results/objective-grid-v1/before.json) · [Changed run](../results/objective-grid-v1/after.json) · [Original-model check](../results/objective-grid-v1/after-check.json)

## Validation

All 56 Release suites and all 56 AddressSanitizer and UndefinedBehaviorSanitizer suites passed. The native-only build passed all 52 applicable suites. Tests include exact enumeration, both objective senses, nonuniform bin costs and coefficients below double precision, continuous costs, and exhausted search budgets. Branch-order, node-limit, diving and cut tests retain their historical assertions with the objective-lattice optimization disabled where necessary to observe branching. [Release](../results/packing-pattern-v1/release.log) · [Sanitizers](../results/packing-pattern-v1/sanitizer.log) · [Native build](../results/packing-pattern-v1/native-release.log)

## Retaining LP bases during branch search

On neos5 the development variant reused a basis on 3,681 of 3,683 LP calls at four seconds. It still hit the 60-second limit, as did the shipped control. Across neos5, beasleyC3 and mad, neither version recorded a correct completion at 60 seconds. The retained-basis variant remains experimental. [Variant](../results/basis-reuse-v1/retained-long/records.jsonl) · [Control](../results/basis-reuse-v1/shipped-long/records.jsonl)

## Routing experiment

A separate prototype ruled out conflicting route assignments on drayage-100-23, computed an exact minimum-cost assignment bound, and recovered start times in source-decimal arithmetic. The original-model checker accepted its witness with zero tolerance. The ordinary release has no routing path for this case yet, so it contributes zero to the reported completion and speed counts. On drayage-25-23 the same search found no feasible assignment at its lower bound in 1,000 attempts. [Development evidence and limitations](../results/assignment-probe-v1/development.json)

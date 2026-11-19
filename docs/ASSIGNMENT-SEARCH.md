# Assignment search with continuous side constraints

The combined C17 release recorded 6 strict paired wins toward the 120-of-240 target. Each solver had one thread, seed 1 and 60 seconds on the same Mac. The runner timed startup through exit. Both solvers had to finish and provide a feasible witness matching the locked objective for a strict win. [Frozen inputs](../results/assignment-release-v2/manifest.json) · [Matched runs](../results/assignment-release-v2/paired-60s/runs.jsonl)

| Model | Candidate seconds | HiGHS seconds | Checked objective |
|---|---:|---:|---:|
| drayage-100-23 | 0.246 | 0.807 | 51666937036514277/500000000000 |
| enlight_hard | 0.024 | 9.983 | 37 |
| neos-3381206-awhea | 0.027 | 1.471 | 453 |
| nw04 | 0.829 | 46.151 | 16862 |
| p200x1188c | 0.028 | 0.738 | 15078 |
| sp150x300d | 0.029 | 0.137 | 69 |

Across both panels the candidate finished 12 selected MILPs correctly at 60 seconds. The separate two-second coverage inventory finished 8 of the 240 locked models. Its other cases included 215 time limits and 17 numerical failures. Coverage used four concurrent lanes. [Additional models](../results/assignment-release-v2/boolean-60s/runs.jsonl) · [Full inventory](../results/assignment-release-v2/gate.json)

## Routing result

On drayage-100-23, 105 binary assignment equations on each side define a matching. Conditional propagation removes a binary value when setting it to one forces its peers to zero and the resulting row activity cannot reach a source bound. The standalone development run fixed 3,341 columns before the matching search.

For a proposed matching, the solver checks exact dual inequalities for all available edges using rational source costs. Difference constraints describe the remaining continuous potentials. Conflicting cycles guide tie breaking. A matching closes the problem only when the exact dual test succeeds and source-feasible potential coordinates can be reconstructed. The output preserves the decimal coordinates so the independent checker can replay the original rows exactly. [Independent zero-tolerance check](../results/assignment-release-v1/dev-exact.json)

The same method leaves drayage-25-23 unresolved. Its assignment relaxation admits a cheaper objective than the locked integer optimum. The solver returns to normal branch and bound when the matching search fails to close the gap.

## Verification

All 57 accelerated Release contracts and all 57 AddressSanitizer and UndefinedBehaviorSanitizer contracts passed. The native-only build passed all 53 applicable contracts. A repeated sparse coefficient regression preserves the unique feasible point in clique probing. [Previous frozen results](../results/assignment-release-v1/summary.json) · [Reproduced regression](../results/assignment-release-v2/duplicate-before.log). New cases compare against exhaustive permutations and test signed costs, differences below binary64 precision, conflicting potential cycles, nonunit coefficients, and decimal witness reconstruction. A checked witness and reference-objective agreement establish the benchmark classifications. The checker does not certify a complete MILP proof. [Release](../results/assignment-release-v2/release.log) · [Sanitizers](../results/assignment-release-v2/sanitizer.log) · [Native build](../results/assignment-release-v2/native-release.log)

## Follow-up experiments

On a separate panel of nine LP relaxations, the QSopt_ex double simplex library returned independently checked feasible points and valid bounds on all nine. The C17 SCIP 9.2.4 build using this LP library then found checked incumbents on 12 of 13 selected MILPs at eight seconds. All 13 timed out. HiGHS completed one. Neither test contributes a production speed win. [Development build](../results/c-simplex-stack-v1/manifest.json) · [Panel results](../results/c-simplex-stack-v1/summary.json)

A development-only reduction grouped interchangeable assignment vertices on drayage-25-23. It shrank the model from 11,090 to 1,106 columns. SCIP reported optimality on the floating-point reduced model in 12.35 seconds. Expanding the group counts and reconstructing the original continuous values produced a witness with zero violation at zero tolerance and an objective matching the locked reference. The original-model checker established feasibility. Source-decimal equivalence of the reduction and transfer of the reported lower bound remain unverified. [Reduction and witness](../results/assignment-symmetry-v2/summary.json)

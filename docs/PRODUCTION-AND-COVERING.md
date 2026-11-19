# Production planning and binary covering

The combined C17 solver recorded **8 strict paired wins** toward 120 on the 240 locked models. Each timed pair used the same limit, one requested thread and seed 1. Each solver had to finish correctly with a checked original-model witness and agreement with the locked reference. The runner measured elapsed time from process startup through exit. [Build manifest](../results/combined-cover-v1/manifest.json) · [Scored results](../results/combined-cover-v1/summary.json)

| Instance | Limit (s) | Candidate (s) | HiGHS (s) | Candidate status | HiGHS status |
|---|---:|---:|---:|---|---|
| decomp2 | 60 | 0.099 | 4.078 | optimal | optimal |
| drayage-100-23 | 60 | 0.236 | 0.808 | optimal | optimal |
| enlight_hard | 60 | 0.024 | 10.299 | optimal | optimal |
| exp-1-500-5-5 | 60 | 0.024 | 4.323 | optimal | optimal |
| neos-3381206-awhea | 60 | 0.025 | 1.450 | optimal | optimal |
| neos-787933 | 180 | 0.540 | 181.034 | optimal | time_limit |
| nw04 | 60 | 0.819 | 47.386 | optimal | optimal |
| p200x1188c | 60 | 0.025 | 0.700 | optimal | optimal |
| sp150x300d | 60 | 0.025 | 0.150 | optimal | optimal |

The two-second sweep produced 13 correct completions across all 240 models. Recorded statuses were 209 time_limit, 13 optimal, 1 infeasible, 17 numerical_failure. Two coverage processes ran concurrently. The selected longer runs produced 15 correct candidate completions. [Full inventory](../results/combined-cover-v1/awake-repeat/gate.json) · [60-second pairs](../results/combined-cover-v1/paired-60s/runs.jsonl) · [180-second pair](../results/combined-cover-v1/awake-repeat/extended-180s/runs.jsonl) · [Other candidate runs](../results/combined-cover-v1/awake-repeat/other-60s/runs.jsonl)

The previous build recorded 8 correct completions at the two-second limit with four concurrent processes. The new sweep used two. Of its 13 correct completions, 13 finished within two seconds measured from process startup to exit. The serial matched pairs retained 6 wins and added decomp2, exp-1-500-5-5. On neos-787933 the candidate reached the reference objective 30 in 0.540 s. HiGHS returned a checked incumbent of 31 when its 180-second budget expired. The strict score requires both solvers to complete. [Preceding release](ASSIGNMENT-SEARCH.md)

## Solver changes

On binary models with many objective terms or fixed columns, the solver can fix each priced variable at the cheapest endpoint of its safely propagated bounds. Kissat searches the remaining variables. The search may omit inequalities it cannot cheaply encode. Before accepting the result, the solver checks both sides of every source row in scaled integer arithmetic. A feasible point on this face attains the objective bound from the box. This closes decomp2 at −160. [Implementation](../src/binary_sat.c) · [Zero-tolerance original-model check](../results/general-search-v1/normal-exact.json)

The production-planning method recognizes paths with nonnegative directional transport costs, positive integer demand, sufficient supply and arc capacity, and binary setup decisions. For a fixed set of open sites, ties in delivered unit cost can be resolved so each supplier serves a contiguous interval of demands. A dynamic program prices each path independently. If two chosen sites compete for an exclusive opening, the solver branches on forbidding either one. Exhausting these branches and reconstructing the original shipments closes exp-1-500-5-5 at 65,887. [Implementation](../src/lot_sizing.c) · [Zero-tolerance original-model check](../results/line-network-v1/checked.json)

On eligible binary covering problems, opening a gate permits every member of a group to take the value one. The solver checks the gate capacity, the absence of auxiliary costs, and nonnegative coverage coefficients. It then substitutes the gate for its members. GLPK proposes row prices and a point on the resulting relaxation. Integer arithmetic checks the reduced-cost contributions of all column bounds, coverage of each demand, and the integer-objective stopping condition. The solution to neos-787933 expands from 1,764 activation decisions to 236,376 original variables and has checked objective 30. [Implementation](../src/activation_cover.c) · [Zero-tolerance original-model check](../results/activation-cover-v1/normal-exact.json)

The independent checker verifies primal feasibility on the original MPS. These measurements rely on solver-reported optimality and reference-objective agreement. No independent MILP optimality proof was submitted.

## Verification and experiments

All 59 accelerated Release contracts, 59 ASan/UBSan contracts and 54 native-only contracts passed. Tests exercise original-model witnesses, exhaustive opening patterns, directed transport, signed objectives, decimal coefficients, insufficient capacity, conflicting openings, omitted Boolean constraints, fractional relaxations and exhausted budgets. Earlier overlapping runs recorded watchdog and test deadline failures. The final validation ran each configuration in sequence. The first 180-second pair and several inventory processes recorded long workstation pauses. The subsequent run repeated the extended pair, all six additional candidates and the entire inventory with an idle-sleep assertion. [Interrupted inventory](../results/combined-cover-v1/gate.json) · [Repeat plan](../results/combined-cover-v1/awake-repeat/plan.json) · [Release tests](../results/combined-cover-v1/release.log) · [Sanitizers](../results/combined-cover-v1/sanitizer.log) · [Native tests](../results/combined-cover-v1/native-release.log)

Enabling GLPK presolve on all transformed simplex proposals lost the checked LP completions on air05 and nw04 at six seconds. The shipped initializer keeps its previous settings. A separate experiment enlarged the bounded exact sparse basis and indexed its singleton rows and columns. It recovered checked LP points and valid bounds on binkar10_1 and cost266-UUE. Neither full MILP completed in the five-second follow-up. [Presolve comparison](../results/lp-presolve-v1/records.jsonl) · [Exact-basis LP checks](../results/sparse-pricing-v2/records.jsonl) · [MILP follow-up](../results/sparse-pricing-v2/mip/records.jsonl)

## Reproducing the release

Use the accelerated C17 Release build with GLPK 5.0 and the patched Kissat 4.0.4 static library. The manifest records source, executable and library hashes. The native-only configuration includes the path production solver. The activated covering and Boolean methods use their respective optional libraries. [Build instructions](NETWORKS-AND-PACKING.md#reproducing-the-measured-build) · [Configuration and source hashes](../results/combined-cover-v1/manifest.json)

From the repository root, repeat the eight 60-second head-to-heads with the frozen binaries. Choose a fresh output directory for each run.

```sh
work/combined-cover-v1/validated/milpbench --mip \
  --suite results/combined-cover-v1/panel.tsv \
  --solvers results/combined-cover-v1/paired.tsv \
  --checker work/combined-cover-v1/validated/milpcheck \
  --limit 60 --grace 15 --seed 1 --out results/combined-cover-repeat-60s
```

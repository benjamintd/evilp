# Native solver performance checkpoint

The current release completed 2 of the 240 locked MIPLIB instances at a two-second limit per model. It returned 227 time limits and 11 numerical failures. Thirteen original-model witnesses passed the independent checker. The 120-win target remains open. The [full inventory gate](../results/solver-iteration-v1/gate.json) records all 240 attempts.

## Selected MILP comparisons

The C17 runner executed the current release and HiGHS 1.15.1 serially on this arm64 Mac with a 60-second limit, one requested thread and seed 1. Times include process startup. The candidate solver is deterministic. The checker used the original MPS files.

| Model | Solver | Status | Wall seconds | Checked objective | Locked optimum |
|---|---|---|---:|---:|---:|
| enlight_hard | lab | optimal | 0.030 | 37 | 37 |
| enlight_hard | highs | optimal | 9.751 | 37 | 37 |
| markshare_4_0 | highs | time_limit | 60.020 | 14999999999999/5000000000000 | 1 |
| markshare_4_0 | lab | optimal | 0.060 | 1 | 1 |

Both candidate objectives agree with their locked references. Optimality is the solver’s reported status. The checker independently establishes witness feasibility. The paired-speed gate requires correct completion on both sides of a comparison and correct completion across the inventory. [Raw runs](../results/solver-iteration-v1/paired-60s/runs.jsonl)

## LP comparisons

Median solver-process CPU across three serial pairs on the final build. Each process had ten seconds. The checker evaluated primal points and lower bounds against the original LP.

| Model | Solver | Outcomes | CPU seconds | Valid points | Valid bounds |
|---|---|---|---:|---:|---:|
| air05 | lab | time_limit | 9.835 | 0/3 | 0/3 |
| air05 | highs | Optimal | 0.325 | 3/3 | 3/3 |
| nw04 | lab | optimal | 1.374 | 3/3 | 3/3 |
| nw04 | highs | Optimal | 1.377 | 3/3 | 3/3 |

On the separate six-model LP panel, the candidate produced checked primal points and valid lower bounds on all six. HiGHS had a lower median CPU on all six. Its exported point on glass4 exceeded the original-row tolerance with a violation of 3.6e-7. [Six-model runs](../results/indexed-factor-v1/paired-lp/records.json)

## Engine changes

- Adaptive simplex pricing now also runs during MILP search. With an eight-second limit on the first node of ran14x18-disj-8, the old MILP default stopped after 5,024 pivots. Adaptive pricing processed the node in 1,847 pivots.
- Presolve derives parities from source-exact integer equalities. On enlight_hard this fixes all 200 variables before any LP pivot.
- A meet-in-the-middle search handles small 0–1 models with a bounded number of rows. Eligible equality slacks are substituted exactly. Limit and unsupported cases fall back to ordinary search.
- For models with at least 1,024 columns, the CLI checks for upper bounds implied by packing rows and for duplicate set-partitioning equalities. Runs that read or write a basis require explicit opt-in with `--packing-bounds`. The compact formulation reduced the nw04 LP from 87,518 rows to 36. The older formulation hit its ten-second limit in all three development repetitions. Use `--no-packing-bounds` to disable the reduction.
- Sparse LU indexes remaining rows by their next nonzero column. The median CPU for eight repetitions of the 4,096-row banded contract fell from 0.353 to 0.010 seconds. Both builds used about 0.745 seconds on the 1,000-pivot ran14x18-disj-8 LP prefix. [Differential factor checks](../results/indexed-factor-v1/differential.json)
- Objective evaluation can establish exact zero contributions from source enclosures. In the preceding full inventory, 15 models stopped at the initial objective evaluation. One stopped there in the current inventory.

## Validation and scope

All 48 suites passed in the normal release and ASan/UBSan builds. The release includes the earlier agent’s working-tree changes. Each experiment has separate model, source and binary hashes. The two-second inventories overlapped correctness testing and support coverage comparisons only. Of the 227 current time limits, 166 stopped with one node and one LP call. The locked best-known reference for supportcase22 remains unresolved. [Release tests](../results/solver-iteration-v1/release.log) · [Sanitizer tests](../results/solver-iteration-v1/sanitizer.log) · [Manifest](../results/solver-iteration-v1/manifest.json) · [Measurements](../results/solver-iteration-v1/summary.json)

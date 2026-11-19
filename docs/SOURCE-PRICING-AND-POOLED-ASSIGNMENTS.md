# Source pricing, exact covers and pooled assignments

The combined release recorded 11 strict paired wins on the locked benchmark of 240 MILPs. The target is 120. The eleven selected pairs ran serially with a 90-second budget, one thread and seed 1. The C17 runner measured solver wall time from process startup through solution output. Both completions had to pass original-model witness and reference-objective checks. [Measured pairs](../results/pooled-assignment-v1/final/paired-90s/runs.jsonl) · [Locked build and campaign](../results/pooled-assignment-v1/final/manifest.json)

| MILP | Candidate wall seconds | HiGHS wall seconds |
|---|---:|---:|
| decomp2 | 0.096 | 3.928 |
| drayage-100-23 | 0.241 | 0.789 |
| enlight_hard | 0.025 | 10.103 |
| exp-1-500-5-5 | 0.022 | 4.042 |
| irp | 0.374 | 28.008 |
| neos-2978193-inde | 0.274 | 52.978 |
| neos-3083819-nubu | 13.099 | 48.618 |
| neos-3381206-awhea | 0.025 | 1.658 |
| nw04 | 0.827 | 47.350 |
| p200x1188c | 0.044 | 0.749 |
| sp150x300d | 0.025 | 0.149 |

Across the separate two-second sweep, 15 models completed correctly and 15 exited within two wall seconds. 26 original-model feasible witnesses were checked. The sweep used two concurrent processes and recorded 1 infeasible, 12 numerical_failure, 15 optimal, 212 time_limit. The two longer selected panels together yielded 18 correct candidate completions. [All 240 outcomes](../results/pooled-assignment-v1/final/gate.json) · [Classification and times](../results/pooled-assignment-v1/final/summary.json)

## How the two new searches work

- On `irp`, bounded source-objective checks made the root LP usable. The solver converts exact decimal cover costs to integer units, repairs the root row prices and searches compatible covers. It prunes when the remaining rows' minimum available costs exceed the budget after sharing each column's cost among at most its number of rows. Enlarged accumulators handle the six-decimal objective. The original file's checked objective is 12,159.492836. [Cover search](../src/partition.c) · [Zero-tolerance original-model check](../results/pooled-assignment-v1/final/exact-irp.json)
- On `neos-2978193-inde`, 64 facilities form eight cost and capacity classes across eight opening-quota groups. The C engine solves a network flow relaxation using GLPK, validates integral flows and complementary slackness in integer arithmetic, and branches on the number of open sites per class. Alternating paths realize individual openings and the solver assigns each client to an open site. The checked original objective is −2.388061686. The method recognizes row structure and compares every original assignment cost. The independent checker verifies the witness. The internal search supplied the optimality claim. [Network solver](../src/pooled_assignment.c) · [Zero-tolerance original-model check](../results/pooled-assignment-v1/final/exact-neos-2978193-inde.json) · [Opening-pattern oracle](../results/pooled-assignment-v1/matching-oracle.json)

## LP relaxations

The previous release stopped on source-objective pricing on all five of these models. The candidate can now retain isolated non-binary coefficients even when the full rational snapshot exceeds its storage budget. It also checks omitted column entries without re-reading represented entries and reconstructs larger coupled basis cores after singleton elimination. Source values, bit limits and the independent checker retain their validation roles. [Prior and candidate diagnostics](../results/source-coverage-v2/development.jsonl)

The LP measurements below use median solver-process CPU seconds over three serial pairs per model. Each solver had ten seconds. The last two columns count independently checked primal points and valid dual bounds. A valid dual bound can differ from the primal objective. [Paired LP measurements](../results/pooled-assignment-v1/final/lp/records.json)

| LP relaxation | Candidate CPU seconds | HiGHS CPU seconds | Candidate points, bounds | HiGHS points, bounds |
|---|---:|---:|---|---|
| irp | 0.213 | 0.230 | 3/3, 3/3 | 3/3, 3/3 |
| neos-2978193-inde | 0.532 | 0.131 | 3/3, 3/3 | 3/3, 3/3 |
| sct2 | 0.418 | 0.252 | 3/3, 3/3 | 3/3, 0/3 |
| sp97ar | 1.308 | 0.714 | 3/3, 3/3 | 0/3, 0/3 |
| sp98ar | 1.246 | 0.774 | 3/3, 3/3 | 0/3, 0/3 |

## A third model in development

On `nu25-pr12`, an experimental C17 presolver enumerated all 2,545 simple routes across 81 independent flow components. The objective and capacity inequalities permit deleting cycles. The presolver fixed 4,554 arc variables to zero and 156 to one. Two in-process SCIP searches produced the checked original objective 53,905 in 2.084 and 1.733 seconds. These development runs used a separate executable and contribute zero wins to this release. [Exact replay and development records](../results/path-pruning-v1/summary.json)

## Reproduction and limits

The build uses the preceding release's locked GLPK, Kissat, SCIP and QSopt_ex dependencies. The new pooled network solver requires GLPK. All 61 accelerated release contracts, 61 accelerated ASan/UBSan contracts and 54 native-only contracts passed. Sanitizers instrumented the project code and used release builds of the imported libraries. The tests exercise equal binary64 coefficients with distinct source decimals, duplicate records, copied models, coupled rational systems, exact decimal costs, permuted facilities, fractional source controls and exhausted budgets. [Build instructions](INTEGER-SEARCH-AND-SPARSE-LP.md#build-and-verification) · [Source hashes](../results/pooled-assignment-v1/validation/source.json) · [Release tests](../results/pooled-assignment-v1/validation/release-contracts.log) · [Sanitizers](../results/pooled-assignment-v1/validation/san-contracts.log) · [Native tests](../results/pooled-assignment-v1/validation/native-contracts.log)

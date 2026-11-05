# Sparse basis updates and complete corpus coverage

The evaluated 0.7 engine is committed at `82c1355`. Documentation and the
reproducible inventory driver bring the frozen campaign revision to `61a738c`.
All development stays on the local `main` branch. The independent checker,
original-model acceptance policies and numerical residual guards are unchanged.

## Implemented foundation

Sparse factorization peels structural singleton rows and columns and factors
only the remaining coupled core. A generated 100,000-dimensional basis has a
three-dimensional core and passes normal, transpose and in-place solves.
Guarded eta updates reuse that factorization, with full-refactor fallback,
32 retained updates and at most four million eta entries. Basis membership
lookups use a bitmap. The generated 2,500-variable LP contract verifies its
known original point, a three-variable core, independent checking and fewer
full refactorizations than `--no-updates`.

Source-aware root propagation returns conservative bounds and redundant-row
masks without mutating the original model. Search uses deterministic
reliability-based pseudocost ranking learned from accepted solved children.
`--no-presolve` and `--first-fractional` retain explicit controls.
[The engine contract](LP-ENGINE.md) describes resource limits and ablations.

The sparse standard-form and basis ceilings are 16 million and eight million.
All 240 recorded original-model shapes fit conservative dimension and
transformation-entry bounds. This does not guarantee memory availability or
that sparse elimination remains within its separate 64-million-entry fill
budget. The dense reference ceiling remains 4096.

## Retained validation

- [Full Release and sanitizer contracts](../results/milestone-2g-contracts-v1/summary.json):
  all 22 suites pass at `82c1355`; matching engine/test hashes are recorded.
- [Native MILP crosscheck](../results/mip-native-v6/summary.json): eight models
  across HiGHS, SCIP, Cbc and warm/cold lab search give 40 matching outcomes;
  all 35 feasible candidates pass independent checking.
- [Full v2 corpus inventory](../results/miplib-candidate-v2/report.html): all
  240 locked instances are attempted using evaluators built directly from the
  archived `61a738c` source. Commands, compiler identity, source hashes,
  evaluator hashes, witnesses and failures remain available.

The corpus campaign uses a two-second CLI budget, fifteen-second watchdog
grace, seed 1 and serial solver processes. Engineering continued concurrently,
so the campaign supplies development coverage evidence. It cannot establish a
performance ranking or a paired HiGHS win.

| Retained outcome | Frozen v1 | Frozen v2 |
| --- | ---: | ---: |
| Unsupported | 209 | 0 |
| Time limit | 31 | 235 |
| Numerical stop | 0 | 5 |
| Checked incumbent | 7 | 9 |
| Correct completed outcome | 0 | 0 |

Both inventories retain all 240 models. Checked incumbents are a subset of
limited runs, not additional completed cases. The v2 campaign has no missing
or duplicate attempts, malformed logs, crashes, watchdog kills, invalid
candidates, checker errors or reference contradictions. Archived source and
evaluator hashes remain unchanged at its end.

Four numerical stops (`fhnw-binpack4-4`, `glass4`, `graphdraw-domain`,
`neos859080`) report an uncertain Phase I dual sign where exact reconstruction
is unavailable. They preserve the existing refusal to infer infeasibility
from an uncertain sign. `mad` fails original-source primal reconstruction with
root presolve; disabling presolve permits root reconstruction in the retained
diagnostic. These remain unresolved outcomes.

## Follow-up checkpoints outside the frozen campaign

`303d9bf` caches structural column norms after scaling and refreshes them when
rows are deleted. This avoids repeated coefficient scans while retaining the
same uncertainty guard. Five focused contracts pass in Release and sanitizer
builds; logs are in `work/norm-cache/`.

`4c02d7e` retains transactional replacement-column overlays instead of copying
the complete current CSC after each accepted eta. Explicit accessors expose
the current matrix; the immutable base cannot be mistaken for it. Iteration
order, rollback, input aliasing, residual thresholds and storage caps remain
unchanged. Expanded standalone linear contracts pass in Release and
ASAN/UBSAN. The combined
[Release run](../results/basis-overlays-contracts-v1/release-ctest.log) passes
all 22 suites, and the
[targeted sanitizer run](../results/basis-overlays-contracts-v1/sanitizer-ctest.log)
passes reoptimization, linear refinement and scalability. A complete sanitizer
run of this later checkpoint is not claimed here.

The [completion gate](AUTORESEARCH.md) remains closed: zero of the 240 corpus
problems completes correctly at the declared budget. Root-LP progress,
source-safe infeasibility validation and the `mad` reconstruction stop are
current engineering work. No paper has been selected and no paper-driven
autoresearch cycle has started. The later speed target still requires at least
120 strict matched HiGHS wins over the unchanged 240-instance denominator.

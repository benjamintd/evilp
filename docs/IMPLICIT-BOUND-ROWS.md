# Implicit bound-row experiment

The control is the sparse-update solver `c4df495`, with report checkpoint
`2be211e`. The [predeclared experiment](../results/implicit-bound-rows-v1/plan.json)
retains the original models and control binary before implementation.

For a lower-anchored variable `x = lower + scale*z`, nonnegativity of `z`
already supplies the lower bound. For an upper-anchored variable
`x = upper - scale*z`, it supplies the upper bound. Omit that explicit row only
when checked original-source data proves that the anchor is the exact bound.
Keep uncertain bounds explicit; retain the opposite finite bound.

Those rows previously supported RHS-only warm updates. To preserve reuse after
their removal, move an implicit anchor when its bound changes and recompute
the affected source-row RHS and explicit bound widths. Keep coefficients,
basis state, source ownership and all numerical guards intact. Test tightening,
fixing and restoration against independent optima and cold solves. The
`--keep-bound-rows` control preserves the earlier representation and retained
basis-file regressions.

Measure actual standard-form size, pivots and CPU on the original
`sp150x300d` and `ran14x18-disj-8` relaxations, then compare ten-second MILP
runs on the same six-model panel as the sparse-update experiment. Source
checking and the full 240-instance gate remain unchanged. Limits, numerical
failures and unchecked solutions never become completed problems.

The implementation passes all 45 Release and ASAN/UBSAN suites. The new
contract checks 308 session points independently against exact Fraction
optima, covering cold/warm agreement, lower/upper anchors, fixing, restoration,
scaling, ranged rows, uncertain source bounds and large source snapshots.
The unchanged original 2,500-column scalability suite also passes in both
builds under `--keep-bound-rows`; default-form capacity tests use 5,000 columns.
Handcrafted historical basis regressions explicitly retain their original
representation. Their expected outcomes and tolerances are unchanged.

Size probes on the original LPs confirm 1,650 to 1,050 rows for sp150x300d,
and 1,455 to 951 for ran14x18-disj-8. sp150x300d reaches the same checked LP
objective and valid lower bound in 447 pivots with median CPU time decreasing
from 0.101 to 0.069 seconds across three alternating pairs. ran14x18-disj-8
still reaches ten seconds in its first LP. The declared ten-second target is
therefore not met by this change alone.

The completed MILP panel uses three serial alternating pairs per model, with
ten seconds per attempt. Median nodes are:

| Model | Control | Compact form | Outcome |
| --- | ---: | ---: | --- |
| sp150x300d | 394 | 505 | time limit |
| ran14x18-disj-8 | 1 | 1 | time limit in first LP |
| gen-ip002 | 2,055 | 2,353 | time limit |
| mas76 | 1,617 | 1,995 | time limit |
| pk1 | 497 | 713 | time limit |
| glass4 | 441 | 534 | numerical failure |

All 30 incumbent files from the 36 runs pass the independent checker. The
best sp150x300d incumbent improves from 97 to 96; its locked optimum is 69.
Node throughput improves by 14–43% on four models that reach the full time
budget. Because the formulation can change search paths, those node counts
are not fixed-work speed ratios. External host activity was uncontrolled.
The two variants have zero completed MILP runs and no watchdogs or signals.

Both forms still stop numerically on glass4. The control stops on an uncertain
Phase I source dual sign at node 441 / 15,655 pivots. The compact form stops
during fixed-integer repair at node 534 / 22,033 pivots. A separate candidate
run with `--keep-bound-rows` reproduces the original control stop. Both
failures and all emitted points are retained. The compact incumbent passes
the checker; that does not validate its unresolved node or complete the solve.
A separate 32-node mas76 trace emits four root cuts, all independently
accepted from original-source certificates. This is a local cut validation,
not a proof of the entire search.

The decision is to keep the reduction at `eda3561`. The frozen
[report](../results/implicit-bound-rows-v1/REPORT.md),
[summary](../results/implicit-bound-rows-v1/summary.json) and
[verification](../results/implicit-bound-rows-v1/verification.json) bind the
129 tested source files, binaries, settings, source models and terminal
process records. Historical results are unchanged. Native finite-bound
handling and pricing are outside this experiment.

Warm-start reuse holds in the targeted contracts, but cold fallbacks remain
substantial in real search. Degenerate pivots also remain: the first compact
ran LP attempt records 3,775 zero-step pivots among 7,070 accepted primal
pivots. The next LP investigation should measure remaining bound-row costs,
zero-step pivots and failed reuse before selecting the next change. No full
240-model campaign was rerun; its completion gate and paper-driven
autoresearch remain closed.

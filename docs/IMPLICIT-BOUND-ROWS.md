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
therefore not met by this change alone. The MILP comparison is recorded after
all panel runs terminate.

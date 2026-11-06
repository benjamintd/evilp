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

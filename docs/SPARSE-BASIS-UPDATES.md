# Sparse basis-update experiment

The control is `6be37e0`, whose 43 Release and sanitizer suites pass. The
[predeclared experiment](../results/sparse-basis-updates-v1/plan.json) freezes
the original MIPLIB models, control executable, fixed-work prefixes and
alternating ten-second MILP comparisons before implementation.

Following the [HiGHS audit](HIGHS-LP-AUDIT.md), retain every exact nonzero
update coefficient in ascending index order. Use packed off-diagonal entries
at density at most one quarter and retain a dense path otherwise. Keep the
current pivot policy, 32-update/four-million-entry capacity limits, original-unit
residual checks, refinement, weak-pivot and underflow guards.

Compensated addition can renormalize even on zero terms. Replay a skipped zero
run until the accumulator reaches a numerical fixed point, rather than assume
that omitting zero operations is harmless. The supported comparison treats
positive and negative zero as the same numerical value; nonzero arithmetic,
finite-input rejection and observed solver outcomes must agree with the dense
control. Differential tests and real fixed-work prefixes test that contract.

Primary LP cases are `sp150x300d` and `ran14x18-disj-8`; a `mas76` MILP prefix
also exercises repeated reoptimization. The six declared MILP models add
`gen-ip002`, `pk1` and `glass4`. They remain a diagnostic subset of the required
240-instance corpus. Limits and numerical failures never count as completed
problems. Retain failed experiments and independently check original-model
witnesses before accepting any speed result.

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

The implementation passes all 44 Release and ASAN/UBSAN suites. Its new test
compares 30,000 sparse/dense eta calculations numerically, checks an exact
dyadic oracle, and retains public finite-input, subnormal-underflow, dense
fallback and dense-equivalent storage-limit controls. Nonzero coefficients
are never thresholded. Only zero signs are outside the equality contract.

The first [fixed-work comparison](../results/sparse-basis-updates-v1/fixed-work-analysis.json)
uses three alternating control/candidate pairs per case:

| Original MIPLIB workload | Control median CPU | Sparse median CPU | Control / sparse |
| --- | ---: | ---: | ---: |
| sp150x300d, 400 LP pivots | 0.234 s | 0.089 s | 2.64 |
| ran14x18-disj-8, 3,000 LP pivots | 3.171 s | 2.975 s | 1.07 |
| mas76, 200 MILP nodes, cuts/diving disabled | 0.672 s | 0.563 s | 1.19 |

Every non-timing JSON result field agrees across each pair. These runs end at
their declared iteration/node limits and are not completed benchmark solves.
All emitted witnesses pass independent original-model checking. The improvement
is workload-dependent; external host activity was uncontrolled, and three short
pairs do not establish a general performance ranking. The ten-second MILP
campaign is retained in the [final experiment report](../results/sparse-basis-updates-v1/REPORT.md).

The six-model MILP campaign uses three alternating pairs per model. Median
nodes in ten seconds are:

| Model | Control | Sparse | Outcome |
| --- | ---: | ---: | --- |
| sp150x300d | 143 | 381 | time limit |
| ran14x18-disj-8 | 1 | 1 | time limit in first LP |
| gen-ip002 | 2,050 | 2,034 | time limit; no clear throughput change |
| mas76 | 1,233 | 1,534 | time limit |
| pk1 | 462 | 493 | time limit |
| glass4 | 441 | 441 | same numerical failure |

Sparse updates find a checked sp150x300d incumbent of 97 in every repetition,
where the control finds none; the locked optimum is 69. The best mas76
incumbent also improves. All 27 emitted MILP Lab witnesses pass original-source
checking, with no watchdogs. Neither variant completes a panel problem.

A posthoc 1,500-node gen-ip002 comparison leaves the implementation unchanged
and confirms neutral CPU cost: 6.735 s versus 6.733 s, with all non-timing
fields equal. The HiGHS reference completes sp150x300d at 69. Its other five
runs hit limits. Its original glass4 file fails the objective-header check;
the unchanged coordinate values pass a separate point-only check. Both
records are preserved, and no tolerance is loosened.

The decision is to keep this change. It improves the measured sparse workload
and preserves the numerical contract, while providing little benefit on some
models. The next LP investigation should address the remaining formulation
and degeneracy costs; this experiment does not establish a competitive MILP
solver or open the full-corpus gate.

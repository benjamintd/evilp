# Correctness and first benchmark completion

The latest full inventory at `95b7450` completes zero of 240 MIPLIB models
within ten seconds. It records 239 time limits, one numerical failure and
15 independently checked incumbents. Of those time limits, 199 stop in the
first LP. Component-test coverage and incumbent counts have improved, but
neither establishes a competitive solver. All 240 models remain required.

Source-decimal cancellation is a correctness blocker. Duplicate input terms
can have a nonzero exact sum even when their binary representatives cancel.
Separate objective coefficients can also cancel along a feasible direction:
maximizing `1.000000000000000001 X - Y`, with `X = Y` and both variables in
`[0, 10^15]`, has optimum `1/1000` at the upper bounds. Rounding each cost to
binary64 before optimization loses that improvement. Fixing the printed
objective alone cannot repair the solution. The new source-accuracy contract
requires actual correct optimal points and independent original-source checks,
including larger models, LP/MIP, fixed/free choices and solver-option controls.
Numerical abstention does not satisfy these regressions.

Feature expansion is deferred until that correctness work passes. The next
performance change must address a measured bottleneck on a predeclared real
model. The frozen [LP diagnostic](../results/root-lp-profile-v1/REPORT-v2.md)
uses `sp150x300d`, a historical HiGHS-easy model, and `ran14x18-disj-8`, the
smallest model by rows plus columns among the 199 first-LP timeouts. Our
original sp relaxation has an independently checked primal and LP lower bound;
its MILP still times out. The ran relaxation times out in our solver while
HiGHS finishes it. These diagnostic runs had concurrent engineering activity
and do not establish fair speed ratios.

Basis-update application and compensated residuals dominated the initial CPU
samples. Stored update directions were frequently sparse, but application
scanned every entry. The first experiment used a sparse representation
using exact zeros only, with a dense fallback. It retained original-unit
residuals, refinement, finite-input rejection, underflow and weak-pivot guards,
and resource limits. Ascending nonzero order alone did not prove arithmetic
equivalence: compensated addition renormalizes even on zero terms, and signed
zero required explicit examination. The experiment validated these semantics
and compared fixed-work prefixes against a frozen version with the same
correctness fixes, retaining all changed outcomes and failed runs.
The [completed sparse-update experiment](SPARSE-BASIS-UPDATES.md) passed this
contract and reduced fixed-work CPU with a workload-dependent benefit.

The subsequent [bound-row experiment](IMPLICIT-BOUND-ROWS.md) removes only
source-proved redundant transform bounds and reanchors reusable LP sessions.
It passes all 45 Release and sanitizer suites. The original sp relaxation
shrinks from 1,650 to 1,050 rows and uses 0.069 versus 0.101 seconds median
CPU, with the same 447 pivots and independently checked point/lower bound.
The six-model search panel processes more nodes on four time-limited models,
but still completes none. ran14x18-disj-8 remains in its first LP; glass4 has
unresolved numerical stops in both forms. The full inventory above was not
rerun and does not describe the new engine's measured coverage.

The next controlled LP investigation should distinguish cost from remaining
explicit finite-bound rows, degenerate zero-step pivots and cold fallbacks.
The compact ran relaxation's first attempt has 3,775 zero steps among 7,070
accepted primal pivots. Removing redundant rows therefore addresses only
part of its difficulty. Native finite-bound handling and stronger pricing
need separate experiments, source-unit checks and preserved controls.

Faster kernels will not by themselves supply HiGHS's cuts, presolve and search
strength. First establish a correct completed real benchmark, then measure the
next limiting operation. Additional cuts, branching policies and heuristics
require evidence from that sequence. Paper-driven autoresearch remains closed
until the [full completion gate](AUTORESEARCH.md) passes; the later speed target
remains at least 120 strict matched wins over the fixed 240-model denominator.

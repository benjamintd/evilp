# LP comparison and optimization

The user requested a direct LP benchmark against HiGHS and continued
optimization toward matching it. The frozen starting engine is `eda3561`,
with report checkpoint `4a3f27e`. The original six MIPLIB models and binaries
are locked before changes in `results/lp-parity-v1/plan.json`.

Compare original LP relaxations, three serial alternating pairs, ten seconds,
one thread, and original-source primal/bound checks at the existing policies.
HiGHS 1.15.1 uses its ordinary solver choice with presolve enabled; this is
the practical target, including its presolve advantage. Record process CPU,
wall time, pivots, memory and every failed or limited attempt. Builds and
tests must finish before timing; external host activity is uncontrolled.

Parity requires both solvers to complete correctly with independently valid
primal and LP-bound records and agreeing original objectives, and our median
process CPU no greater than HiGHS on each model. A limited subset cannot
establish global parity; all 240 models remain in the MILP completion gate.
No limit or numerical failure counts as a solved LP. The first investigation
isolates the permanent switch to Bland pricing on a single zero step, with
cycle protection, legacy controls and numerical/source guards retained.

Each additional cohort or budget is declared before its runs. Improvements
require exact/adversarial contracts and sanitizers, then a fresh frozen
candidate comparison. This engineering is separate from paper-driven
autoresearch, whose full MILP gate remains closed.

# Autoresearch gate and future protocol

Autoresearch starts only after every MILP problem in the existing benchmarks
passes. Until then, work remains benchmark completion and solver correctness.
Do not start paper selection or an autoresearch implementation cycle early.
Once the gate is met, each cycle finds an interesting primary research paper
and implements one clearly identified idea from it.

## Completion gate

Freeze a versioned inventory of problem membership, expected outcomes,
configurations, reference versions, tolerances and budgets. The current inventory
includes all current CTest suites and their generated/adversarial cases, the
eight native MILP fixture models in the
[crosscheck](../results/mip-native-v5/summary.json), the selected development
MILP cases, and all 240 MIPLIB benchmark instances identified by the
[dataset lock](../config/miplib.lock.json) and imported suite. Newly added
regressions remain part of the gate; no failing case may be dropped or skipped.
The original 18 suites remain in this inventory. Presolve, branch selection,
benchmark-gate, scalability, source-candidate and source-sign repair contracts
are additive requirements. The 24-suite source-safety checkpoint `ade4924`
passes Release and sanitizer checks, with
[recorded logs and source hashes](../results/source-safety-contracts-v1/summary.json).
The later `e02a45e` candidate retry passes all 24 Release suites and its focused
certificate sanitizer checks. Later additions remain required even when an
older frozen checkpoint has fewer suites. Earlier
[18-suite Release and sanitizer evidence](../results/primal-regressions-v1/final-validation.json)
remains a frozen record of that checkpoint, not evidence for the added suites.

For the candidate solver and its declared required configurations:

- Every contract passes, including independent-checker adversarial controls.
- Every MILP problem completes with its correct declared outcome. A known
  optimal reference requires correct completion and an independently checked
  feasible witness whose objective matches the locked reference policy.
- No candidate witness is invalid and no checker error is ignored. Known
  infeasible or unbounded cases must meet their declared validation contract;
  an unbounded LP relaxation is not proof of MILP unboundedness.
- `unsupported`, `numerical_failure`, `time_limit`, `node_limit` and
  `iteration_limit` are unresolved outcomes, never problem-level passes.
  Unknown references remain unresolved until their expected outcome and
  validation basis are established; do not invent proof or optimality status.

A passing campaign means its assertions passed. It may intentionally retain
old failures or rejected witnesses as negative controls. Reference agreement,
feasible witnesses, numerical statuses and independently checked proofs remain
distinct; individual LP certificates do not constitute a complete MILP tree proof.
Any budget change creates a fresh declared campaign, preserving earlier results.

The gate is currently **not complete**. The latest frozen
[ten-second full inventory](../results/miplib-candidate-10s-v1/inventory.json)
contains one attempt for every locked instance: 232 `time_limit`, four
`numerical_failure`, four `watchdog_timeout`, 13 independently checked
incumbents and zero correctly completed problems. Its
[gate record](../results/miplib-candidate-10s-v1/gate.json) has no missing
instances and leaves both completion and speed gates closed. It evaluates
`332afbd`, before the later rational reconstruction and root cuts, while
engineering ran concurrently. The earlier
[two-second v3 inventory](../results/miplib-candidate-v3/inventory.json) retains
its 238 time limits, two numerical failures and nine checked incumbents.
These are coverage and failure evidence, explicitly not comparable timing
campaigns. The `=best=` reference for `supportcase22`
also remains an unresolved expected-outcome policy; it cannot be promoted to a
known optimum. Passing unit tests or finding incumbents does not open the gate.

Version 0.7 adds singleton peeling, guarded basis updates, larger finite sparse
capacity guards, source-aware root propagation and learned branch ranking.
The [earlier v1 inventory](../results/miplib-candidate-v1/inventory.json) retains
its 209 size-cap stops and 31 time limits. The v2 inventory retains its 235 time limits and five numerical stops. Each
archive and evaluator remains unchanged. Later pricing, scaling shortcuts and
candidate constructors are outside the v3 evaluated source; every subsequent
candidate requires fresh evidence.
The capacity guard and the cost of factoring a large remaining coupled core
are separate blockers. Larger allowed dimensions do not establish completion,
memory feasibility or a speed advantage. The earlier
[smoke](../results/miplib-smoke-2s-v4/runs/summary.html) and
[selected comparison](../results/primal-recovery-v1/report.html) remain available
as historical development evidence.

The [source-safety follow-up](SOURCE-SAFETY.md) gets past three original
numerical stops in selected development reruns. Those runs still hit limits;
they do not replace a complete versioned 240-instance campaign or open the gate.

After correct completion, the requested performance target is strictly faster
than HiGHS on at least 120 of the fixed 240 MIPLIB instances. Use matched
budgets, hardware, threads, seeds, solver identities and timing accounting, with
correct completion and independently checked candidate witnesses. The
denominator remains 240: unsupported cases, limits, invalid candidates and
ties cannot become strict wins or be removed to improve the fraction. An
unmatched or concurrently disturbed development run cannot satisfy this paired
speed gate. No speed claim is made by either inventory or the 0.7 contract tests.

Pricing and pricing-integration tests are also additive prerequisites. All 26
suites pass in both builds at `af61788` and the combined `332afbd` checkpoint;
[recorded evidence](../results/simplified-source-contracts-v1/summary.json) binds
tested source files to the commit. The ten-second campaign preserves every
member and earlier result. A longer declared budget cannot turn a limit into a
pass. Native tableau, GMI arithmetic/oracle, independent cut checker,
augmentation, root-cut pipeline/search and cut-aware trace contracts are
additional required prerequisites once registered. Root cuts retain their
explicit local MIP-box proof scope; they do not open this gate by themselves.

## One-paper, one-idea cycle after the gate

1. Identify an observed solver bottleneck and find a relevant primary paper.
   Record its citation, accessible source, precise algorithm or theorem, and
   any code/data license and provenance relevant to reuse.
2. Isolate one implementable idea. Write the hypothesis, assumptions,
   correctness obligations, expected benefit, possible regressions and
   predeclared keep/reject metrics before changing the solver.
3. Implement a reversible, isolated change in C. Add meaningful exact or
   adversarial regressions. Keep the independent checker separate, preserve
   tolerances and resource guards, and never specialize on filenames, hashes
   or stored answers.
4. Compare against frozen controls with identical budgets, hardware, compiler,
   solver settings and accounting. Record threads, seeds, memory limits and
   model identities. Run timed solver processes serially unless the experiment
   explicitly measures parallel execution under controlled conditions.
5. Evaluate both development cases and predeclared held-out model families or
   seeds. Preserve raw positive and negative results, including every failed
   run. Do not tune on a holdout and continue describing it as held out.
6. Accept only after unchanged correctness checks pass and the declared metrics
   support the hypothesis. Include limits and failures in the analysis; report
   solved counts, time or gap metrics and memory according to the stated plan.
   Otherwise reject or revise the idea in a new recorded experiment. Retain
   rejected evidence and explain the decision.

No paper has been selected and no autoresearch cycle has started for this
milestone. The immediate work is satisfying the completion gate.

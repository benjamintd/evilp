# Autoresearch gate and future protocol

Autoresearch starts only after every MILP problem in the existing benchmarks
passes. Until then, work remains benchmark completion and solver correctness.
Do not start paper selection or an autoresearch implementation cycle early.
Once the gate is met, each cycle finds an interesting primary research paper
and implements one clearly identified idea from it.

## Completion gate

Freeze a versioned inventory of problem membership, expected outcomes,
configurations, reference versions, tolerances and budgets. The current inventory
includes all eighteen CTest suites and their generated/adversarial cases, the
eight native MILP fixture models in the
[crosscheck](../results/mip-native-v5/summary.json), the selected development
MILP cases, and all 240 MIPLIB benchmark instances identified by the
[dataset lock](../config/miplib.lock.json) and imported suite. Newly added
regressions remain part of the gate; no failing case may be dropped or skipped.

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

The gate is currently **not complete**. The
[latest smoke](../results/miplib-smoke-2s-v4/runs/summary.html) has six size-cap
failures and six time limits; none reaches a reference optimum. The two full
MILP development runs in the [selected comparison](../results/primal-recovery-v1/report.html)
stop at 1,000 nodes. The complete 240-instance candidate campaign still needs
successful completion and validation. Passing unit tests or finding four
incumbents does not open the gate.

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

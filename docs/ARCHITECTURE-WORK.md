# Evaluator and persistent LP sessions

**Coordination, November 22, 2026.** Implementation and verification are complete. The evaluator lives in `scripts/run_evaluation.py`. The native session API lives in `include/lp.h`, with its implementation in `src/lp.c`, sparse row construction in `src/augment.c` and branch-and-bound integration in `src/bb.c`. Agents can build on these interfaces. The v22 dependency lock identifies the feature-enabled verification build.

## Run evidence

The C runner owns the timed observation and records the solver claim, model and binary identities, effective settings, checker invocation and hashes of retained evidence. The independent C checker evaluates the original model. The Lab LP adapter retains Farkas certificates alongside its existing support for solutions and dual bounds.

Use `scripts/run_evaluation.py` to interpret a checker report. `assess_evidence(run, reference, relaxation=...)` returns the reported status and bound, checked objective and scope, reference agreement, certified bound and independent proof flags. Use `classify(model, run)` for the additional locked MILP completion policy. Campaign identity and strict paired timing belong to `benchmark_gate.py`. The gate, inventory, HTML/JSON report and replay consume the shared evaluator. Replay preserves the original solver claim and wall time and records the evaluator hash with each fresh assessment. Gate outputs, reports and inventory provenance identify the same policy by its SHA-256.

## Native LP sessions

`LabLpSession` borrows the immutable input and owns the transformed sparse matrix, numerical basis and bound overlays. A structural edit gives the session ownership of the augmented source matrix, including decimal provenance. The original model stays available to the checker.

- `lab_lp_session_begin_probe` and `lab_lp_session_rollback_probe` bracket tentative bound changes. The native fixed-integer repair uses these operations. Rollback restores the column bounds and compatible basis. A probe that required a different standard form forces a fresh solve when revisiting the parent. The original deadline applies throughout.
- `lab_lp_session_append_rows` accepts a sparse batch with sorted column indices. The dense-input convenience entry point serves the existing GMI cut producer. Both use one column allocation plan for the entire batch. Compatible inequality batches extend the basis with slacks and reoptimize by dual simplex. Changed transforms, removed basis rows or a new equality can require Phase I. Root cuts retain their existing box and proof checks.
- `lab_lp_session_set_tolerances` accepts tighter original-unit primal residual and relative reduced-cost tolerances. Zero selects the previous defaults. Reduced-cost comparisons retain the column-dependent floating-point roundoff guard. These settings apply to the native numerical engine. The independent checker records its own tolerance policy.
- General SCIP branch-and-cut continues to use its locked QSopt_ex LP interface. The native session does not implement the complete SCIP LPI. Experiments that swap that backend need a separate dependency identity and paired measurements.

## Verification

Release and sanitizer contracts cover dense and sparse reoptimization, rollback after infeasibility, owned batch storage, warm cut insertion, atomic rejection, exact source residuals, and relocation of primal, dual and Farkas evidence. A replay of the eight saved audit-panel attempts reproduced every checker report, solver status and wall time. Seven witnesses passed the checker and none had an independent optimality proof. This replay evaluates historical evidence and supplies no new timing result.

On a synthetic maximization LP with 32 independent pairs, each pair satisfies `x[i] + y[i] = 1` on the unit box and contributes `x[i]` to the objective. The root took 64 pivots and returned 32. A batch adding `x[0] <= 0.5` and `x[1] <= 0.5` brought the objective to 31. The session took two dual pivots. A cold native solve of the identical augmented model took 64 pivots.

The default native Release build passed all 56 contracts. A frozen source snapshot with SCIP/QSopt_ex v22, GLPK and Kissat enabled passed all 69 contracts in Release and all 69 with project code instrumented by ASan/UBSan. A separate overlapping three-row CSC regression passed against both Release and sanitized libraries. Final report, gate and replay checks also cover the shared policy digest. The solver dependency archives used for these tests came from the existing locked build.

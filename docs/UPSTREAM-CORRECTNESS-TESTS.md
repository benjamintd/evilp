# Tests for MILP terminal claims

Reviewed November 22, 2026 against the shared working tree, HiGHS v1.15.1 test sources and the SCIP source revision used by Lab. Recommendations below concern evidence for numerical claims. The frequency of incorrect claims across the benchmark remains unmeasured.

## A one-row counterexample with an exact optimum

Minimize $0.0220000004a + 0.022b$, subject to $a+b\ge 10^6$, with nonnegative continuous $a$ and nonnegative integer $b$. Thirty zero-cost columns put this through Lab's general SCIP dispatch. Every feasible objective is at least $0.022(a+b)\ge 22000$, and $b=10^6$, $a=0$ attains it.

| Run on the saved fixture | Claimed optimum in original objective units |
|---|---:|
| Lab v23, general SCIP | 22000.0004 |
| HiGHS 1.15.1, both MIP gap targets zero | 22000.0004 |
| Lab v23, native search | 22000 |
| Lab v23, general SCIP, costs multiplied by 1,000 and result divided by 1,000 | 22000 |

The independent checker accepts both the worse point and the exact minimizer at zero tolerance. The excess is 0.0004. The project's reference comparison allowance at 22000 is 0.000122. The observed dependency on objective scaling makes a regression easy to write. This experiment does not locate the internal cutoff or establish that fleet has the same root cause.

The reported general SCIP result and numerical bound are both 22000.0004. At that objective, the adapter's comparison of its backend's primal and dual values with the returned original objective allows approximately 2.2001 in objective units. Tightening that comparison alone cannot detect a backend whose own values agree on the worse point. [Result acceptance](../src/scip_mip.c) · [Saved run evidence](../results/upstream-correctness-review-v1/verification.json)

Reproduce using a fresh temporary directory with the feature-enabled binary. These commands exercise the saved fixture. These commands reproduce the historical behavior with the v23 binary. The current feature-enabled CTest suite exercises the retained model and its equivalent formulations.

```sh
out=$(mktemp -d)
model=tests/fixtures/near_tied_costs_order.mps
build-v23/milpsolve "$model" --mip --solution "$out/lab.sol"
build-v23/milpcheck --model "$model" --solution "$out/lab.sol" --exact
build-v23/milpcheck --model "$model" --solution tests/fixtures/near_tied_costs_order.sol --exact
build-v23/milpsolve "$model" --mip --no-scip
```

## Upstream patterns to adopt

| Priority | Concrete upstream example | What to require in Lab |
|---|---|---|
| 1 | SCIP's test-set evaluator distinguishes best-known feasible values and known optima, and can flag a claimed dual bound that cuts off a known solution. | Index retained checker-valid original-MIP witnesses by model hash, objective sense and checking policy. Across versions, seeds and backends, flag any claimed bound excluding the best witness beyond the declared objective tolerance, an optimal claim worse than that witness, or infeasibility despite that witness. For maximization reverse the inequality. Include censored runs that returned witnesses. Keep the flag distinct from independent proof, particularly with positive feasibility tolerances. |
| 2 | HiGHS `MIP-max-offset-test` solves a model, adds an offset, and then negates objective and sense. Its issue-2173 and issue-2585 tests request zero MIP gap and assert a stored objective as well as optimal status. | Build a small exact-oracle corpus around close costs times large activities. Transform objective units, offsets, sense and column order. Map results back to original units before applying the tolerance. Exercise the same family on native, contracted and general SCIP routes where applicable, and assert the route for each intended case. |
| 3 | HiGHS `presolve-solve-postsolve-mip` solves the exported presolved model in a second solver object. Its MIP postsolve checks primal feasibility and objective reconstruction and expects unknown model status when only the solution is supplied. OR-Tools CP-SAT has explicit presolve-on/off feasibility regressions and debug checks for a presolver breaking a hint. | Add a round trip at the Lab-to-SCIP boundary. Re-read the surviving coefficients, objective, row sides, bounds and integrality; reconstruct the original witness and objective. For custom reductions require the presolver to fire in a positive case and abstain in a mathematically adjacent negative case. A surviving witness alone cannot establish an optimal status. |
| 4 | HiGHS `HighsDebugSol` checks cuts, bound changes and pruned nodes against a debug solution. SCIP provides its own debug-solution mode. | In a diagnostic run, supply a known optimum for small cases or a better retained fleet witness. Follow only its compatible subtree. Record the first reduction, cut or bound that excludes it without sufficient objective justification. The reference point is diagnostic data and must stay outside production search. |
| 5 | The pinned SCIP source ships LP-interface tests for feasibility flags, ray availability, basis operations, and small bound changes. | Run the applicable upstream assertions against each patched QSopt_ex build. Cover primal and dual simplex after edits, finite limits, singular bases and fresh versus restored solves. Explicitly record unsupported operations. Current dependency builds set `BUILD_TESTING=OFF`; local QSopt_ex tests exercise selected regressions. |
| 6 | MIP-DD reduces a model and its solver settings using a fixed reference witness. SCIP 10 also supports an exact-solving and independently checkable certificate path. | Reduce fleet 23:00 while preserving both the better witness and the contradictory claim on the reduced model. For small held-out MILPs, use enumeration, rational bounds, or a separately checked proof. SCIP 10's exact mode is an external oracle option. Its certificate scope depends on presolve and cut settings. |

Sources for the comparison

- [SCIP automated testing](https://scipopt.org/doc-10.0.0/html/TEST.php) and [debug solution](https://www.scipopt.org/doc/html/DEBUG.php).
- HiGHS v1.15.1 [MIP tests](https://github.com/ERGO-Code/HiGHS/blob/v1.15.1/check/TestMipSolver.cpp), [presolve tests](https://github.com/ERGO-Code/HiGHS/blob/v1.15.1/check/TestPresolve.cpp) and [debug solution](https://github.com/ERGO-Code/HiGHS/blob/v1.15.1/highs/mip/HighsDebugSol.cpp).
- [OR-Tools CP-SAT solver regressions](https://github.com/google/or-tools/blob/stable/ortools/sat/cp_model_solver_test.cc), particularly `PresolveChangesFeasibility` and `PresolveChangesFeasibility2`.
- [Pinned SCIP LP-interface suite](https://github.com/scipopt/scip/tree/a94f778ac864f6c488e64f10f683c9d05253aa79/tests/src/lpi) and [Lab's dependency builder](../scripts/build_scip_integer.py).
- [MIP-DD reduction procedure](https://github.com/scipopt/MIP-DD), [SCIP exact mode](https://scipopt.org/doc/html/EXACT.php) and [certificate scope](https://www.scipopt.org/doc/html/FAQ.php).

## What Lab already covers

The native presolver has 1,100 enumerated MILP cases. The solver has 300 tiny exact MILP oracles, additional scaled LP/MILP tests, LP/Farkas certificate checks and GMI cut replay. The SCIP suite checks an exact oracle on 16 padded generated cases. Each case now names its expected final route. Fourteen require general SCIP and two require a presolve contradiction. The fleet release check compares claims and raw bounds against retained checked schedules. The [correctness plan](SOLVER-CORRECTNESS-PLAN.md) records the candidate, the pricing regression it exposed and the original-model replay results. [Presolve contracts](../tests/test_presolve.c) · [MILP oracle](../tests/test_mip_solver.py) · [Scaling contracts](../tests/test_numerics.py) · [SCIP contracts](../tests/test_scip_mip.py)

Before adding comparisons across runs, replaying the saved fleet 23:00 Lab run gave `reference_objective_agreement=false`, `reference_contradiction=false`, and `reported_bound=null`. The full benchmark classifier refuses completion for this best-known reference. The more informative failure is that the better checked HiGHS witness beats Lab's optimal claim by 0.000400022916529. Both schedules passed the original-model checker with row tolerance 1e-10 and zero integrality tolerance. [Shared evaluator](../scripts/run_evaluation.py) · [Fleet replay](../results/fleet-repair-60s-v1/fleet-2300-tolerance-check.json)

As a spot check, the local v23 CLI matched HiGHS's expected objectives on issue-2173 and issue-2204, with valid witnesses. It returned the expected infeasible status on issue-2402. The separate flugpl run exhausted its five-second budget. These observations cover four imported MPS inputs, not execution of the entire upstream test harness.


## Implemented checks

The general SCIP adapter tests whether its configured comparison epsilon merges distinct objective coefficients before importing the model. It repeats a lossy import using the existing strict threshold. If that threshold still merges distinct numerical prices, it returns control to the fallback path. This comparison covers prices representable as distinct doubles. It cannot certify a branch-and-bound tree.

The SCIP contract now solves 32 variants of the retained one-row problem across the general and native routes. It changes column order, objective units, offset and minimization/maximization, checks the actual dispatch, checks the solution at zero tolerance, and converts both the reported objective and bound back to the original units. All must equal the analytical optimum of 22000. [Production-route regression](../tests/test_scip_mip.py)

The report and completion gate now index policy-checked original-MIP incumbents across the supplied runs. A time-limited solver can supply the best witness. A claimed optimal objective or numerical dual bound that excludes a better witness beyond the reference allowance produces a conflict. An infeasible claim on a model with a checked witness also produces a conflict. The gate excludes conflicted attempts from completion and timing credit. Each conflict identifies the witness, tolerances and available hashes. Pass additional archived JSONL files with repeated `--witness-runs PATH` arguments to `report.py`, `benchmark_gate.py` or `revalidate.py`. Archived witnesses do not add timed attempts to the campaign. New Lab runner records retain the reported numerical bound from the solver output. [Shared policy](../scripts/run_evaluation.py)


The new policy finds two contradictions in retained fleet evidence. At 23:00 the archived v19 optimum of 22021.524649683968 is worse than the checked HiGHS schedule at 22021.52424966105. At 23:30 the archived v20 optimum and reported lower bound of 20259.870667388142 both exclude a checked Lab schedule at 20259.870267366074. The second case is visible even though the published fleet reference is weaker than both schedules.

In a fresh 60-second Lab run per snapshot with the v25 dependency lock, 23:00 reached 22021.52424966103. The other five runs reached the time limit. All six schedules passed the standard independent checker. Four also passed at a linear tolerance of 1e-10 and zero integrality tolerance. The 22:30 and 22:45 schedules had maximum row residuals of 1.6255e-10. None of the six runs supplied an independent optimality certificate. These runs overlap with contract tests on the same host. [Fleet assessments and input hashes](../results/checked-witness-fleet-v1/fleet-audit.json) · [Stricter witness checks](../results/checked-witness-fleet-v1/fleet-strict-checks.json)

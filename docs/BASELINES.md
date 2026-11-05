# Milestone 1 baseline — 2026-11-05

The C17 benchmark runner and independent rational checker are implemented.
All 240 official MIPLIB 2017 benchmark-v2 instances are imported and SHA-256
locked, with reference metadata v37. The final checker accepts the model grammar
of all 240 instances. The project does not yet implement a new optimizer.

## Completed experiments

| Solver | Version | Full-suite attempts | Passed checker | Failed checker policy | Reached reference |
| --- | --- | ---: | ---: | ---: | ---: |
| HiGHS | 1.15.1 | 240 | 110 | 0 | 12 |
| SCIP | 10.1.0 / SoPlex 8.1.0 | 240 | 109 | 1 | 8 |
| Cbc | 2.10.13, local Apple Clang build | 240 | 78 | 6 | 6 |

The [full report](../results/miplib-2s-seed1/validated.html) contains all 720
solver–instance pairs: no missing or duplicate runs. Of these, 297 pass the
checker, seven fail its fixed acceptance policy, 369 supply no candidate, and
47 have no retained candidate available for validation. There are no final
checker/parser errors and no contradictions against the pinned reference data.
Reference agreement compares objectives; it does not certify optimality.

Each run requested one thread, seed 1, a two-second internal solver limit and
zero absolute/relative MIP gap. The runner executes serially and rotates solver
order by instance. Its watchdog allows 15 seconds beyond the internal limit,
then sends SIGTERM and, if needed, SIGKILL one second later. Internal solver
clocks account for setup differently; wall times include the whole process.

This is a pipeline calibration baseline, not a solver performance ranking.
Runs used a shared Apple M2 Pro workstation (10 logical CPUs, 16 GiB RAM,
Darwin 24.6.0). Brief build/test jobs overlapped parts of the campaign. A single
seed, short limits, parsing costs, watchdog outcomes and shared-machine load
prevent a fair speed comparison. Longer dedicated multi-seed runs are still
needed before making performance claims.

A separate [smoke campaign](../results/smoke-5s-seed1/validated-final.html)
completed 12 instances × three solvers at five-second internal limits: 29
candidates pass, one fails the checker policy, five have no solution, and one
is unavailable. No outcomes were removed from either campaign.

## Independent acceptance policy

The checker uses GMP exact rationals for model coefficients, supplied decimal
values and all residual calculations. It reads the original MPS model and does
not link a solver or reuse a solver's feasibility API. Acceptance uses absolute
row/bound tolerance 1e-4, integrality tolerance 1e-5, and claimed-objective
consistency 1e-4. `--exact` sets all three tolerances to zero.

The seven full-suite rejections are retained with their raw solutions and
violation evidence:

| Instance | Solver | Reason |
| --- | --- | --- |
| cbs-cta | Cbc | Claimed-objective consistency only; no row/bound/integer violation |
| germanrr | Cbc | Row residual and claimed-objective consistency |
| leo1 | Cbc | Row residual and claimed-objective consistency |
| leo2 | Cbc | Row residual and claimed-objective consistency |
| mas76 | SCIP | One row residual |
| neos-1122047 | Cbc | Row residuals |
| sct2 | Cbc | Row residual and claimed-objective consistency |

These are failures of the recorded common policy, not established solver bugs.
Solver scaling, native tolerance rules and solution-file precision differ.
For example, the mas76 residual is approximately 0.1075 on a right-hand side
of approximately 157345 (relative residual approximately 6.8e-7); the independent
absolute policy rejects it. An additional Python Fraction calculation confirms
this residual. The checker establishes primal feasibility and objective
consistency under its policy; it never certifies optimality or infeasibility.

## Provenance and validation

Original timings use Git tag `evaluator-v1` (commit `2a86842`). Their original
`run.json`, logs, solutions and timings remain untouched. The experiment folders
contain an archive and file-hash lock of that source revision. Final validation
uses the hardened C evaluator from `milestone-1-evaluator` (`1a40b18`) and the
versioned validation driver from `4e9c025`. All 720 retained run records were
rechecked with the same final checker identity:

```text
original runner SHA-256:
736f411821b93e81e120366839e9942349a586a5537a2ce5fe4b03eeec7f31df
final checker SHA-256:
0e5008fc201485219ab9e1d6c6b657ec1ce3dac6273862415a29133dddaa5163
validation driver SHA-256:
0642d7f8fc2b5740758482228fdd16c65a574361131ad6c749ed46bccbe8b8ec
```

Each final report links the original run record/log and the new per-run
validation record. `final-validation.jsonl` and its metadata file preserve the
new results without replacing earlier validation versions. Coverage and import
audits are saved under `results/miplib-2s-seed1/coverage.json` and
`results/import-audit-final/summary.json`.

Release and AddressSanitizer/UndefinedBehaviorSanitizer CTest suites pass,
including 120 randomized exact-arithmetic oracle cases, malformed-input and
format regressions, watchdog/cancellation behavior and reporting contracts.
Six additional runs with the real native engines pass known minimization and
maximization adapter checks. Cbc's source build also completed upstream tests.
The repository includes a CI workflow; no remote CI execution is claimed.

## Reproducing and extending the baseline

Follow [the README](../README.md) to build, fetch the locked tools/dataset and
run tests. Use a fresh results directory for each experiment. Preserve solver
settings, model/tool/checker hashes, failures and retained solutions. The
[evaluation contract](PROTOCOL.md) documents supported formats and limitations.

The source archive contains project code, configuration, locks, tests and these
documents. Imported MIPLIB data, downloaded binaries, native build trees and
generated result folders stay outside Git and the source archive. Upstream
licenses apply to those artifacts; unresolved per-instance redistribution
rights are recorded in the dataset manifest.

For the next milestone, build the initial C LP/MILP implementation behind the
same runner interface. Keep the evaluator in a separately controlled checkout
before autonomous optimization experiments, and confirm any improvement on
longer budgets, several seeds and prospective held-out model families.

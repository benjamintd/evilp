# SCIP v9: QSopt_ex tolerances, root budgets and a component profile

This checkpoint fixes a root-LP rejection in the SCIP/QSopt_ex engine, caps
the binary-SAT root budget, and records where SCIP spends its time across all
240 locked MIPLIB instances. Evidence is in
[results/scip-v9-60s-v1](../results/scip-v9-60s-v1/).

## Method

All runs use a 60-second limit, one thread and seed 1 on an arm64 Mac, five
instances in parallel. A Lab *win* is a reported optimal or infeasible outcome
with a checker-valid original-model witness (for optimal), no reference
contradiction, and a wall time below HiGHS 1.15.1's — or HiGHS not finishing
within the same limit. HiGHS statuses are taken as reported; the scorer does
not check them against references (see `binkar10_1` below). Runners:
[screen_v6.py](../results/scip-v9-60s-v1/screen_v6.py) and
[fleet_compare.py](../results/scip-v9-60s-v1/fleet_compare.py).

Run-to-run variability is large near the time limit: the same build moved
individual instances between solved and unsolved across repeated screens.
Single-instance changes are not evidence; only panel-level counts are.

## Root LP rejected after an optimal QSopt_ex solve

On the fleet models, SCIP reported `unresolved numerical troubles in LP 1`
and fell back to the pseudo bound (−2.13e7). QSopt_ex had in fact solved the
root LP to optimality, but its built-in dual feasibility tolerance (~1e-6) is
looser than SCIP's `numerics/dualfeastol` (1e-7). SCIP's post-solve check
flagged the reduced costs (`pfeas=1, dfeas=0`) and tried to resolve with a
tighter tolerance; `lpi_qso` returned `SCIP_PARAMETERUNKNOWN` for
`FEASTOL`/`DUALFEASTOL`, so the root relaxation was discarded.

The scaled-copy behaviour of QSopt_ex (scale a temporary copy, keep only the
basis, re-solve unscaled) was suspected but is not the cause. It also costs
little: the root LP takes 0.6–1.1 s on the fleet models and turning scaling
off entirely changes root time by less than a second on the tested models.
A persistent-scaling rewrite of the LP interface was therefore not pursued.

Fix: [qsopt-ex-c17-v9.patch](../config/patches/qsopt-ex-c17-v9.patch) adds
per-LP primal/dual feasibility tolerance parameters;
[scip-9.2.4-c17-v6.patch](../config/patches/scip-9.2.4-c17-v6.patch) maps
SCIP's `FEASTOL`/`DUALFEASTOL` onto them (never looser than QSopt's 1e-6).
Lock: [scip-integer-v9](../config/scip-integer-v9.lock.json). All 62 contract
tests pass.

## Binary SAT root budget

SAT closes MIPLIB instances quickly (all twelve SAT closures finish in under
ten seconds at a 60-second limit) or not at all, but could consume 85% of the
budget before SCIP started (`neos-1582420`: 35 s). The budget is now
`min(0.85, max(25 s, 0.25 × remaining))`; a 10–12 s cap lost `wachplan`
because SAT's internal descent steps are fractions of its deadline. On the 45
SAT-active instances: 15 solved (was 14), no lost solves.

## Results

| Screen (240 instances, 60 s) | Correct solves | Wins vs HiGHS | of which HiGHS unsolved |
|---|---:|---:|---:|
| HiGHS 1.15.1 | 43 reported | — | — |
| Lab, previous full screen (`lab-screen-60s-v6`) | 40 | 33 | — |
| Lab v9 | 43 | 36 | 16 |

Fleet models (user data, no locked reference; every Lab witness checker-valid):
Lab proves optimality on 2 of 6 (2300 in 43.2 s, 2330 in 36.1 s; HiGHS 1 of 6,
2330 in 57.5 s). Before the fix, fleet 2215 ended with no incumbent and the
pseudo bound. HiGHS keeps a stronger dual bound on 2215 and a better
incumbent on 2230. These pairs ran while another benchmark loaded the
machine; one rerun of 2215 alone changed Lab's incumbent from 823,953 to
22,531, so the fleet table needs an unloaded repeat.

`binkar10_1`: the locked v37 reference is 6741.38, but both HiGHS and Lab
report 6742.20 as optimal. Lab's claim is a recorded contradiction; HiGHS's
equivalent claim is not checked by the scorer. The reference is unchanged.

## Component profile

[scip-component-profile-240.txt](../results/scip-v9-60s-v1/scip-component-profile-240.txt)
sums SCIP statistics over all 240 instances (60 s each, SCIP called
directly, about 13,800 solver-seconds). Times nest (heuristics include their
own LP time), so shares are indicative:

- node dual LP solves ~6,000 s; strong branching ~2,200 s; diving/probing
  LPs ~1,200 s — LP speed is the dominant cost;
- sub-MIP heuristics stall on individual instances for most of the budget:
  `vbounds` on both drayage instances (59 s), `rens` on `glass-sc` (58 s),
  `clique` on `neos-5114902-kasavu` (54 s), `locks` on `dano3_3`/`dano3_5`
  (one cold probing LP of ~48 s);
- the aggregation separator takes 21 of 46 s on `h80x6320d`.

The 191 Lab time limits include only 11 with a final gap below 1% and 30
below 5%; 52 have no incumbent. More than 120 wins requires substantially
faster LP re-solves and stronger root bounds, not budget tuning alone.

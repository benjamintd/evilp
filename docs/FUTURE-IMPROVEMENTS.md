# MILP Lab: Future Improvements Roadmap (Performance, Stability & Capabilities)

This document synthesizes measured profiling data, full-suite MIPLIB 2017 screens,
`fleet` real-world diagnostics, and lessons from Waves 1–6 into a prioritized
engineering roadmap for future agents.

### Implementation progress

- **§6.1, elapsed time.** Lab, the C process runner, SCIP wall clocks and
  QSopt_ex simplex deadlines now use awake monotonic wall time on macOS
  with the [v14 backend](../config/scip-integer-v14.lock.json). CPU timers
  remain available for profiling. On macOS and Linux the runner records
  host suspension separately. More than five seconds invalidates the paired
  timing campaign while preserving the independently checked result.
- **§4.1, heuristic budgets.** ALNS, the scheduler, local branching and trust
  region now honor the optional child time fraction. Locks has a separate
  `heuristics/locks/maxlpiters` setting and returns control to the main
  search on iteration exhaustion. Both caps are opt-in. The feasibility pump,
  presolve polling and watchdog work in §6.2 remain open.
- **§2.1, factor storage.** The eta index and coefficient arrays now grow on
  first use and retain capacity across basis refactorizations. The existing
  logical update limit, pivot order and numeric operations are preserved.
  Caching the other LU arrays remains open.

Every idea below respects [AGENTS.md](../AGENTS.md): pure C17 numerical engine,
checker, and runner; unpresolved GMP rational witness verification; strictly
general-purpose or mathematically general structure-detecting algorithms (never
specialized on filenames, hashes, or stored answers); and fresh versioned result
directories.

---

## 1. Current Baseline & Ground Truth (as of `be556ab` + WIP `scip-v10`)

From [SCIP-V9-TOLERANCES-AND-BUDGETS.md](SCIP-V9-TOLERANCES-AND-BUDGETS.md) and
[results/scip-v9-60s-v1/](../results/scip-v9-60s-v1/) (awake 60s serial runs, 1
thread, seed 1, gap = 0, arm64 Mac, 5-way parallel):

- **240-instance MIPLIB 2017 screen (60s)**:
  - **HiGHS 1.15.1**: `43` reported solves (`42` optimal, `1` infeasible), `8` watchdogs.
  - **MILP Lab (`scip-v9`)**: `43` checker-valid non-contradictory solves (`41` optimal, `3` infeasible, minus `binkar10_1` pinned-reference discrepancy), **`36` wins vs HiGHS** (`16` where HiGHS timed out, `20` both-solved faster), `5` watchdogs.
  - **Both-solved breakdown (`27` instances)**: Lab is faster on **`20 / 27`**, slower on `7 / 27` (`neos-860300` 3.20x, `neos859080` 2.62x, `mik-250-20-75-4` 2.21x, `h80x6320d` 1.73x, `piperout-27` 1.34x, `neos17` 1.30x, `app1-1` 1.22x).
  - **Unsolved breakdown (`196` instances)**: `57` have **no incumbent**, `11` have `<1%` gap (`neos-957323` 0.00%, `neos-4647030-tutaki` 0.02%, `gmu-35-40` 0.05%, `neos-3627168-kasai` 0.11%, `gmu-35-50` 0.11%, `tr12-30` 0.13%, `unitcal_7` 0.17%, `pg5_34` 0.35%, `gen-ip002` 0.60%, `seymour1` 0.66%, `neos-827175` 0.89%), `19` have `1–5%` gap, `24` have `5–20%` gap, `85` have `>20%` gap.
  - **Zero-LP / Presolve stalls (`8` instances)**: `5` watchdogs (`neos-2075418-temuka`, `neos-3402454-bohle`, `neos-5114902-kasavu`, `square41`, `square47`) plus `k1mushroom` (73.9s), `splice1k1` (68.1s), `supportcase19` (64.3s) all end with `0` LP iterations.
- **240-instance SCIP component profile** (~13,800 total solver-seconds in [scip-component-profile-240.txt](../results/scip-v9-60s-v1/scip-component-profile-240.txt)):
  - `LP/dual LP`: **6,028.4 s** (197 instances >3s)
  - `Branching Rules/relpscost`: **2,315.3 s** (104 >3s) — of which `LP/strong branching` is **2,184.0 s** (99 >3s)
  - `LP/diving/probing LP`: **1,198.5 s** (103 >3s)
  - Sub-MIP & probing heuristics: `feaspump` **423.5 s**, `rens` **375.0 s**, `locks` **366.9 s**, `vbounds` **350.5 s**, `clique` **282.1 s**
  - Presolvers: `probing` **281.3 s**, `linear` **231.7 s**, `logicor` **159.7 s**, `knapsack` **153.0 s**, `dualsparsify` **93.4 s**
  - Separators: `gomory` **258.9 s**, `aggregation` **248.9 s**, `clique` **93.7 s**, `cut pool` **64.5 s**
  - Conflict analysis: `propagation` **194.8 s**
- **Uncommitted WIP already in working tree (audit & finish first)**:
  - `config/scip-integer-v10.lock.json` + `config/patches/scip-9.2.4-c17-v7.patch` + `src/scip_mip.c` (`heuristics/subsciptimefrac` to cap heuristic sub-SCIP stalls like `vbounds`, `rens`, `clique`).
  - `src/pair_assignment.c` / `include/pair_assignment.h` / `tests/test_pair_assignment.py` (audit carefully for generality vs `AGENTS.md` before keeping).
  - `src/model.c`, `src/exact.c`, `tests/test_exact_model.c`.

> **Strategic reality check on `>120 / 240` wins vs HiGHS**: At 60 seconds serial,
> HiGHS 1.15.1 itself only solves 43/240 instances. Crossing `>120 / 240` wins
> requires **both** (a) major general-purpose gains in LP re-solve speed,
> strong branching, presolve, root cuts, heuristics, and SAT/PB optimization, and
> (b) evaluating at a budget (e.g. 180s–600s) where >120 MIPLIB instances are
> actually tractable in serial — after first establishing clear dominance at 60s.

---

## 2. Track 1 — LP Engine Performance (`QSopt_ex` & Native Dual Simplex)

**~68% of all time across the 240 MIPLIB instances (`9,411 s / 13,800 s`) is inside LP solves** (`dual LP` + `strong branching` + `diving/probing LP`). Even where Lab already needs *fewer* LP iterations than HiGHS (`app1-1`: 83k vs 123k; `h80x6320d`: 36.5k vs 39.9k), Lab loses on wall time because `QSopt_ex` is ~1.8–2.5x slower per iteration on medium/large bases and drops to `168–258 it/s` on huge/dense bases (`neos-631710`, `supportcase10`, `seymour`, `opm2-z10-s4`, `s100`).

### 2.1 Capacity-Cached `factor_work` Across Refactorizations (`QSopt_ex`)
- **Evidence**: In `qsopt_ex/basis.c:ILLbasis_factor`, every ~100 pivots (every refactorization) calls `ILLfactor_free_factor_work(lp->f)` followed immediately by `ILLfactor_create_factor_work(lp->f, lp->O->nrows)`. On `binkar10_1` and `mik-250-20-75-4`, `dbl_ILLfactor` + `free_factor_work` (`madvise`) + `create_factor_work` (`__bzero`) is the #1–#2 CPU hotspot even after `scip-v8` cached `internal_lpinfo`.
- **Idea**: Add `dim_cap` and `*_space_cap` to `factor_work` (`factor.h`). In `ILLbasis_factor`, if `lp->f != NULL && lp->f->dim_cap >= dim`, reset metadata (`nzcnt`, `etacnt`, `rperm`/`cperm`, `work_indx`, `ur_inf`/`uc_inf` headers) in-place and reuse the existing `urindx`/`urcoef`/`ucindx`/`uccoef`/`lcindx`/`lccoef`/`lrindx`/`lrcoef`/`erindx`/`ercoef` buffers when they fit, only freeing on `ILLbasis_free_basisinfo`.
- **Expected impact**: Zero `madvise`/`__bzero` during every 100-pivot refactorization; 5–15% pure LP speedup with **bit-for-bit identical numerical trajectory**.

### 2.2 Cache-Compact, Bit-for-Bit Order-Preserving `btranu3` / `ftranu3` (`QSopt_ex`)
- **Evidence**: On `h80x6320d` (`6,332` rows), `btranu3_process2` + `btranu3_delay2` + `ftranu3_process2` + `ftranu3_delay2` consume **~38% of total CPU**. Why:
  1. `ur_info` is 32 bytes/row (`EGLPNUM_TYPE max` + 6 ints = 202 KB for 6,332 rows, thrashing L1 cache), and `uc_info` is 20 bytes/col.
  2. Every node performs **double indirection** `f->rperm[f->crank[r]]` or `f->cperm[f->rrank[c]]` in *both* `delay2` and `process2` (`4x` dependent loads per node). Yet `shift_permutations()` in `factor.c` shifts `rperm` and `cperm` in exact lockstep, so `row_to_col[r] = cperm[rrank[r]]` and `col_to_row[c] = rperm[crank[c]]` are **invariant across `ILLfactor_update`**!
- **Idea**:
  - Maintain compact flat arrays built once per `ILLfactor`: `int *u_delay` (4 bytes/entry — **8x denser in L1 than `ur_inf`**), `int *row_to_col`, `int *col_to_row`, and compact 8-byte `(beg, nzcnt)` pairs for U rows/cols.
  - Add a **1-nonzero trivial diagonal fast path** (`anzcnt == 1 && nzcnt == 1`: single division, zero `delay` or `work` touch).
  - **Keep exact Kahn topological order** (or verify any dense switch only triggers at high fill-in `> 40–50%`, NOT `12%` — see Section 9 dead ends!).
  - In `dbl` mode, test storing reciprocal U diagonals `inv_u_diag[i] = 1.0 / coef[0]` (updated in `eliminate_row`/`update_column`) to replace 10-cycle ARM64 `fdiv` with 3-cycle `fmul`.

### 2.3 True Zero-Copy Strong Branching Snapshots & 1-Pivot Lookahead (`2,184 s` Hotspot!)
- **Evidence**: `LP/strong branching` alone burns **2,184 s across 99 instances** (`mushroom-best`: 54s, `supportcase40`: 51s, `neos-3024952-loue`: 47s).
  - Currently in `qsopt_ex/lib.c:ILLlib_strongbranch_cached`, `base_ready` is **invalidated after any single pivot** (`lp->basisid != bid_before`), forcing a full basis reload + LU refactorization (`ILLbasis_load` + `ILLfactor`) on the next branch probe!
  - Furthermore, `lpi_qso.c` calls `SCIPlpiStrongbranchFrac` one candidate at a time in many SCIP code paths and re-solves the parent LP at `EndStrongbranch`.
- **Idea**:
  1. **Dual-bound 1-pivot analytical bound (Driebeek–Tomlin / BFRT dual step)** inside `ILLlib_strongbranch_cached` before invoking full dual simplex: for a basic variable `x_j = v`, compute the exact bound change from the first dual ratio test (one BTRAN + one `compute_zA` + BFRT scan) with **zero basis mutation**. For `itlim == 1` or infeasible/pruned branches, return immediately without touching `lp->f`!
  2. **Rollback-able Forrest-Tomlin eta stack** during strong branching: save `etacnt`, `er_freebeg`, `baz/nbaz/vstat/vindex` delta list, and `rownorms` deltas when a branch takes `<= 10` pivots; undo those `<= 10` eta updates (or restore saved `factor_work` snapshot!) instead of running a full `O(nnz)` Markowitz `ILLfactor` refactorization between every strong-branch candidate!
  3. **Adaptive strong-branching budget**: test tuning SCIP's `branching/relpscost/maxlookahead`, `sbiterquot`, `sbiterofs`, and `inititer` specifically for `QSopt_ex`'s iteration cost (e.g., capping `inititer` and `maxlookahead` on large/slow LPs where `iter/s < 1000`).

### 2.4 Dense-Basis Pricing & Hyper-Sparse `ftrane2` / `btran` Improvements
- **Evidence**: The 10 slowest LP instances (`168–258 it/s` on `neos-631710`, `supportcase10`, `seymour`, `opm2-z10-s4`, `s100`, `roi5alpha10n8`) have large, relatively dense bases where:
  - DSE (`QS_PRICE_DSTEEP`) FTRAN/BTRAN updates are expensive,
  - `ILLfactor_ftrane2` loops over all `etacnt` eta rows scanning every nonzero even when `work_coef` is sparse (128 samples on `h80x6320d`),
  - `ILLfct_compute_zA` scatters through dense rows.
- **Ideas**:
  - **Adaptive DSE → Devex switch** when average BTRAN/FTRAN density exceeds ~20% (or basis dimension > 15,000 with slow iter/s): Devex skips the extra FTRAN per pivot needed to maintain exact steepest-edge row norms, cutting per-pivot triangular solves by ~33%!
  - **Sparse-aware `ILLfactor_ftrane2` / `btrane2`**: maintain a tiny `uint64_t` bitmask bloom filter (64-bit hash of touched rows in `work_coef` / eta row support) to skip eta rows with guaranteed zero inner product in 2 CPU cycles!
  - **Adaptive `etamax`**: `QSopt_ex` hardcodes `f->etamax = 100`. On sparse/small bases (`dim < 500`), `etamax = 150–200` cuts refactorization overhead (`dbl_ILLfactor` was #1 on `binkar10_1`); on dense bases (`nnz_factor > 5 * nnz_orig`), `etamax = 50–70` keeps `ftrane2`/`btrane2` fast.

### 2.5 Longer-Term Architectural Leap: Native High-Speed C17 Dual Simplex Backend
- **Why**: `QSopt_ex` was designed in 2001–2006 primarily for exact rational verification, not modern MIP re-optimization. Even heavily patched, it lacks:
  - **SIMD-friendly column/row-wise pricing** (HiGHS `HDualRow` / `HMatrix` split),
  - **Multi-choice / minor-iteration PAMI (Parallel/Suboptimization Dual Simplex)** (Huangfu & Hall 2015),
  - **Modern Forrest-Tomlin / Bartels-Golub / Product-Form-of-Inverse `HFactor`** with adaptive hyper-sparsity from the ground up,
  - **Integrated persistent Ruiz + geometric scaling** baked directly into the solver with clean primal/dual residual polishing.
- **Idea**: Evolve MILP Lab's own pure-C17 LP engine (`src/lp.c`, `src/factor.c`, `src/pricing.c`, `src/pivot.c`) into a full **SCIP LPI plugin (`lpi_lab.c`)** implementing Hall–McKinnon hyper-sparse dual simplex + DSE + BFRT + bound-shifting perturbation + Forrest-Tomlin updates. That would remove the `QSopt_ex` ceiling permanently *and* unify the codebase under MIT pure C17.

---

## 3. Track 2 — Presolve, Cuts & Root Bound Strengthening

In `lab-v9-full.jsonl`, HiGHS closes **13 instances at Node 1** (`drayage-25-23` 1.24s, `neos-860300` 15.3s, `piperout-08` 28.9s, `piperout-27` 40.4s, `neos-4413714-turia` 34.9s, `physiciansched6-2` 26.6s, `neos-960392` 26.7s, `neos-950242` 30.9s, `mzzv42z` 52.8s, `neos-933966` 54.2s, etc.) where SCIP enters deep trees or stalls. Because SCIP builds with `PAPILO=OFF` (PaPILO is C++), closing this presolve + root-loop gap is high-leverage.

### 3.1 Pure-C17 PaPILO-Class Presolver Plugins (Extending `src/scip_presol_parallel.c`)
Following the clean SCIP presolver plugin pattern in [src/scip_presol_parallel.c](../src/scip_presol_parallel.c):
1. **Parallel-Column & Column-Dominance Merging (`scip_presol_parcol.c`)**:
   - Hash columns by normalized row-index support + quantized coefficient ratios (exact mirror of `scip_presol_parallel.c`!).
   - For two parallel columns `x_1, x_2` (`A_2 = alpha * A_1`, same integrality/sign structure):
     - If continuous (or integer with `alpha = ±1`), merge into a single aggregated variable `z = x_1 + alpha * x_2` when `c_2 = alpha * c_1`, or tighten/fix the dominated column via dual dominance when `c_2 >= alpha * c_1`!
   - Especially impactful on network-design, scheduling, and fleet/routing formulations with duplicate arcs/shifts.
2. **Bounded-Fill-In Doubleton & Implied-Free Substitution (`scip_presol_subst.c`)**:
   - PaPILO's single highest-reduction pass on models like `drayage-25-23` (`4,630r → 732r` in 0.07s!).
   - Detect equality rows `a_1 x_1 + a_2 x_2 = b` (or short equations `len <= 5`) where a continuous or implied-free variable `x_1` has `col_nnz <= 4` and Markowitz fill-in `(row_len - 1) * (col_nnz - 1) <= 12` and coefficient dynamism `<= 1e3`; substitute `x_1` via `SCIPaggregateVars` / `SCIPmultiaggregateVar` in strict work-bounded passes!
3. **Fast Work-Bounded Coefficient Strengthening & Implied-Bound Lifting (`scip_presol_coefftight.c`)**:
   - Using `SCIPmatrixCreate`, tighten big-M coefficients `a_{ij}` on binary variables `x_j` against row min/max activity (PaPILO's `CoeffTightening`) in a single linear `O(nnz)` scan — zero quadratic blowup, huge tightening of LP relaxations on big-M fleet/scheduling/network MIPs!
4. **Fast 1-Pass Binary Clique / Implication Probing**:
   - SCIP's `presol_probing` burns `281 s` across 20 instances (33s on `tbfp-network`, 33s on `s100`). Implementing a lightweight, strictly operation-counted (`max_ops = 5e6`) implication-graph transitive closure + 2-literal fixing pass gives most of PaPILO's probing benefit in `<0.2 s` without stalls.

### 3.2 Root Restarts + Reduced-Cost Fixing + Root LNS Synergy
- **Evidence**: On `drayage-25-23`, HiGHS solves in `1.24s` at Node 1 via: `Feasibility Jump` → `Root Cuts` → `Sub-MIP (L)` (1.95% gap) → **Restart #1** (18.3% reduced-cost fixings!) → `Root Cuts` → `Sub-MIP (L)` (0.60% gap) → **Restart #2** (11.1% fixings!) → `Root Cuts` → `Sub-MIP (L)` → **0.00% gap at Node 1!**
- **Ideas**:
  - Enable `propagating/redcost/useimplics = TRUE` (currently `FALSE` by default in SCIP's `prop_redcost.c:82`) — strengthens root reduced-cost fixing using binary implications, triggering more root fixings and earlier restarts!
  - Enable `propagating/probing/maxruns = -1` (currently `DEFAULT_MAXRUNS = 1`, so probing **never** runs again after a root restart!) with a tight work cap so post-restart presolving cascades!
  - Allow fast capped root LNS (`rens` / `dins` / `alns` with `subsciptimefrac` from `scip-v10`!) to run on **subsequent root restarts** once a tight primal bound is found, closing the `cuts → primal → reduced-cost fixing → restart → tighter cuts` loop!

### 3.3 Cutting-Plane Speedups Without Losing Cut Strength
- **Evidence**: `sepa_gomory` (`258.9 s`, 22 instances >3s) and `sepa_aggregation` (`248.9 s`, 23 instances >3s; `21 s` of `46 s` on `h80x6320d`!) are the two heaviest separators.
- **Crucial lesson**: Naively capping `separating/aggregation/maxrounds = 1` regressed `binkar10_1` by 2.2x (`19.5s → 42.9s`) and `h80x6320d` by 1.64x — the cuts themselves are essential!
- **Ideas (make separation *faster* and *stronger*, not weaker)**:
  1. **Profile-driven C17 optimization of `sepa_aggregation.c` & `cuts.c` inside SCIP patch**:
     - `sample` showed `computeMIREfficacy`, `SCIPcalcFlowCover`, `SCIPcalcIntegralScalar`, and `SCIProwGetLPEfficacy` dominate.
     - Cache row Euclidean norms and LP-solution dense array `solval[col_idx]` once per separation round instead of calling `SCIPvarGetSol` (`36` samples on `h80x6320d`) inside inner loops!
     - Skip rows with density `> maxrowdensity` or slack `> maxslack` *before* allocating aggregation buffers; early-abort `SCIPcalcFlowCover` when row has no continuous/variable-bound structure.
  2. **Sparse `SCIPlpiGetBInvARow` in `lpi_qso.c` (`sepa_gomory`: `258.9 s`)**:
     - Right now `lpi_qso.c:3156` forces `*ninds = -1` (dense!) and allocates/copies `ncols` doubles on every Gomory tableau row, even when the tableau row has only 15 nonzeros (`sorrell3`: 27s in `sepa_gomory`!).
     - Exposing QSopt_ex's sparse `svector` from `ILLbasis_tableau_row` directly through `SCIPlpiGetBInvARow` with real `inds` + `ninds` will speed up `sepa_gomory` by **3–10x** on sparse models!
  3. **New general cut separator plugins**:
     - **Lifted Knapsack Cover / Extended Cover** cuts with fast exact DP lifting on LP-tight binary knapsack rows (`scip_sepa_liftedcover.c`) — complements SCIP's default knapsack separator when `knapsack/disaggregation = FALSE` on large models.
     - **Conflict-Graph / Odd-Cycle & Clique-Cut Separator** using fast incremental BFS on the fractional LP support graph (`0 < x_j < 1`), much faster than `sepa_clique`'s 30s stall on `neos-4722843-widden`.

---

## 4. Track 3 — Primal Heuristics & Breaking the 57-Instance "No Incumbent" Wall

**57 / 240 MIPLIB instances (`24%`) hit 60s with NO incumbent at all**, while SCIP's built-in sub-MIP and probing heuristics burn **`1,798 s`** (`feaspump` 424s, `rens` 375s, `locks` 367s, `vbounds` 351s, `clique` 282s) stalling for 30–60s on single instances.

### 4.1 Hard Time & Effort Caps on Sub-MIP and Probing Heuristics
- **Finish & validate `scip-v10` (`heuristics/subsciptimefrac`)**:
  - Cap each sub-SCIP heuristic (`vbounds`, `rens`, `clique`, `rins`, `alns`, `crossover`, `dins`) to e.g. `min(3–5 s, 0.12 * remaining)` (with adaptive scaling: tighter on root when LP hasn't even solved yet, slightly larger once a gap exists).
  - Directly cures: `vbounds` 59s stall on `drayage-25-23` & `drayage-100-23`, `rens` 58s on `glass-sc`, `clique` 54s on `neos-5114902-kasavu` & 50s on `neos-4532248-waihi`!
- **Cap `heur_locks` and `heur_feaspump` probing/LP work**:
  - `heur_locks` burns **60 s on `dano3_3`**, **59 s on `dano3_5`**, and **45 s on `proteindesign122trx11p8`** inside a single uncapped probing LP!
  - Set `heuristics/locks/minfixingrate` / LP iteration limit inside `heur_locks.c` (or time-cap probing LP calls in `lpi_qso.c` / `heur_locks.c` to `min(2 s, 0.08 * remaining)`).
  - Cap `heuristics/feaspump/maxlpiterquot` / per-call time (`423.5 s` total across 34 instances >3s!).

### 4.2 Feasibility Jump v2 (`src/scip_feasjump.c` — Currently Only `78.5 s` Across All 240 Instances!)
Our pure-C17 Feasibility Jump (`src/scip_feasjump.c`) is extremely cheap (only 1 instance >3s across the entire 240-instance suite!) and already finds early incumbents (e.g. `drayage-25-23` in 0.15s). Upgrading it to state-of-the-art **FJ 2.0 (Luteberget & Sartor 2023 / HiGHS-style)** is one of the highest-ROI ways to break the 57-instance no-incumbent wall:
1. **LP-Assisted Continuous Repair ("FJ + LP Polish")**:
   - Currently FJ tries to jump continuous variables individually, which struggles on coupled continuous balance equations (flow conservation, fleet/energy balance, multi-period inventory — exactly where many of the 57 no-incumbent instances live!).
   - **Upgrade**: Run FJ primarily on **integer/binary variables** (using constraint slack projected onto integers or relaxed continuous bounds); once integer violations drop to zero (or low), fix integers in a tiny SCIP probing LP (with a 0.2s / 2000-iter cap) to solve all continuous variables *exactly*!
2. **2-Variable Swap / Compound Moves for Equality Constraints**:
   - On set-partitioning (`sum x_i = 1`), assignment, and flow equalities, flipping a single binary variable ALWAYS creates a new violation.
   - When stuck at a local minimum on a violated equality row `r`, evaluate **pair swaps** `(j_1 -> v_1, j_2 -> v_2)` among variables sharing `r` (especially `0 <-> 1` swaps in `sum x = k` rows!).
3. **Periodic Post-Cut & In-Tree Execution (`freq = 15`, light 0.15s budget)**:
   - Currently `freq = 0` (root only). Running a fast 0.15s LP-seeded FJ warm-started from the current node LP solution (or perturbed incumbent) whenever `SCIPgetNSols(scip) == 0` or gap hasn't improved for 50 nodes dramatically improves time-to-first-incumbent and primal-dual integral.
4. **Weight Decay + Dynamic Tabu Tenure + Parallel-Restart Portfolio**:
   - Add smooth multiplicative weight decay (`w[r] *= 0.99` periodically, or PAWS-style decay on feasible hits) instead of only halving at `1e5`.
   - Seed diverse initial points: (a) LP relaxation rounded, (b) LP relaxation with fractional vars biased toward objective, (c) trivial/zero/midpoint, (d) constraint-lock-guided initialization.

### 4.3 Dedicated Polishing & Objective-Guided Local Search
- For the `11` near-closing (`<1%` gap) and `19` (`1–5%` gap) instances, finding the optimal incumbent 5 seconds earlier lets reduced-cost fixing and branch-and-bound pruning collapse the remaining tree.
- Add a lightweight **1-opt / 2-opt + RINS-lite polishing pass** right after any new incumbent is accepted at depth 0–5 (bounded to `<0.3 s`), using `feasjump`'s exact objective-cutoff row `c^T x <= z* - eps`!

---

## 5. Track 4 — SAT / Pseudo-Boolean, Symmetry & General Combinatorial Search

Currently `src/binary_sat.c` closes ~12–15 MIPLIB instances rapidly (all under 10s!), often where HiGHS times out at 60s. Expanding its mathematical reach is a direct multiplier on wins where LP relaxations are weak.

### 5.1 Incremental C17 CDCL + Core-Guided MaxSAT (`OLL` / `RC2`)
- **Current bottleneck**: `src/binary_sat.c` links Kissat 4.0.4 (non-incremental, no assumptions API), shells out via `fork()` per bound probe, builds monolithic adder/totalizer encodings over the entire objective (`MAX_OBJ_TERMS = 1024`), and descends top-down (`SAT-UNSAT`) from a loose upper bound.
- **Idea**: Integrate or build a pure-C17 **incremental assumption-capable CDCL solver** (IPASIR-style `add`/`assume`/`solve`/`failed`, e.g., a clean C17 CDCL engine with 2-watched literals, 1-UIP learning, LBD clause management, VSIDS, restarts, and assumption core extraction) and implement **Stratified OLL / RC2 Core-Guided MaxSAT** (Morgado et al. 2014; Ignatiev et al. 2019):
  - Assume all objective literals at their optimal value (`lb = 0`).
  - On each `UNSAT` core `C`, raise the **exact dual lower bound** `lb += w_min`, subtract `w_min` from literals in `C`, and add a **tiny incremental totalizer** over *only* `C` (typically `|C| = 2..15` literals instead of 5,000!).
  - **Raises `MAX_OBJ_TERMS` from `1024` to `50,000+`!**
  - **Dual-bound handoff to SCIP**: Even if RC2 hits its time cap before closing, its proven exact lower bound `lb` can be passed into SCIP as a **valid global objective cutoff/lower bound** (`c^T x >= lb`) AND its best model passed as an incumbent!
  - Target classes: weighted set covering/packing/partitioning, scheduling, timetabling, coloring, graph/clique models (`cod105`, `sorrell3`, `graph20-20-1rand`, `neos-631710`, `supportcase10`, `peg-solitaire-a3`, `highschool1-aigio`, `academictimetablesmall`, `neos-950242`).

### 5.2 LP-Bound-Seeded Lower Bound (`pb_lo`) for Existing PB Descent
- **Quick win even before incremental CDCL**: In `src/binary_sat.c` (`pb_mode`), `pb_lo` currently starts at `0` (or `1` if the 0-cost face is UNSAT).
- Running a fast capped (<0.3 s) root LP relaxation before entering PB descent gives an immediate valid lower bound `pb_lo = ceil((z_LP - offset - tol) / scaled_gcd)`!
- On binary optimization models where the LP gap is 0–3 units, `pb_hi` jumps straight to `pb_lo` and closes in 1–2 SAT calls instead of binary-searching across a `[0, 500]` range!

### 5.3 Expanding General Cardinality / Pseudo-Boolean Encodings in `src/binary_sat.c`
- Support **Generalized Totalizer (GTE)** (Joshi, Martins, Manquinho 2015) for weighted PB constraints with few distinct weights (common in scheduling/resource-allocation MIPs where weights are `{1, 2, 4}` or small integers), avoiding dense adder circuits.
- Support **At-Most-One (AMO) Commander / Bimander / Product encodings** on large clique/set-packing rows (`m > 16`) for fewer auxiliary variables and stronger unit propagation.
- **Exact PB Constraint Propagation / Probing** before CNF emission: simplify fixed literals, detect implied cliques, and eliminate redundant constraints so `base_cnf` stays well below `MAX_CLAUSES = 3.5M`.

### 5.4 Upgrading Pure-C17 Symmetry Detection (`src/scip_symmetry.c`)
- Currently `src/scip_symmetry.c` caps graphs at `SYM_MAX_VERTICES = 250,000`, `SYM_MAX_EDGES = 600,000`, `SYM_MAX_TREE_NODES = 2,000`.
- **Ideas**:
  - **Fast 1-WL color-refinement early exit**: if after equitable refinement all variable cells are singletons (as on `fleet-*`), abort in `<10 ms` without constructing full automorphism search state.
  - **Value / signed symmetry support** (already wired to `SYM_SYMTYPE_SIGNPERM` — audit coverage and orbitopal fixing compatibility in SCIP 9.2.4).
  - **Component-wise symmetry detection**: split disconnected symmetry graph components before McKay/Bliss search so one large asymmetric block doesn't exhaust `SYM_MAX_TREE_NODES` before symmetric blocks are found.

### 5.5 Generality Audit for Specialized Root Detectors (`src/bb.c`)
- Before merging `src/pair_assignment.c` (or keeping existing specialized detectors in `src/bb.c:1050-1140`), audit each detector against `AGENTS.md`:
  - Ensure recognition is purely algebraic/structural (never tied to benchmark instance artifacts) and has a **<5 ms rejection fast-path** on non-matching models so overhead across general MIPs is strictly zero.
  - Prefer general mathematical abstractions (e.g., general network-flow/TU-submatrix extraction, general conflict/SAT/PB engines, general symmetry/presolve) over narrow one-family detectors.

---

## 6. Track 5 — Stability, Numerical Robustness & Watchdog Elimination

Zero watchdogs, zero contradictions, and rock-solid numerical behavior across both MIPLIB 2017 and messy real-world models (`fleet-*`) are non-negotiable.

### 6.1 Fix macOS Sleep Clock Drift in `src/common.c:lab_now()`
- **Root cause documented in `SCIP-V9-TOLERANCES-AND-BUDGETS.md`**: `src/common.c:39-44` uses `CLOCK_MONOTONIC`, which advances during macOS system sleep, whereas Python's `time.monotonic()` uses `CLOCK_UPTIME_RAW` (pauses during sleep). Overnight sleep corrupted 160/240 rows of an earlier screen.
- **Fix**: On `#ifdef __APPLE__`, use `CLOCK_UPTIME_RAW` in `lab_now()` (and ensure SCIP's `clocktype = 2` wall clock or `QSopt_ex`'s `zeit.c` uses the same clock source!).

### 6.2 Eliminate All 5 Watchdogs (`>75 s` on a 60 s Limit)
All 5 watchdogs in `lab-v9-full.jsonl` (`neos-2075418-temuka`, `neos-3402454-bohle`, `neos-5114902-kasavu`, `square41`, `square47`) plus 3 late finishers (`k1mushroom` 73.9s, `splice1k1` 68.1s, `supportcase19` 64.3s) hung with **`0` LP iterations** inside presolving or root heuristics:
1. **`presol_logicor` (`neos-3402454-bohle`: 59s, `neos-2075418-temuka`: 27s)**: add `SCIPisStopped(scip)` checks inside pairwise subsumption / backward-strengthening inner loops (every 1,024 iterations) and cap pairwise comparisons on huge logicor sets (`nconss > 20,000`).
2. **`presol_knapsack` (`splice1k1`: 59s, `k1mushroom`: 59s)**: `simplifyinequalities` (`cons_knapsack.c:8047`) extracts massive clique/setpack sets without checking `SCIPisStopped()` frequently enough. Add periodic `SCIPisStopped(scip)` checks and cap clique extraction per round.
3. **`presol_linear` (`square47`: 41s, `rd-rplusc-21`: 40s, `square41`: 18s)**: add `SCIPisStopped(scip)` inside pairwise constraint comparison / dual-presolving inner loops and bound `constraints/linear/nmincomparisons` / time per round.
4. **`QSopt_ex` time-limit polling**: `QSopt_ex` currently checks `lp->maxtime` only in the main simplex outer loop, not inside a huge Markowitz `ILLfactor` factorization or dense `compute_zA`. Adding an every-4096-ops deadline check inside `ILLfactor` prevents single-call LP overruns (`supportcase10`: 74s, `satellites2-40`: 73s, `buildingenergy`: 73s).
5. **Reclaim 4% lost wall time in `src/scip_mip.c:246`**: Once watchdogs are fixed, change `limits/time = remaining * 0.96` (which throws away **2.4 s** of every 60s budget, leaving `neos-957323` at `0.00%` gap right at 57.6s!) to `fmax(0.01, remaining - 0.4)` (or `fmin(remaining - 0.3, remaining * 0.99)`).

### 6.3 Complete the `lpi_qso.c` + `QSopt_ex` Numerical Stability Contract
Building on `scip-v9`'s `FEASTOL`/`DUALFEASTOL` fix:
1. **Wire `SCIP_LPPAR_FROMSCRATCH`** in `lpi_qso.c` + `QSclear_basis()` in `QSopt_ex`: currently `lpi_qso.c` returns `SCIP_PARAMETERUNKNOWN`, preventing SCIP's `lpSolveStable` from clearing a corrupted/singular basis during numerical recovery.
2. **Prevent Silent Tolerance Drift in `qsopt_ex/simplex.c`**:
   - Thirteen recovery branches in `simplex.c` multiply `lp->tol->pfeas_tol` and `lp->tol->dfeas_tol` by `SIMPLEX_FACTOR = 5U` without ever restoring them! After 2–3 numerical restarts, `QSopt_ex` silently declares `QS_LP_OPTIMAL` at `1e-4` violation, causing SCIP to reject the LP!
   - **Fix**: Cap relaxation at `<= 100 * initial_tol` AND if relaxed tolerances were used to finish, restore `initial_tol` and run a bounded polish phase before returning `QS_LP_OPTIMAL` (or return `QS_LP_NUMERR` if polishing fails, so SCIP knows honestly instead of trusting a drifted basis!).
3. **Optional Persistent Power-of-2 Ruiz Scaling for Extreme-Dynamism LPs (`> 1e8` coeff range)**:
   - While `scip-v9` proved `FEASTOL`/`DUALFEASTOL` alone fixed the `fleet` root LP rejection (`[0.2, 2e4]` range), models with `1e9–1e12` matrix/RHS dynamism still suffer unscaled basis instability across B&B nodes because `qsopt.c:290` discards `p2`'s scaling after the first solve.
   - If pursued later, implement power-of-2 Ruiz scaling cleanly inside `lpi_qso.c` **only when matrix dynamism exceeds `1e6`**, including exact `QSopt_strongbranch_bds` (since `floor(psol/cscale) != floor(psol)/cscale`!) and exact `BInvRow`/`BInvARow` unscaling.
4. **Resolve `binkar10_1` Exact Reference Discrepancy**:
   - Both HiGHS 1.15.1 and SCIP/QSopt_ex prove `6742.200024` vs locked MIPLIB v37 `.solu` `6741.38`. Use MILP Lab's own exact rational LP (`QSopt_ex` `mpq`) or check MIPLIB 2017 v38+ errata to document mathematically whether `6741.38` was a known upstream `.solu` bug (never silently edit locked references without versioning a new suite lock!).

---

## 7. Track 6 — New Capabilities & Independent Certification

### 7.1 Exact Safe Dual Bounds (Neumaier–Shcherbina) & Rational LP Polishing
- `QSopt_ex` **already compiles and links its `mpq` exact rational simplex engine** inside `libqsopt.a`!
- **Ideas**:
  1. **Neumaier–Shcherbina safe dual bounding** at every B&B leaf/pruning decision (or root LP): compute a rigorous floating-point-directed-rounding lower bound on the LP dual objective using bounded variable intervals `min/max(d_j * l_j, d_j * u_j)` with `FE_DOWNWARD`/`FE_UPWARD` (already checked by SCIP's CMake!). This guarantees zero false bound-pruning on ill-conditioned MIPs.
  2. **On-demand `mpq` exact LP fallback** when `dbl` `QSopt_ex` returns `QS_LP_NUMERR` or contradicts SCIP's feasibility checks at the root node — completely eliminates numerical-failure aborts on pathological exact/ill-conditioned models (`AGENTS.md` exact/adversarial cases!).

### 7.2 Independent Optimality & Infeasibility Proofs (VIPR Certificates + Pure-C17 GMP Verifier)
- Today `milpcheck` certifies **primal feasibility** in GMP rational arithmetic (`checked_feasible: true`), while `certifies_optimality` and `certifies_infeasibility` are `false` for MIPs.
- **Capability breakthrough**:
  1. Add an optional `--vipr <path>` flag to `milpsolve` (supported directly by SCIP's exact/certificate hooks or a custom SCIP event/proof handler recording original-space cuts, branching disjunctions, Farkas infeasibility rays, and leaf dual multipliers).
  2. Implement a **pure-C17 GMP VIPR proof checker** inside `milpcheck` (linking zero optimization solvers, per `AGENTS.md`!).
  3. This transforms MILP Lab from a fast numerical MILP solver into a **formally certifying MILP solver** — a rare capability even among commercial solvers!

### 7.3 Deterministic Multi-Threaded Mode (`--threads N`)
- Keep `--threads 1` as the default serial benchmark contract, and add optional `--threads N` supporting:
  1. **Concurrent Portfolio Root** (Thread 0: SCIP/QSopt_ex B&C; Thread 1: Core-Guided SAT/PB / RC2; Thread 2: Aggressive Feasibility Jump + LNS primal finder; Thread 3: Alternative branching/cut profile) sharing incumbents and exact dual bounds via lock-free deterministic sync points!
  2. **Parallel B&B Subtree Racing** (ramp-up from open SCIP frontier nodes using independent C17 worker contexts sharing global cutoff bounds).

### 7.4 Broader Modeling & Input Capabilities
- **Streaming compressed MPS / LP reader** (`src/model.c`) with zero fixed-size caps (continuing the `savsched1` `524,288`-column fix toward fully dynamic growth for `1M+` column models).
- **SOS1 / SOS2 and Indicator Constraints** (standard MPS `SOS` / `INDICATORS` sections, currently marked `unsupported`) mapped cleanly into SCIP's native `cons_sos1`, `cons_sos2`, `cons_indicator` handlers and verified directly in `milpcheck`!
- **Warm-start & Solution Polishing API**: accept `--start-solution <sol>` to inject a user incumbent into both `feasjump`/LNS and SCIP for iterative real-world re-optimization workflows (like the 15-minute rolling-horizon `fleet-20261117T2215..2330` instances!).

---

## 8. Track 7 — Benchmarking & Verification Discipline for Future Agents

Hard-won operational rules from Waves 1–6:

1. **Always run screens under `caffeinate -i`** and flag `sleep_affected` runs (`solve_wall_s - wall_s > 5s`).
2. **Never trust cross-hour comparisons without a back-to-back control**: macOS thermal/battery throttling shifted `pg` on identical `scip-v7` code from `9.48 s` to `13.74 s` (1.45x!). Always run A and B **back-to-back** (or interleaved) in the same thermal state.
3. **Use a fast 20-instance stratified A/B micro-panel before 240-instance screens**, covering all 6 bottlenecks:
   - *LP-throughput*: `pg`, `binkar10_1`, `app1-1`, `gmu-35-40`, `tr12-30`, `mik-250-20-75-4`
   - *Strong branching*: `mushroom-best`, `supportcase40`
   - *Root cuts & restarts*: `drayage-25-23`, `h80x6320d`, `neos-860300`, `piperout-27`
   - *Sub-MIP/probing stalls & watchdogs*: `glass-sc`, `dano3_3`, `neos-3402454-bohle`, `square41`
   - *No-incumbent / FJ*: `seymour`, `neos-950242`, `physiciansched6-2`
   - *Real-world ill-conditioned*: `fleet-20261117T2215`, `fleet-20261117T2300`
4. **Never judge by 1 instance** — SCIP's tree has high performance variability; require **panel-level median/geomean speedup AND win-count gain** with **zero new contradictions or watchdogs**.
5. **Always verify Release (62/62, including `harness_contract` outside sandbox) + ASan/UBSan + Native-only (`LAB_SCIP_INTEGER=OFF`)** before committing and pushing to `origin/main`.

---

## 9. Dead Ends / Already-Disproven Ideas (Do Not Repeat)

| Idea Tested | Outcome | Why It Failed |
|---|---|---|
| `PARAM_MAX_NOPROG = 40` in `QSopt_ex` | **Regressed** medium LPs (`gmu-35-40`, `tr12-30`, `binkar10_1`, `neos17`) | Premature perturbation triggered constant unshift + re-solve + refactor loops. Keep `300`. |
| Quick-Start DSE threshold `32` | **Regressed** 400–800 row LPs (+60% iters/node) | Approximate norms degraded pivot quality on medium bases. Sweet spot is `1024`. |
| Iterative DFS + 12% dense bailout in `factor.c` (`btranu3`/`ftranu3`) | **Regressed** `binkar10_1`, `pg`, `h80x6320d`, `drayage-25-23` | 12% threshold is far too low (dense scan touches 88% zeros); reloading node metadata per edge broke L1 locality; changed output index order. Keep Kahn order or only switch to dense at `>40%` fill-in. |
| Capping `separating/aggregation/maxrounds = 1` + `maxtries = 50` in tree | **Regressed** `binkar10_1` (2.2x slower), `h80x6320d` (1.64x), `pg` (1.52x) | In-tree CMIR/flow-cover cuts are essential for tree pruning; speed up the separator's code instead of starving it. |
| Capping `separating/gomory/maxsepacutsroot = 60` + `maxstallroundsroot = 5` | **Regressed** `nursesched-sprint02` and `n5-3` (solved → timeout) | Needed full root Gomory rounds to close in ~3 nodes. Keep SCIP defaults (`200`, `10`). |
| Aggressive `rins` (`nwaitingnodes=0`, `nodesofs=300`, `minfixingrate=0.25`) | **Regressed** `binkar10_1` (`20.4s → >24s`) | Too many early RINS sub-MIP calls. Keep `rins/freq = 15`. |
| 10–12 s flat SAT root budget cap | **Lost** `wachplan` | SAT descent steps scale as fractions of deadline. Keep `min(0.85, max(25s, 0.25 * rem))` (or decouple step timing from total cap). |
| Rewriting `lpi_qso.c` for persistent scaling *just* to fix `fleet` root LP | **Unnecessary for `fleet`** | Root cause was missing `FEASTOL`/`DUALFEASTOL` mapping (`scip-v9` fixed it cleanly). Persistent scaling only makes sense for extreme (`>1e8`) dynamism or as part of a native dual simplex backend. |

---

## 10. Suggested Execution Order for Next Agents

1. **Stabilize & Clean Up (Quick, High-Confidence)**:
   - Finish/validate uncommitted `scip-v10` (`heuristics/subsciptimefrac`), audit `src/pair_assignment.c` for generality, and commit cleanly.
   - Fix `src/common.c:lab_now()` (`CLOCK_UPTIME_RAW` on macOS) + add `SCIPisStopped()` guards in `presol_logicor`/`presol_knapsack`/`presol_linear` + cap `heur_locks` probing LP → **0/240 watchdogs**.
   - Wire `SCIP_LPPAR_FROMSCRATCH` + cap/polish `SIMPLEX_FACTOR` tolerance drift in `QSopt_ex`, and reclaim the 4% wall-time margin (`limits/time`).
2. **LP & Strong-Branching Core Speedups (`68%` of total runtime)**:
   - Capacity-cached `factor_work` in `ILLbasis_factor` (zero-drift!).
   - Cache-compact, bit-for-bit order-preserving `btranu3`/`ftranu3` (`u_delay[]`, `row_to_col[]`, `col_to_row[]`, 1-nz fast path).
   - Sparse `SCIPlpiGetBInvARow` (`3–10x` faster `sepa_gomory`) + 1-pivot dual lookahead & rollback-able eta stack for strong branching.
3. **Presolve, Primal & Root Loop**:
   - Pure-C17 PaPILO-class presolvers: `scip_presol_coefftight.c` (O(nnz) big-M tightening), `scip_presol_parcol.c`, `scip_presol_subst.c`.
   - `prop_redcost/useimplics = TRUE` + post-restart probing + capped root LNS restarts.
   - Feasibility Jump v2 (`src/scip_feasjump.c`): integer-FJ + continuous LP repair, 2-variable equality swaps, periodic in-tree runs when no incumbent exists.
4. **SAT/PB & Exact Certification**:
   - LP-bound-seeded `pb_lo` in `src/binary_sat.c`, followed by pure-C17 incremental CDCL + Stratified OLL/RC2 Core-Guided MaxSAT with exact dual-bound handoff to SCIP.
   - Neumaier–Shcherbina safe dual bounds, `mpq` root-LP numerical fallback, and VIPR proof emission + GMP checker in `milpcheck`.
   - Re-evaluate full 240-instance MIPLIB 2017 at 60s and 180–300s (awake, quiescent) to push decisively past `>120 / 240` wins vs HiGHS.

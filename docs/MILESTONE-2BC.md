# Sparse LP engine and first branch-and-bound solver

The C17 optimizer now supports both continuous LP relaxation and numerical
mixed-integer search. The original checker remains a separate process with its
own reader and exact rational arithmetic. Existing milestone-1 and milestone-2A
Git tags and experiment records are preserved.

## Sparse LP storage and factorization

The optimizer keeps both the original model and transformed standard form
sparse. Columns store only nonzero entries. Basis assembly uses CSC data for
the default sparse factorization. No dense matrix is allocated in that path.
`--backend dense` selects the retained reference factorization for comparison.

Sparse LU uses sorted sparse rows, deterministic partial row pivoting, sparse
elimination merges and a sparse transpose index for transpose solves. Exact
floating zeros are dropped; there is no tolerance-based coefficient dropping.
Factorization rejects nonfinite/unstable pivots and has a 20-million-entry fill
budget. Pivot search still scans rows and no fill-reducing ordering or factor
update scheme is implemented. The change establishes sparse memory use; it
does not establish a general performance advantage.

The dimension limit applies to the expanded standard form, including slacks
and artificial columns: default 4096, sparse ceiling 32768, dense ceiling 4096.
An explicit transformation-entry budget is checked before construction. Large
or badly scaled models may return unsupported or numerical failure. A read or
factorization in progress is not preempted by the internal time check; the
benchmark watchdog remains responsible for process-level deadlines.

## Basis sessions and dual repair

`LabLpSession` owns bound overlays and a saved basis while sharing immutable
original model data. Bound updates can tighten, fix and restore values within
the original domain. Stable transformations retain bound equations, so these
updates change their right-hand sides without changing the matrix. Dual revised
simplex repairs the resulting primal infeasibility while preserving dual
feasibility. Phase I is skipped on a successful warm solve.

Changes to finite/free/fixed-variable topology, incompatible transformations
and unstable warm solves trigger an explicit cold rebuild. A session has one
cumulative deadline; each solve receives the remaining budget. Removed-row
identities survive updates, and a session never changes the caller's original
bound arrays. File basis persistence still validates model identity and canonical
row IDs; it is separate from in-memory bound reoptimization.

Diagnostics include warm starts, dual iterations, cold fallbacks, active
standard-form nonzeros and peak logical factor entries. Factor entries are not
resident bytes: row metadata and the transpose index also occupy memory.
Refactorization still occurs at each pivot; factor updates are future work.

## Numerical branch-and-bound

```sh
mkdir -p work
build/milpsolve tests/fixtures/tiny.mps --mip \
  --solution work/tiny-mip.sol --trace work/tiny-mip.jsonl
build/milpcheck --model tests/fixtures/tiny.mps --solution work/tiny-mip.sol --exact
build/milpsolve tests/fixtures/tiny.mps --mip --cold --backend dense
```

The search uses deterministic depth-first exploration and the first fractional
integer variable, branching at floor/ceil. Original integer bounds are rounded
inward. Every near-integer candidate gets a fixed-integer LP repair, followed
by original-model row/bound/integrality checks. Simply rounding a coupled
continuous/integer solution is insufficient.

Objective comparisons use direct cost-weighted differences and a fixed
reference vector. Objective offsets and fixed-variable constants cannot mask a
worse repaired assignment. Closing a repaired node uses an absolute 1e-7
objective gap. Bounds used for numerical pruning include a conservative margin
and remain explicitly uncertified. Integer bounds or inferred integer values
beyond the supported double integer spacing (absolute value greater than
9007199254740991) return unsupported.

`--time-limit` and `--node-limit` retain an incumbent if one exists. A node-limit
counts search nodes; fixed-integer repair can add LP solves. `--iteration-limit`
limits each LP call. Pending/current node bounds are retained for the reported
global bound on interrupted searches. An unresolved node stops the search with
an uncertainty/limit status; it is not discarded as proved infeasible. LP
unboundedness returns `lp_unbounded`, since it alone does not establish MILP
unboundedness.

Traces record deterministic node, branch, incumbent and prune events without
wall-clock fields. Output and trace files must be fresh. `--mip` and `--relax`
are exclusive and required explicitly. LP-only dual/basis file options are
rejected in MIP mode. The benchmark adapter supports explicit MIP mode with
native solvers and keeps integer checking enabled; LP mode keeps separate
reference objectives.

The solver's `optimal`, `infeasible` and global-bound outputs are numerical
claims. `certifies_optimality` and `bound_certified` remain false for search.
The independent checker proves only the scope its witness supports: original
primal acceptance, an exact rational LP bound, or an exactly matching feasible
primal/LP-bound certificate. There is no full search-tree proof, Farkas/ray
verification, presolve, cutting-plane separation or parallel search yet.

## Correctness gates

The test set includes:

- 1,000 generated LPs checked by an independent exact vertex oracle.
- 300 tiny MILPs checked by integer enumeration and exact residual LP solves.
- Sparse/dense factorization, normal and transpose solves, permutation, fill,
  singular inputs and a 4096-row banded storage case.
- 200 bound-overlay sequences across both backends, checking a closed-form
  optimum plus cold/warm agreement, fixes, sibling restores and infeasibility.
- Original-model primal/dual certificates, parser adversarial cases, cancellation,
  tiny coefficients/costs, objective offsets, fixed-variable constants and safe
  integer-range boundaries.
- Limits with retained incumbents, deterministic trace replay, incompatible
  basis rejection, cumulative deadlines, exclusive output and benchmark modes.

The optional `tests/test_live_mip.py` generates eight independently solved tiny
fixtures and compares HiGHS, SCIP, Cbc, warm lab search and cold lab search in
fresh directories. This checks 40 reference outcomes, including an infeasible
integer model. It does not constitute a performance ranking or general proof.

## Next optimization work

Profile LP iterations, factorization fill/time and branch counts on a declared
development corpus. The next likely components are factor updates and ordering,
then stronger presolve, branching and primal heuristics. Before autoresearch,
control the evaluator separately from optimizer changes, retain all failures,
and confirm improvements on held-out generated/model families and longer,
multi-seed runs. GPU or parallel work should follow a measured bottleneck.

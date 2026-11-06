# Bounded root primal diving

The MILP solver tries one bounded primal dive after solving the root relaxation
and its available root cuts. `--no-root-diving` disables it for comparisons.
The heuristic receives at most one second and five percent of the remaining
search budget, 128 fixing steps and eight sibling retries. Separate size caps
bound its temporary model/state.

The dive first rounds integer coordinates and evaluates the complete point
against the original model and immutable root domain. If that fails, it fixes
one eligible integer coordinate to its closest permitted integer and solves
the remaining LP. Objective direction and column order break distance ties.
An infeasible tentative path can retry the sibling integer or backtrack within
the fixed limits. Nearly integral points can receive an all-integers-fixed LP
solve so coupled continuous coordinates can move before final evaluation.
Values outside the exactly spaced integer range cause abstention.

The heuristic owns its copied bounds, points and LP session. It does not change
the main search session, factorization, basis, node bounds or original model.
Every active cut's root-box premise remains in force. The BB hook repeats the
ordinary original-model candidate check and objective comparison before copying
an improved candidate into the incumbent. It emits the existing `incumbent`
trace event; separate root-dive metrics identify the heuristic work.

Failed, limited, numerically unresolved or locally infeasible dives establish
no original-model infeasibility or objective bound and cannot prune a node.
The MILP search continues through its normal branch policy. All heuristic
time counts against the search budget. Extra LP solves/iterations, steps,
backtracks, candidates, accepted incumbents and elapsed time are reported
separately from the main search counters.

## Evidence

The focused frozen record in `work/root-diving-next/handoff.json` covers 107
calls and 100 independently checked candidates per Release/sanitizer profile,
including 48 exhaustively enumerated MIP families in warm and cold modes,
source feasibility, bound containment, sibling retry, continuous repair,
large integer rejection and main-session/basis immutability.

The final selected first-node diagnostic finds independently checked incumbents
on six of ten MIPLIB cases, versus none with diving disabled. The sanitizer
variant finds two; limited discovery varies with the time cap under concurrent
engineering. Every run remains `node_limit`. These are development discovery
results, not problem completions, a speed ranking or a HiGHS comparison.

The additive [combined contract record](../results/primal-diving-contracts-v1/summary.json)
tests the heuristic together with native exact/point source cuts, source-unit
pivots and the independent checker. The eight original native fixtures remain
separate completion prerequisites. The full 240-instance MIPLIB denominator and
the autoresearch completion gate are unchanged.

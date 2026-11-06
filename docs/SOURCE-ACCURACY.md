# Original-decimal optimization safeguards

The source-accuracy regressions distinguish a correct feasible point from a
correct optimum. For example, maximizing
`1.000000000000000001 X - Y` subject to `X = Y` and
`0 <= X, Y <= 10^15` has optimum `1/1000` at the upper bounds.
The zero point is feasible but suboptimal. Rounding the two costs to binary64
first hides the improving direction. Duplicate row or objective records can
lose their entire nonzero sum through the same cancellation.

The reader coalesces retained source rationals before rounding the resulting
coefficient. Objective-only source data has independent ownership and limits;
it does not turn on the small-model exact-basis path. Candidate evaluation uses
original-source enclosures, with bounded exact reconstruction of ambiguous row
activities. Objective values and differences accumulate source terms before
rounding their final scalar. Branch-and-bound compares those differences and
keeps fixed contributions and offsets out of improvement tests.

Solution files accumulate the source objective at the shortest decimal
coordinates actually exported. This preserves small improvements beside large
offsets and makes the header agree with independent rational scoring. If the
bounded serializer cannot establish a decimal objective claim, it emits the
point without that optional claim. The checker still scores the original model.

When rounded objective costs could hide a material improvement over the current
box, the LP stopping path reconstructs an exact source dual and reduced costs.
An improving source sign selects an ordinary guarded floating-point pivot.
This preserves the existing pivot, residual and feasibility checks. A bounded
exact source dual can also close a support-gap check that fails solely because
the recovered binary dual cannot represent the necessary multiplier. Missing
or exhausted source checks stop conservatively. Proved binary source data
retains the prior fast support check.

When exact support reconstruction is unavailable, complete source enclosures
can instead establish a rigorous support-gap upper bound. Acceptance requires
that upper bound to be at most `1e-7`, which is stricter than the existing
variation-scaled allowance. Ambiguous signs require every relevant finite
support side. Wide intervals, missing provenance and limits never establish
optimality. A source-feasible point accompanying an improving-ray claim does
not require an optimal support bound; it retains the original primal checks.

These are optimizer safeguards. The reported LP/MIP status is still numerical;
internal source-dual closure is not an independently checked optimality
certificate. The independent checker continues to parse and validate original
models separately. Safe numerical failures are unresolved problems, not
completed solves.

Storage and reconstruction remain bounded. The feasibility-only source cache
allows 2,048 rows, 4,096 columns and 100,000 numeric records. Objective-only
storage allows 65,536 columns and 100,000 captured numeric records. The source
pricing path allows 4,096 original rows, 16,384 original columns, 100,000 source
nonzeros and the existing 1,024-row exact sparse basis with a coupled core of
at most 128. Rational operations have bit and time guards; source-guided pivots
also have a 32-pivot cap. Exact scalar activity reconstruction has separate
limits. Reaching a limit never establishes optimality or infeasibility.

The [frozen combined validation](../results/source-accuracy-contracts-v2/summary.json)
passes all 43 suites in Release and ASAN/UBSAN, with no skips. It includes 29
unchanged source-accuracy variants, 76 source-pricing/session controls, eight
large-offset optimal solves, four serialization-only point controls, and 264
interval models checked against 64 exact correlated corners each. Supporting
unit controls also exercise unknown provenance, bounds, masks and resource
limits. The initial combined v1 run and all six failures are retained; the
repairs add exact-dual retry and rigorous interval support, without changing
the original expectations. Earlier informal reports counted 45 source variants;
the actual unchanged total is 29, as the archived correction explains.

This repairs concrete wrong answers and invalid objective headers. It does not
establish universal LP completion: the retained proportional/free-variable LP
diagnostics and decimal fixed-point examples still contain numerical failures.
Those remain unresolved and are not problem-level passes. No new full-corpus
completion or competitive speed result is claimed.

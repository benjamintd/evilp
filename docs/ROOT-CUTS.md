# Exact root GMI cuts

The numerical search can derive up to four Gomory mixed-integer cuts from its
first optimal fractional root LP. This is one bounded round, with at most one
second and five percent of the remaining search budget available for exact
derivation. `--no-root-cuts` preserves the previous search for controlled
comparisons. LP-relaxation mode never activates integer cuts. `--no-exact`
causes cut construction to abstain.

This implements ordinary solver foundations. Paper selection and the
paper-driven autoresearch cycle remain behind the unchanged completion gate.

## Derivation and publication

`native_tableau.inc` reconstructs source-arithmetic basis values and inverse
rows, then checks the canonical basis identity and its identity in original
variables. Original equality rows and fixed box coordinates can occur as
explicit identity premises. Integer distance coordinates require an original
integer variable and an exactly integral anchor; row and bound slacks remain
continuous. A scaled coordinate's integrality is never inferred from its
original variable's flag. Free split coordinates cause abstention.

`gmi.c` applies the exact rational GMI formula to a tableau row. `root_cut.c`
maps its nonnegative native distances back to original coordinates and
converts the cut to a finite binary64 inequality. Every nonzero coefficient
keeps its sign and support. The right-hand side is weakened using the exact
lower support of the coefficient conversion error over the recorded box.
Required missing box sides, nonfinite conversion, bit/storage limits and
expired budgets cause abstention. No proof or feasibility tolerance is added.

`MILPLAB_GMI_1` records the complete box, native identity, integrality premises
and exact dyadic published inequality. The independent checker uses its own
original-source reader and GMI implementation. It does not link optimizer
objects or call engine feasibility functions. See
[the full format and checking rules](GMI-CERTIFICATES.md).

## Search lifecycle and scope

`augment.c` constructs an owned CSC model with appended cut rows. It copies
original names, exact decimal coefficients and outward enclosures without
changing the source model. New cut rows retain their exact binary64 rational
meaning in the exact and enclosure arrays. Existing bound overlays and
presolve row masks remain in the LP view, while original source provenance
continues to identify unchanged source sides.

Search creates the new model and session before swapping its active view.
Sessions are freed before the models they borrow. Cold search also derives
cuts from the initial root session, then uses fresh LP solves as before.
Candidate solutions and objectives are always checked against the original
model. Each node and fixed-integer repair box must be contained in every
active cut's recorded root box, checked with exact source-aware arithmetic.

A cut is valid on **the original MIP intersected with its recorded box**. Its
certificate does not prove that presolve implies that box. A cut-dependent
contradiction is consequently a MIP-in-box contradiction, not an original LP
infeasibility proof. Traces use separate `root_gmi_cut` and
`mip_in_box_infeasibility_certificate` events, retain cut IDs and complete
derivations, and never index appended row names through the original model.
Neither these certificates nor the numerical tree establish an independent
full-tree or global optimality proof.

## Current limits and evidence

The source snapshot is still limited to 128 original rows, 256 columns and
20,000 coefficient records. Native reconstruction additionally limits basis
dimension to 1,024 and basis/canonical nonzeros to 100,000, with a bounded exact
coupled core and rational bit budgets. Only ten locked MIPLIB models meet the
original source-shape limits. Large or unsupported exact reconstructions
continue ordinary search without cuts; they do not count as benchmark passes.

Focused Release and ASAN/UBSAN evidence is retained under
`work/root-cut-contracts-v1` and `work/gmi-checker-integration-v1`. The engine
pipeline independently checks 53 generated certificates per build; native
algebra replay covers 34 exported rows, and Fraction oracles cover the GMI
formula and binary weakening. Integrated search tests compare generated
knapsacks against exhaustive integer enumeration and independently check
returned original-model witnesses and every exported cut.
The [combined 35-suite record](../results/root-cut-combined-contracts-v1/summary.json)
additionally validates search and source-unit pivot recovery together. Its two
initial full runs each retain a Python 3.9 oracle portability failure; the other
34 suites pass, and the repaired oracle passes both unchanged engine binaries
in a targeted rerun. Every tested source file is bound to the final commit.

Controlled packing and parity examples reduce node counts from five to three
and three to one respectively. Those counts are not MIPLIB speed measurements.
The frozen ten-second 240-instance campaign predates these changes and
completes zero problems. All 240 remain in the completion and paired-speed
denominators.

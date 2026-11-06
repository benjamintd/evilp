# Independent GMI derivation certificates

This C17 checker verifies a conditional cut derivation against the **original MIP
intersected with a fully declared rational box**. It does not prove that box
follows from presolve, that all descendants remain in it, or that a full search
tree is complete. No solver or GMI engine code/API is linked. GMI checking is a
separate mode of `milpcheck`:

```sh
milpcheck --model original.mps.gz --gmi cut.gmi
```

The mode rejects primal/dual/Farkas/inspect and tolerance options. No proof
tolerances are accepted.

## Certificate grammar

One record per line; tokens are whitespace separated. Exact numbers are decimal,
scientific-decimal or numerator/positive-denominator rational strings. The first
noncomment record is the header. Blank lines and lines beginning with `#` are
allowed. Names use the exact original MPS namespace. Term IDs are consecutive
decimal integers beginning with zero. `*` is permitted only for infinite box
sides. Explicit rational zeros are allowed; duplicate records are rejected.

```
MILPLAB_GMI_1
SCOPE ORIGINAL_MIP_IN_BOX
BASIC original_column
BETA exact_native_rhs
BOX original_column lower_or_* upper_or_*
... exactly every original column once ...
TERM 0 exact_tableau_coefficient I_or_C VAR_or_ROW_or_BOUND original_name L_or_U anchor
TERM 1 ...
EQ original_equality_row exact_multiplier
FIX fixed_box_column exact_multiplier
PUB_RHS exact_finite_binary64_rational
PUB original_column exact_finite_binary64_rational
... exactly every original column once, including zeros ...
END
```

The single BASIC/BETA/SCOPE/PUB_RHS records and complete BOX/PUB vectors are
required. TERM, EQ and FIX may be empty. Except for the header and END, record
groups may be interleaved. EQ/FIX are unique by source name. END is required and
noncomment data after it is invalid. Unknown columns/rows, N rows as constraint
descriptors, bad fractions, NUL bytes, and malformed/excess fields are rejected.

The native exact tableau convention is

    x_basic + sum(a_k * t_k(x)) = beta.

The normalized coefficient of the **original** basic variable is exactly one.
The term descriptors are native affine forms without a scale parameter:

| Descriptor | L side | U side |
| --- | --- | --- |
| VAR or BOUND column j | `x_j-anchor` | `anchor-x_j` |
| ROW constraint i | `A_i*x-anchor` | `anchor-A_i*x` |

VAR/BOUND anchors must equal the selected finite original source bound **or**
the selected finite recorded-box side. Either premise establishes nonnegativity
on original-model feasible points in that box. ROW anchors must equal the exact
selected finite source row side, after original MPS RANGES interpretation.
The checker never rounds those anchors. A wider recorded box is allowed; the
scope remains an intersection with the original model. A narrower box is not
claimed globally implied.

BASIC must be an original integer column. I is accepted only for VAR with an
original integer column and an integral exact anchor. C is always conservative
once nonnegativity is established. ROW and BOUND must be C in this version.
An original integer column flag does not establish integrality of a scaled
canonical coordinate; an unbounded split-free variable is not a native affine
distance of original x. Such unsupported coordinates must be omitted by the
exporter, transformed with a separate exact justification, or cause abstention.

## Identity modulo original equalities and fixed box coordinates

The checker verifies the exact polynomial identity in original variables:

    E(x) = x_basic + sum(a_k*t_k(x)) - beta
         = sum(mu_i * (A_i*x-h_i)) + sum(nu_j * (x_j-fixed_j)).

Every EQ row must have exactly equal finite original lower/upper row sides,
including valid zero-width ranged rows. Every FIX column must have exactly
equal finite sides in the recorded box. The signs of mu/nu are unrestricted.
Both the original-column coefficients and the affine constant must be exactly
zero after subtracting the right-hand expression. No tolerance is accepted.

This mode reuses only `src/check.c`'s independent exact decimal/rational
parser, original row-bound interpretation, model and sparse matrix streaming.
Its Farkas-mode weighted `A^T*y` accumulation is used as an algebraic matrix
product, without calling a feasibility routine. The two passes independently
assemble the identity residual and the mapped exact cut. There is no dense
matrix allocation, basis inverse, optimizer parser, or optimizer feasibility API.
`src/check_gmi.inc` shares the independent parser within the checker translation
unit; the checker executable links no optimizer objects or APIs.

Input model bytes must stay immutable across these passes, as in the existing
independent checker. GMI mode hashes both inputs before and after checking,
rejects a changed input and records the initial SHA-256 identities in its JSON.
Callers must also preserve the model and certificate alongside the result. The certificate itself is checked by exact
names/algebra, with no trusted stored answers or model-specific behavior.

## Independently derived cut and binary weakening

The checker independently evaluates the standard exact GMI formula from beta,
the tableau coefficients and the checked I/C premises. It never includes,
links or calls `src/gmi.c`. It substitutes every native affine term
back into original coordinates, obtaining the exact original cut `c*x>=d`.
Integral beta returns unavailable; a declared false identity is invalid first.

Every PUB value, including PUB_RHS, must be exactly equal to an actual finite
binary64 double. `1/10` and `0.10000000000000001` are not exact binary64 values,
even when the latter round-trips to the nearest `0.1` double. A fraction or full
exact finite decimal of that dyadic value is accepted. A future replay MPS must
serialize that full exact decimal, not `%.17g`.

The conservative acceptance rules preserve each mapped coefficient's sign and
zero/nonzero support, matching the work-only engine converter. Thus dropping,
introducing, or reversing a coefficient is rejected even if a different, fully
compensated weakening could mathematically be valid. The checker accepts other
finite same-sign roundings, not only the builder's nearest adjacent choice.

For published coefficients p, it computes the exact lower support

    L = inf_{x in recorded box} sum((p_j-c_j)*x_j).

A positive error needs the finite recorded lower side; a negative error needs
the finite recorded upper side; exact zero needs neither. A missing required
side is unavailable, without silently using zero or relying on unrecorded
source bounds. The published rhs must satisfy `published_rhs <= d+L` exactly.
This certifies the rational meaning of the binary inequality on the stated
scope. Floating dot-product accuracy, recovered point feasibility and any
subsequent LP dual bound still require their own unchanged guards.

The GMI JSON reports the exact cut, exact L, maximum safe published rhs, and exact
nonnegative rounding slack for a valid certificate. A cut whose exact maximum
over the recorded box is strictly below its rhs additionally reports
`certifies_mip_in_box_infeasibility`. This remains a scoped MIP claim. All flags
for original LP infeasibility, original global MIP infeasibility, optimality,
checked primal witness, globally implied box and full-tree proof remain false.
The control `integer x, 2x=1, 0<=x<=1` gives the valid constant cut `0>=1/2`,
while the separate original checker accepts x=1/2 in LP relaxation mode.

Exit codes are 0 valid, 1 invalid, 3 unsupported source format, 4 unavailable,
and 64 CLI misuse. There are no proof tolerances or a flag that can weaken them.

## Limits and tests

Limits are 4,096 terms, 16,384 source columns/rows, 1,048,576 source
coefficient records, 16 MiB certificate text, 65,536 bits per rational integer,
an aggregate linear prepass cap and a 20-second cooperative deadline. Exact
operations are checked between stages/terms. The existing MPS parser's full
streaming pass is not asynchronously interrupted and GMP retains normal fatal
system-OOM allocation semantics; this mode is not a hardened hostile-input
resource boundary. Resource/deadline limits cause abstention, not invalid proof
acceptance or an approximate fallback.

The Python contract methods include pure exact source scaling, original decimal
duplicates, native upper/lower distances, continuous and integer nonbasics,
ranged-row and bound slacks, equalities plus fixed-box multipliers, exact binary
coefficient-error compensation, valid alternate weakenings, and explicit local
scope. They reject tampered basic integrality, false integer distances,
anchors/scales, tableau coefficients/constants, equality/fixed multipliers,
boxes, support deletion/sign reversal, unsafe rhs, nonbinary serialization,
duplicate/unknown/missing/malformed records and partial oversized certificates.
Seventy-two generated tableau families are checked against independent Fraction
MIP points. The original primal checker is separately compiled and used for
three witness/LP-typing controls. No optimizer is run.

Normal and ASAN/UBSAN contract runs exercise the same checker executable for
GMI proofs and original primal controls. The existing independent parser
contracts remain required alongside the GMI contracts.

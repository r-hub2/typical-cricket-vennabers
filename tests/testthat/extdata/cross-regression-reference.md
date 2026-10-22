# Cross regression reference evidence

CVAR here is the manual R contract from the
[regression resolution](https://github.com/RMKruse/vennabers/issues/5#issuecomment-5561777290).
Python venn-abers 1.5.4 has no corresponding public manual CVAR wrapper.
The fixture generator uses the same pinned **IVAR components** documented in
[regression-reference.md](regression-reference.md), followed by independent
`math.fsum` of points, lower bounds and upper bounds, each divided by fold count.
No fold-size weighting, pooling or second point merger enters the expectations.

`tools/generate-regression-fixtures.py` regenerates both IVAR and CVAR evidence
with the pinned checkout's `src` on `PYTHONPATH`. It verifies upstream commit
`0e7b4bc575ce230134576ba1c76688d46bf99ac7`, version 1.5.4, and source module
SHA-256 `ea7816e8832eef3f0fa129d2bb80dffae4efcb785f37f083fc087136edb720a5`.
The Cross environment JSON records actual Python/dependency versions, seed,
reference method, counts, and SHA-256 hashes of the four CSV files.
Normal tests use only the installed R package and saved data, without Python
or networking.

Generator version: `cbd2729:tools/generate-regression-fixtures.py` in this
repository's Git history. The Cross CSVs are unchanged from
`cbd2729:inst/extdata`; their exact hashes are recorded in
`cross-regression-environment.json`. The generator is currently git-ignored;
restoring reproducible release tooling belongs to ticket 19.

There are 80 cases: four deterministic and 16 random cases (seed 20260908),
each in ordinary, small, large and shifted outcome scales. They contain 220
folds, 1,540 complete fold predictions and 560 aggregate predictions. Fold
counts vary from two to four, with unequal sizes and noncontiguous calibration
membership. Options include the implicit default, explicit m=1 or m=2, and
epsilon=0.5 with fold-dependent effective moderation. Cases include minimum
size, constant truth, ties, exterior scores, nonconstant zero merger weights
and very small positive weights.

The positive unit scale `s` is set before reference generation: 1, 1e-12,
or 1e120; shifted cases use `s=1` and an additive shift of 1e5. Every numerical
comparison uses `abs(actual-expected) <= 1e-10*s + 1e-10*abs(expected)`.
These tolerances are not enlarged for reference defects. Schema, double types,
row/fold identity and finiteness are checked separately. Ordering is checked
within eight machine epsilons of endpoint magnitude, independently of reference
tolerances, to account for IVAR endpoint interpolation roundoff.

The calibration CSV stores row order, fold membership, actual truth/model
predictions and both reference-moderated truth vectors. The options CSV stores
requested access, fold size, effective m/epsilon, unit/shift and outer bounds.
The predictions CSV stores each fold's named new inputs and full expected
outputs. The aggregate CSV stores the independent equal-weight means.

As with IVAR, the valid lower reference is a reflected upper isotonic component;
the direct Python lower scan has known tie, ordering and exterior-score defects.
Direct Python points/bounds are retained in the predictions CSV: 1,243 valid
rows are compared directly, and 297 defective rows remain visible. Every row
is checked against the valid IVAR component result. Undefined zero-weight
reference mergers use the agreed finite midpoint, never a reference NaN.

Installed R contract tests independently distinguish equal weighting from size
weighting and pooling, and averaging completed points from merging averaged
bounds or taking midpoints. They also cover validation, all column-name errors,
fold-local invalid moderation, metadata, positional alignment, empty/singleton
output, constant subnormal values, overflow-scale means, and cancellation of
opposite-signed folds at magnitudes 1e16 and 1e308 across all six permutations
(independent expected mean 1/3). They also check fresh-process RDS
restoration. The installed fold-model example checks training exclusion,
reconstruction and RMSE/MAE on separate assessment data. Its help example and
the package acceptance gate execute it against the installed artifact.
Independent tied-score calculations additionally verify epsilon=.5 on five
and eight rows (effective m=1 and m=2), complete fold results and their separate
means. Numeric 0/1 outcomes are checked through actual regression points and
bounds, including their unequal-fold aggregate.

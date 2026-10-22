# Inductive regression reference evidence

Source: [Python venn-abers 1.5.4](https://github.com/ip200/venn-abers/tree/0e7b4bc575ce230134576ba1c76688d46bf99ac7),
commit `0e7b4bc575ce230134576ba1c76688d46bf99ac7`, source module SHA-256
`ea7816e8832eef3f0fa129d2bb80dffae4efcb785f37f083fc087136edb720a5`.
The [agreed regression contract](https://github.com/RMKruse/vennabers/issues/5#issuecomment-5561777290)
defines the R interface and numerical exceptions.

`tools/generate-regression-fixtures.py` verifies the source checkout commit,
version and module hash. Run it with that checkout's `src` on `PYTHONPATH` and
NumPy, SciPy, pandas and scikit-learn installed. The environment JSON records
actual versions, seed, options, coverage and SHA-256 hashes of every CSV.
Normal installed-package tests need neither Python nor a network connection.

Generator version: `cbd2729:tools/generate-regression-fixtures.py` in this
repository's Git history. Its IVAR fixtures are unchanged from
`1cc4f79:inst/extdata`; the CSV hashes in `regression-environment.json`
identify the exact stored inputs and outputs. The generator is currently
git-ignored; restoring reproducible release tooling belongs to ticket 19.

The 116 cases comprise five deterministic and 24 seeded random cases, each in
four target transformations: ordinary (`s=1`), small (`s=1e-12`), large
(`s=1e120`) and shifted (`s=1`, additive shift `1e5`). Here `s` is the positive
unit scale fixed before generation, not estimated from errors. Every comparison
uses `abs(actual-expected) <= 1e-10*s + 1e-10*abs(expected)`. Separate checks
cover shape, double types, finiteness, ordering and positional identity.

The options CSV stores requested and effective moderation, unit scale and outer
bounds. The calibration CSV stores actual truth/prediction pairs and the two
moderated target vectors produced by upstream `VennAbersCV.fit`. An identity
estimator supplies exact predictions; a nonshuffled split places the actual
calibration rows last, after three unused proper-training rows. This harness
exercises model-backed Python components; it does not imply a public manual
Python regression wrapper API or package-managed R model training.

## Valid components and retained defects

`VennAbersCV.predict_interval` supplies direct endpoints and Formula 6 points.
Its upper regression endpoint is valid. Its lower scan has defects, including
an incorrect lower endpoint outside the calibration score range, tie handling,
and dependence on equal-score row ordering. Instead, the lower reference is
the negated **upper** regression component on negated predictions and negated
lower-moderated truths. Reflection exchanges isotonic lower and upper insertion
problems. All ordinary-scale endpoints are additionally checked against the
independent max-min characterization of isotonic regression in R, with the
hypothetical observation explicitly appended and equal scores pooled.

For example, truth `(-10,0,2,6,20)`, all scores zero, `m=1`, and new score `-1`
give moderated lower truth `(-10,0,2,6,6)` and `a=-10`. The correct lower
endpoint is `-10`; the direct Python lower endpoint is `-1`. The valid upper
endpoint is `8`. At the tied score zero, lower `-1` and upper `8` merge to
`29/7`; these values are independently hand-checked in the R contract tests.

The predictions CSV preserves direct Python points and both direct bounds,
including nonfinite outputs. `direct_valid` identifies all three direct outputs
agreeing with the component result within the fixed tolerance. There are 1,099
valid rows, compared directly in R, and 337 defective rows retained visibly.
Every row, including the latter, is compared against the valid component
endpoints and their Formula 6 merger. These counts are asserted in tests.
This is component parity with explicit defect repairs, not blanket wrapper parity.

When both merger weights vanish, upstream divides by zero. R uses the stable
midpoint under the agreed contract. Constant truth, nonconstant zero weights,
very small positive weights, coincident bounds, extreme finite values, epsilon
rounding/domain errors and fresh-process RDS restoration have separate public
contract tests. Reference NaNs are never treated as desired predictions.
Independent hand calculations also check maximum moderation (`m=2` on five
tied scores: lower 1, upper 4, point 3) and numeric 0/1 regression outcomes.

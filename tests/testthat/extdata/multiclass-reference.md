# Manual multiclass numerical fixtures

`ovo-calibration.csv` and `ovo-input.csv` are the shared, identical inputs for
OvO and OvR. They store hexadecimal doubles, including recurring fixed-case
values, rather than relying on equal random seeds in Python and R.
`ovo-pairs.csv` / `ovr-classes.csv` contain binary endpoints and points;
`ovo-output.csv` / `ovr-output.csv` contain the finished class probabilities.
There are 33 cases: one fixed and ten random cases for each of two, three,
and four classes, with precision `NULL` (`none`), `0`, `2`, and `-1`, under
both Log and Brier regularization. Random inputs use seed 20260907 and
positive multiples of 1/128; fixed cases include uniform predictions and ties.

The oracle is Python `venn-abers` **1.5.4**, commit
`0e7b4bc575ce230134576ba1c76688d46bf99ac7`. Numerical module SHA-256:
`ea7816e8832eef3f0fa129d2bb80dffae4efcb785f37f083fc087136edb720a5`;
release sdist SHA-256:
`a6bf5383e618ce42d0ecb2335bba0d6b104ecbe979041f10019b57572364bc69`.
`ovo-environment.json` and `ovr-environment.json` record the actual generating
Python/platform/dependency versions, reference identities, and all CSV hashes.

The combined generator is preserved at Git commit `02380de` as
`tools/generate-multiclass-fixtures.py`; the earlier OvO-only generator is at
`1341198` as `tools/generate-ovo-fixtures.py`. Generation checks the installed
module hash and package version before running. Both originally wrote to
`inst/extdata`; the unchanged fixtures now live in `tests/testthat/extdata`.
Restoring generators for a fresh checkout belongs to v1.0 ticket 19.

## Reference components and R extensions

For Log, the generator checks reciprocal OvO coupling against the actual
manual `VennAbersCalibrator.predict_proba` path. OvR is checked against
`predict_proba_prefitted_va(..., "one_vs_rest")`, and also the facade for
three/four classes. The facade dispatches two classes to its binary path,
so the explicit multiclass component is used to preserve two-class OvR
rounding. The R two-class factor contract retains multiclass normalization,
the first-level tie rule, and the unhalved multiclass Brier score.

Binary endpoints and Log/Brier points come from the actual model-backed
`VennAbersCV.predict_proba` merger with table predictors, a one-fold harness,
and explicit precision/loss. Brier multiclass outputs then use independent
coupling/normalization. This harness is reference tooling: it is not a public
manual Python Brier API, and does not change R Cross's minimum of two folds.
Manual Brier and manual shared-fold Cross multiclass are intentional R
extensions. Cross is checked by public binary components followed by
independent equal-weight aggregation and coupling/normalization, without
claiming a public Python manual Cross-multiclass whole-procedure API.

## Separate parity boundaries

- Python divides by zero OvO pair mass without defining that case; NaN
  scores can even lead to finite outputs through lookup placement. Those
  outputs are not oracles. R explicitly uses binary event probability .5
  for zero pairs at fit and prediction; fit rounding still applies, so
  `precision = 0` rounds the replacement to zero.
- Python's manual `predict(..., va_type="one_vs_rest")` does not forward
  the method and can return OvO classes. R classes come from the selected
  method's final matrix. The `random_3_02` fixture distinguishes the two
  methods at row 10 for both losses without rounding.
- R uses native `round()`, not NumPy emulation. At decimal boundaries such
  as .545, the pinned NumPy and local R can differ. Dedicated R checks
  compare public binary fits after native rounding; new scores stay
  unrounded. OvO pair normalization precedes fit rounding.
- Python's Cross product merger can underflow. Stable R aggregation is
  checked separately against exact reductions; see [cross-reference.md](cross-reference.md).
  The known `All-NaN slice encountered` intermediate warning is suppressed
  only during reference generation; every stored result must be finite.

Installed-package tests require maximum absolute error `1e-10`, with no
relative allowance, for endpoints, points, and final matrices. The `1e-8`
input row-sum tolerance is a separate validation rule: accepted matrices
are normalized before calibration. Shapes, finite values, names, class/fold/
row identities, ties, empty/singleton subproblems, and zero/one-row outputs
are checked independently. Public diagnostics reconstruct finished matrices;
unequal folds distinguish equal weights from pooling, size weights, or
averaging finished multiclass fold predictions.

Scores and reliability groups use final matrices. Hand-worked checks protect
unhalved Brier, unclipped Log Loss (including `Inf`), and literal plot-group
means, rates, and counts, including an unobserved assessment class.
Normal installation, tests, and examples use R and stored fixtures only;
they neither execute Python nor download reference data.

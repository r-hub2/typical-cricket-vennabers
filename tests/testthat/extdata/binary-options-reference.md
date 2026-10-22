# Binary option numerical fixtures

`binary-options-folds.csv` stores Inductive endpoints and log/Brier point
probabilities by precision, case, fold, and new row. `binary-options-combined.csv`
stores the corresponding Cross probabilities. Inputs are reused exactly from
`cross-calibration.csv` and `cross-predictions.csv`: all five deterministic
cases and the first twenty seeded random cases described in
[cross-reference.md](cross-reference.md). Each runs with precision `NULL`
(serialized as `none`), `0`, `2`, and `-1`, under both losses. Stored inputs,
rather than matching seeds in different languages, define parity.

The oracle is Python `venn-abers` **1.5.4**, commit
`0e7b4bc575ce230134576ba1c76688d46bf99ac7`. Module/sdist identities, actual
dependency versions, input and fixture SHA-256 hashes, and case coverage are
in `binary-options-environment.json`. Generation verifies the module hash.
Regenerate in the pinned environment from the checkout:

```sh
python tools/generate-binary-options-fixtures.py
```

The original generator is preserved at Git commit `37f3def` as
`tools/generate-binary-options-fixtures.py`; its original output directory
was `inst/extdata`. The unchanged fixtures now live in `tests/testthat/extdata`.
Generator/tool restoration for a fresh checkout belongs to v1.0 ticket 19.

The generator drives the actual `VennAbersCV.predict_proba` model-backed
classification merger with table predictors, matched fold models,
`cv_ensemble=True`, and explicit precision/loss. One-fold reductions of that
same merger supply the Inductive Brier oracle. This is a version-specific
research harness, not a public manual Python Brier API or public one-fold R
Cross API. Manual Brier selection in R is an intentional extension.

Installed-package R tests compare both public fits/predict methods, class
decisions, and fold details against these stored results at maximum absolute
error `1e-10`, without relative tolerance. They also reconstruct both Cross
aggregates from the public table, verify unchanged endpoints, match reversed
fold columns, and check finite probabilities and ordered endpoints. Separate
hand-worked cases check Brier points `3/8`, `1/2`, `5/8`, and equal-weight
aggregation `79/144` from unequal calibration folds. Existing v0.1/v0.2
reference cases and log-space underflow regressions remain in the suite.

The public binary-options tests additionally reverse truth and probabilities
to select either original class as the event, for both losses and all four
precision settings. Unrounded endpoints use the pinned example; Brier points
and the single balanced knots after rounding are independently hand-worked.
At negative precision both orientations round to zero, so calibrated results
need not be complements. Classes are checked against the expected points,
including an exact 0.5. A 1,024-fold public Cross check uses constant nine-row
folds: `(9/10, 1)` gives Log `10/11` and Brier `181/200`; `(0, 1/10)` gives
Log `1/11` and Brier `19/200`. These exact reductions check underflow stability
through the installed public API without using the defective Python products.

## Deliberate parity boundaries

Native R `round()` controls calibration ties. The pinned NumPy environment
rounds `.545` to `.55` and `.575` to `.57` at two digits; local R 4.6.1 on
macOS ARM rounds them to `.54` and `.58`. R/platform versions may differ at
such decimal boundaries. Tests compare calibration with native R rounding
directly, separately from the comparable Python fixtures; no global tolerance
is relaxed. Exactly representable halfway values `.125` and `.375` also
exercise R's ties-to-even behavior. New predictions remain unrounded.

The pinned Python product merger's underflow defects are documented in
`cross-reference.md`. Its unstable values are not oracles; production R
log-space aggregation is checked against algebraic identical-fold reductions.
The reference's known `All-NaN slice encountered` warning is suppressed only
during fixture generation; every stored result must be finite and in range.

Normal installation, tests, help examples, and standalone examples use only R
and the stored data. They neither invoke Python nor download references.

## Local acceptance evidence

On 2026-09-07, R 4.6.1 on macOS 26.6.2 ARM
(`aarch64-apple-darwin25.4.0`) passed the complete installed-artifact gate:
all six package test files, help examples, four README R snippets, and both
standalone examples. `R CMD check --as-cran --no-manual` finished with
`Status: OK` (zero errors, warnings, or NOTEs). The new binary option fixtures
had maximum absolute error `3.330669e-16`. R syntax and Rd parsing also passed.

The gate ran inside a sandbox denying networking and Python execution:

```sh
ACCEPTANCE_OUTPUT_DIR=.scratch/binary-options-acceptance \
  /usr/bin/sandbox-exec \
  -p '(version 1)(allow default)(deny network*)(deny process-exec (regex #".*/python[^/]*$"))' \
  /opt/homebrew/opt/r/bin/Rscript --vanilla tools/acceptance.R
```

Local logs and environment records are retained in that ignored directory;
the console log is `.scratch/binary-options-acceptance.log`. This is local
ARM evidence, not a claim that the hosted OS/R-version matrix has run.

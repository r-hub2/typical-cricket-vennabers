# vennabers parity and capability index

This index lists every capability vennabers supports, what kind of statement
it makes about the Python reference, and the executable evidence behind it.
README, DESCRIPTION, help pages and NEWS refer to this file for parity claims.
"Complete" refers only to the examined and agreed capability scope below. It
does not cover arbitrary sklearn estimators, Tailor step sequences or other
engines' raw interfaces. A newly found reference capability or gap is added
here and re-evaluated visibly.

## Reference and tolerances

The reference is Python `venn-abers` 1.5.4, commit
`0e7b4bc575ce230134576ba1c76688d46bf99ac7`; numerical module SHA-256
`ea7816e8832eef3f0fa129d2bb80dffae4efcb785f37f083fc087136edb720a5`, release
sdist SHA-256 `a6bf5383e618ce42d0ecb2335bba0d6b104ecbe979041f10019b57572364bc69`.
Stored fixtures (identical inputs, not matching seeds) define parity. The
reference notes (`cross-reference.md`, which also covers the randomized
fixtures, `binary-options-reference.md`, `multiclass-reference.md` for OvO and
OvR, `regression-reference.md`, `cross-regression-reference.md`), the
`*-environment.json` records and the Sigmoid generator
`generate-sigmoid-reference.py` in `tests/testthat/extdata` describe generators,
Python environments and hashes. Tests and examples need neither Python nor a
network connection.

- Classification cores: maximum absolute difference `1e-10`, no relative part.
  The `1e-8` input row-sum check is a separate validation rule.
- Regression cores: `abs(actual - expected) <= 1e-10 * s + 1e-10 * abs(expected)`,
  with the positive unit scale `s` fixed before comparison.
- Direct Sigmoid outputs: `1e-7` absolute. Direct isotonic values: `1e-10`
  absolute plus `1e-10` relative.
- Shape, types, finiteness, names and class, fold and row identities are
  checked separately from every numerical tolerance.

## Kinds of statement

- **N numerical parity**: the same valid inputs, calibration/test predictions,
  options and fold assignments give the reference results within the
  tolerances above.
- **T training and data parity**: source models, pair/fold splits,
  preprocessing and calibration roles follow the reference procedure. Equal
  arrays alone do not establish this.
- **R native use**: R functions, workflows, tables and outputs meet the R
  contract. Python class names, tuple shapes and wrapper inheritance are not
  part of the R API.
- **E R extension**: a deliberate capability without a Python counterpart. Its
  evidence is the R contract, not the reference.
- **D reference deviation**: a documented reference defect or defined
  difference, replaced by its own contract tests. Tolerances are never widened
  for it.
- **L accepted limit**: a documented limit of the 1.x contract. The composition
  fails completely with context instead of approximating.

## Capability index

Tests are in `tests/testthat` of the source package (installed with
`R CMD INSTALL --install-tests`); examples are installed under
`system.file("examples", package = "vennabers")`.

| Capability | Statement | Evidence |
| --- | --- | --- |
| Binary Inductive VA endpoints and Log points | N, D | `test-venn-abers.R`, `test-binary-options-parity.R`, `test-randomized-parity.R`; `extdata/binary-options-*`, `extdata/randomized-*`; `examples/binary-options.R`. Native R `round()` controls ties at decimal boundaries, where NumPy can differ (D). |
| Binary Cross VA aggregation | N, D | `test-cross-parity.R`, `test-cross-venn-abers.R`; `extdata/cross-*`; `examples/cross-venn-abers.R`. Stable log-space aggregation replaces the reference product merger, which can underflow; checked against exact reductions (D). |
| Brier regularization, precision, class predictions and fold details | N, E, D | `test-binary-options.R`, `test-binary-options-parity.R`, `test-cross-venn-abers.R`; `extdata/binary-options-*`. Manual Brier selection is an R extension; its reference comes from the model-backed Python merger (E). |
| Manual multiclass one-vs-one | N, D | `test-inductive-multiclass.R`; `extdata/ovo-*`; `examples/inductive-ovo.R`. A zero pair mass uses event probability 0.5 instead of dividing by zero (D). |
| Manual multiclass one-vs-rest | N, D | `test-inductive-ovr.R`, `test-inductive-multiclass.R`; `extdata/ovr-*`; `examples/inductive-ovr.R`. Classes come from the selected method's final matrix; the reference manual OvR predict can return OvO classes (D). |
| Manual Brier multiclass and Cross multiclass with shared folds | E | `test-inductive-multiclass.R`, `test-cross-multiclass.R`; `examples/cross-multiclass.R`. Binary components are compared numerically; the whole procedure has no public Python counterpart. |
| Two-class factor truth with the multiclass contract | D | `test-inductive-multiclass.R`. The reference facade dispatches two classes to its binary path; R keeps multiclass normalization, tie rule and Brier scale. |
| One-hot class decisions | R, D | `test-one-hot.R`; `examples/one-hot-classes.R`. `venn_abers_one_hot()` represents the existing R class decision as a double 0/1 matrix with every class as a column (R). A binary 0.5 tie sets the event column, as in Python `VennAbers.predict` (`>= 0.5`); the Python facade's `argmax` over sorted classes puts it in the first class instead (D). Columns and multiclass ties follow the stored level order, which matches the facade's sorted labels for default `factor()` levels. |
| Brier score, log loss and calibration plots | R, D | `test-venn-abers.R`, `test-inductive-multiclass.R`, `test-cross-multiclass.R`; `examples/inductive-ovo.R`. Unhalved multiclass Brier, unclipped log loss and literal plot groups are hand-checked; no Python scoring export is assumed. |
| Fit diagnoses of fitted calibrators | E | `test-fit-plots.R`; `examples/calibration-graphics.R`. `plot()` draws the stored calibration map and the Venn-Abers width `upper - lower` of binary, multiclass and regression Inductive/Cross calibrators and of workflow fits, without new truth values; the invisible table equals `predict()` at the stored knots. Python venn-abers has no plotting interface. |
| CORP reliability, Brier decomposition and regression assessment graphics | E | `test-corp.R`, `test-regression-plot.R`; `examples/calibration-graphics.R`. `corp_plot()`, `corp_decomposition_plot()` and `regression_plot()` assess held-out pairs; CORP follows Dimitriadis, Gneiting and Jordan (2021). With the same seed, both resampling bands were checked numerically against the R package reliabilitydiag 0.2.1 by `tools/corp-parity-reliabilitydiag.R` in the source repository (maximum difference 0). Two deviations are by design: the cases are sized so that a replicate drawing one unique probability, which the reference mishandles, is practically impossible, and with one outcome only vennabers omits the confidence band that the reference still computes. That check is developer evidence outside the package tests, not the Python reference. Band limits are not Venn-Abers bounds. |
| Moderated Inductive regression (IVAR) with m/epsilon | N, E, D | `test-inductive-regression.R`, `test-regression-parity.R`; `extdata/regression-*`; `examples/inductive-regression.R`. Components come from the model-backed Python path; the manual regression interface itself is an R extension without a public Python wrapper (E). 1,099 valid reference rows compare directly; the defective reference lower scan is replaced by the reflected valid upper component, and vanishing merger weights use the stable midpoint (D). |
| Cross regression (CVAR) with fold results | E, N | `test-cross-regression.R`, `test-cross-regression-parity.R`; `extdata/cross-regression-*`; `examples/cross-regression.R`. No public Python manual CVAR exists; fold components are pinned IVAR components followed by an equal-weight mean. Bounds are descriptive, not intervals. |
| Shared trained workflows: binary, OvO, OvR and regression, Inductive or Cross | T, R, E | `test-calibrated-workflow.R`, `test-cross-workflow.R`, `test-multiclass-workflow.R`, `test-cross-multiclass-workflow.R`, `test-regression-workflow.R`, `test-local-plans.R`, `test-integration-processes.R`; `examples/calibrated-workflow.R`. One shared multiclass model and trained OvR are R extensions, not Python training parity (E). |
| Own models through `model_adapter()` callbacks | R, L | `test-model-adapter.R`, `test-model-adapter-compositions.R`, `test-model-adapter-processes.R`; `examples/model-adapter.R`. R-native counterpart to passing an arbitrary Python estimator: user `fit(x, y)` and `predict(object, new_data)` functions form an unfitted parsnip specification inside the existing workflow lifecycle (Inductive/Cross, matched/full, outer resampling, frozen adoption, saved fits and workers), with probability columns aligned by original class name (R). The same callbacks carry pairwise OvO on each pair's rows and two original levels, Sigmoid probability, Tailor logistic/multinomial/linear and numeric isotonic postprocessing on B with holdout or inner ensembles, and Grid/Bayes searches over calibration and recipe parameters, each matching the native engine; a tuned adapter fit replays on workers and in a fresh process (R). No T: learning only from the local `x`/`y` is a user contract the package cannot check. `sigmoid_input = "glm_link"` requires the native GLM engine, callback state with its own save/restore procedure (external handles) is not covered, and callbacks expose no tunable or replaceable parameters: grid columns or `update()` arguments for the adapter are rejected (L). |
| Fixed engine options and model-internal validation | R, L | `test-engine-options.R`; `examples/engine-options.R`. R-native counterpart to Python fit kwargs that configure the estimator: fixed options in the model specification and `parsnip::set_engine()` reach every local fit (Inductive/Cross sources, class pairs, Full1 and postprocessing B, also in inner ensembles), each on its own rows; native XGBoost early stopping (`stop_iter` with `validation`) translates to a validation share split from each local fit's prepared rows (R). xgboost is not a dependency; `tools/xgboost-early-stopping.R` in the source repository runs the real fits and checks local rows, early stopping and the fold, pair or B/P context of an invalid setting. No T: internal validation follows the local data roles, not a reference procedure. Extra fit arguments to `calibrated_workflow()` or `fit()` are rejected; general Python fit kwargs, per-row fit data (weights, offsets) and an external `eval_set` are not supported (L). |
| Complete outer resampling, last_fit and native scoring | R | `test-outer-resampling.R`, `test-outer-identities.R`, `test-matched-last-fit.R`; `examples/calibrated-resampling.R`, `examples/check-encoded-evaluation.R`. Exact `(id, id2, .row, .config)` coverage is checked, not only a metric table. |
| Automatic, explicit and ordered reference splits | T | `test-reference-splits.R`, `test-local-plans.R`; `examples/reference-splits.R`. Paper-style ordered holdouts illustrate data roles; they do not reproduce paper benchmarks. |
| Independent Inductive training and calibration sizes | T, R, D, L | `test-inductive-sizes.R`; `test-inductive-sizes-processes.R`; `examples/inductive-sizes.R`. `training_size`/`calibration_size` follow sklearn's `train_size`/`test_size` rounding on the local rows (training floor, calibration ceiling, whole numbers exact, missing size as complement). Unassigned rows stay out of source training, B/P and calibration and join the F of Full1, as in the reference; pairwise multiclass resolves the sizes on each class pair, as the reference does per pair (T). R-native arguments of `calibrated_workflow()`, stored `unassigned` positions (per pair) and resolved counts; calibration selection resolves on the final-fit F (R). Stratified draws allocate exact totals per class by largest remainders, rows are not drawn identically to sklearn, and Full1 exists only with `prediction_model = "full"` (D). Absolute row counts are rejected in searches with Venn-Abers candidates (calibration and joint model searches), whose candidates resolve the sizes on resample rows; proportions are accepted there and two-stage model searches accept both (L). |
| Reversible schema encoding of special names | R | `test-schema.R`, `test-schema-workflow.R`, `test-schema-resampling.R`, `test-schema-multiclass-regression.R`, `test-schema-processes.R`; `examples/encoded-binary.R`, `examples/encoded-multiclass-regression.R`. |
| Calibration of existing predictions and frozen workflows | N, E | Manual cores above (N); `test-frozen-workflow.R`, `test-frozen-evaluation.R`, `test-frozen-schema.R`, `test-frozen-multiclass-regression.R`, `test-frozen-processes.R`; `examples/frozen-binary.R`, `examples/frozen-multiclass-regression.R`, `examples/frozen-schema.R`. Direct adoption of a fitted workflow is an R extension (E). |
| Python-style pairwise OvO training | T | `test-pairwise-workflow.R`, `test-pairwise-cross-workflow.R`, `test-pairwise-processes.R`; `examples/pairwise-inductive.R`, `examples/pairwise-cross.R`. Each class pair has its own models and local plans. |
| Inductive and Cross full-model transfer | T, N | `test-full-model-workflow.R`, `test-full-pairwise-workflow.R`, `test-full-multiclass-workflow.R`, `test-full-regression-workflow.R`, `test-full-multiclass-processes.R`, `test-full-regression-processes.R`; `extdata/full-model-reference.csv`; `examples/full-model-workflow.R`, `examples/full-multiclass-workflow.R`, `examples/full-regression-workflow.R`. Full models do not inherit split-calibration guarantees. |
| Joint and two-stage Grid selection without postprocessing | E | `test-binary-grid.R`, `test-multiclass-grid.R`, `test-pairwise-grid.R`, `test-pairwise-grid-cache.R`, `test-pairwise-grid-outer.R`, `test-regression-grid.R`, `test-tuning-plan.R`, `test-tuning-score.R`, `test-tuning-replay-grid-binary.R`; `examples/binary-grid.R`, `examples/shared-multiclass-grid.R`, `examples/pairwise-grid.R`, `examples/regression-grid.R`. All 40 cells of the training matrix. |
| Joint and two-stage Bayes selection without postprocessing | E | `test-binary-bayes.R`, `test-multiclass-bayes.R`, `test-pairwise-bayes.R`, `test-regression-bayes.R`, `test-tuning-replay-bayes-binary.R`; `examples/shared-multiclass-bayes.R`, `examples/pairwise-bayes.R`, `examples/regression-bayes.R`. All 40 cells of the training matrix. |
| Sigmoid holdouts and inner ensembles | N, E, D | `test-sigmoid-holdout.R`, `test-sigmoid-processes.R`, `test-binary-ensemble.R`; `extdata/sigmoid-reference.csv`, `extdata/sigmoid-line-search.csv`; `examples/sigmoid-holdout.R`, `examples/binary-ensemble.R`. Direct outputs match sklearn's smoothed-target fit within `1e-7` on identical inputs; constant P keeps the exact zero-slope solution (D). The composition itself is an R extension (E). |
| Tailor logistic holdouts and ensembles | E | `test-tailor-holdout.R`, `test-tailor-logistic-contract.R`, `test-tailor-processes.R`, `test-tailor-boundary.R`, `test-multiclass-ensemble.R`; `examples/tailor-holdout.R`, `examples/multiclass-ensemble.R`. Python has no tidymodels postprocessing. |
| Tailor multinomial holdouts and ensembles | E, L | `test-tailor-multiclass.R`, `test-tailor-multinomial-compositions.R`, `test-tailor-multiclass-processes.R`; `examples/tailor-multiclass-holdout.R`, `examples/multiclass-ensemble.R`. The P fit uses probably's fixed 100 `nnet::multinom()` iterations; weakly informative inputs can fail to converge, and the composition then fails completely with context (L, maintainer decision of 2026-10-04). |
| Pairwise OvO postprocessing | E | `test-pairwise-postprocessing.R`, `test-pairwise-postprocessing-contract.R`, `test-pairwise-postprocessing-processes.R`; `examples/pairwise-postprocessing.R`. |
| Regression Tailor linear and reference isotonic postprocessing | E, N | `test-tailor-regression.R`, `test-tailor-regression-processes.R`, `test-isotonic-regression.R`, `test-isotonic-regression-processes.R`, `test-regression-ensemble.R`, `test-regression-postprocessing-compositions.R`; `examples/tailor-regression-holdout.R`, `examples/isotonic-regression-holdout.R`, `examples/regression-ensemble.R`. Isotonic values follow the semantics of sklearn's `IsotonicRegression(out_of_bounds = "clip")` and are checked against hand-worked values; there is no stored Python fixture. |
| The 80 fixed postprocessing compositions | E | `test-full-postprocessing.R`, `test-full-postprocessing-processes.R`; `examples/full-postprocessing.R`. Method x holdout/ensemble x Inductive/Cross x matched/full, each with native outer evaluation, last_fit and fresh RDS prediction. |
| Grid and Bayes selection inside postprocessed B | E | `test-postprocessing-grid.R`, `test-postprocessing-grid-processes.R`, `test-postprocessing-bayes.R`, `test-postprocessing-bayes-processes.R`; `examples/postprocessing-grid.R`, `examples/postprocessing-bayes.R`. A targeted 20 + 20 case matrix (five methods x joint/two-stage), not every Cartesian combination of further options. |
| Saved fits, fresh processes and workers | R | `test-integration-processes.R`, `test-schema-processes.R`, `test-binary-ensemble-processes.R`, `test-regression-ensemble-processes.R`, `test-multiclass-ensemble-processes.R`, `test-frozen-multiclass-processes.R`, `test-tuning-replay-grid-ovo.R`, `test-tuning-replay-bayes-ovr.R`; `examples/full-postprocessing.R`. Predict never refits; complete pre- and postprocessing state travels with the fit. |

## Outside the agreed scope

- Out-of-fold postprocessing followed by a complete engine refit, further
  postprocessing methods, foreign engines beyond `model_adapter()` callbacks,
  Tailor step sequences, and other new public constructors.
- General group or time-series procedures (group IDs, rolling windows,
  horizons, embargoes, online updates) beyond the reference splits; the
  manual cores remain available for external splits.
- Standalone multicalibration, quantile or conformal extensions.
- Reproducing every paper dataset or benchmark table, and any blanket promise
  of better accuracy, log loss, Brier score or RMSE, or of interval coverage.

## Open gaps

No known gap remains inside the agreed scope. Platform checks and benchmark
runs are release evidence, not capabilities; this index makes no claim for
any platform.

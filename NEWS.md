# vennabers 1.0.0

Initial CRAN release.

## Manual calibration

- `venn_abers()` fits an Inductive Venn-Abers predictor to calibration
  probabilities and labels from any binary model; `cross_venn_abers()` combines
  the predictors of matched cross-fitted models. `predict()` returns the lower
  and upper Venn-Abers probabilities, a regularized probability chosen by log
  loss or Brier score, and class decisions. Calibration probabilities can be
  rounded with `precision`.
- `brier_score()`, `log_loss()` and `calibration_plot()` assess held-out
  predictions.

## Multiclass

- `venn_abers()` and `cross_venn_abers()` accept factor truth and named
  probability matrices, with one-vs-one coupling (default) or one-vs-rest
  normalization. Class diagnostics are reconstructible from the stored binary
  subproblems.
- `venn_abers_one_hot()` represents existing class decisions as a 0/1 matrix.

## Regression

- `venn_abers_regression()` and `cross_venn_abers_regression()` calibrate
  real-valued predictions into a point with lower and upper bounds in target
  units, with moderation set by `m` or `epsilon`.

## Workflows

- With the optional 'tidymodels' packages installed, `calibrated_workflow()`
  trains a model and its Venn-Abers calibrator together: Inductive or Cross,
  binary, shared one-vs-one or one-vs-rest, pairwise one-vs-one and regression.
  Calibration splits are created locally (`calibration_prop`, or independent
  `training_size` and `calibration_size`), or given explicitly through
  `local_calibration_splits()`. `prediction_model = "full"` refits the model on
  all local rows while keeping the held-out calibrators.
- `calibrate_workflow()` adds a calibrator to an already fitted workflow
  without refitting it.
- `model_adapter()` turns your own fit and predict functions into a model for
  workflows. `venn_abers_schema()`, `venn_abers_encode()` and
  `venn_abers_decode()` keep original column names and class order.
  `venn_abers_details()` returns bounds, fold identities and source-model
  diagnostics.
- `calibrated_workflow()` specifications can be assessed by native outer
  resampling with `tune::fit_resamples()` and `tune::last_fit()`.

## Tuning

- `calibrated_workflow(tuning = ...)` selects model, recipe and calibration
  parameters by Grid or Bayes search, jointly or in two stages, inside each
  local fit. `venn_abers_tuning()` returns the stored search diagnostics.

## Postprocessing

- `calibrated_workflow(postprocessing = ...)` learns a postprocessor on rows of
  its own before the Venn-Abers calibrator: Sigmoid (Platt) scaling, Tailor
  logistic, multinomial or linear calibration, or numeric isotonic regression,
  as a single holdout or an inner ensemble.

## Graphics

- `plot()` methods draw the fitted calibration map or the Venn-Abers width of
  manual calibrators and workflow fits.
- `corp_plot()` draws CORP reliability diagrams with resampling bands,
  `corp_decomposition_plot()` splits the Brier score into miscalibration,
  discrimination and uncertainty, and `regression_plot()` shows held-out
  regression predictions. All graphics use base R.

## Documentation

- Guides on the method and its guarantees, manual calibration, workflows,
  tuning, diagnostics and a comparison with Platt and isotonic calibration,
  plus four collections of the installed examples (`browseVignettes("vennabers")`).
- The installed `PARITY.md` (`system.file("PARITY.md", package = "vennabers")`)
  separates numerical parity with the 'Python' package 'venn-abers' 1.5.4 from
  R extensions, documented deviations and accepted limits.
- `citation("vennabers")` lists the package and the method papers.

## Compatibility promise

- The user-facing functions listed in `help(package = "vennabers")`, their
  arguments and defaults, and the output formats their help pages document
  stay compatible throughout 1.x. Incompatible changes wait for 2.0.0; new
  options keep the previous behavior unless enabled explicitly. Documented bug
  fixes may still change numerical results, and NEWS names each one.
- The parsnip engine entry points `fit_vennabers()`, `fit_vennabers_proxy()`
  and `predict_vennabers_*()` are exported only for model registration; use
  `calibrated_workflow()` instead.
- The minimum R version is not part of this promise. Each release supports at
  least the previous R minor series and may raise the minimum, noted here.

## Development history

- Versions 0.1.0 to 0.9.0 were internal development stages and were never
  published. Their change log is kept in `tools/NEWS-development.md` in the
  source repository.

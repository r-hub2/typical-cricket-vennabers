# B trains the full recipe/engine. Disjoint P trains numeric isotonic
# calibration; C trains IVAR/CVAR. Outer assessment A is evaluation only.
packages <- c("vennabers", "workflows", "recipes", "rsample", "tune", "yardstick")
cases <- list()
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  source(system.file("examples", "assertions.R", package = "vennabers"), local = TRUE)
  original <- data.frame("input value" = sin(1:120),
    "target y" = 10 + 3 * sin(1:120) + cos(1:120) + sin((1:120) / 4), check.names = FALSE)
  schema <- vennabers::venn_abers_schema(original, "target y")
  d <- vennabers::venn_abers_encode(original, schema)
  outcome <- schema$encoded_outcome
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(reformulate(setdiff(names(d), outcome), outcome), d) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  metrics <- yardstick::metric_set(yardstick::rmse, yardstick::mae)
  post <- list(method = "isotonic_numeric", strategy = "holdout", prop = .3)
  outer <- withr::with_seed(74, rsample::vfold_cv(d, v = 2))
  split <- rsample::make_splits(list(analysis = 1:90, assessment = 91:120), d)
  for (method in c("inductive", "cross")) {
    options <- if (method == "inductive") list(m = 1) else list(epsilon = .5)
    spec <- do.call(vennabers::calibrated_workflow, c(list(workflow = base,
      outcome = outcome, schema = schema, method = method, seed = 42,
      postprocessing = post), options))
    resampled <- tune::fit_resamples(spec, outer, metrics = metrics,
      control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, allow_par = FALSE))
    assert_conditions(nrow(tune::collect_notes(resampled)) == 0L,
      nrow(tune::collect_predictions(resampled)) == 120L)
    for (i in seq_along(outer$splits)) {
      state <- tune::collect_extracts(resampled)$.extracts[[i]]
      Q <- as.integer(outer$splits[[i]], "analysis")
      A <- as.integer(outer$splits[[i]], "assessment")
      for (plan in state$postprocessing) assert_conditions(
        !any(Q[c(plan$B, plan$P, plan$C)] %in% A),
        !length(intersect(plan$B, plan$P)), !length(intersect(plan$Q, plan$C)),
        plan$observed_base_fits == 1L, plan$observed_postprocessor_fits == 1L)
    }
    last <- tune::last_fit(spec, split, metrics = metrics)
    assert_conditions(nrow(tune::collect_notes(last)) == 0L)
    fit <- last$.workflow[[1L]]
    new <- d[91:120, setdiff(names(d), outcome), drop = FALSE]
    native <- vennabers::venn_abers_details(fit, new)
    decoded <- vennabers::venn_abers_decode(native, schema)
    collected <- tune::collect_predictions(last)
    assert_conditions(identical(collected$.row, 91:120),
      identical(collected$.pred, native$prediction$.pred),
      all(native$details$postprocessing_method == "isotonic_numeric"),
      all(native$details$postprocessing_input == "numeric"),
      all(is.finite(decoded$prediction)))
    before <- workflows::extract_fit_engine(fit)
    for (column in c(names(d)[1L], outcome)) {
      changed <- d
      changed[91:120, column] <- changed[91:120, column] + 100
      other <- tune::last_fit(spec, rsample::make_splits(list(analysis = 1:90, assessment = 91:120), changed),
        metrics = metrics)$.workflow[[1L]]
      after <- workflows::extract_fit_engine(other)
      assert_conditions(identical(before$postprocessing, after$postprocessing),
        identical(before$calibrator, after$calibrator),
        identical(predict(fit, new), predict(other, new)))
      for (id in names(before$models)) assert_conditions(
        identical(before$models[[id]]$isotonic, after$models[[id]]$isotonic),
        identical(coef(workflows::extract_fit_engine(before$models[[id]]$workflow)),
                  coef(workflows::extract_fit_engine(after$models[[id]]$workflow))))
    }
    cases[[method]] <- list(fit = fit, new = new, schema = schema)
  }
}

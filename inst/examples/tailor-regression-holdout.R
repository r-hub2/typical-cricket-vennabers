# Installed example: B learns the recipe/engine, disjoint P learns Tailor's
# linear adjustment, and held-out C learns IVAR/CVAR. Outer A learns nothing.
packages <- c("vennabers", "workflows", "tailor", "probably", "recipes", "tune", "yardstick")
cases <- list()
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  source(system.file("examples", "assertions.R", package = "vennabers"), local = TRUE)
  audit <- new.env()
  sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), audit)
  original <- data.frame("input value" = sin(1:240), "input.value" = cos(1:240),
    ".pred" = 10 + 3 * sin(1:240) + cos(1:240) + sin((1:240) / 3), check.names = FALSE)
  schema <- vennabers::venn_abers_schema(original, ".pred")
  d <- vennabers::venn_abers_encode(original, schema)
  outcome <- schema$encoded_outcome
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(reformulate(setdiff(names(d), outcome), outcome), d) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  # prop reserves P's fraction of Q. P needs at least three rows and an actual
  # full-rank intercept/prediction design. smooth=FALSE is fixed; finite input,
  # learned method and convergence are checked. No fallback or refit at predict.
  # Automatic B/P plans are unstratified, including numeric 0/1 outcomes.
  post <- list(method = "tailor_linear", strategy = "holdout", prop = .3)
  metrics <- yardstick::metric_set(yardstick::rmse, yardstick::mae)
  for (method in c("inductive", "cross")) {
    spec <- vennabers::calibrated_workflow(base, outcome, schema = schema,
      method = method, m = 1, seed = 42, postprocessing = post)
    outer <- withr::with_seed(74, rsample::vfold_cv(d, v = 2, repeats = 2))
    result <- tune::fit_resamples(spec, outer, metrics = metrics,
      control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, allow_par = FALSE))
    audit$check_encoded_evaluation(result, outer, schema)
    states <- tune::collect_extracts(result)$.extracts
    for (i in seq_along(states)) for (plan in states[[i]]$postprocessing) {
      Q <- as.integer(outer$splits[[i]], "analysis")
      A <- as.integer(outer$splits[[i]], "assessment")
      assert_conditions(!any(Q[c(plan$B, plan$P, plan$C)] %in% A),
        !length(intersect(plan$B, plan$P)), !length(intersect(plan$Q, plan$C)),
        plan$observed_base_fits == 1L, plan$observed_postprocessor_fits == 1L)
    }
    split <- rsample::make_splits(list(analysis = 1:180, assessment = 181:240), d)
    last <- tune::last_fit(spec, split, metrics = metrics)
    fit <- last$.workflow[[1L]]
    native <- vennabers::venn_abers_details(fit, d[181:240, setdiff(names(d), outcome)])
    decoded <- vennabers::venn_abers_decode(native, schema)
    predictions <- tune::collect_predictions(last)
    scores <- tune::collect_metrics(last)
    truth <- original[[schema$outcome]][181:240]
    reference <- c(rmse = yardstick::rmse_vec(truth, decoded$prediction),
      mae = yardstick::mae_vec(truth, decoded$prediction))
    assert_conditions(nrow(tune::collect_notes(last)) == 0L, identical(predictions$.row, 181:240),
      !anyDuplicated(predictions$.row), all(predictions$.config == "pre0_mod0_post0"),
      nrow(scores) == 2L, setequal(scores$.metric, names(reference)), all(is.finite(scores$.estimate)),
      all(abs(scores$.estimate - reference[scores$.metric]) <= 1e-10),
      identical(predictions$.pred, native$prediction$.pred),
      identical(decoded$postprocessing[[1L]]$columns, names(original)),
      all(native$details$postprocessing_method == "tailor_linear"),
      all(native$details$postprocessing_input == "numeric"))
    # Predictor-only and truth-only outer A mutations leave all learned state fixed.
    for (column in c(names(d)[1L], outcome)) {
      changed <- d; changed[181:240, column] <- changed[181:240, column] + 1000
      other <- tune::last_fit(spec, rsample::make_splits(list(analysis = 1:180, assessment = 181:240), changed), metrics = metrics)
      before <- workflows::extract_fit_engine(fit)
      after <- workflows::extract_fit_engine(other$.workflow[[1L]])
      assert_conditions(identical(before$postprocessing, after$postprocessing),
        identical(before$calibrator, after$calibrator),
        identical(predict(fit, d), predict(other$.workflow[[1L]], d)))
      for (id in names(before$models)) assert_conditions(
        identical(workflows::extract_recipe(before$models[[id]])$steps, workflows::extract_recipe(after$models[[id]])$steps),
        identical(coef(workflows::extract_fit_engine(before$models[[id]])), coef(workflows::extract_fit_engine(after$models[[id]]))))
    }
    cases[[method]] <- list(fit = fit, schema = schema, new = d[setdiff(names(d), outcome)])
  }
}

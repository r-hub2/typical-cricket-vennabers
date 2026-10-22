testthat::test_that("regression Bayes selects moderation and preserves its fitted prediction", {
  skip_if_no_integration("tune", "dials", "yardstick")
  data <- data.frame(x = sin(seq_len(100) / 7))
  data$truth <- 3 * data$x + cos(seq_len(nrow(data))) / 5
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  params <- dials::parameters(dials::new_quant_param(
    type = "integer", range = c(1L, 2L), inclusive = c(TRUE, TRUE),
    label = c(moderation = "Moderation")))
  config <- list(strategy = "joint", selection_prop = .5,
    calibration = list(method = "bayes", initial = 2L, iter = 1L,
      param_info = params, resamples = function(d) rsample::vfold_cv(d, v = 2),
      control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 2L),
      report_metrics = list(yardstick::new_numeric_metric(
        function(...) stop("report unavailable"), direction = "minimize"))))
  fitted <- generics::fit(calibrated_workflow(base, "truth", m = tune::tune("moderation"),
    tuning = config, seed = 7), data)
  selected <- venn_abers_tuning(fitted)$calibration
  testthat::expect_true(selected$chosen %in% selected$candidates$id)
  testthat::expect_true(all(selected$candidates$status %in% c("eligible", "excluded")))
  testthat::expect_true(is.finite(selected$candidates$estimate[
    selected$candidates$id == selected$chosen]))
  testthat::expect_true("m" %in% names(tune::collect_predictions(selected$native)))
  testthat::expect_match(selected$report_metrics[[1L]][[1L]][[1L]], "report unavailable")
  testthat::expect_identical(selected$stop_reason, "repeated_candidate")
  testthat::expect_true(length(selected$steps$step_1$repeated) > 0L)
  state <- workflows::extract_fit_engine(fitted)
  core <- venn_abers_regression(state$calibration$Inductive1$truth,
    state$calibration$Inductive1$prediction, m = state$parameters$requested$m)
  expected <- predict(core, predict(state$models$Inductive1, data)$.pred)$prediction
  testthat::expect_true(all(abs(predict(fitted, data)$.pred - expected) <=
    1e-10 + 1e-10 * abs(expected)))
  saved <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  log <- tempfile(fileext = ".log")
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  saveRDS(list(fitted = fitted, data = data, prediction = predict(fitted, data)), saved)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "library(vennabers)", paste0("saved <- readRDS(", deparse(saved), ")"),
    'trace("tune_bayes", where = asNamespace("tune"), tracer = quote(stop("Unexpected search")), print = FALSE)',
    'trace("fit.workflow", where = asNamespace("workflows"), tracer = quote(stop("Unexpected fit")), print = FALSE)',
    'stopifnot(identical(predict(saved$fitted, saved$data), saved$prediction))'
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
    stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

testthat::test_that("regression Bayes does not replace a failed final fit", {
  skip_if_no_integration("tune", "dials")
  audit <- new.env(parent = emptyenv())
  audit$fits <- list()
  final_fit <- function(x, y, variant) {
    audit$fits[[length(audit$fits) + 1L]] <- list(rows = nrow(x), variant = variant)
    if (variant == 1L && nrow(x) > 60L) stop("planned final fit failure")
    list(variant = variant, mean = mean(y))
  }
  final_predict <- function(object, new_data) {
    if (object$variant == 1L) new_data$x else rep(object$mean, nrow(new_data))
  }
  register_reference_engine("linear_reg", "regression", "final_fit", "final_predict", "numeric",
    engine = "bayes_final_reference", new_model = FALSE)
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg() |>
      parsnip::set_engine("bayes_final_reference", variant = tune::tune("variant")))
  data <- data.frame(x = seq(-1, 1, length.out = 100))
  data$truth <- data$x
  params <- dials::parameters(dials::new_quant_param(
    type = "integer", range = c(0L, 1L), inclusive = c(TRUE, TRUE),
    label = c(variant = "Variant")))
  config <- list(strategy = "two_stage", model = list(method = "bayes",
    initial = 2L, iter = 1L, param_info = params,
    resamples = function(d) rsample::vfold_cv(d, v = 2),
    control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 2L)))
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth",
    tuning = config, seed = 7), data), "planned final fit failure")
  testthat::expect_equal(sum(vapply(audit$fits, function(f)
    f$variant == 0L && f$rows > 60L, logical(1))), 0L)
})

testthat::test_that("regression Bayes rejects locally invalid moderation and records time stops", {
  skip_if_no_integration("tune", "dials")
  data <- data.frame(x = sin(seq_len(100) / 7))
  data$truth <- 3 * data$x + cos(seq_len(nrow(data))) / 5
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  params <- function(range, type) dials::parameters(dials::new_quant_param(
    type = type, range = range, inclusive = c(TRUE, TRUE),
    label = c(moderation = "Moderation")))
  config <- list(strategy = "joint", selection_prop = .5,
    calibration = list(method = "bayes", initial = 2L, iter = 2L,
      param_info = params(c(1L, 100L), "integer"),
      resamples = function(d) rsample::vfold_cv(d, v = 2),
      control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 2L)))
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth",
    m = tune::tune("moderation"), tuning = config, seed = 7), data),
    "complete initial calibration candidates")
  config$calibration$param_info <- params(c(.001, .02), "double")
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth",
    epsilon = tune::tune("moderation"), tuning = config, seed = 7), data),
    "complete initial calibration candidates")
  config$calibration$param_info <- params(c(1L, 2L), "integer")
  config$calibration$control <- tune::control_bayes(seed = 4, save_pred = TRUE,
    no_improve = 2L, time_limit = 1e-9)
  fitted <- generics::fit(calibrated_workflow(base, "truth",
    m = tune::tune("moderation"), tuning = config, seed = 7), data)
  selected <- venn_abers_tuning(fitted)$calibration
  testthat::expect_identical(selected$stop_reason, "time_limit")
  testthat::expect_identical(selected$attempts, 2L)
  testthat::expect_true(is.finite(selected$elapsed_seconds))
})

testthat::test_that("Cross regression Bayes keeps scale and zero-target parity", {
  skip_if_no_integration("tune", "dials")
  data <- data.frame(x = sin(seq_len(120) / 7))
  unit_truth <- 3 * data$x + cos(seq_len(nrow(data))) / 5
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  params <- dials::parameters(dials::new_quant_param(
    type = "double", range = c(.3, .4), inclusive = c(TRUE, TRUE),
    label = c(moderation = "Moderation")))
  config <- list(strategy = "two_stage", selection_prop = .5,
    calibration = list(method = "bayes", initial = 3L, iter = 1L,
      param_info = params, resamples = function(d) rsample::vfold_cv(d, v = 2),
      control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 2L)))
  for (scale in c(0, 1e-9, 1e9)) {
    data$truth <- scale * unit_truth
    fitted <- generics::fit(calibrated_workflow(base, "truth", method = "cross", v = 2,
      epsilon = tune::tune("moderation"), tuning = config, seed = 7), data)
    state <- workflows::extract_fit_engine(fitted)
    held <- state$calibration
    core <- cross_venn_abers_regression(
      truth = unlist(lapply(held, `[[`, "truth")),
      prediction = unlist(lapply(held, `[[`, "prediction")),
      fold = rep(names(held), lengths(lapply(held, `[[`, "truth"))),
      epsilon = state$parameters$requested$epsilon)
    source <- lapply(state$sources, function(id) predict(state$models[[id]], data)$.pred)
    expected <- predict(core, do.call(cbind, source))
    actual <- venn_abers_details(fitted, data)$aggregate
    for (column in c("prediction", "lower", "upper")) {
      testthat::expect_true(all(abs(actual[[column]] - expected[[column]]) <=
        1e-10 * max(scale, 1e-9) + 1e-10 * abs(expected[[column]])))
    }
  }
})

testthat::test_that("regression Bayes spends continuation budget on new candidates", {
  skip_if_no_integration("tune", "dials")
  data <- data.frame(x = sin(seq_len(150) / 7))
  data$truth <- 3 * data$x + cos(seq_len(nrow(data))) / 5
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  params <- dials::parameters(dials::new_quant_param(
    type = "double", range = c(.25, .45), inclusive = c(TRUE, TRUE),
    label = c(moderation = "Moderation")))
  config <- list(strategy = "joint", selection_prop = .5,
    calibration = list(method = "bayes", initial = 3L, iter = 2L,
      param_info = params, resamples = function(d) rsample::vfold_cv(d, v = 2),
      control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 3L)))
  fitted <- generics::fit(calibrated_workflow(base, "truth", calibration_prop = .5,
    epsilon = tune::tune("moderation"), tuning = config, seed = 7), data)
  selected <- venn_abers_tuning(fitted)$calibration
  testthat::expect_gt(selected$attempts, 3L)
  testthat::expect_true(any(vapply(selected$steps[-1L], function(step)
    is.data.frame(step) && nrow(step) > 0L, logical(1))))
  eligible <- subset(selected$candidates, status == "eligible")
  testthat::expect_identical(selected$chosen, eligible$id[which.min(eligible$estimate)])
  testthat::expect_equal(selected$cost$candidate_evaluations,
    selected$attempts * length(selected$resamples))
  testthat::expect_true("epsilon" %in% names(tune::collect_predictions(selected$native)))
})

testthat::test_that("regression Bayes fits local model and moderation choices in both strategies", {
  skip_if_no_integration("tune", "dials")
  data <- data.frame(x = sin(seq_len(200) / 7))
  data$truth <- data$x^2 + cos(seq_len(nrow(data))) / 5
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x, data) |>
      recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  model_params <- dials::parameters(dials::new_quant_param(
    type = "integer", range = c(1L, 3L), inclusive = c(TRUE, TRUE),
    label = c(degree = "Degree")))
  for (strategy in c("joint", "two_stage")) for (method in c("inductive", "cross")) {
    for (prediction_model in c("matched", "full")) {
      option <- if (prediction_model == "matched") "m" else "epsilon"
      moderation_params <- dials::parameters(dials::new_quant_param(
        type = if (option == "m") "integer" else "double",
        range = if (option == "m") c(1L, 2L) else c(.3, .4),
        inclusive = c(TRUE, TRUE), label = c(moderation = "Moderation")))
      search <- function(params, seed) list(method = "bayes", initial = 3L, iter = 1L,
        param_info = params, resamples = function(d) rsample::vfold_cv(d, v = 2),
        control = tune::control_bayes(seed = seed, save_pred = TRUE, no_improve = 2L))
      config <- list(strategy = strategy, selection_prop = .5,
        model = search(model_params, 3), calibration = search(moderation_params, 4))
      args <- list(workflow = base, outcome = "truth", method = method,
        prediction_model = prediction_model,
        tuning = config, seed = 31)
      args[[option]] <- tune::tune("moderation")
      if (method == "cross") args$v <- 2L else args$calibration_prop <- .5
      fitted <- generics::fit(do.call(calibrated_workflow, args, quote = TRUE), data)
      selection <- venn_abers_tuning(fitted)
      state <- workflows::extract_fit_engine(fitted)
      testthat::expect_setequal(names(selection$models),
        c(names(state$folds), if (prediction_model == "full") "Full1"))
      testthat::expect_length(intersect(selection$calibration$selection_rows,
        selection$calibration$final_rows), 0L)
      testthat::expect_true(all(is.finite(predict(fitted, data)$.pred)))
      testthat::expect_true(all(vapply(selection$models, function(x)
        x$chosen %in% x$candidates$id, logical(1))))
      final <- data[selection$calibration$final_rows, , drop = FALSE]
      for (id in names(selection$models)) {
        local <- if (id == "Full1") final else
          final[state$folds[[id]]$analysis, , drop = FALSE]
        testthat::expect_equal(as.data.frame(selection$models[[id]]$native$splits[[1L]]$data),
          local, ignore_attr = TRUE)
      }
      restored <- unserialize(serialize(fitted, NULL))
      testthat::expect_equal(predict(restored, data), predict(fitted, data))
      held <- state$calibration
      core_args <- list(truth = unlist(lapply(held, `[[`, "truth")),
        prediction = unlist(lapply(held, `[[`, "prediction")))
      core_args[[option]] <- state$parameters$requested[[option]]
      if (method == "cross") core_args$fold <- rep(names(held),
        lengths(lapply(held, `[[`, "truth")))
      core <- do.call(if (method == "cross") cross_venn_abers_regression else
        venn_abers_regression, core_args)
      source <- lapply(state$sources, function(id) predict(state$models[[id]], data)$.pred)
      expected <- predict(core, if (method == "cross") do.call(cbind, source) else source[[1L]])
      actual <- venn_abers_details(fitted, data)$prediction$.pred
      testthat::expect_true(all(abs(actual - expected$prediction) <=
        1e-10 + 1e-10 * abs(expected$prediction)))
    }
  }
})

testthat::test_that("Bayes accounts for discarded and repeated proposals without postprocessing", {
  skip_if_no_integration("tune", "dials")
  audit <- new.env(parent = emptyenv())
  strict_fit <- function(x, y, variant) {
    audit$fits <- c(audit$fits, variant)
    if (!variant %in% c(1L, 3L, 6L)) stop("planned continuation failure")
    list(variant = variant)
  }
  strict_predict <- function(object, new_data) new_data$x / object$variant
  register_reference_engine("linear_reg", "regression", "strict_fit", "strict_predict", "numeric",
    engine = "strict_bayes_reference", new_model = FALSE)
  data <- data.frame(x = seq(-1, 1, length.out = 120))
  data$truth <- data$x
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg() |>
      parsnip::set_engine("strict_bayes_reference", variant = tune::tune("variant")))
  params <- dials::parameters(dials::new_quant_param(type = "integer", range = c(1L, 6L),
    inclusive = c(TRUE, TRUE), label = c(variant = "Variant")))
  run <- function(no_improve = 10L, time_limit = NA_real_) {
    audit$fits <- integer()
    config <- list(strategy = "two_stage", model = list(method = "bayes", initial = 3L, iter = 3L,
      param_info = params, resamples = function(d) rsample::vfold_cv(d, v = 2),
      control = tune::control_bayes(seed = 4, save_pred = TRUE,
        no_improve = no_improve, time_limit = time_limit)))
    fitted <- generics::fit(calibrated_workflow(base, "truth", seed = 2, tuning = config), data)
    selection <- venn_abers_tuning(fitted)
    testthat::expect_identical(selection$cost$observed_base_fits, length(audit$fits))
    selection$models$Inductive1
  }
  selected <- run()
  testthat::expect_identical(selected$steps$initial$variant, c(6L, 1L, 3L))
  testthat::expect_identical(selected$candidates$id, sprintf("model_%03d", 1:4))
  testthat::expect_identical(selected$steps$step_1$variant, 2L)
  testthat::expect_identical(selected$steps$step_1$status, "excluded")
  testthat::expect_match(selected$steps$step_1$reason, "planned continuation failure")
  testthat::expect_true(any(grepl("planned continuation failure", selected$notes$note)))
  testthat::expect_setequal(tune::collect_metrics(selected$native)$variant, c(1L, 3L, 6L))
  testthat::expect_identical(selected$steps$step_2$repeated$variant, 3L)
  testthat::expect_identical(selected$stop_reason, "repeated_candidate")
  testthat::expect_identical(selected$attempts, 5L)
  testthat::expect_identical(selected$no_improve, 2L)
  testthat::expect_identical(selected$cost$candidate_evaluations, 10L)
  testthat::expect_identical(selected$cost$observed_base_fits, 10L)
  testthat::expect_identical(selected$chosen, "model_002")
  stale <- run(no_improve = 1L)
  testthat::expect_identical(stale$stop_reason, "no_improve")
  testthat::expect_identical(stale$attempts, 4L)
  testthat::expect_identical(stale$cost$observed_base_fits, 8L)
  timed <- run(time_limit = 1e-9)
  testthat::expect_identical(timed$stop_reason, "time_limit")
  testthat::expect_named(timed$steps, "initial")
  testthat::expect_identical(timed$attempts, 3L)
  testthat::expect_identical(timed$cost$observed_base_fits, 6L)
})

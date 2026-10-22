testthat::test_that("regression model grids use local numeric predictions and RMSE", {
  testthat::skip_if_not_installed("tune")
  data <- data.frame(x = seq(-2, 2, length.out = 80))
  data$truth <- data$x^2 + sin(seq_len(nrow(data))) / 10
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x, data) |>
      recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  config <- list(strategy = "two_stage", model = list(
    method = "grid", grid = data.frame(degree = c(1L, 2L)),
    resamples = function(d) rsample::vfold_cv(d, v = 2)))
  outer <- calibrated_workflow(base, "truth", tuning = config, m = 1, seed = 42)
  testthat::expect_equal(nrow(generics::tune_args(outer)), 0L)
  fitted <- generics::fit(outer, data)
  selected <- venn_abers_tuning(fitted)$models$Inductive1
  testthat::expect_identical(selected$selected$degree, 2L)
  testthat::expect_identical(selected$primary_metric, yardstick::rmse)
  testthat::expect_identical(selected$candidates$status, rep("eligible", 2L))
  testthat::expect_identical(selected$cost$observed_base_fits, 4L)
  predictions <- tune::collect_predictions(selected$native)
  testthat::expect_true(all(is.finite(predictions$.pred)))
  first <- subset(predictions, candidate_id == "model_001" & id == "Fold1")
  testthat::expect_equal(selected$scores[[1L]][[1L]],
    sqrt(mean((first$truth - first$.pred)^2)))
  testthat::expect_true(all(is.finite(predict(fitted, data)$.pred)))
  config$model$grid <- 2L
  config$model$param_info <- dials::parameters(dials::new_quant_param(
    type = "integer", range = c(1L, 2L), inclusive = c(TRUE, TRUE),
    label = c(degree = "Degree")))
  generated <- generics::fit(calibrated_workflow(base, "truth", tuning = config, m = 1, seed = 42), data)
  testthat::expect_setequal(venn_abers_tuning(generated)$models$Inductive1$grid$degree, c(1L, 2L))
})

testthat::test_that("regression moderation grids validate options before fitting", {
  testthat::skip_if_not_installed("tune")
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  config <- list(strategy = "joint", selection_prop = .5,
    calibration = list(method = "grid", grid = data.frame(moderation = 1),
      resamples = function(d) rsample::vfold_cv(d, v = 2)))
  for (option in c("m", "epsilon")) {
    args <- list(workflow = base, outcome = "truth", tuning = config)
    args[[option]] <- tune::tune("moderation")
    bad_values <- if (option == "m") list(0, -1, 1.5, Inf, NA, "1") else
      list(0, -1, 1, Inf, NA, ".5")
    for (bad in bad_values) {
      args$tuning$calibration$grid <- data.frame(moderation = bad)
      testthat::expect_error(do.call(calibrated_workflow, args, quote = TRUE),
        "m|epsilon|missing")
    }
    args$tuning$calibration$grid <- data.frame(moderation = if (option == "m") 1 else .5)
    args[[if (option == "m") "epsilon" else "m"]] <- .5
    testthat::expect_error(do.call(calibrated_workflow, args, quote = TRUE), "only one")
  }
  testthat::expect_error(calibrated_workflow(base, "truth", m = tune::tune("moderation")), "require.*tuning")
  testthat::expect_error(calibrated_workflow(base, "truth", m = tune::tune("moderation"),
    epsilon = tune::tune("other"), tuning = config), "only one")
  for (option in c("loss", "precision", "classification", "event_level")) {
    args <- list(workflow = base, outcome = "truth", m = tune::tune("moderation"), tuning = config)
    args[[option]] <- if (option == "classification") "binary" else
      if (option == "event_level") "first" else tune::tune("classification_option")
    testthat::expect_error(do.call(calibrated_workflow, args, quote = TRUE), "Regression does not accept")
  }
  config$calibration$control <- tune::control_grid(save_pred = FALSE)
  testthat::expect_error(calibrated_workflow(base, "truth", m = tune::tune("moderation"),
    tuning = config), "save_pred")
  config$calibration$control <- NULL
  config$selection_prop <- NULL
  testthat::expect_error(calibrated_workflow(base, "truth", m = tune::tune("moderation"),
    tuning = config), "exactly one")
})

testthat::test_that("regression grids preserve scales and zero targets with explicit unequal folds", {
  testthat::skip_if_not_installed("tune")
  data <- data.frame(x = sin(seq_len(80) / 7))
  unit_truth <- 3 * data$x + cos(seq_len(nrow(data)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  folds <- function(d) {
    held <- list(seq_len(floor(nrow(d) / 3)), seq.int(floor(nrow(d) / 3) + 1L, nrow(d)))
    rsample::manual_rset(lapply(held, function(rows) rsample::make_splits(list(
      analysis = setdiff(seq_len(nrow(d)), rows), assessment = rows), d)), c("small", "large"))
  }
  for (scale in c(1, 1e-9, 1e9, 0)) {
    data$truth <- scale * unit_truth
    selection_split <- rsample::make_splits(list(analysis = 1:40, assessment = 41:80), data)
    config <- list(strategy = "two_stage", selection_split = selection_split,
      calibration = list(method = "grid", grid = data.frame(level = c(.4, .6)),
        resamples = withr::with_seed(7, rsample::vfold_cv(data[41:80, ], v = 2))))
    fitted <- generics::fit(calibrated_workflow(base, "truth", method = "cross", splits = folds,
      epsilon = tune::tune("level"), tuning = config, seed = 8), data)
    # Training rows are calibration scores, where the last bit of a prediction can depend
    # on the batch size (OpenBLAS) and differently at each scale; predict unseen rows.
    result <- venn_abers_details(fitted, data.frame(x = sin((seq_len(80) + .5) / 7)))
    values <- cbind(result$prediction$.pred, result$aggregate$lower, result$aggregate$upper)
    if (scale == 1) reference <- values
    expected <- scale * reference
    # Unit scale is transformed with the targets, never fitted to the error.
    testthat::expect_true(all(abs(values - expected) <=
      1e-10 * (if (scale == 0) 1 else scale) + 1e-10 * abs(expected)))
    state <- workflows::extract_fit_engine(fitted)
    testthat::expect_false(identical(unname(state$parameters$effective$m[1L]),
      unname(state$parameters$effective$m[2L])))
    details <- split(result$details, result$details$fold)
    for (i in seq_along(c("prediction", "lower", "upper"))) {
      column <- c("prediction", "lower", "upper")[i]
      mean_folds <- (details$small[[column]] + details$large[[column]]) / 2
      testthat::expect_true(all(abs(values[, i] - mean_folds) <=
        1e-10 * max(scale, 1e-9) + 1e-10 * abs(mean_folds)))
    }
  }
})

testthat::test_that("regression grids reconstruct IVAR/CVAR across strategies and model paths", {
  testthat::skip_if_not_installed("tune")
  data <- data.frame(x = sin(seq_len(240) / 7))
  data$truth <- data$x^2 + cos(seq_len(nrow(data))) / 5
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x, data) |>
      recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  for (strategy in c("joint", "two_stage")) for (method in c("inductive", "cross")) {
    for (mode in c("matched", "full")) {
      option <- if (mode == "matched") "m" else "epsilon"
      config <- list(strategy = strategy, selection_prop = .5,
        model = list(method = "grid", grid = data.frame(degree = c(1L, 2L)),
          resamples = function(d) rsample::vfold_cv(d, v = 2)),
        calibration = list(method = "grid", grid = data.frame(moderation =
          if (option == "m") c(1, 2) else c(.4, .6)),
          resamples = function(d) rsample::vfold_cv(d, v = 2)))
      args <- list(workflow = base, outcome = "truth", method = method,
        prediction_model = mode, tuning = config, seed = 31)
      args[[option]] <- tune::tune("moderation")
      if (method == "cross") args$v <- 2L
      fitted <- generics::fit(do.call(calibrated_workflow, args, quote = TRUE), data)
      state <- workflows::extract_fit_engine(fitted)
      selection <- venn_abers_tuning(fitted)
      final <- data[selection$calibration$final_rows, , drop = FALSE]
      testthat::expect_setequal(names(selection$models),
        c(names(state$folds), if (mode == "full") "Full1"))
      for (id in names(selection$models)) {
        search <- selection$models[[id]]
        local <- if (id == "Full1") final else final[state$folds[[id]]$analysis, , drop = FALSE]
        testthat::expect_equal(as.data.frame(search$native$splits[[1L]]$data), local,
          ignore_attr = TRUE)
        testthat::expect_true(search$chosen %in% search$candidates$id)
        testthat::expect_identical(search$seed, 31)
      }
      held <- state$calibration
      options <- setNames(list(state$parameters$requested[[option]]), option)
      core_args <- c(list(truth = unlist(lapply(held, `[[`, "truth")),
        prediction = unlist(lapply(held, `[[`, "prediction"))), options)
      cross <- method == "cross"
      if (cross) core_args$fold <- rep(names(held), lengths(lapply(held, `[[`, "truth")))
      core <- do.call(if (cross) cross_venn_abers_regression else venn_abers_regression, core_args)
      source <- lapply(state$sources, function(id) predict(state$models[[id]], data)$.pred)
      expected <- predict(core, if (cross) do.call(cbind, source) else source[[1L]])
      actual <- venn_abers_details(fitted, data)
      testthat::expect_true(all(abs(actual$prediction$.pred - expected$prediction) <=
        1e-10 + 1e-10 * abs(expected$prediction)))
      bounds <- if (cross) actual$aggregate else actual$details
      for (column in c("lower", "upper")) testthat::expect_true(all(
        abs(bounds[[column]] - expected[[column]]) <= 1e-10 + 1e-10 * abs(expected[[column]])))
      if (strategy == "two_stage") for (id in names(selection$calibration$resamples)) {
        choices <- lapply(selection$calibration$decisions,
          function(candidate) lapply(candidate[[id]], `[[`, "selected"))
        testthat::expect_identical(choices[[1L]], choices[[2L]])
      }
      # Rebuild one model-selection candidate at the public seam: joint scores
      # the calibrated predictor; stage one scores the uncalibrated workflow.
      search <- selection$models[[1L]]
      split <- search$native$splits[[1L]]
      candidate <- tune::finalize_workflow(base, list(degree = 1L))
      if (strategy == "joint") {
        candidate_args <- args
        candidate_args$workflow <- candidate
        candidate_args$tuning <- NULL
        candidate_args[[option]] <- options[[option]]
        candidate <- do.call(calibrated_workflow, candidate_args)
      }
      # Native tune workers use L'Ecuyer-CMRG, including serial workers.
      rebuilt <- withr::with_seed(31, generics::fit(candidate, rsample::analysis(split)),
        .rng_kind = "L'Ecuyer-CMRG")
      assessment <- rsample::assessment(split)
      expected_score <- sqrt(mean((assessment$truth - predict(rebuilt, assessment)$.pred)^2))
      testthat::expect_equal(search$scores[[1L]][[1L]], expected_score, tolerance = 1e-10)
      before <- selection$cost
      predict(fitted, data)
      testthat::expect_identical(venn_abers_tuning(fitted)$cost, before)
    }
  }
})

testthat::test_that("moderation grids retain S/F roles and reject locally invalid candidates", {
  testthat::skip_if_not_installed("tune")
  data <- data.frame(x = sin(seq_len(120) / 7))
  data$truth <- 3 * data$x + cos(seq_len(nrow(data)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  for (option in c("m", "epsilon")) {
    grid <- if (option == "m") c(1, 100) else c(.5, .001)
    config <- list(strategy = "joint", selection_prop = .5,
      calibration = list(method = "grid", grid = data.frame(moderation = grid),
        resamples = function(d) rsample::vfold_cv(d, v = 2)))
    args <- list(workflow = base, outcome = "truth", seed = 12, tuning = config)
    args[[option]] <- tune::tune("moderation")
    fitted <- generics::fit(do.call(calibrated_workflow, args, quote = TRUE), data)
    state <- workflows::extract_fit_engine(fitted)
    selected <- venn_abers_tuning(fitted)$calibration
    testthat::expect_equal(selected$grid$moderation, grid)
    saved <- tune::collect_predictions(selected$native)
    testthat::expect_true(option %in% names(saved))
    testthat::expect_true(all(saved[[option]] == grid[1L]))
    testthat::expect_identical(selected$chosen, "calibration_001")
    testthat::expect_identical(selected$candidates$status, c("eligible", "excluded"))
    testthat::expect_true(nrow(tune::collect_notes(selected$native)) > 0L)
    testthat::expect_match(selected$candidates$reason[2L], "Fold")
    testthat::expect_length(intersect(selected$selection_rows, selected$final_rows), 0L)
    testthat::expect_equal(sort(c(selected$selection_rows, selected$final_rows)), seq_len(nrow(data)))
    testthat::expect_equal(state$parameters$requested[[option]], grid[1L])
    held <- state$calibration$Inductive1
    core <- do.call(venn_abers_regression, c(list(truth = held$truth,
      prediction = held$prediction), setNames(list(grid[1L]), option)))
    expected <- predict(core, predict(state$models$Inductive1, data)$.pred)$prediction
    # The target's declared unit scale is 1, independent of observed error.
    testthat::expect_true(all(abs(predict(fitted, data)$.pred - expected) <=
      1e-10 + 1e-10 * abs(expected)))
    args$tuning$calibration$grid <- data.frame(moderation = grid[2L])
    testthat::expect_error(generics::fit(do.call(calibrated_workflow, args, quote = TRUE), data),
      "calibration size")
  }
})

testthat::test_that("regression candidate failures and report errors remain isolated", {
  testthat::skip_if_not_installed("tune")
  audit <- new.env()
  audit$fits <- list()
  numeric_fit <- function(x, y, variant) {
    audit$fits[[length(audit$fits) + 1L]] <- list(x = x, y = y, variant = variant)
    if (variant == 1 && nrow(x) < 40) stop("one failed local fold")
    if (variant == 5 && nrow(x) > 70) stop("selected final fit failed")
    list(mean = mean(y), variant = variant)
  }
  numeric_predict <- function(object, new_data) {
    if (object$variant == 2) return(rep(NA_real_, nrow(new_data)))
    if (object$variant == 3) return(rep(object$mean, nrow(new_data) - 1L))
    if (object$variant == 4) return(rep(Inf, nrow(new_data)))
    rep(object$mean, nrow(new_data))
  }
  register_reference_engine("linear_reg", "regression", "numeric_fit", "numeric_predict", "numeric",
    engine = "grid_reference", new_model = FALSE)
  spec <- parsnip::linear_reg() |>
    parsnip::set_engine("grid_reference", variant = tune::tune("variant"))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |> workflows::add_model(spec)
  data <- data.frame(x = seq_len(100), truth = rep(c(0, 1), 50))
  plans <- function(d) {
    held <- list(seq_len(floor(nrow(d) / 3)), seq.int(floor(nrow(d) / 3) + 1L, nrow(d)))
    rsample::manual_rset(lapply(held, function(rows) rsample::make_splits(list(
      analysis = setdiff(seq_len(nrow(d)), rows), assessment = rows), d)), c("small", "large"))
  }
  broken <- yardstick::new_numeric_metric(function(...) stop("report unavailable"), direction = "minimize")
  config <- list(strategy = "two_stage", model = list(method = "grid",
    grid = data.frame(variant = 0:4), resamples = plans, report_metrics = list(broken)))
  fitted <- generics::fit(calibrated_workflow(base, "truth", tuning = config, seed = 8), data)
  selection <- venn_abers_tuning(fitted)$models$Inductive1
  testthat::expect_identical(selection$candidates$status, c("eligible", rep("excluded", 4)))
  testthat::expect_identical(selection$chosen, "model_001")
  testthat::expect_match(selection$report_metrics[[1L]][[1L]][[1L]], "report unavailable")
  testthat::expect_match(paste(tune::collect_notes(selection$native)$note, collapse = " "), "one failed local fold")
  testthat::expect_true(all(is.finite(predict(fitted, data)$.pred)))
  fits_before <- length(audit$fits)
  predict(fitted, data)
  venn_abers_details(fitted, data)
  testthat::expect_equal(length(audit$fits), fits_before)
  config$model$grid <- data.frame(variant = c(0, 0))
  config$model$primary_metric <- yardstick::mae
  repeated <- generics::fit(calibrated_workflow(base, "truth", tuning = config, seed = 8), data)
  testthat::expect_identical(venn_abers_tuning(repeated)$models$Inductive1$chosen, "model_001")
  config$model$grid <- data.frame(variant = 1L)
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", tuning = config, seed = 8), data),
    "No eligible model")
  config$model$grid <- data.frame(variant = c(5L, 0L))
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", tuning = config, seed = 8), data),
    "selected final fit failed")
  testthat::expect_equal(sum(vapply(audit$fits, function(f) f$variant == 5 && nrow(f$x) > 70, logical(1))), 1L)
  config$model$grid <- data.frame(variant = 0L)
  config$model$primary_metric <- yardstick::new_numeric_metric(
    function(...) stop("primary unavailable"), direction = "minimize")
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", tuning = config, seed = 8), data),
    "primary unavailable")
})

testthat::test_that("regression selection stays inside outer training and survives native evaluation", {
  testthat::skip_if_not_installed("tune")
  withr::local_seed(9)
  data <- data.frame(x = sin(seq_len(160) / 7))
  data$truth <- 3 * data$x + cos(seq_len(nrow(data)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  config <- list(strategy = "two_stage", selection_prop = .5,
    calibration = list(method = "grid", grid = data.frame(moderation = c(1, 2)),
      resamples = function(d) rsample::vfold_cv(d, v = 2)))
  outer <- calibrated_workflow(base, "truth", method = "cross", v = 2,
    prediction_model = "full", m = tune::tune("moderation"), tuning = config, seed = 17)
  fitted <- generics::fit(outer, data)
  repeated <- generics::fit(outer, data)
  testthat::expect_identical(venn_abers_tuning(fitted)$calibration$chosen,
    venn_abers_tuning(repeated)$calibration$chosen)
  testthat::expect_equal(predict(fitted, data), predict(repeated, data))
  restored <- unserialize(serialize(fitted, NULL))
  testthat::expect_equal(predict(restored, data), predict(fitted, data))
  testthat::expect_equal(venn_abers_tuning(restored)$calibration$candidates,
    venn_abers_tuning(fitted)$calibration$candidates)
  saved <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  log <- tempfile(fileext = ".log")
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  saveRDS(list(fitted = fitted, data = data, prediction = predict(fitted, data),
    candidates = venn_abers_tuning(fitted)$calibration$candidates), saved)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "library(vennabers)", paste0("saved <- readRDS(", deparse(saved), ")"),
    'trace("tune_grid", where = asNamespace("tune"), tracer = quote(stop("Unexpected search")), print = FALSE)',
    'trace("fit.workflow", where = asNamespace("workflows"), tracer = quote(stop("Unexpected fit")), print = FALSE)',
    'stopifnot(identical(predict(saved$fitted, saved$data), saved$prediction),',
    '  identical(venn_abers_tuning(saved$fitted)$calibration$candidates, saved$candidates))'
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
  plan <- withr::with_seed(9, rsample::vfold_cv(data, v = 2, repeats = 2))
  result <- tune::fit_resamples(outer, resamples = plan, metrics = yardstick::metric_set(yardstick::rmse),
    control = tune::control_resamples(save_pred = TRUE, allow_par = FALSE,
      extract = function(w) venn_abers_tuning(w)))
  predictions <- tune::collect_predictions(result)
  testthat::expect_equal(nrow(predictions), 2L * nrow(data))
  testthat::expect_length(unique(predictions$.config), 1L)
  testthat::expect_setequal(unique(paste(predictions$id, predictions$id2)), paste(plan$id, plan$id2))
  testthat::expect_equal(anyDuplicated(predictions[c("id", "id2", ".config", ".row")]), 0L)
  testthat::expect_true(all(is.finite(predictions$.pred)))
  for (i in seq_len(nrow(plan))) {
    rows <- predictions$.row[predictions$id == plan$id[i] & predictions$id2 == plan$id2[i]]
    expected <- match(rsample::assessment(plan$splits[[i]])$x, data$x)
    testthat::expect_equal(sort(rows), sort(expected))
  }
  testthat::expect_equal(predictions$truth, data$truth[predictions$.row])
  metrics <- tune::collect_metrics(result, summarize = FALSE)
  testthat::expect_equal(nrow(metrics), nrow(plan))
  testthat::expect_true(all(is.finite(metrics$.estimate)))
  testthat::expect_equal(nrow(tune::collect_notes(result)), 0L)
  testthat::expect_equal(nrow(tune::collect_extracts(result)), 4L)
  train <- 1:120
  test <- 121:160
  final <- tune::last_fit(outer, rsample::make_splits(list(analysis = train, assessment = test), data),
    metrics = yardstick::metric_set(yardstick::rmse), control = tune::control_last_fit(allow_par = FALSE))
  changed <- data
  changed$truth[test] <- 1e6
  mutated <- tune::last_fit(outer, rsample::make_splits(list(analysis = train, assessment = test), changed),
    metrics = yardstick::metric_set(yardstick::rmse), control = tune::control_last_fit(allow_par = FALSE))
  testthat::expect_equal(tune::collect_predictions(final)$.pred, tune::collect_predictions(mutated)$.pred)
  testthat::expect_identical(venn_abers_tuning(tune::extract_workflow(final))$calibration$chosen,
    venn_abers_tuning(tune::extract_workflow(mutated))$calibration$chosen)
  # Changing S must not change fitted coefficients when the only allowed option
  # is fixed by the one-row grid. F is explicitly the same data in both runs.
  config$calibration$grid <- data.frame(moderation = 1)
  config$selection_prop <- NULL
  config$selection_split <- function(d) rsample::make_splits(
    list(analysis = 1:80, assessment = 81:160), d)
  outer <- calibrated_workflow(base, "truth", method = "cross", v = 2, m = tune::tune("moderation"),
    prediction_model = "full", tuning = config, seed = 17)
  fitted <- generics::fit(outer, data)
  changed <- data
  changed$truth[81:160] <- changed$truth[81:160] + 1000
  changed_fit <- generics::fit(outer, changed)
  testthat::expect_equal(predict(fitted, data), predict(changed_fit, data))
  # A foreign plan that smuggles in calibration rows is rejected at binding.
  recipe <- recipes::recipe(truth ~ x, data) |> recipes::step_poly(x, degree = tune::tune("degree"))
  tuned <- workflows::workflow() |> workflows::add_recipe(recipe) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  bad <- list(strategy = "two_stage", model = list(method = "grid", grid = data.frame(degree = 1L),
    resamples = rsample::vfold_cv(data, v = 2)))
  testthat::expect_error(generics::fit(calibrated_workflow(tuned, "truth", tuning = bad, seed = 17), data), "binding")
})

testthat::test_that("Grid and Bayes searches exclude S and each fold's held calibration rows", {
  testthat::skip_if_not_installed("tune")
  data <- withr::with_seed(84, data.frame(x = stats::rnorm(480), truth = stats::rnorm(480)))
  base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, data) |>
    recipes::step_normalize(x) |> recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  for (search in c("grid", "bayes")) for (strategy in c("joint", "two_stage")) {
    tuning <- list(strategy = strategy,
      selection_split = function(d) rsample::make_splits(list(analysis = 1:240, assessment = 241:480), d),
      model = list(method = "grid", grid = data.frame(degree = 1:2),
        resamples = function(d) rsample::vfold_cv(d, v = 2)),
      calibration = list(method = "grid", grid = data.frame(moderation = .3),
        resamples = function(d) rsample::vfold_cv(d, v = 2)))
    if (search == "bayes") {
      params <- function(id, range, type) dials::parameters(dials::new_quant_param(
        type = type, range = range, inclusive = c(TRUE, TRUE), label = stats::setNames(id, id)))
      tuning$model <- list(method = "bayes", initial = 2L, iter = 1L,
        param_info = params("degree", c(1L, 2L), "integer"),
        resamples = function(d) rsample::vfold_cv(d, v = 2),
        control = tune::control_bayes(seed = 4, save_pred = TRUE))
      # A tied selection metric fixes the chosen option while S is mutated.
      fixed <- yardstick::new_numeric_metric(function(data, truth, estimate, ...) {
        out <- yardstick::rmse(data, {{ truth }}, {{ estimate }}, ...)
        out$.metric <- "fixed"
        out$.estimate[] <- 0
        out
      }, direction = "minimize")
      tuning$calibration <- list(method = "bayes", initial = 2L, iter = 1L,
        param_info = params("moderation", c(.3, .4), "double"), primary_metric = fixed,
        resamples = function(d) rsample::vfold_cv(d, v = 2),
        control = tune::control_bayes(seed = 5, save_pred = TRUE, time_limit = 1e-9))
    }
    # Bind the same row roles in every local candidate, regardless of its size.
    splits <- function(d) {
      half <- seq_len(nrow(d) %/% 2L)
      rsample::manual_rset(lapply(list(half, setdiff(seq_len(nrow(d)), half)), function(rows)
        rsample::make_splits(list(analysis = setdiff(seq_len(nrow(d)), rows), assessment = rows), d)),
        c("left", "right"))
    }
    model <- calibrated_workflow(base, "truth", method = "cross", splits = splits,
      prediction_model = "full", epsilon = tune::tune("moderation"), tuning = tuning, seed = 31)
    fit <- generics::fit(model, data)
    diagnostic <- venn_abers_tuning(fit)
    if (search == "bayes") {
      testthat::expect_identical(diagnostic$calibration$chosen, "calibration_001")
      testthat::expect_identical(diagnostic$calibration$candidates$estimate, c(0, 0))
    }
    changed <- data
    changed$x[241:480] <- changed$x[241:480] + 100
    changed$truth[241:480] <- changed$truth[241:480] - 100
    other <- generics::fit(model, changed)
    testthat::expect_identical(venn_abers_details(fit, data[1:7, ]), venn_abers_details(other, data[1:7, ]))
    for (id in names(diagnostic$models)) {
      original <- diagnostic$models[[id]]
      mutated <- venn_abers_tuning(other)$models[[id]]
      testthat::expect_identical(original$candidates, mutated$candidates)
      testthat::expect_identical(original$selected, mutated$selected)
      training_rows <- if (id == "Full1") 1:240 else diagnostic$fit_roles[[id]]$analysis
      testthat::expect_equal(as.data.frame(original$native$splits[[1L]]$data), data[training_rows, ], ignore_attr = TRUE)
    }
    # C changes the calibrator and the sibling/full searches, but its own B
    # search (including learned preprocessing and nested candidates) is fixed.
    changed <- data
    changed$x[1:120] <- changed$x[1:120] + 20
    changed$truth[1:120] <- changed$truth[1:120] + 30
    other <- generics::fit(model, changed)
    original <- diagnostic$models$left
    mutated <- venn_abers_tuning(other)$models$left
    testthat::expect_identical(original$candidates, mutated$candidates)
    testthat::expect_identical(tune::collect_predictions(original$native), tune::collect_predictions(mutated$native))
    testthat::expect_false(identical(workflows::extract_fit_engine(fit)$calibration$left,
      workflows::extract_fit_engine(other)$calibration$left))
    testthat::expect_false(identical(diagnostic$models$Full1$candidates,
      venn_abers_tuning(other)$models$Full1$candidates))
  }
})

testthat::test_that("moderation metadata cannot overwrite regression outcome names", {
  testthat::skip_if_not_installed("tune")
  for (search in c("grid", "bayes")) for (option in c("m", "epsilon")) {
    data <- data.frame(x = sin(seq_len(120) / 7))
    data[[option]] <- 3 * data$x + cos(seq_len(nrow(data)))
    base <- workflows::workflow() |> workflows::add_formula(stats::reformulate("x", option)) |>
      workflows::add_model(parsnip::linear_reg(engine = "lm"))
    args <- list(workflow = base, outcome = option, seed = 19,
      tuning = list(strategy = "joint", selection_prop = .5,
        calibration = list(method = "grid", grid = data.frame(status = if (option == "m") 1 else .5),
          resamples = function(d) rsample::vfold_cv(d, v = 2))))
    args[[option]] <- tune::tune("status")
    if (search == "bayes") args$tuning$calibration <- list(method = "bayes", initial = 2L, iter = 2L,
      param_info = dials::parameters(dials::new_quant_param(
        type = if (option == "m") "integer" else "double",
        range = if (option == "m") c(1L, 2L) else c(.3, .4), inclusive = c(TRUE, TRUE),
        label = c(status = "Moderation"))), resamples = function(d) rsample::vfold_cv(d, v = 2),
      control = tune::control_bayes(seed = 4, save_pred = TRUE))
    fitted <- generics::fit(do.call(calibrated_workflow, args, quote = TRUE), data)
    selected <- venn_abers_tuning(fitted)$calibration
    predictions <- tune::collect_predictions(selected$native)
    testthat::expect_true(all(selected$candidates$status == "eligible"))
    if (search == "bayes") {
      testthat::expect_true(is.numeric(selected$selected$status))
      testthat::expect_true(selected$selected$status %in% selected$candidates$status_1)
      testthat::expect_false(anyDuplicated(selected$candidates$id) > 0L)
      testthat::expect_setequal(selected$candidates$status_1, tune::collect_metrics(selected$native)$status)
    }
    testthat::expect_equal(predictions[[option]], data[[option]][selected$selection_rows][predictions$.row])
    testthat::expect_true(paste0("va_", option) %in% names(predictions))
  }
})

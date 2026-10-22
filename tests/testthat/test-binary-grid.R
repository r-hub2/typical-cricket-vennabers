testthat::test_that("binary grid resolves a local model choice and predicts", {
  skip_if_no_integration("tune")
  data <- data.frame(x = rep(seq(-2, 2, length.out = 20), 3),
                     truth = factor(rep(c("yes", "no", "yes"), each = 20)))
  recipe <- recipes::recipe(truth ~ x, data = data) |>
    recipes::step_poly(x, degree = tune::tune("degree"))
  base <- workflows::workflow() |>
    workflows::add_recipe(recipe) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  config <- list(strategy = "two_stage", model = list(
    method = "grid", grid = data.frame(degree = c(1L, 2L)),
    resamples = function(d) rsample::vfold_cv(d, v = 2)))
  fitted <- generics::fit(calibrated_workflow(base, "truth", seed = 42,
                                              tuning = config), data)
  selected <- venn_abers_tuning(fitted)$models$Inductive1
  testthat::expect_true(selected$chosen %in% selected$candidates$id)
  testthat::expect_identical(venn_abers_tuning(fitted)$cost$observed_base_fits, 5L)
  testthat::expect_true(inherits(selected$native, "tune_results"))
  testthat::expect_setequal(tune::collect_metrics(selected$native)$candidate_id,
                           selected$candidates$id)
  testthat::expect_identical(selected$budget$gross_base_fit_estimate, 4)
  testthat::expect_identical(selected$cost$observed_base_fits, 4L)
  testthat::expect_setequal(tune::collect_predictions(selected$native)$.row,
                           unlist(lapply(selected$resamples, `[[`, "assessment")))
  testthat::expect_equal(nrow(predict(fitted, data, type = "prob")), nrow(data))
  repeated <- generics::fit(calibrated_workflow(base, "truth", seed = 42,
                                                tuning = config), data)
  testthat::expect_identical(venn_abers_tuning(repeated)$models$Inductive1$chosen,
                             selected$chosen)
  testthat::expect_equal(predict(repeated, data, type = "prob"),
                         predict(fitted, data, type = "prob"))
  skip_if_no_integration("dials")
  config$model$grid <- 2L
  config$model$param_info <- dials::parameters(dials::new_quant_param(
    type = "integer", range = c(1L, 2L), inclusive = c(TRUE, TRUE),
    label = c(degree = "Degree")))
  fitted <- generics::fit(calibrated_workflow(base, "truth", seed = 42,
                                              tuning = config), data)
  testthat::expect_equal(nrow(venn_abers_tuning(fitted)$models$Inductive1$candidates), 2L)
})

testthat::test_that("binary calibration grid keeps S out of F and matches the core", {
  skip_if_no_integration("tune")
  data <- data.frame(x = rep(seq(-1, 1, length.out = 20), 4),
                     truth = factor(rep(c("yes", "no", "yes", "no"), each = 20)))
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (event_level in c("first", "second")) {
    config <- list(strategy = "joint", selection_prop = .3, calibration = list(
      method = "grid", grid = data.frame(va_loss = c("log", "brier", "log", "brier"),
                                         va_rounding = c("none", "0", "2", "-2")),
      resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
    fitted <- generics::fit(calibrated_workflow(base, "truth", event_level = event_level,
      loss = tune::tune("va_loss"), precision = tune::tune("va_rounding"),
      seed = 23, tuning = config), data)
    state <- workflows::extract_fit_engine(fitted)
    diagnostic <- venn_abers_tuning(fitted)$calibration
    testthat::expect_length(intersect(diagnostic$selection_rows, diagnostic$final_rows), 0L)
    testthat::expect_equal(sort(c(diagnostic$selection_rows, diagnostic$final_rows)), seq_len(nrow(data)))
    testthat::expect_true(diagnostic$chosen %in% diagnostic$candidates$id)
    testthat::expect_identical(venn_abers_tuning(fitted)$cost$observed_base_fits, 9L)
    testthat::expect_null(state$parameters$requested$.fit_counter)
    source <- state$models$Inductive1
    calibration <- state$calibration$Inductive1
    core <- venn_abers(calibration$truth, calibration$probability,
                       precision = state$parameters$requested$precision)
    p <- predict(source, data, type = "prob")[[paste0(".pred_", state$event_class)]]
    expected <- predict(core, p, loss = state$loss)$probability
    actual <- predict(fitted, data, type = "prob")[[paste0(".pred_", state$event_class)]]
    testthat::expect_lte(max(abs(expected - actual)), 1e-10)
  }
})

testthat::test_that("joint and two-stage grids select Cross and full models locally", {
  skip_if_no_integration("tune")
  data <- data.frame(x = sin(seq_len(240) / 7),
                     truth = factor(rep(c("yes", "no", "no", "yes"), 60)))
  recipe <- recipes::recipe(truth ~ x, data = data) |>
    recipes::step_poly(x, degree = tune::tune("degree"))
  base <- workflows::workflow() |>
    workflows::add_recipe(recipe) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (strategy in c("joint", "two_stage")) {
    config <- list(strategy = strategy, selection_prop = .4,
      model = list(method = "grid", grid = data.frame(degree = c(1L, 2L)),
                   resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)),
      calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")),
                         resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
    fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth", method = "cross", v = 2,
      prediction_model = "full", loss = tune::tune("va_loss"), seed = 41,
      tuning = config), data))
    state <- workflows::extract_fit_engine(fitted)
    selection <- venn_abers_tuning(fitted)
    testthat::expect_setequal(names(selection$models), c(names(state$folds), "Full1"))
    testthat::expect_true(all(vapply(selection$models, function(x) x$chosen %in% x$candidates$id, logical(1))))
    testthat::expect_equal(length(selection$models$Full1$resamples[[1L]]$analysis) +
                           length(selection$models$Full1$resamples[[1L]]$assessment),
                           length(selection$calibration$final_rows))
    testthat::expect_true(all(state$sources == "Full1"))
    testthat::expect_true(all(is.finite(as.matrix(predict(fitted, data, type = "prob")))))
    calibration <- state$calibration
    truth <- unlist(lapply(calibration, `[[`, "truth"), use.names = FALSE)
    held_p <- unlist(lapply(calibration, `[[`, "probability"), use.names = FALSE)
    fold <- rep(names(calibration), lengths(lapply(calibration, `[[`, "truth")))
    core <- cross_venn_abers(truth, held_p, fold,
                             precision = state$parameters$requested$precision)
    full_p <- predict(state$models$Full1, data, type = "prob")[[paste0(".pred_", state$event_class)]]
    new_p <- matrix(rep(full_p, length(calibration)), nrow(data), length(calibration),
                    dimnames = list(NULL, names(calibration)))
    expected <- predict(core, new_p, loss = state$loss, type = "details")$prediction$probability
    actual <- predict(fitted, data, type = "prob")[[paste0(".pred_", state$event_class)]]
    testthat::expect_lte(max(abs(expected - actual)), 1e-10)
    if (strategy == "two_stage") for (id in names(selection$calibration$resamples)) {
      choices <- lapply(selection$calibration$decisions,
                        function(candidate) vapply(candidate[[id]], `[[`, "", "chosen"))
      testthat::expect_identical(choices[[1L]], choices[[2L]])
    }
  }
})

testthat::test_that("grid input and candidate failures are strict", {
  skip_if_no_integration("tune", "yardstick")
  data <- data.frame(x = sin(seq_len(80) / 4),
                     truth = factor(rep(c("yes", "no", "no", "yes"), 20)))
  recipe <- recipes::recipe(truth ~ x, data = data) |>
    recipes::step_poly(x, degree = tune::tune("degree"))
  base <- workflows::workflow() |>
    workflows::add_recipe(recipe) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  model <- list(method = "grid", grid = data.frame(degree = c(1L, 100L)),
                resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0))
  fit <- generics::fit(calibrated_workflow(base, "truth", seed = 8,
    tuning = list(strategy = "two_stage", model = model)), data)
  choices <- venn_abers_tuning(fit)$models$Inductive1
  testthat::expect_identical(choices$chosen, "model_001")
  testthat::expect_identical(choices$candidates$status, c("eligible", "excluded"))
  testthat::expect_identical(choices$cost$observed_base_fits, 4L)
  testthat::expect_match(choices$candidates$reason[2L], "Fold")
  broken_report <- yardstick::new_prob_metric(function(...) stop("report unavailable"),
                                               direction = "minimize")
  model$grid <- data.frame(degree = 1L)
  model$report_metrics <- list(broken_report)
  fit <- generics::fit(calibrated_workflow(base, "truth", seed = 8,
    tuning = list(strategy = "two_stage", model = model)), data)
  choices <- venn_abers_tuning(fit)$models$Inductive1
  testthat::expect_identical(choices$candidates$status, "eligible")
  testthat::expect_match(choices$report_metrics[[1L]][[1L]][[1L]], "report unavailable")
  model$report_metrics <- NULL
  model$grid <- data.frame(degree = c(1L, 1L))
  fit <- generics::fit(calibrated_workflow(base, "truth", seed = 8,
    tuning = list(strategy = "joint", model = model)), data)
  testthat::expect_identical(venn_abers_tuning(fit)$models$Inductive1$chosen, "model_001")
  model$control <- tune::control_grid(save_pred = FALSE)
  testthat::expect_error(calibrated_workflow(base, "truth",
    tuning = list(strategy = "joint", model = model)), "save_pred")
  model$control <- tune::control_grid(save_pred = TRUE, allow_par = FALSE)
  testthat::expect_s3_class(calibrated_workflow(base, "truth",
    tuning = list(strategy = "joint", model = model)), "workflow")
  model$control <- NULL
  testthat::expect_error(calibrated_workflow(base, "truth", loss = tune::tune("degree"),
    tuning = list(strategy = "joint", model = model, selection_prop = .2,
      calibration = list(method = "grid", grid = data.frame(degree = "log"),
        resamples = function(d) rsample::vfold_cv(d, v = 2)))), "unique")
  fixed <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  testthat::expect_error(calibrated_workflow(fixed, "truth",
    precision = tune::tune("rounding"),
    tuning = list(strategy = "joint", selection_prop = .2,
      calibration = list(method = "grid", grid = data.frame(rounding = "1"),
        resamples = function(d) rsample::vfold_cv(d, v = 2)))), "precision")
  bad_split <- rsample::initial_split(data[-1L, ], prop = .7)
  config <- list(strategy = "joint", selection_split = bad_split,
    calibration = list(method = "grid", grid = data.frame(va_loss = "log"),
      resamples = function(d) rsample::vfold_cv(d, v = 2)))
  testthat::expect_error(generics::fit(calibrated_workflow(fixed, "truth",
    loss = tune::tune("va_loss"), tuning = config), data), "binding")
})

testthat::test_that("a failed selected final fit has no fallback candidate", {
  skip_if_no_integration("tune")
  fit_sizes <- integer()
  grid_final_fit <- function(x, y) {
    fit_sizes <<- c(fit_sizes, nrow(x))
    if (nrow(x) > 50L) stop("planned final fit failure")
    structure(list(levels = levels(y)), class = "va_grid_final")
  }
  grid_final_prob <- function(object, new_data) {
    p <- data.frame(a = rep(.5, nrow(new_data)), b = rep(.5, nrow(new_data)))
    names(p) <- paste0(".pred_", object$levels)
    p
  }
  register_reference_engine("va_grid_final", "classification", "grid_final_fit", "grid_final_prob", "prob")
  spec <- parsnip::new_model_spec("va_grid_final", args = list(), eng_args = NULL,
                                  mode = "classification", method = NULL, engine = "reference")
  data <- data.frame(x = sin(seq_len(100) / 8), truth = factor(rep(c("yes", "no"), 50)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |> workflows::add_model(spec)
  calibration <- list(method = "grid", grid = data.frame(va_loss = c("log", "brier")),
                resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0))
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", seed = 3,
    loss = tune::tune("va_loss"), tuning = list(strategy = "joint", selection_prop = .3,
      calibration = calibration)), data), "planned final fit failure")
  testthat::expect_equal(sum(fit_sizes > 50L), 1L)
})

testthat::test_that("binary grid survives native resampling and serialization", {
  skip_if_no_integration("tune", "yardstick")
  data <- data.frame(x = sin(seq_len(96) / 6),
                     truth = factor(rep(c("yes", "no", "no", "yes"), 24)))
  recipe <- recipes::recipe(truth ~ x, data) |>
    recipes::step_poly(x, degree = tune::tune("degree"))
  base <- workflows::workflow() |> workflows::add_recipe(recipe) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  config <- list(strategy = "two_stage", model = list(
    method = "grid", grid = data.frame(degree = c(1L, 2L)),
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
  outer <- calibrated_workflow(base, "truth", tuning = config, seed = 16)
  fitted <- generics::fit(outer, data)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fitted, path)
  restored <- readRDS(path)
  testthat::expect_identical(venn_abers_tuning(restored)$models$Inductive1$chosen,
                             venn_abers_tuning(fitted)$models$Inductive1$chosen)
  testthat::expect_equal(predict(restored, data, type = "prob"),
                         predict(fitted, data, type = "prob"))
  splits <- rsample::vfold_cv(data, v = 2, strata = truth, pool = 0)
  result <- tune::fit_resamples(outer, resamples = splits,
    metrics = yardstick::metric_set(yardstick::mn_log_loss),
    control = tune::control_resamples(save_pred = TRUE,
      extract = function(w) venn_abers_tuning(w), allow_par = FALSE))
  testthat::expect_equal(nrow(tune::collect_predictions(result)), nrow(data))
  testthat::expect_equal(nrow(tune::collect_extracts(result)), 2L)
  split <- rsample::initial_split(data, prop = .8, strata = truth, pool = 0)
  final <- tune::last_fit(outer, split,
    metrics = yardstick::metric_set(yardstick::mn_log_loss),
    control = tune::control_last_fit(allow_par = FALSE))
  testthat::expect_equal(nrow(tune::collect_predictions(final)), nrow(rsample::testing(split)))
})

testthat::test_that("a class metric uses the second event on tied probabilities", {
  skip_if_no_integration("tune", "yardstick")
  tie_fit <- function(x, y) structure(list(levels = levels(y)), class = "va_grid_tie")
  tie_prob <- function(object, new_data) {
    p <- data.frame(a = rep(.5, nrow(new_data)), b = rep(.5, nrow(new_data)))
    names(p) <- paste0(".pred_", object$levels)
    p
  }
  register_reference_engine("va_grid_tie", "classification", "tie_fit", "tie_prob", "prob")
  spec <- parsnip::new_model_spec("va_grid_tie", args = list(), eng_args = NULL,
                                  mode = "classification", method = NULL, engine = "reference")
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |> workflows::add_model(spec)
  data <- data.frame(x = seq_len(100),
                     truth = factor(c(rep("yes", 20), rep("no", 20),
                                      rep("yes", 20), rep("no", 40)), levels = c("yes", "no")))
  selection_split <- rsample::make_splits(list(analysis = 1:40, assessment = 41:100), data)
  local_resamples <- function(d) {
    train <- c(which(d$truth == "yes")[1:10], which(d$truth == "no")[1:10])
    rsample::manual_rset(list(rsample::make_splits(list(
      analysis = train, assessment = setdiff(seq_len(nrow(d)), train)), d)), "Holdout1")
  }
  config <- list(strategy = "joint", selection_split = selection_split,
    calibration = list(method = "grid", grid = data.frame(va_loss = "log"),
                       resamples = local_resamples, primary_metric = yardstick::accuracy))
  fitted <- generics::fit(calibrated_workflow(base, "truth", event_level = "second",
    loss = tune::tune("va_loss"), tuning = config, seed = 2), data)
  testthat::expect_equal(venn_abers_tuning(fitted)$calibration$scores[[1L]][[1L]], .75)
  testthat::expect_true(all(predict(fitted, data, type = "class")$.pred_class == "no"))
})

testthat::test_that("shared multiclass calibration grid scores complete final matrices", {
  skip_if_no_integration("tune", "yardstick", "nnet")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ .) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  grid <- data.frame(va_loss = c("log", "brier", "log", "brier"),
                     rounding = c("none", "0", "2", "-2"))
  broken_report <- yardstick::new_prob_metric(function(...) stop("report unavailable"),
                                               direction = "minimize")
  config <- list(strategy = "joint", selection_prop = .4, calibration = list(
    method = "grid", grid = grid,
    report_metrics = list(broken_report),
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
  outer <- calibrated_workflow(base, "truth", classification = "ovo",
    loss = tune::tune("va_loss"), precision = tune::tune("rounding"),
    tuning = config, seed = 21)
  fitted <- generics::fit(outer, data)
  selection <- venn_abers_tuning(fitted)$calibration
  testthat::expect_identical(selection$primary_metric, yardstick::mn_log_loss)
  testthat::expect_identical(selection$candidates$status, rep("eligible", 4L))
  testthat::expect_match(selection$report_metrics[[1L]][[1L]][[1L]], "report unavailable")
  predicted <- predict(fitted, data, type = "prob")
  testthat::expect_identical(names(predicted), paste0(".pred_", levels(data$truth)))
  testthat::expect_lte(max(abs(rowSums(predicted) - 1)), 1e-8)
  native <- tune::collect_predictions(selection$native)
  for (id in selection$candidates$id) for (fold in names(selection$resamples)) {
    rows <- selection$resamples[[fold]]$assessment
    candidate <- native[native$candidate_id == id & native$id == fold, ]
    candidate <- candidate[match(rows, candidate$.row), ]
    matrix <- as.matrix(candidate[paste0(".pred_", levels(data$truth))])
    colnames(matrix) <- levels(data$truth)
    manual <- log_loss(candidate$truth, matrix)
    testthat::expect_equal(selection$scores[[id]][[fold]], manual, tolerance = 1e-10)
    native_brier <- rlang::inject(yardstick::brier_class(candidate, truth = truth,
      !!!rlang::syms(paste0(".pred_", levels(data$truth)))))$.estimate
    testthat::expect_equal(native_brier, brier_score(candidate$truth, matrix) / 2,
                           tolerance = 1e-10)
  }
  resampled <- tune::fit_resamples(outer,
    resamples = rsample::vfold_cv(data, v = 2, strata = truth, pool = 0),
    metrics = yardstick::metric_set(yardstick::mn_log_loss),
    control = tune::control_resamples(save_pred = TRUE, allow_par = FALSE,
      extract = function(w) venn_abers_tuning(w)))
  testthat::expect_equal(nrow(tune::collect_predictions(resampled)), nrow(data))
  testthat::expect_equal(nrow(tune::collect_extracts(resampled)), 2L)
})

testthat::test_that("shared multiclass model grids choose local Cross and full settings", {
  skip_if_no_integration("tune", "nnet")
  data <- iris
  names(data)[5L] <- "truth"
  recipe <- recipes::recipe(truth ~ ., data) |>
    recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))
  base <- workflows::workflow() |>
    workflows::add_recipe(recipe) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  for (strategy in c("joint", "two_stage")) {
    config <- list(strategy = strategy, selection_prop = .4,
      model = list(method = "grid", grid = data.frame(degree = c(1L, 2L)),
        resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)),
      calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")),
        resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
    fitted <- generics::fit(calibrated_workflow(base, "truth", classification = "ovr",
      method = "cross", v = 2L, prediction_model = "full", seed = 19,
      loss = tune::tune("va_loss"), tuning = config), data)
    state <- workflows::extract_fit_engine(fitted)
    selection <- venn_abers_tuning(fitted)
    testthat::expect_identical(selection$classes, levels(data$truth))
    testthat::expect_identical(selection$classification, "ovr")
    testthat::expect_identical(selection$sources, state$sources)
    testthat::expect_identical(selection$fit_roles, state$folds)
    testthat::expect_identical(selection$models$Full1$tie_break,
                               "first candidate in grid order")
    testthat::expect_setequal(names(selection$models), c(names(state$folds), "Full1"))
    testthat::expect_true(all(vapply(selection$models,
      function(x) x$chosen %in% x$candidates$id, logical(1))))
    testthat::expect_identical(unname(state$sources), rep("Full1", length(state$folds)))
    testthat::expect_equal(length(selection$models$Full1$resamples[[1L]]$analysis) +
      length(selection$models$Full1$resamples[[1L]]$assessment), length(selection$calibration$final_rows))
    predicted <- predict(fitted, data, type = "prob")
    testthat::expect_identical(names(predicted), paste0(".pred_", levels(data$truth)))
    testthat::expect_lte(max(abs(rowSums(predicted) - 1)), 1e-8)
    calibration <- state$calibration
    truth <- do.call(c, lapply(calibration, `[[`, "truth"))
    held <- do.call(rbind, lapply(calibration, `[[`, "probability"))
    fold <- rep(names(calibration), lengths(lapply(calibration, `[[`, "truth")))
    core <- cross_venn_abers(truth, held, fold, method = "ovr",
      precision = state$parameters$requested$precision)
    full <- predict(state$models$Full1, data, type = "prob")
    matrix <- as.matrix(full)
    colnames(matrix) <- levels(data$truth)
    expected <- predict(core, stats::setNames(rep(list(matrix), length(state$folds)),
                                               names(state$folds)), loss = state$loss)
    testthat::expect_lte(max(abs(as.matrix(predicted) - expected)), 1e-10)
  }
})

testthat::test_that("pairwise tuning remains explicit in multiclass search", {
  skip_if_no_integration("tune", "nnet")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ .) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  pairwise <- calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", loss = tune::tune("va_loss"),
    tuning = list(strategy = "joint", selection_prop = .3,
      calibration = list(method = "grid", grid = data.frame(va_loss = "log"),
        resamples = function(d) rsample::vfold_cv(d, v = 2))))
  requested <- rlang::eval_tidy(workflows::extract_spec_parsnip(pairwise)$eng_args$parameters)$requested
  testthat::expect_identical(requested$model_training, "pairwise")
})

testthat::test_that("two and four class grids retain the multiclass core contract", {
  skip_if_no_integration("tune", "yardstick", "nnet")
  fourth <- iris[iris$Species == "setosa", ]
  fourth$Sepal.Length <- fourth$Sepal.Length + 1.2
  fourth$Species <- "fourth"
  for (k in c(2L, 4L)) {
    data <- if (k == 2L) iris[iris$Species != "virginica", ] else rbind(iris, fourth)
    if (k == 4L) data <- rbind(data, data)
    data$truth <- factor(as.character(data$Species), levels =
      if (k == 2L) c("setosa", "versicolor") else c(levels(iris$Species), "fourth"))
    data$Species <- NULL
    base <- workflows::workflow() |>
      workflows::add_formula(truth ~ Sepal.Length) |>
      workflows::add_model(if (k == 2L) parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(engine = "nnet"))
    tuned <- workflows::workflow() |>
      workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
        recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
      workflows::add_model(if (k == 2L) parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(engine = "nnet"))
    for (classification in c("ovo", "ovr")) {
      cross <- k == 4L
      mode <- if (classification == "ovo") "matched" else "full"
      config <- list(strategy = if (classification == "ovo") "joint" else "two_stage",
        selection_prop = .3, calibration = list(method = "grid",
          grid = data.frame(va_loss = c("log", "brier")),
          resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
      if (classification == "ovo") config$model <- list(method = "grid",
        grid = data.frame(degree = c(1L, 2L)),
        resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0))
      if (classification == "ovr") config$calibration$primary_metric <- yardstick::brier_class
      precision <- if (classification == "ovo") 0 else -2
      candidate <- if (classification == "ovo") tuned else base
      outer <- if (cross) calibrated_workflow(candidate, "truth", classification = classification,
        prediction_model = mode, loss = tune::tune("va_loss"), precision = precision,
        method = "cross", v = 2L, tuning = config, seed = 31) else
        calibrated_workflow(candidate, "truth", classification = classification,
          prediction_model = mode, loss = tune::tune("va_loss"), precision = precision,
          tuning = config, seed = 31)
      fitted <- generics::fit(outer, data)
      state <- workflows::extract_fit_engine(fitted)
      held <- state$calibration
      truth <- do.call(c, lapply(held, `[[`, "truth"))
      probability <- do.call(rbind, lapply(held, `[[`, "probability"))
      core <- if (cross) cross_venn_abers(truth, probability,
        rep(names(held), lengths(lapply(held, `[[`, "truth"))),
        method = classification, precision = state$parameters$requested$precision) else
        venn_abers(truth, probability, method = classification,
          precision = state$parameters$requested$precision)
      source <- lapply(state$sources, function(id) {
        matrix <- as.matrix(predict(state$models[[id]], data, type = "prob"))
        colnames(matrix) <- levels(data$truth)
        matrix
      })
      expected <- predict(core, if (cross) source else source[[1L]], loss = state$loss)
      actual <- as.matrix(predict(fitted, data, type = "prob"))
      testthat::expect_lte(max(abs(actual - expected)), 1e-10)
      testthat::expect_identical(colnames(actual), paste0(".pred_", levels(data$truth)))
      testthat::expect_lte(max(abs(rowSums(actual) - 1)), 1e-8)
      selection <- venn_abers_tuning(fitted)$calibration
      testthat::expect_true(selection$chosen %in% selection$candidates$id)
      if (classification == "ovo") testthat::expect_true(all(vapply(
        venn_abers_tuning(fitted)$models,
        function(x) x$chosen %in% x$candidates$id, logical(1))))
      if (classification == "ovr") {
        native <- tune::collect_predictions(selection$native)
        rows <- selection$resamples[[1L]]$assessment
        candidate <- native[native$candidate_id == selection$candidates$id[1L] &
          native$id == names(selection$resamples)[1L], ]
        candidate <- candidate[match(rows, candidate$.row), ]
        probability <- as.matrix(candidate[paste0(".pred_", levels(data$truth))])
        colnames(probability) <- levels(data$truth)
        testthat::expect_equal(selection$scores[[1L]][[1L]],
          brier_score(candidate$truth, probability) / 2, tolerance = 1e-10)
      }
    }
  }
})

testthat::test_that("regression-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")
if (all(vapply(c("parsnip", "workflows", "rsample", "recipes", "withr"),
               requireNamespace, logical(1), quietly = TRUE))) {
  # Fixed before comparison: one outcome unit; multiply s when changing units.
  close <- function(actual, expected, s = 1) {
    assert_conditions(identical(dim(actual), dim(expected)), length(actual) == length(expected),
              all(is.finite(actual)), all(is.finite(expected)),
              all(abs(actual - expected) <= 1e-10 * s + 1e-10 * abs(expected)))
  }
  data <- data.frame(x = sin(1:80), truth = 10 + cos(1:80) + 3 * sin(1:80))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  plan <- withr::with_seed(42, rsample::initial_split(data, .8))
  held <- as.integer(plan, data = "assessment")
  for (mode in c("matched", "full")) {
    fitted <- generics::fit(calibrated_workflow(base, "truth", prediction_model = mode, seed = 42), data)
    state <- workflows::extract_fit_engine(fitted)
    actual <- venn_abers_details(fitted, data)
    stopifnot(identical(state$folds$Inductive1$assessment, held),
      all(actual$details$prediction_model == mode),
      identical(names(state$models), if (mode == "full") "Full1" else "Inductive1"))
    if (mode == "full") next
    expected <- predict(venn_abers_regression(data$truth[held],
      predict(state$models$Inductive1, data[held, ])$.pred),
      predict(state$models$Inductive1, data)$.pred)
    close(actual$prediction$.pred, expected$prediction)
    close(as.matrix(actual$details[c("prediction", "lower", "upper")]), as.matrix(expected))
    stopifnot(identical(actual$prediction, predict(fitted, data, type = "numeric")),
      is.double(actual$prediction$.pred), identical(names(actual$prediction), ".pred"))
  }
}

if (all(vapply(c("parsnip", "workflows", "rsample", "recipes", "withr"),
               requireNamespace, logical(1), quietly = TRUE))) {
  error <- function(expr, text = "") {
    e <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(e, "error"), grepl(text, conditionMessage(e), fixed = TRUE))
  }
  counts <- new.env()
  counts$fit <- counts$predict <- 0L
  regression_fit <- function(x, y, output = "normal") {
    counts$fit <- counts$fit + 1L
    if (output == "fail") stop("instrumented training failure")
    structure(list(x = x, y = y, output = output, offset = mean(y), scale = sd(y)),
              class = "regression_reference")
  }
  regression_predict <- function(object, new_data) {
    counts$predict <- counts$predict + 1L
    stopifnot(nrow(new_data) > 0L)
    p <- if (object$output == "raw") new_data$id else object$offset + object$scale * new_data$x
    if (object$output == "nonfinite") p[1L] <- Inf
    p
  }
  register_reference_engine("regression_reference", "regression", "regression_fit", "regression_predict", "numeric")
  spec <- parsnip::new_model_spec("regression_reference", args = list(), eng_args = NULL,
                                 mode = "regression", method = NULL, engine = "reference")
  make_base <- function(data, output = "normal") workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x + id, data) |>
      recipes::step_impute_median(x) |> recipes::step_normalize(x)) |>
    workflows::add_model(parsnip::set_engine(spec, "reference", output = output))
  make_plan <- function(data, held = list(1:3, 4:11, 12:30)) {
    rsample::manual_rset(lapply(held, function(h) rsample::make_splits(
      list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data)),
      paste0("fold", seq_along(held)))
  }
  # Unit scale is fixed at s=1 for the reference truth, then scaled explicitly.
  original <- data.frame(x = c(100, NA, sin(3:30)), id = 1:30, truth = 2 + cos(1:30) + sin(1:30))
  for (scale in c(1e-9, 1, 1e9)) for (shift in c(0, 100)) {
    data <- original
    data$truth <- scale * (data$truth + shift)
    base <- make_base(data)
    for (method in c("inductive", "cross")) for (option in list(list(), list(m = 1), list(epsilon = .5))) {
      plan <- make_plan(data)
      if (method == "inductive") plan <- plan$splits[[2L]]
      before <- counts$fit
      fitted <- generics::fit(do.call(calibrated_workflow,
        c(list(workflow = base, outcome = "truth", method = method, splits = plan), option)), data)
      state <- workflows::extract_fit_engine(fitted)
      assert_conditions(counts$fit == before + length(state$folds),
                identical(names(state$models), names(state$folds)),
                state$parameters$requested$mode == "regression",
                is.null(state$class_levels), is.null(state$event_class), is.null(state$loss),
                is.null(state$parameters$effective$classification), is.null(state$parameters$effective$loss))
      calibration <- truth <- fold <- list()
      new_p <- matrix(NA_real_, nrow(data), length(state$folds), dimnames = list(NULL, names(state$folds)))
      for (id in names(state$folds)) {
        f <- state$folds[[id]]
        trained <- workflows::extract_fit_engine(state$models[[id]])
        recipe <- workflows::extract_recipe(state$models[[id]])
        x <- data$x[f$analysis]
        x[is.na(x)] <- median(x, na.rm = TRUE)
        close(recipe$steps[[1L]]$medians[["x"]], median(x))
        close(recipe$steps[[2L]]$means[["x"]], mean(x))
        close(recipe$steps[[2L]]$sds[["x"]], sd(x))
        stopifnot(identical(trained$x$id, data$id[f$analysis]),
                  identical(as.vector(trained$y), data$truth[f$analysis]),
                  f$training_size == length(f$analysis), f$calibration_size == length(f$assessment))
        calibration[[id]] <- predict(state$models[[id]], data[f$assessment, ])$.pred
        truth[[id]] <- data$truth[f$assessment]
        fold[[id]] <- rep(id, f$calibration_size)
        new_p[, id] <- predict(state$models[[id]], data)$.pred
        expected_m <- if (length(option$epsilon)) floor(option$epsilon * (f$calibration_size + 1) / 2) else 1
        stopifnot(f$moderation$m == expected_m,
                  state$parameters$effective$m[[id]] == expected_m,
                  f$moderation$epsilon == 2 * expected_m / (f$calibration_size + 1),
                  identical(f$moderation$requested, option))
      }
      if (method == "cross") {
        core <- do.call(cross_venn_abers_regression,
          c(list(truth = unlist(truth), prediction = unlist(calibration), fold = unlist(fold)), option))
        expected <- predict(core, new_p, type = "details")
        points <- expected$prediction
      } else {
        core <- do.call(venn_abers_regression, c(list(truth = truth[[1L]], prediction = calibration[[1L]]), option))
        points <- predict(core, new_p[, 1L])
      }
      actual <- venn_abers_details(fitted, data)
      close(actual$prediction$.pred, points$prediction, scale)
      assert_conditions(identical(actual$prediction, predict(fitted, data)),
                identical(actual$prediction, predict(fitted, data, type = "numeric")),
                identical(names(actual$details), c("row", "fold", "prediction", "lower", "upper", "prediction_model", "model_source")),
                all(actual$details$prediction_model == "matched"),
                identical(actual$details$model_source, actual$details$fold))
      if (method == "cross") {
        for (id in names(state$folds)) close(as.matrix(actual$details[actual$details$fold == id,
          c("prediction", "lower", "upper")]), as.matrix(expected$folds[[id]]), scale)
        close(as.matrix(actual$aggregate[c("lower", "upper")]), as.matrix(points[c("lower", "upper")]), scale)
        # Detect wrong model assignment and fold-size weighting separately.
        if (scale == 1 && shift == 0 && !length(option)) {
          wrong <- new_p[, 3:1]; colnames(wrong) <- colnames(new_p)
          stopifnot(any(abs(predict(core, wrong)$prediction - points$prediction) > 1e-5),
            any(abs(as.vector(matrix(actual$details$prediction, nrow(data)) %*% c(3, 8, 19) / 30) - points$prediction) > 1e-5))
        }
      } else close(as.matrix(actual$details[c("prediction", "lower", "upper")]), as.matrix(points), scale)
      fit_count <- counts$fit
      pred_count <- counts$predict
      for (n in 0:1) {
        small <- data[seq_len(n), ]
        details <- venn_abers_details(fitted, small)
        assert_conditions(identical(details$prediction, predict(fitted, small)),
                  nrow(details$details) == n * length(state$folds),
                  is.integer(details$details$row), is.character(details$details$fold),
                  all(vapply(details$details[c("prediction", "lower", "upper")], is.double, logical(1))))
        if (!n) stopifnot(counts$predict == pred_count)
        if (method == "cross") stopifnot(identical(names(details$aggregate),
          c("row", "lower", "upper", "prediction_model", "model_source")), nrow(details$aggregate) == n,
          identical(details$aggregate$model_source, rep(list(names(state$models)), n)))
      }
      batched <- do.call(rbind, lapply(unname(split(seq_len(nrow(data)), rep(1:3, each = 10))),
        function(rows) predict(fitted, data[rows, ])))
      stopifnot(identical(batched, actual$prediction), counts$fit == fit_count)
      for (type in c("prob", "class", "conf_int", "pred_int")) error(predict(fitted, data, type = type), "type")
      error(predict(fitted, data, loss = "log"), "Regression")
      error(venn_abers_details(fitted, data, loss = NULL), "Regression")
    }
  }

  data <- original
  base <- make_base(data)
  for (method in c("inductive", "cross")) {
    auto <- generics::fit(calibrated_workflow(base, "truth", method = method, seed = 42), data)
    automatic <- withr::with_seed(42, if (method == "cross") rsample::vfold_cv(data, v = 5) else rsample::initial_split(data, .8))
    explicit <- generics::fit(calibrated_workflow(base, "truth", method = method,
      splits = local_calibration_splits(list(automatic), names(data))), data)
    stopifnot(identical(venn_abers_details(auto, data), venn_abers_details(explicit, data)),
              identical(workflows::extract_fit_engine(auto)$parameters$effective$strata, FALSE))
    # Constant outcomes exercise the zero-weight midpoint rule; 0/1 stays numeric.
    for (truth in list(rep(7, 30), rep(c(0, 1), 15))) {
      d <- data; d$truth <- truth
      f <- generics::fit(calibrated_workflow(make_base(d), "truth", method = method,
                                             m = NULL, epsilon = NULL, seed = 42), d)
      stopifnot(is.double(predict(f, d)$.pred), is.null(workflows::extract_fit_engine(f)$class_levels))
      if (all(truth == 7)) stopifnot(all(predict(f, d)$.pred == 7))
    }
  }
  # Nonconstant zero-weight intervals have an independently known midpoint.
  zero <- data.frame(x = sin(1:8), id = rep(as.double(1:4), 2), truth = rep(c(0, 0, 2, 2), 2))
  for (method in c("inductive", "cross")) {
    plan <- make_plan(zero, list(1:4, 5:8))
    if (method == "inductive") plan <- plan$splits[[1L]]
    fitted <- generics::fit(calibrated_workflow(make_base(zero, "raw"), "truth",
                                               method = method, splits = plan), zero)
    details <- venn_abers_details(fitted, data.frame(x = 0, id = 2.5))
    stopifnot(all(details$details$lower == 0), all(details$details$upper == 2),
              all(details$details$prediction == 1), identical(details$prediction$.pred, 1))
  }
  # Explicit m above the default is forwarded unchanged when every fold permits it.
  moderated <- generics::fit(calibrated_workflow(base, "truth", method = "cross",
    splits = make_plan(data, list(1:5, 6:30)), m = 2), data)
  stopifnot(all(workflows::extract_fit_engine(moderated)$parameters$effective$m == 2))
  before <- counts$fit
  for (option in list(list(classification = "binary"), list(event_level = "first"), list(precision = NULL),
                      list(loss = "log"), list(m = 1, epsilon = .5), list(m = TRUE), list(m = 1.5),
                      list(m = numeric()), list(epsilon = 0), list(epsilon = 1), list(epsilon = NA_real_))) {
    error(do.call(calibrated_workflow, c(list(workflow = base, outcome = "truth"), option)))
  }
  for (truth in list(rep(NA_real_, 30), rep(Inf, 30), rep(TRUE, 30), rep("1", 30), factor(rep(1:2, 15)))) {
    d <- data; d$truth <- truth
    error(generics::fit(calibrated_workflow(base, "truth"), d), "regression truth")
  }
  for (option in list(list(m = 2), list(epsilon = .5 - 1e-12), list(epsilon = .999))) {
    error(generics::fit(do.call(calibrated_workflow,
      c(list(workflow = base, outcome = "truth", method = "cross", splits = make_plan(data)), option)), data), if (identical(option, list(epsilon = .999))) "fold2 calibration size 8" else "fold1 calibration size 3")
  }
  # An invalid last fold must prevent training of the valid first fold too.
  error(generics::fit(calibrated_workflow(base, "truth", method = "cross",
    splits = make_plan(data, list(1:28, 29:30))), data), "fold2 calibration size 2")
  error(generics::fit(calibrated_workflow(base, "truth", calibration_prop = .01), data), "calibration size")
  stopifnot(counts$fit == before)
  error(generics::fit(calibrated_workflow(make_base(data, "fail"), "truth", seed = 42), data), "Inductive1 base workflow training")
  error(generics::fit(calibrated_workflow(make_base(data, "nonfinite"), "truth", seed = 42), data), "base regression prediction")
}

if (all(vapply(c("tune", "yardstick", "rsample", "workflows", "recipes"),
               requireNamespace, logical(1), quietly = TRUE))) {
  example <- new.env()
  sys.source(system.file("examples", "calibrated-resampling.R", package = "vennabers"), example)
  data <- data.frame(x = c(100, NA, sin(3:100)), id = 1:100, truth = 2 + cos(1:100) + sin(1:100))
  base <- make_base(data)
  outer <- withr::with_seed(71, rsample::vfold_cv(data, v = 5))
  scoring <- yardstick::metric_set(yardstick::rmse, yardstick::mae)
  run <- function(model, outer, metrics = scoring) withr::with_seed(72, tune::fit_resamples(
    model, outer, metrics = metrics, control = tune::control_resamples(save_pred = TRUE,
      extract = workflows::extract_fit_engine, allow_par = FALSE)))
  for (method in c("inductive", "cross")) {
    prepared <- withr::with_seed(73, lapply(outer$splits, function(s) {
      local <- rsample::analysis(s)
      if (method == "cross") rsample::vfold_cv(local, v = 5) else rsample::initial_split(local, .8)
    }))
    for (splits in list(NULL, local_calibration_splits(prepared, names(data)))) {
      before <- counts$fit
      model <- calibrated_workflow(base, "truth", method = method, splits = splits, seed = 73)
      result <- run(model, outer)
      checked <- example$check_evaluation(result, outer, classification = "regression")
      states <- tune::collect_extracts(result)$.extracts
      stopifnot(counts$fit == before + if (method == "cross") 25L else 5L, length(states) == 5L)
      for (i in seq_along(states)) {
        local <- rsample::analysis(outer$splits[[i]])
        assessment <- rsample::assessment(outer$splits[[i]])
        state <- states[[i]]
        for (id in names(state$folds)) {
          f <- state$folds[[id]]
          engine <- workflows::extract_fit_engine(state$models[[id]])
          recipe <- workflows::extract_recipe(state$models[[id]])
          stopifnot(identical(engine$x$id, local$id[f$analysis]),
                    identical(as.vector(engine$y), local$truth[f$analysis]),
                    !length(intersect(engine$x$id, assessment$id)),
                    !length(intersect(local$id[f$assessment], assessment$id)))
          x <- local$x[f$analysis]; x[is.na(x)] <- median(x, na.rm = TRUE)
          close(recipe$steps[[1L]]$medians[["x"]], median(x))
          close(recipe$steps[[2L]]$means[["x"]], mean(x))
          close(recipe$steps[[2L]]$sds[["x"]], sd(x))
        }
        # Reconstruct the manual core from actual fold models and local splits.
        calibration <- truth <- fold <- list()
        p <- checked$predictions[checked$predictions$id == outer$id[i], ]
        new_p <- matrix(NA_real_, nrow(p), length(state$folds), dimnames = list(NULL, names(state$folds)))
        for (id in names(state$folds)) {
          held <- state$folds[[id]]$assessment
          calibration[[id]] <- predict(state$models[[id]], local[held, ])$.pred
          truth[[id]] <- local$truth[held]
          fold[[id]] <- rep(id, length(held))
          new_p[, id] <- predict(state$models[[id]], data[p$.row, ])$.pred
        }
        expected <- if (method == "cross") {
          predict(cross_venn_abers_regression(unlist(truth), unlist(calibration), unlist(fold)), new_p)$prediction
        } else predict(venn_abers_regression(truth[[1L]], calibration[[1L]]), new_p[, 1L])$prediction
        close(p$.pred, expected)
      }
    }
    # Change only this outer holdout's truth: fitted states/predictions are fixed,
    # while both scores change. No calibration or preprocessing sees that truth.
    split <- rsample::make_splits(list(analysis = 1:80, assessment = 81:100), data)
    heldout <- rsample::manual_rset(list(split), "Test")
    altered <- data; altered$truth[81:100] <- altered$truth[81:100] + 1000
    altered_outer <- rsample::manual_rset(list(rsample::make_splits(
      list(analysis = 1:80, assessment = 81:100), altered)), "Test")
    model <- calibrated_workflow(base, "truth", method = method, seed = 73)
    one <- run(model, heldout); two <- run(model, altered_outer)
    example$check_evaluation(one, heldout, classification = "regression")
    example$check_evaluation(two, altered_outer, classification = "regression")
    stopifnot(identical(tune::collect_predictions(one)$.pred, tune::collect_predictions(two)$.pred),
              all(abs(tune::collect_metrics(one)$mean - tune::collect_metrics(two)$mean) > 1))
    a <- tune::collect_extracts(one)$.extracts[[1L]]
    b <- tune::collect_extracts(two)$.extracts[[1L]]
    stopifnot(identical(a$folds, b$folds), identical(a$parameters, b$parameters),
              identical(a$calibrator, b$calibrator))
    for (id in names(a$models)) {
      stopifnot(identical(workflows::extract_fit_engine(a$models[[id]]),
                          workflows::extract_fit_engine(b$models[[id]])),
                identical(workflows::extract_recipe(a$models[[id]])$steps,
                          workflows::extract_recipe(b$models[[id]])$steps))
    }
    # Partial outer failures and missing metrics are incomplete evaluation.
    partly_bad <- rsample::manual_rset(list(split, rsample::make_splits(
      list(analysis = 1:10, assessment = 11:100), data)), c("good", "bad"))
    failed <- suppressWarnings(run(model, partly_bad))
    stopifnot(nrow(tune::collect_notes(failed)) > 0L)
    error(example$check_evaluation(failed, partly_bad, classification = "regression"))
    incomplete <- run(model, heldout, yardstick::metric_set(yardstick::rmse))
    error(example$check_evaluation(incomplete, heldout, classification = "regression"))
  }
  close(yardstick::rmse_vec(c(0, 3), c(0, 1)), sqrt(2))
  close(yardstick::mae_vec(c(0, 3), c(0, 1)), 1)
}

  testthat::expect_true(TRUE)
})

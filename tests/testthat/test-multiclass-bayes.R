testthat::test_that("shared multiclass Bayes selects complete calibrated matrices", {
  skip_if_no_integration("tune", "dials", "nnet")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  params <- dials::parameters(dials::new_quant_param(type = "integer", range = c(1L, 4L),
    inclusive = c(TRUE, TRUE), label = c(degree = "Degree")))
  search <- list(method = "bayes", initial = 2L, iter = 1L, param_info = params,
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
    control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 2L))
  fitted <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
    tuning = list(strategy = "joint", model = search), seed = 7), data)
  selection <- venn_abers_tuning(fitted)$models$Inductive1
  testthat::expect_true(selection$chosen %in% selection$candidates$id)
  testthat::expect_true(length(selection$steps) >= 1L)
  testthat::expect_identical(selection$budget$candidates, 3L)
  testthat::expect_identical(selection$cost$candidate_evaluations,
    2L * (2L + length(selection$steps) - 1L))
  testthat::expect_identical(selection$attempts, 2L + length(selection$steps) - 1L)
  testthat::expect_true(all(selection$steps$initial$id %in% selection$candidates$id))
  prediction <- predict(fitted, data, type = "prob")
  testthat::expect_identical(names(prediction), paste0(".pred_", levels(data$truth)))
  testthat::expect_lte(max(abs(rowSums(prediction) - 1)), 1e-8)
})

testthat::test_that("two and four classes retain shared Bayes matrix paths", {
  skip_if_no_integration("tune", "dials", "nnet")
  fourth <- iris[iris$Species == "setosa", ]
  fourth$Sepal.Length <- fourth$Sepal.Length + 1.2
  fourth$Species <- "fourth"
  for (classes in c(2L, 4L)) {
    data <- if (classes == 2L) iris[iris$Species != "virginica", ] else rbind(iris, fourth)
    data$truth <- factor(as.character(data$Species), levels =
      if (classes == 2L) c("setosa", "versicolor") else c(levels(iris$Species), "fourth"))
    data$Species <- NULL
    base <- workflows::workflow() |>
      workflows::add_formula(truth ~ Sepal.Length) |>
      workflows::add_model(if (classes == 2L) parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(engine = "nnet"))
    classification <- if (classes == 2L) "ovo" else "ovr"
    strategy <- if (classes == 2L) "joint" else "two_stage"
    search <- list(method = "bayes", initial = 2L, iter = 1L,
      param_info = dials::parameters(dials::new_qual_param(type = "character",
        values = c("log", "brier"), label = c(va_loss = "Loss"))),
      resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
      control = tune::control_bayes(seed = 16, save_pred = TRUE, no_improve = 2L))
    fitted <- generics::fit(calibrated_workflow(base, "truth",
      classification = classification, loss = tune::tune("va_loss"),
      tuning = list(strategy = strategy, selection_prop = .4, calibration = search),
      seed = 17), data)
    state <- workflows::extract_fit_engine(fitted)
    held <- state$calibration$Inductive1
    core <- venn_abers(held$truth, held$probability, method = classification,
      precision = state$parameters$requested$precision)
    raw <- as.matrix(predict(state$models$Inductive1, data, type = "prob"))
    colnames(raw) <- levels(data$truth)
    expected <- predict(core, raw, loss = state$loss)
    actual <- as.matrix(predict(fitted, data, type = "prob"))
    testthat::expect_lte(max(abs(actual - expected)), 1e-10)
    testthat::expect_identical(colnames(actual), paste0(".pred_", levels(data$truth)))
    testthat::expect_identical(venn_abers_tuning(fitted)$classification, classification)
  }
})

testthat::test_that("shared OvR Bayes calibration keeps matrix scores and saved predictions", {
  skip_if_no_integration("tune", "dials", "yardstick", "nnet")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ Sepal.Length + Petal.Length) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  params <- dials::parameters(
    dials::new_qual_param(type = "character", values = c("log", "brier"),
      label = c(va_loss = "Loss")),
    dials::new_qual_param(type = "character", values = c("none", "0", "2", "-2"),
      label = c(rounding = "Rounding")))
  search <- list(method = "bayes", initial = 4L, iter = 1L, param_info = params,
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
    primary_metric = yardstick::brier_class,
    control = tune::control_bayes(seed = 8, save_pred = TRUE, no_improve = 2L))
  fitted <- generics::fit(calibrated_workflow(base, "truth", classification = "ovr",
    loss = tune::tune("va_loss"), precision = tune::tune("rounding"),
    tuning = list(strategy = "joint", selection_prop = .4, calibration = search), seed = 9), data)
  selection <- venn_abers_tuning(fitted)$calibration
  testthat::expect_true(selection$chosen %in% selection$candidates$id)
  testthat::expect_length(intersect(selection$selection_rows, selection$final_rows), 0L)
  testthat::expect_equal(sort(c(selection$selection_rows, selection$final_rows)), seq_len(nrow(data)))
  testthat::expect_identical(selection$primary_metric, yardstick::brier_class)
  testthat::expect_true(all(selection$candidates$status %in% c("eligible", "excluded")))
  native <- tune::collect_predictions(selection$native)
  for (id in selection$candidates$id[selection$candidates$status == "eligible"]) {
    for (fold in names(selection$resamples)) {
      rows <- selection$resamples[[fold]]$assessment
      candidate <- native[native$id == fold & native$va_loss ==
        selection$candidates$va_loss[match(id, selection$candidates$id)] &
        native$rounding == selection$candidates$rounding[match(id, selection$candidates$id)], ]
      candidate <- candidate[match(rows, candidate$.row), ]
      matrix <- as.matrix(candidate[paste0(".pred_", levels(data$truth))])
      colnames(matrix) <- levels(data$truth)
      testthat::expect_equal(selection$scores[[id]][[fold]],
        brier_score(candidate$truth, matrix) / 2, tolerance = 1e-10)
    }
  }
  saved <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  log <- tempfile(fileext = ".log")
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  expected <- predict(fitted, data[1:3, ], type = "prob")
  saveRDS(list(fit = fitted, new = data[1:3, ], expected = expected), saved)
  testthat::expect_identical(predict(readRDS(saved)$fit, data[1:3, ], type = "prob"), expected)
  writeLines(c(installed_vennabers_lines(),
    "library(workflows)",
    paste0("saved <- readRDS(", deparse(saved), ")"),
    "stopifnot(identical(predict(saved$fit, saved$new, type = 'prob'), saved$expected),",
    "  venn_abers_tuning(saved$fit)$calibration$chosen %in%",
    "  venn_abers_tuning(saved$fit)$calibration$candidates$id)"), script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_identical(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

testthat::test_that("two-stage Cross Bayes chooses its full model independently", {
  skip_if_no_integration("tune", "dials", "nnet")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  model <- list(method = "bayes", initial = 2L, iter = 1L,
    param_info = dials::parameters(dials::new_quant_param(type = "integer",
      range = c(1L, 4L), inclusive = c(TRUE, TRUE), label = c(degree = "Degree"))),
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
    control = tune::control_bayes(seed = 10, save_pred = TRUE, no_improve = 2L))
  calibration <- list(method = "bayes", initial = 2L, iter = 1L,
    param_info = dials::parameters(dials::new_qual_param(type = "character",
      values = c("log", "brier"), label = c(va_loss = "Loss"))),
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
    control = tune::control_bayes(seed = 11, save_pred = TRUE, no_improve = 2L))
  fitted <- generics::fit(calibrated_workflow(base, "truth", classification = "ovr",
    method = "cross", v = 2L, prediction_model = "full", loss = tune::tune("va_loss"),
    tuning = list(strategy = "two_stage", selection_prop = .4,
      model = model, calibration = calibration), seed = 12), data)
  selection <- venn_abers_tuning(fitted)
  state <- workflows::extract_fit_engine(fitted)
  testthat::expect_setequal(names(selection$models), c(names(state$folds), "Full1"))
  testthat::expect_identical(unname(state$sources), rep("Full1", length(state$folds)))
  testthat::expect_true(all(vapply(selection$models,
    function(x) x$chosen %in% x$candidates$id, logical(1))))
  testthat::expect_true(selection$calibration$chosen %in% selection$calibration$candidates$id)
  truth <- do.call(c, lapply(state$calibration, `[[`, "truth"))
  held <- do.call(rbind, lapply(state$calibration, `[[`, "probability"))
  fold <- rep(names(state$calibration), lengths(lapply(state$calibration, `[[`, "truth")))
  core <- cross_venn_abers(truth, held, fold, method = "ovr",
    precision = state$parameters$requested$precision)
  full <- as.matrix(predict(state$models$Full1, data, type = "prob"))
  colnames(full) <- levels(data$truth)
  expected <- predict(core, stats::setNames(rep(list(full), length(state$folds)),
    names(state$folds)), loss = state$loss)
  actual <- as.matrix(predict(fitted, data, type = "prob"))
  testthat::expect_lte(max(abs(actual - expected)), 1e-10)
})

testthat::test_that("shared Bayes rejects a search with no complete initial matrices", {
  skip_if_no_integration("tune", "dials")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  search <- list(method = "bayes", initial = 2L, iter = 1L,
    param_info = dials::parameters(dials::new_quant_param(type = "integer",
      range = c(1L, 4L), inclusive = c(TRUE, TRUE), label = c(degree = "Degree"))),
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
    control = tune::control_bayes(seed = 4, save_pred = TRUE))
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", tuning = list(strategy = "two_stage", model = search)), data),
    "Bayes needs at least 2 complete initial model candidates; found 0")
})

testthat::test_that("shared Bayes excludes incomplete fold and row candidates", {
  skip_if_no_integration("tune", "dials")
  data <- data.frame(x = seq_len(90), truth = factor(rep(c("a", "b", "c"), 30)))
  multi_fit <- function(x, y, variant) {
    if (variant == 1L && nrow(x) < 30L) stop("planned fold failure")
    list(classes = levels(y), variant = variant)
  }
  multi_prob <- function(object, new_data) {
    raw <- vapply(seq_along(object$classes), function(i)
      exp(i * new_data$x / 100 + i * object$variant / 20), numeric(nrow(new_data)))
    out <- as.data.frame(raw / rowSums(raw))
    names(out) <- paste0(".pred_", object$classes)
    if (object$variant == 2L && nrow(new_data) < 30L) out <- out[-1L, ]
    if (object$variant == 3L && nrow(new_data) < 30L) out <- rbind(out, out[1L, ])
    out
  }
  register_reference_engine("multinom_reg", "classification", "multi_fit",
    "multi_prob", "prob", engine = "bayes_reference", new_model = FALSE)
  spec <- parsnip::multinom_reg() |>
    parsnip::set_engine("bayes_reference", variant = tune::tune("variant"))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(spec)
  plans <- function(d) {
    held <- list(seq_len(floor(nrow(d) / 4)), seq.int(floor(nrow(d) / 4) + 1L, nrow(d)))
    rsample::manual_rset(lapply(held, function(rows) rsample::make_splits(list(
      analysis = setdiff(seq_len(nrow(d)), rows), assessment = rows), d)), c("small", "large"))
  }
  search <- list(method = "bayes", initial = 6L, iter = 1L,
    param_info = dials::parameters(dials::new_quant_param(type = "integer",
      range = c(0L, 5L), inclusive = c(TRUE, TRUE), label = c(variant = "Variant"))),
    resamples = plans, control = tune::control_bayes(seed = 20, save_pred = TRUE))
  fitted <- generics::fit(calibrated_workflow(base, "truth", classification = "ovr",
    tuning = list(strategy = "two_stage", model = search), seed = 21), data)
  selection <- venn_abers_tuning(fitted)$models$Inductive1
  testthat::expect_true(selection$chosen %in% selection$candidates$id)
  testthat::expect_gte(length(selection$steps), 2L)
  testthat::expect_setequal(selection$candidates$variant[selection$candidates$status == "excluded"], 1:3)
  testthat::expect_match(paste(tune::collect_notes(selection$native)$note, collapse = " "),
    "planned fold failure")
  testthat::expect_true(all(selection$steps$initial$id %in% selection$candidates$id))
  testthat::expect_identical(selection$cost$candidate_evaluations,
    selection$attempts * length(selection$resamples))
  testthat::expect_true(selection$stop_reason %in%
    c("budget", "no_improve", "time_limit", "repeated_candidate"))
})

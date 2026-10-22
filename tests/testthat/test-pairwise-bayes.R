testthat::test_that("pairwise joint Bayes scores complete coupled predictions", {
  skip_if_no_integration("tune", "dials")
  data <- iris[iris$Species != "virginica", c("Sepal.Length", "Species")]
  names(data)[2L] <- "truth"
  data$truth <- droplevels(data$truth)
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  search <- list(method = "bayes", initial = 2L, iter = 1L,
    param_info = dials::parameters(dials::new_quant_param(type = "integer",
      range = c(1L, 4L), inclusive = c(TRUE, TRUE), label = c(degree = "Degree"))),
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
    control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 2L))
  fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", tuning = list(strategy = "joint", model = search,
      joint_resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)),
    seed = 7), data))
  selection <- venn_abers_tuning(fitted)
  testthat::expect_identical(selection$joint$pair_ids, "pair_1")
  testthat::expect_identical(selection$joint$candidates$status, "eligible")
  testthat::expect_true(selection$models$pair_1$Inductive1$chosen %in%
    selection$models$pair_1$Inductive1$candidates$id)
  prediction <- predict(fitted, data, type = "prob")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(fit = fitted, new = data, expected = prediction), path)
  testthat::expect_identical(predict(readRDS(path)$fit, data, type = "prob"), prediction)
  testthat::expect_lte(max(abs(rowSums(prediction) - 1)), 1e-8)
  script <- tempfile(fileext = ".R")
  log <- tempfile(fileext = ".log")
  on.exit(unlink(c(script, log)), add = TRUE)
  writeLines(c(installed_vennabers_lines(),
    "library(workflows)",
    paste0("saved <- readRDS(", deparse(path), ")"),
    "stopifnot(venn_abers_tuning(saved$fit)$joint$candidates$status == 'eligible',",
    "  identical(predict(saved$fit, saved$new, type = 'prob'), saved$expected))"), script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_identical(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

testthat::test_that("pairwise Bayes keeps local searches separate from global calibration", {
  skip_if_no_integration("tune", "dials")
  data <- iris[c("Sepal.Length", "Species")]
  names(data)[2L] <- "truth"
  rownames(data) <- NULL
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  model <- list(method = "bayes", initial = 2L, iter = 1L,
    param_info = dials::parameters(dials::new_quant_param(type = "integer",
      range = c(1L, 2L), inclusive = c(TRUE, TRUE), label = c(degree = "Degree"))),
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
    control = tune::control_bayes(seed = 9, save_pred = TRUE, no_improve = 2L))
  calibration <- list(method = "bayes", initial = 2L, iter = 1L,
    param_info = dials::parameters(dials::new_qual_param(type = "character",
      values = c("log", "brier"), label = c(va_loss = "Loss"))),
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
    control = tune::control_bayes(seed = 10, save_pred = TRUE, no_improve = 2L))
  for (strategy in c("joint", "two_stage")) {
    fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
      classification = "ovo", model_training = "pairwise", prediction_model = "full",
      loss = tune::tune("va_loss"), tuning = list(strategy = strategy, selection_prop = .6,
        model = model, calibration = calibration), seed = 11), data))
    selection <- venn_abers_tuning(fitted)
    testthat::expect_identical(names(selection$models), paste0("pair_", 1:3))
    testthat::expect_true(selection$calibration$chosen %in% selection$calibration$candidates$id)
    testthat::expect_true(all(selection$calibration$candidates$status == "eligible"))
    testthat::expect_length(intersect(selection$calibration$selection_rows,
      selection$calibration$final_rows), 0L)
    for (id in names(selection$models)) {
      testthat::expect_true("Full1" %in% names(selection$models[[id]]))
      for (path in names(selection$models[[id]])) {
        decision <- selection$models[[id]][[path]]
        testthat::expect_identical(decision$pair_id, id)
        testthat::expect_identical(decision$fold_id, path)
        testthat::expect_true(decision$chosen %in% decision$candidates$id)
        testthat::expect_identical(decision$budget$candidates, 3L)
        testthat::expect_true(length(decision$steps) >= 1L)
        testthat::expect_true(decision$stop_reason %in%
          c("budget", "no_improve", "time_limit", "repeated_candidate"))
      }
    }
    testthat::expect_true(all(vapply(selection$sources,
      function(sources) all(sources == "Full1"), logical(1))))
    prediction <- predict(fitted, data, type = "prob")
    testthat::expect_identical(names(prediction), paste0(".pred_", levels(data$truth)))
    testthat::expect_lte(max(abs(rowSums(prediction) - 1)), 1e-8)
  }
  grid_model <- list(method = "grid", grid = data.frame(degree = 1:2),
    resamples = model$resamples)
  mixed <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", loss = tune::tune("va_loss"),
    tuning = list(strategy = "joint", selection_prop = .6,
      model = grid_model, calibration = calibration), seed = 11), data))
  testthat::expect_true(venn_abers_tuning(mixed)$calibration$chosen %in%
    venn_abers_tuning(mixed)$calibration$candidates$id)
  testthat::expect_identical(names(venn_abers_tuning(mixed)$models), paste0("pair_", 1:3))
  grid_calibration <- list(method = "grid", grid = data.frame(va_loss = c("log", "brier")),
    resamples = calibration$resamples)
  mixed <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", loss = tune::tune("va_loss"),
    tuning = list(strategy = "joint", selection_prop = .6,
      model = model, calibration = grid_calibration), seed = 11), data))
  testthat::expect_true(venn_abers_tuning(mixed)$calibration$chosen %in%
    venn_abers_tuning(mixed)$calibration$candidates$id)
})

testthat::test_that("four-class pairwise Bayes uses unequal local plans and a complete global score", {
  skip_if_no_integration("tune", "dials", "yardstick")
  fourth <- head(iris[iris$Species == "setosa", ], 25L)
  fourth$Sepal.Length <- fourth$Sepal.Length + 1
  fourth$Species <- "fourth"
  data <- rbind(iris, fourth)
  data$truth <- factor(as.character(data$Species), levels = c(levels(iris$Species), "fourth"))
  data <- data[c("Sepal.Length", "truth")]
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ Sepal.Length) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  plans <- function(d) rsample::vfold_cv(d,
    v = if (identical(levels(d$truth), c("setosa", "versicolor"))) 2L else 3L,
    strata = truth, pool = 0)
  calibration <- list(method = "bayes", initial = 2L, iter = 1L,
    param_info = dials::parameters(dials::new_qual_param(type = "character",
      values = c("log", "brier"), label = c(va_loss = "Loss"))),
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
    primary_metric = yardstick::brier_class,
    control = tune::control_bayes(seed = 12, save_pred = TRUE, no_improve = 2L))
  fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", method = "cross",
    prediction_model = "full", splits = plans, loss = tune::tune("va_loss"),
    tuning = list(strategy = "joint", selection_prop = .6, calibration = calibration),
    seed = 13), data))
  selection <- venn_abers_tuning(fitted)
  testthat::expect_identical(names(selection$models), paste0("pair_", 1:6))
  testthat::expect_gt(length(unique(vapply(workflows::extract_fit_engine(fitted)$pair_models,
    function(pair) length(pair$rows), 0L))), 1L)
  testthat::expect_identical(lengths(selection$fit_roles), c(pair_1 = 2L, pair_2 = 3L,
    pair_3 = 3L, pair_4 = 3L, pair_5 = 3L, pair_6 = 3L))
  testthat::expect_true(selection$calibration$chosen %in% selection$calibration$candidates$id)
  testthat::expect_true(all(selection$calibration$candidates$status == "eligible"))
  testthat::expect_identical(selection$calibration$primary_metric, yardstick::brier_class)
  rows <- selection$calibration$resamples[[1L]]$assessment
  native <- tune::collect_predictions(selection$calibration$native)
  candidate <- native[native$id == names(selection$calibration$resamples)[1L] &
    native$va_loss == selection$calibration$candidates$va_loss[1L], ]
  candidate <- candidate[match(rows, candidate$.row), ]
  probability <- as.matrix(candidate[paste0(".pred_", levels(data$truth))])
  colnames(probability) <- levels(data$truth)
  testthat::expect_equal(selection$calibration$scores[[1L]][[1L]],
    brier_score(candidate$truth, probability) / 2, tolerance = 1e-10)
  testthat::expect_lte(max(abs(rowSums(predict(fitted, data, type = "prob")) - 1)), 1e-8)
})

testthat::test_that("four-class model Bayes keeps pair budgets and resamples local", {
  skip_if_no_integration("tune", "dials", "rpart")
  fourth <- head(iris[iris$Species == "setosa", ], 20L)
  fourth$Sepal.Length <- fourth$Sepal.Length + 1
  fourth$Species <- "fourth"
  data <- rbind(head(iris[iris$Species == "setosa", ], 30L),
    head(iris[iris$Species == "versicolor", ], 30L),
    head(iris[iris$Species == "virginica", ], 30L), fourth)
  data$truth <- factor(as.character(data$Species), levels = c(levels(iris$Species), "fourth"))
  data <- data[c("Sepal.Length", "truth")]
  rownames(data) <- NULL
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ Sepal.Length) |>
    workflows::add_model(parsnip::decision_tree(tree_depth = tune::tune("depth"),
      mode = "classification", engine = "rpart"))
  model <- list(method = "bayes", initial = 2L, iter = 1L,
    param_info = dials::parameters(dials::new_quant_param(type = "integer",
      range = c(1L, 2L), inclusive = c(TRUE, TRUE), label = c(depth = "Depth"))),
    resamples = function(d) rsample::vfold_cv(d,
      v = if (identical(levels(d$truth), c("setosa", "versicolor"))) 2L else 3L,
      strata = truth, pool = 0),
    control = tune::control_bayes(seed = 14, save_pred = TRUE, no_improve = 2L))
  fitted <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", tuning = list(strategy = "joint", model = model,
      joint_resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)),
    seed = 15), data)
  selection <- venn_abers_tuning(fitted)
  testthat::expect_identical(selection$joint$pair_ids, paste0("pair_", 1:6))
  testthat::expect_identical(selection$joint$candidates$status, "eligible")
  testthat::expect_identical(vapply(selection$models, function(pair)
    length(pair$Inductive1$resamples), 0L),
    c(pair_1 = 2L, pair_2 = 3L, pair_3 = 3L, pair_4 = 3L, pair_5 = 3L, pair_6 = 3L))
  testthat::expect_gt(length(unique(vapply(workflows::extract_fit_engine(fitted)$pair_models,
    function(pair) length(pair$rows), 0L))), 1L)
  testthat::expect_true(all(vapply(selection$models, function(pair)
    pair$Inductive1$chosen %in% pair$Inductive1$candidates$id, logical(1))))
  testthat::expect_true(any(vapply(selection$models, function(pair)
    identical(pair$Inductive1$stop_reason, "repeated_candidate"), logical(1))))
  testthat::expect_lte(max(abs(rowSums(predict(fitted, data, type = "prob")) - 1)), 1e-8)
})

testthat::test_that("a missing pair excludes every global Bayes candidate", {
  skip_if_no_integration("tune", "dials")
  data <- iris[c("Sepal.Length", "Species")]
  names(data)[2L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ Sepal.Length) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  split <- rsample::make_splits(list(analysis = c(1:25, 51:75, 101:150),
    assessment = c(26:50, 76:100)), data)
  calibration <- list(method = "bayes", initial = 2L, iter = 1L,
    param_info = dials::parameters(dials::new_qual_param(type = "character",
      values = c("log", "brier"), label = c(va_loss = "Loss"))),
    resamples = function(d) rsample::vfold_cv(d, v = 2),
    control = tune::control_bayes(seed = 3, save_pred = TRUE))
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", loss = tune::tune("va_loss"),
    tuning = list(strategy = "joint", selection_split = split,
      calibration = calibration), seed = 4), data),
    "Bayes needs at least 2 complete initial calibration candidates; found 0")
})

testthat::test_that("one partial pair candidate cannot win a complete coupled search", {
  skip_if_no_integration("tune", "dials")
  pair_bayes_fit <- function(x, y, variant) list(classes = levels(y), variant = variant)
  pair_bayes_prob <- function(object, new_data) {
    p <- stats::plogis(new_data$Sepal.Length / 4 + object$variant / 10)
    if (object$variant == 3L && identical(object$classes, c("setosa", "versicolor")))
      p <- p[-1L]
    out <- data.frame(1 - p, p)
    names(out) <- paste0(".pred_", object$classes)
    out
  }
  register_reference_engine("logistic_reg", "classification", "pair_bayes_fit",
    "pair_bayes_prob", "prob", engine = "pair_bayes_reference", new_model = FALSE)
  data <- iris[c("Sepal.Length", "Species")]
  names(data)[2L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ Sepal.Length) |>
    workflows::add_model(parsnip::logistic_reg() |>
      parsnip::set_engine("pair_bayes_reference", variant = tune::tune("variant")))
  model <- list(method = "bayes", initial = 3L, iter = 1L,
    param_info = dials::parameters(dials::new_quant_param(type = "integer",
      range = c(1L, 3L), inclusive = c(TRUE, TRUE), label = c(variant = "Variant"))),
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
    control = tune::control_bayes(seed = 22, save_pred = TRUE, no_improve = 2L))
  fitted <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", tuning = list(strategy = "joint", model = model,
      joint_resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)),
    seed = 23), data)
  selection <- venn_abers_tuning(fitted)
  first <- selection$models$pair_1$Inductive1
  testthat::expect_identical(first$candidates$status[first$candidates$variant == 3L], "excluded")
  testthat::expect_match(first$candidates$reason[first$candidates$variant == 3L],
    "Native predictions omit")
  testthat::expect_true(first$selected$variant %in% 1:2)
  testthat::expect_identical(selection$joint$candidates$status, "eligible")
  testthat::expect_true(all(vapply(selection$joint$decisions[[1L]], function(fold)
    identical(names(fold), paste0("pair_", 1:3)), logical(1))))
  grid <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", tuning = list(strategy = "joint",
      model = list(method = "grid", grid = data.frame(variant = 1:3)),
      joint_resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)),
    seed = 23), data)
  global <- venn_abers_tuning(grid)$joint
  partial <- global$grid$`..pair_1_source` == 3L
  testthat::expect_true(all(global$candidates$status[partial] == "excluded"))
  testthat::expect_true(any(global$candidates$status == "eligible"))
  testthat::expect_identical(global$candidates$status[
    match(global$chosen, global$candidates$id)], "eligible")
  testthat::expect_true(all(is.finite(global$candidates$estimate[
    global$candidates$status == "eligible"])))
})

testthat::test_that("pure Grid joint still rejects local model resamples", {
  skip_if_no_integration("tune")
  data <- iris[c("Sepal.Length", "Species")]
  names(data)[2L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  folds <- function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)
  config <- list(strategy = "joint",
    model = list(method = "grid", grid = data.frame(degree = 1:2), resamples = folds),
    joint_resamples = folds)
  testthat::expect_error(calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", tuning = config), "omit `tuning\\$model\\$resamples`")
  attr(config, "vennabers_local_pair_joint") <- TRUE
  testthat::expect_error(calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", tuning = config), "omit `tuning\\$model\\$resamples`")
})

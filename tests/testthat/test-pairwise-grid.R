testthat::test_that("pairwise calibration grid scores the coupled complete matrix", {
  testthat::skip_if_not_installed("tune")
  data <- rbind(iris, iris[1:10, ])
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ Sepal.Length) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  plans <- function(d) rsample::vfold_cv(d,
    v = if (identical(levels(d$truth), c("setosa", "versicolor"))) 2L else 3L,
    strata = truth, pool = 0)
  broken_report <- yardstick::new_prob_metric(function(...) stop("report unavailable"),
                                               direction = "minimize")
  config <- list(strategy = "joint", selection_prop = .5,
    calibration = list(method = "grid",
      grid = data.frame(va_loss = c("log", "brier", "log", "brier"),
                        rounding = c("none", "0", "2", "-2")),
      report_metrics = list(broken_report),
      resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
  fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", method = "cross", prediction_model = "full",
    splits = plans, loss = tune::tune("va_loss"), precision = tune::tune("rounding"),
    tuning = config, seed = 17), data))
  state <- workflows::extract_fit_engine(fitted)
  selection <- venn_abers_tuning(fitted)
  testthat::expect_identical(selection$classification, "ovo")
  testthat::expect_identical(selection$model_training, "pairwise")
  testthat::expect_identical(selection$calibration$budget$pair_count, 3)
  testthat::expect_true(is.na(selection$calibration$budget$paths))
  testthat::expect_identical(names(selection$models), paste0("pair_", 1:3))
  testthat::expect_identical(lengths(selection$fit_roles), c(pair_1 = 2L, pair_2 = 3L, pair_3 = 3L))
  testthat::expect_true(length(unique(vapply(state$pair_models, function(x) length(x$rows), 0L))) > 1L)
  testthat::expect_true(all(vapply(selection$sources,
    function(source) all(source == "Full1"), TRUE)))
  testthat::expect_length(intersect(selection$calibration$selection_rows,
                                    selection$calibration$final_rows), 0L)
  testthat::expect_identical(selection$calibration$candidates$status, rep("eligible", 4L))
  testthat::expect_match(selection$calibration$report_metrics[[1L]][[1L]][[1L]],
                         "report unavailable")
  native <- tune::collect_predictions(selection$calibration$native)
  for (candidate_id in selection$calibration$candidates$id) {
    for (fold in names(selection$calibration$resamples)) {
      rows <- selection$calibration$resamples[[fold]]$assessment
      candidate <- native[native$candidate_id == candidate_id & native$id == fold, ]
      candidate <- candidate[match(rows, candidate$.row), ]
      probability <- as.matrix(candidate[paste0(".pred_", levels(data$truth))])
      colnames(probability) <- levels(data$truth)
      testthat::expect_equal(selection$calibration$scores[[candidate_id]][[fold]],
        log_loss(candidate$truth, probability), tolerance = 1e-10)
    }
  }
  predicted <- as.matrix(predict(fitted, data, type = "prob"))
  pairs <- utils::combn(levels(data$truth), 2L, simplify = FALSE)
  points <- sapply(seq_along(pairs), function(i) {
    component <- state$pair_models[[i]]
    core <- cross_venn_abers(
      unlist(lapply(component$calibration, `[[`, "truth"), use.names = FALSE),
      unlist(lapply(component$calibration, `[[`, "probability"), use.names = FALSE),
      rep(names(component$calibration), lengths(lapply(component$calibration, `[[`, "truth"))),
      precision = state$parameters$requested$precision)
    model <- component$models$Full1
    p <- predict(model, data, type = "prob")[[paste0(".pred_", pairs[[i]][2L])]]
    input <- matrix(rep(p, length(component$folds)), nrow(data), length(component$folds),
                    dimnames = list(NULL, names(component$folds)))
    predict(core, input, loss = state$loss)$probability
  })
  testthat::expect_lte(max(abs(predicted - vennabers:::.couple_ovo(points, pairs,
                                                         levels(data$truth)))), 1e-10)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fitted, path)
  restored <- readRDS(path)
  testthat::expect_equal(predict(restored, data, type = "prob"),
                         predict(fitted, data, type = "prob"))
})

testthat::test_that("pairwise model grids choose each fold and full path locally", {
  testthat::skip_if_not_installed("tune")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  config <- list(strategy = "two_stage", model = list(method = "grid",
    grid = data.frame(degree = c(1L, 2L)),
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
  fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", method = "cross", v = 2L,
    prediction_model = "full", tuning = config, seed = 2), data))
  selected <- venn_abers_tuning(fitted)$models
  testthat::expect_identical(names(selected), paste0("pair_", 1:3))
  for (pair_id in names(selected)) for (fold_id in names(selected[[pair_id]])) {
    decision <- selected[[pair_id]][[fold_id]]
    testthat::expect_identical(decision$pair_id, pair_id)
    testthat::expect_identical(decision$fold_id, fold_id)
    testthat::expect_true(decision$chosen %in% decision$candidates$id)
    testthat::expect_identical(decision$primary_metric, yardstick::mn_log_loss)
  }
  testthat::expect_true(all(vapply(selected, function(pair) "Full1" %in% names(pair), TRUE)))
  testthat::expect_lte(max(abs(rowSums(predict(fitted, data, type = "prob")) - 1)), 1e-8)
})

testthat::test_that("two and four class pairwise grids retain their global class contract", {
  testthat::skip_if_not_installed("tune")
  fourth <- iris[iris$Species == "setosa", ]
  fourth$Sepal.Length <- fourth$Sepal.Length + 1
  fourth$Species <- "fourth"
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ Sepal.Length) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (k in c(2L, 4L)) {
    data <- if (k == 2L) iris[iris$Species != "virginica", ] else rbind(iris, fourth)
    data$truth <- factor(as.character(data$Species), levels =
      if (k == 2L) c("setosa", "versicolor") else c(levels(iris$Species), "fourth"))
    data$Species <- NULL
    config <- list(strategy = "joint", selection_prop = .5,
      calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")),
        primary_metric = if (k == 4L) yardstick::brier_class else yardstick::mn_log_loss,
        resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
    fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
      classification = "ovo", model_training = "pairwise", loss = tune::tune("va_loss"),
      tuning = config, seed = 4), data))
    probability <- predict(fitted, data, type = "prob")
    testthat::expect_identical(names(probability), paste0(".pred_", levels(data$truth)))
    testthat::expect_lte(max(abs(rowSums(probability) - 1)), 1e-8)
    testthat::expect_length(venn_abers_tuning(fitted)$models, choose(k, 2L))
    if (k == 4L) {
      selection <- venn_abers_tuning(fitted)$calibration
      native <- tune::collect_predictions(selection$native)
      fold <- names(selection$resamples)[1L]
      candidate <- native[native$candidate_id == selection$candidates$id[1L] &
                            native$id == fold, ]
      probability <- as.matrix(candidate[paste0(".pred_", levels(data$truth))])
      colnames(probability) <- levels(data$truth)
      testthat::expect_equal(selection$scores[[1L]][[fold]],
        brier_score(candidate$truth, probability) / 2, tolerance = 1e-10)
    }
  }
})

testthat::test_that("joint and two-stage pairwise model and calibration grids fit", {
  testthat::skip_if_not_installed("tune")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (strategy in c("joint", "two_stage")) {
    config <- list(strategy = strategy, selection_prop = .6,
      model = list(method = "grid", grid = data.frame(degree = c(1L, 2L)),
        resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)),
      calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")),
        resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
    if (strategy == "joint") config$model$resamples <- NULL
    fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
      classification = "ovo", model_training = "pairwise",
      loss = tune::tune("va_loss"), tuning = config, seed = 19), data))
    selection <- venn_abers_tuning(fitted)
    testthat::expect_identical(selection$calibration$budget$paths, 3)
    testthat::expect_true(selection$calibration$chosen %in% selection$calibration$candidates$id)
    testthat::expect_identical(names(selection$models), paste0("pair_", 1:3))
    if (strategy == "joint") {
      testthat::expect_equal(nrow(selection$joint$candidates), 16L)
      testthat::expect_identical(selection$joint$chosen, selection$calibration$chosen)
      testthat::expect_true(all(vapply(selection$models,
        function(pair) identical(pair$Inductive1$chosen, selection$joint$chosen), TRUE)))
      native <- tune::collect_predictions(selection$joint$native)
      fold <- names(selection$joint$resamples)[1L]
      testthat::expect_identical(
        selection$joint$decisions[[selection$joint$chosen]][[fold]]$pair_1$Inductive1$chosen,
        selection$joint$chosen)
      candidate <- native[native$candidate_id == selection$joint$chosen & native$id == fold, ]
      probability <- as.matrix(candidate[paste0(".pred_", levels(data$truth))])
      colnames(probability) <- levels(data$truth)
      testthat::expect_equal(selection$joint$scores[[selection$joint$chosen]][[fold]],
        log_loss(candidate$truth, probability), tolerance = 1e-10)
    } else testthat::expect_true(all(vapply(selection$models,
      function(pair) pair$Inductive1$chosen %in% pair$Inductive1$candidates$id, TRUE)))
    testthat::expect_lte(max(abs(rowSums(predict(fitted, data, type = "prob")) - 1)), 1e-8)
  }
})

testthat::test_that("model-only joint selection uses a separate global Cartesian plan", {
  testthat::skip_if_not_installed("tune")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  config <- list(strategy = "joint", model = list(method = "grid",
    grid = data.frame(degree = c(1L, 2L))),
    joint_resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0))
  fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", tuning = config, seed = 19), data))
  selection <- venn_abers_tuning(fitted)
  testthat::expect_null(selection$calibration)
  testthat::expect_equal(nrow(selection$joint$candidates), 8L)
  testthat::expect_identical(selection$joint$pair_ids, paste0("pair_", 1:3))
  testthat::expect_true(selection$joint$chosen %in% selection$joint$candidates$id)
  testthat::expect_lte(max(abs(rowSums(predict(fitted, data, type = "prob")) - 1)), 1e-8)
  two <- data[data$truth != "virginica", ]
  two$truth <- droplevels(two$truth)
  full <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", prediction_model = "full",
    tuning = config, seed = 19), two))
  full_selection <- venn_abers_tuning(full)
  testthat::expect_equal(nrow(full_selection$joint$candidates), 4L)
  testthat::expect_identical(full_selection$models$pair_1$Full1$selected,
                             full_selection$joint$selected$pairs$pair_1$full)
  testthat::expect_equal(length(workflows::extract_fit_engine(full)$pair_models$pair_1$rows),
                         nrow(two))
  config$joint_resamples <- NULL
  testthat::expect_error(calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", tuning = config), "joint_resamples")
})

testthat::test_that("missing pair predictions invalidate the global candidate", {
  testthat::skip_if_not_installed("tune")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ Sepal.Length) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  split <- rsample::make_splits(list(analysis = 1:120, assessment = 121:150), data)
  config <- list(strategy = "joint", selection_split = split,
    calibration = list(method = "grid", grid = data.frame(va_loss = "log"),
      resamples = function(d) rsample::vfold_cv(d, v = 2)))
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise",
    loss = tune::tune("va_loss"), tuning = config, seed = 7), data),
    "No eligible calibration grid candidate")
})

testthat::test_that("explicit global S/F roles stay separate from local pair plans", {
  testthat::skip_if_not_installed("tune")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ Sepal.Length) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  final_rows <- c(1:25, 51:75, 101:125)
  selection_rows <- setdiff(seq_len(nrow(data)), final_rows)
  config <- list(strategy = "joint",
    selection_split = rsample::make_splits(list(analysis = final_rows,
                                                assessment = selection_rows), data),
    calibration = list(method = "grid", grid = data.frame(va_loss = "log"),
      resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
  local_plan <- function(d) rsample::initial_split(d, prop = .7, strata = truth, pool = 0)
  fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", splits = local_plan,
    loss = tune::tune("va_loss"), tuning = config, seed = 8), data))
  selection <- venn_abers_tuning(fitted)$calibration
  testthat::expect_identical(selection$final_rows, final_rows)
  testthat::expect_identical(selection$selection_rows, selection_rows)
  testthat::expect_true(all(vapply(workflows::extract_fit_engine(fitted)$pair_models,
    function(pair) length(pair$rows) == sum(data$truth[final_rows] %in% pair$class_levels), TRUE)))
})

testthat::test_that("pairwise grid diagnostics and schema survive RDS restore", {
  testthat::skip_if_not_installed("tune")
  original <- iris
  names(original)[5L] <- "truth"
  schema <- venn_abers_schema(original, "truth", classification = "ovo")
  data <- venn_abers_encode(original, schema)
  outcome <- schema$encoded_outcome
  base <- workflows::workflow() |>
    workflows::add_formula(stats::reformulate(setdiff(names(data), outcome), outcome)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  config <- list(strategy = "joint", selection_prop = .5,
    calibration = list(method = "grid", grid = data.frame(va_loss = "log"),
      resamples = function(d) rsample::vfold_cv(d, v = 2)))
  fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, outcome,
    schema = schema, model_training = "pairwise", loss = tune::tune("va_loss"),
    tuning = config, seed = 12), data))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fitted, path)
  restored <- readRDS(path)
  testthat::expect_identical(venn_abers_tuning(restored)$schema, schema)
  testthat::expect_identical(names(venn_abers_tuning(restored)$models), paste0("pair_", 1:3))
  testthat::expect_equal(predict(restored, data, type = "prob"),
                         predict(fitted, data, type = "prob"))
  decoded <- venn_abers_decode(predict(restored, data, type = "prob"), schema)
  testthat::expect_identical(colnames(decoded$probability), levels(original$truth))
})

testthat::test_that("joint pairwise grids generated from dials sizes match explicit grids", {
  skip_if_no_integration("dials")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  folds <- function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)
  explicit <- list(strategy = "joint", selection_prop = .6,
    model = list(method = "grid", grid = data.frame(degree = c(1L, 2L))),
    calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")),
      resamples = folds))
  generated <- explicit
  generated$model <- list(method = "grid", grid = 2L, param_info = dials::parameters(dials::new_quant_param(
    type = "integer", range = c(1L, 2L), inclusive = c(TRUE, TRUE), label = c(degree = "Degree"))))
  generated$calibration$grid <- 2L
  generated$calibration$param_info <- dials::parameters(dials::new_qual_param(
    type = "character", values = c("log", "brier"), label = c(va_loss = "Loss")))
  select <- function(config) venn_abers_tuning(suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise",
    loss = tune::tune("va_loss"), tuning = config, seed = 19), data)))$joint
  reference <- select(explicit)
  selection <- select(generated)
  # Two degrees for each of the three pairs, crossed with two losses.
  testthat::expect_equal(nrow(selection$candidates), 16L)
  testthat::expect_equal(sort(selection$candidates$estimate), sort(reference$candidates$estimate))
})

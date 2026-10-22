testthat::test_that("joint Cartesian grids fit each pair candidate once per resample", {
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
  fit <- function() suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", tuning = config, seed = 19), data))
  fitted <- fit()
  selection <- venn_abers_tuning(fitted)
  testthat::expect_equal(nrow(selection$joint$candidates), 8L)
  testthat::expect_equal(selection$joint$cost$observed_base_fits, 2L * 3L * 2L)
  testthat::expect_equal(selection$cost$observed_base_fits, 2L * 3L * 2L + 3L)
  testthat::expect_equal(selection$joint$budget$reused_base_fits, 36L)
  uncached <- testthat::with_mocked_bindings(fit(),
    .fit_pairwise_grid_domains = function(data, domains, task, settings, tuning_plan) {
      vennabers:::.fit_domains(domains, task, settings, tuning_plan)
    }, .package = "vennabers")
  reference <- venn_abers_tuning(uncached)$joint
  testthat::expect_equal(reference$cost$observed_base_fits, 48L)
  testthat::expect_equal(selection$joint$candidates, reference$candidates)
  testthat::expect_identical(selection$joint$chosen, reference$chosen)
  testthat::expect_equal(selection$joint$decisions, reference$decisions)
  testthat::expect_equal(tune::collect_predictions(selection$joint$native),
    tune::collect_predictions(reference$native))
  testthat::expect_equal(predict(fitted, data, type = "prob"), predict(uncached, data, type = "prob"))
  # Memoised grid models never reach the final fit.
  testthat::expect_false(any(vapply(workflows::extract_fit_engine(fitted)$pair_models,
    function(component) any(vapply(component$models, inherits, TRUE, "vennabers_memoised")), TRUE)))
})

testthat::test_that("Cross full grids reuse raw predictions across calibration choices", {
  testthat::skip_if_not_installed("tune")
  data <- droplevels(iris[iris$Species != "virginica", ])
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  config <- list(strategy = "joint", selection_prop = .5,
    model = list(method = "grid", grid = data.frame(degree = c(1L, 2L))),
    calibration = list(method = "grid",
      grid = data.frame(va_loss = c("log", "brier"), rounding = c("none", "0")),
      resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
  fit <- function() suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", method = "cross", v = 2,
    prediction_model = "full", loss = tune::tune("va_loss"), precision = tune::tune("rounding"),
    tuning = config, seed = 19), data))
  predictions <- 0L
  original_predict <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  registerS3method("predict", "workflow", function(object, ...) {
    predictions <<- predictions + 1L
    original_predict(object, ...)
  }, envir = asNamespace("stats"))
  fitted <- tryCatch(fit(), finally = registerS3method("predict", "workflow", original_predict,
    envir = asNamespace("stats")))
  selection <- venn_abers_tuning(fitted)$joint
  testthat::expect_equal(selection$cost$observed_base_fits, 12L)
  testthat::expect_equal(predictions, 19L)
  testthat::expect_equal(selection$budget$reused_base_fits, 36L)
  uncached <- testthat::with_mocked_bindings(fit(),
    .fit_pairwise_grid_domains = function(data, domains, task, settings, tuning_plan) {
      vennabers:::.fit_domains(domains, task, settings, tuning_plan)
    }, .package = "vennabers")
  reference <- venn_abers_tuning(uncached)$joint
  testthat::expect_equal(selection$candidates, reference$candidates)
  testthat::expect_identical(selection$chosen, reference$chosen)
  testthat::expect_equal(selection$decisions, reference$decisions)
  testthat::expect_equal(tune::collect_predictions(selection$native),
    tune::collect_predictions(reference$native))
  testthat::expect_equal(predict(fitted, data, type = "prob"), predict(uncached, data, type = "prob"))
})

testthat::test_that("failed pair fits exclude every dependent combination without retrying", {
  testthat::skip_if_not_installed("tune")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  config <- list(strategy = "joint", model = list(method = "grid",
    grid = data.frame(degree = c(1L, 0L))),
    joint_resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0))
  fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", tuning = config, seed = 19), data))
  selection <- venn_abers_tuning(fitted)$joint
  testthat::expect_identical(selection$candidates$status, c("eligible", rep("excluded", 7L)))
  testthat::expect_identical(selection$chosen, "joint_0001")
  testthat::expect_equal(selection$cost$observed_base_fits, 12L)
})

testthat::test_that("stochastic pair fits are independent of Cartesian candidate order", {
  testthat::skip_if_not_installed("tune")
  fit_base <- vennabers:::.fit_base_model
  base_probability <- vennabers:::.base_probability
  testthat::local_mocked_bindings(
    .fit_base_model = function(workflow, data, counter) {
      fitted <- fit_base(workflow, data, counter)
      # Simulate an engine whose parameter choice changes its RNG consumption.
      attr(fitted, "draw") <- tail(runif(ncol(workflows::extract_mold(fitted)$predictors)), 1L)
      fitted
    },
    .base_probability = function(model, ...) {
      .8 * base_probability(model, ...) + .2 * attr(model, "draw")
    }, .package = "vennabers")
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
      recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  scores <- lapply(list(1:2, 2:1), function(degrees) {
    config <- list(strategy = "joint", model = list(method = "grid", grid = data.frame(degree = degrees)),
      joint_resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0))
    fitted <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
      model_training = "pairwise", tuning = config, seed = 19), data))
    selection <- venn_abers_tuning(fitted)$joint
    names(selection$scores) <- apply(selection$grid[paste0("..pair_", 1:3, "_source")], 1L,
      function(indices) paste(degrees[indices], collapse = "/"))
    selection$scores[sort(names(selection$scores))]
  })
  testthat::expect_equal(scores[[1L]], scores[[2L]])
})

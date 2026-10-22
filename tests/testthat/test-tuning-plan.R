testthat::test_that("nested candidates fit without rebuilding public workflows", {
  testthat::skip_if_not_installed("tune")
  data <- data.frame(x = sin(seq_len(160) / 7),
    truth = factor(rep(c("yes", "no"), 80)))
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x, data) |>
      recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  resamples <- function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)
  outer <- calibrated_workflow(base, "truth", method = "cross", v = 2,
    prediction_model = "full", loss = tune::tune("va_loss"), seed = 41,
    tuning = list(strategy = "joint", selection_prop = .5,
      model = list(method = "grid", grid = data.frame(degree = 1L), resamples = resamples),
      calibration = list(method = "grid", grid = data.frame(va_loss = "log"), resamples = resamples)))
  public_fit <- fit_vennabers
  public_calls <- 0L
  testthat::local_mocked_bindings(
    calibrated_workflow = function(...) stop("Candidates must not rebuild public workflows."),
    fit_vennabers = function(...) {
      public_calls <<- public_calls + 1L
      if (public_calls > 1L) stop("Candidates must use the proxy entry point.")
      public_fit(...)
    }, .package = "vennabers")
  fitted <- generics::fit(outer, data)
  state <- workflows::extract_fit_engine(fitted)
  diagnostic <- venn_abers_tuning(fitted)
  testthat::expect_identical(public_calls, 1L)
  testthat::expect_identical(diagnostic$calibration$candidates$status, "eligible")
  # Two calibration assessment fits and the final fit each need three source/
  # full paths, each with two candidate fits of three paths plus one base fit.
  testthat::expect_identical(diagnostic$cost$observed_base_fits, 63L)
  testthat::expect_true(all(state$sources == "Full1"))
  testthat::expect_length(diagnostic$calibration$decisions[[1L]], 2L)
  testthat::expect_length(intersect(diagnostic$calibration$selection_rows,
    diagnostic$calibration$final_rows), 0L)
  testthat::expect_false(any(startsWith(names(state$parameters$requested), ".")))
  testthat::expect_false(any(startsWith(names(state$parameters$effective), ".")))
  testthat::expect_equal(rowSums(predict(fitted, data, type = "prob")), rep(1, nrow(data)))
})

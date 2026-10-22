testthat::test_that("numeric isotonic holdout preserves tie weights and interpolates the complete workflow", {
  skip_if_no_integration()
  d <- data.frame(x = c(1:6, 0, 0, 1, 2, 3, 4, 1:8),
                  truth = c(1:6, 10, 0, 1, 8, 6, 7, 1:8))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  va <- function(data) rsample::make_splits(list(analysis = 1:12, assessment = 13:20), data)
  post <- list(method = "isotonic_numeric", strategy = "holdout", splits = function(data)
    rsample::make_splits(list(analysis = 1:6, assessment = 7:12), data))
  fit <- generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d)
  state <- workflows::extract_fit_engine(fit)
  new <- data.frame(x = c(-1, 0, .5, 1.5, 2.5, 4, 5))
  expected <- c(11/3, 11/3, 11/3, 16/3, 7, 7, 7)
  actual <- predict(state$models[[1L]], new, type = "numeric")$.pred
  testthat::expect_true(all(abs(actual - expected) <= 1e-10 + 1e-10 * abs(expected)))
  shuffled <- d
  shuffled[7:12, ] <- d[c(10, 7, 12, 9, 8, 11), ]
  shuffled_fit <- generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), shuffled)
  shuffled_actual <- predict(workflows::extract_fit_engine(shuffled_fit)$models[[1L]], new)$.pred
  testthat::expect_true(all(abs(shuffled_actual - expected) <= 1e-10 + 1e-10 * abs(expected)))
  testthat::expect_identical(state$postprocessing[[1L]]$input, "numeric")
  testthat::expect_identical(state$postprocessing[[1L]]$P_size, 6L)
  testthat::expect_identical(state$postprocessing[[1L]]$observed_base_fits, 1L)
  testthat::expect_identical(state$postprocessing[[1L]]$observed_postprocessor_fits, 1L)
  testthat::expect_setequal(c(state$postprocessing[[1L]]$B, state$postprocessing[[1L]]$P), state$postprocessing[[1L]]$Q)
  testthat::expect_length(intersect(state$postprocessing[[1L]]$Q, state$postprocessing[[1L]]$C), 0L)
  core <- venn_abers_regression(state$calibration[[1L]]$truth,
                                state$calibration[[1L]]$prediction)
  expected_core <- predict(core, actual)
  details <- venn_abers_details(fit, new)
  for (column in c("prediction", "lower", "upper")) {
    testthat::expect_true(all(abs(details$details[[column]] - expected_core[[column]]) <=
      1e-10 + 1e-10 * abs(expected_core[[column]])))
  }
  testthat::expect_identical(details$prediction, predict(fit, new))
  testthat::expect_identical(names(venn_abers_details(fit, new[FALSE, , drop = FALSE])$details), names(details$details))
  testthat::expect_equal(predict(state$models[[1L]], new[1L, , drop = FALSE])$.pred, expected[1L])
  testthat::expect_error(predict(fit, new, type = "pred_int"), "type")
  flat <- d; flat$truth[7:12] <- 12
  flat_fit <- generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), flat)
  testthat::expect_equal(predict(workflows::extract_fit_engine(flat_fit)$models[[1L]], new)$.pred,
                         rep(12, nrow(new)), tolerance = 1e-10)
})

testthat::test_that("isotonic accepts two P rows, one support and arbitrary numeric scales", {
  skip_if_no_integration()
  for (scale in c(1e-9, 1, 1e9)) {
    d <- data.frame(x = c(1:6, 0, 0, 1:8), truth = c(1:6, 2, 8, 1:8) * scale)
    base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
      workflows::add_model(parsnip::linear_reg(engine = "lm"))
    va <- function(data) rsample::make_splits(list(analysis = 1:8, assessment = 9:16), data)
    post <- list(method = "isotonic_numeric", strategy = "holdout", splits = function(data)
      rsample::make_splits(list(analysis = 1:6, assessment = 7:8), data))
    fit <- generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d)
    model <- workflows::extract_fit_engine(fit)$models[[1L]]
    actual <- predict(model, data.frame(x = c(-2, 0, 3)))$.pred
    expected <- rep(5 * scale, 3)
    testthat::expect_true(all(abs(actual - expected) <= 1e-10 + 1e-10 * abs(expected)))
  }
})

testthat::test_that("CVAR consumes separately learned isotonic workflow outputs", {
  skip_if_no_integration()
  close <- function(actual, expected, scale) testthat::expect_true(all(
    abs(actual - expected) <= 1e-10 * scale + 1e-10 * abs(expected)))
  for (scale in c(1e-9, 1, 1e9)) for (option in list(list(m = 1), list(epsilon = .5))) {
    d <- data.frame(x = sin(1:80), truth = scale * (10 + 3 * sin(1:80) + cos(1:80)))
    base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
      workflows::add_model(parsnip::linear_reg(engine = "lm"))
    spec <- do.call(calibrated_workflow, c(list(workflow = base, outcome = "truth", method = "cross",
      v = 4, seed = 42, postprocessing = list(method = "isotonic_numeric", strategy = "holdout", prop = .3)), option))
    fit <- generics::fit(spec, d)
    state <- workflows::extract_fit_engine(fit)
    held <- new <- truth <- ids <- list()
    for (id in names(state$models)) {
      plan <- state$postprocessing[[id]]
      base_fit <- generics::fit(base, d[plan$B, , drop = FALSE])
      p <- predict(base_fit, d[plan$P, , drop = FALSE])$.pred
      iso <- stats::isoreg(p, d$truth[plan$P])
      sorted_rows <- order(p)
      distinct_support <- !duplicated(p[sorted_rows])
      x <- p[sorted_rows][distinct_support]; y <- iso$yf[distinct_support]
      reference <- function(data) stats::approx(x, y, xout = predict(base_fit, data)$.pred, rule = 2)$y
      held[[id]] <- reference(d[plan$C, , drop = FALSE])
      new[[id]] <- reference(d["x"])
      truth[[id]] <- d$truth[plan$C]
      ids[[id]] <- rep(id, length(plan$C))
      close(state$calibration[[id]]$prediction, held[[id]], scale)
      close(predict(state$models[[id]], d["x"])$.pred, new[[id]], scale)
    }
    core <- do.call(cross_venn_abers_regression, c(list(truth = unlist(truth),
      prediction = unlist(held), fold = unlist(ids)), option))
    expected <- predict(core, do.call(cbind, new), type = "details")
    actual <- venn_abers_details(fit, d["x"])
    close(actual$prediction$.pred, expected$prediction$prediction, scale)
    for (column in c("lower", "upper")) close(actual$aggregate[[column]], expected$prediction[[column]], scale)
    for (id in names(expected$folds)) for (column in c("prediction", "lower", "upper")) {
      close(actual$details[actual$details$fold == id, column][[1L]], expected$folds[[id]][[column]], scale)
    }
    testthat::expect_identical(actual$prediction, predict(fit, d["x"]))
  }
})

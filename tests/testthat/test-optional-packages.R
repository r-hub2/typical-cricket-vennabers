# Optional integration packages are checked at the call that needs them, with
# the minimum versions declared in DESCRIPTION Suggests.
declared_minimum <- function(package) {
  suggests <- trimws(strsplit(utils::packageDescription("vennabers")$Suggests, ",")[[1L]])
  entry <- suggests[sub("[[:space:]]*[(].*", "", suggests) == package]
  sub(".*>=[[:space:]]*([^)[:space:]]+).*", "\\1", entry)
}

integration <- c("parsnip", "workflows", "hardhat", "rsample", "recipes", "rlang",
                 "tibble", "withr", "generics", "tidyselect")

# Simulates an installed but outdated optional package without touching the library.
local_outdated <- function(package, env = parent.frame()) {
  original <- utils::packageVersion
  testthat::local_mocked_bindings(packageVersion = function(pkg, ...) {
    if (identical(pkg, package)) package_version("0.0.1") else original(pkg, ...)
  }, .package = "utils", .env = env)
}

testthat::test_that("tuning names an outdated tune, dials or yardstick at construction and fit", {
  for (pkg in c(integration, "tune", "dials", "yardstick")) testthat::skip_if_not_installed(pkg)
  data <- data.frame(x = rep(seq(-2, 2, length.out = 20), 3),
                     truth = factor(rep(c("yes", "no", "yes"), each = 20)))
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x, data = data) |>
      recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  config <- list(strategy = "two_stage", model = list(
    method = "grid", grid = data.frame(degree = c(1L, 2L)),
    resamples = function(d) rsample::vfold_cv(d, v = 2)))
  spec <- calibrated_workflow(base, "truth", seed = 42, tuning = config)
  for (package in c("tune", "dials", "yardstick")) local({
    minimum <- declared_minimum(package)
    testthat::expect_length(minimum, 1L)
    expected <- paste0("Tuning requires optional package ", sQuote(package), " >= ", minimum)
    local_outdated(package)
    testthat::expect_error(calibrated_workflow(base, "truth", seed = 42, tuning = config),
                           expected, fixed = TRUE)
    testthat::expect_error(generics::fit(spec, data), expected, fixed = TRUE)
  })
  # Without tuning, an outdated tuning package is irrelevant.
  fixed <- workflows::update_recipe(base, recipes::recipe(truth ~ x, data = data))
  local_outdated("tune")
  testthat::expect_s3_class(generics::fit(calibrated_workflow(fixed, "truth", seed = 42), data), "workflow")
})

testthat::test_that("workflow entry points name an outdated integration package", {
  for (pkg in integration) testthat::skip_if_not_installed(pkg)
  data <- data.frame(x = sin(1:80), truth = factor(rep(c("yes", "no", "no", "yes"), 20)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  earlier <- generics::fit(base, data[1:40, ])
  split <- rsample::initial_split(data)
  expected <- function(package) paste0("Workflow integration requires optional package ", sQuote(package),
                                       " >= ", declared_minimum(package))
  for (package in integration) local({
    testthat::expect_length(declared_minimum(package), 1L)
    local_outdated(package)
    testthat::expect_error(calibrated_workflow(base, "truth"), expected(package), fixed = TRUE)
  })
  # The same check guards the other entry points.
  local_outdated("parsnip")
  testthat::expect_error(calibrate_workflow(earlier, data[41:60, ], "truth"), expected("parsnip"), fixed = TRUE)
  testthat::expect_error(local_calibration_splits(split, c("x", "truth")), expected("parsnip"), fixed = TRUE)
})

testthat::test_that("Tailor postprocessing names an outdated tailor or probably; Sigmoid needs neither", {
  for (pkg in c(integration, "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  data <- data.frame(x = sin(seq_len(160)), truth = factor(rep(c("yes", "no", "no", "yes"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  tailor_post <- list(method = "tailor_logistic", strategy = "holdout", prop = .3)
  sigmoid_post <- list(method = "sigmoid", strategy = "holdout", prop = .3, sigmoid_input = "probability")
  spec <- calibrated_workflow(base, "truth", seed = 1, postprocessing = tailor_post)
  for (package in c("tailor", "probably")) local({
    minimum <- declared_minimum(package)
    testthat::expect_length(minimum, 1L)
    expected <- paste0("Postprocessing requires optional package ", sQuote(package), " >= ", minimum)
    local_outdated(package)
    testthat::expect_error(calibrated_workflow(base, "truth", postprocessing = tailor_post), expected, fixed = TRUE)
    testthat::expect_error(generics::fit(spec, data), expected, fixed = TRUE)
    fitted <- generics::fit(calibrated_workflow(base, "truth", seed = 1, postprocessing = sigmoid_post), data)
    testthat::expect_identical(nrow(predict(fitted, data[1:3, "x", drop = FALSE], type = "prob")), 3L)
  })
})

testthat::test_that("skip_if_no_integration() skips below a declared minimum and runs otherwise", {
  for (pkg in c(integration, "tune")) testthat::skip_if_not_installed(pkg)
  # In the full environment the shared test guard must not skip.
  testthat::expect_no_condition(skip_if_no_integration("tune"), class = "skip")
  for (package in c("recipes", "tidyselect", "tune")) local({
    local_outdated(package)
    testthat::expect_condition(skip_if_no_integration("tune"), declared_minimum(package), class = "skip")
  })
  local_outdated("tune")
  testthat::expect_no_condition(skip_if_no_integration(), class = "skip")
})

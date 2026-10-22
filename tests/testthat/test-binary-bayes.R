testthat::test_that("binary Bayes selection records strict, reproducible steps", {
  testthat::skip_if_not_installed("tune")
  data <- data.frame(x = sin(seq_len(100) / 7),
                     truth = factor(rep(c("yes", "no"), 50)))
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x, data) |>
      recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  params <- dials::parameters(dials::new_quant_param(
    type = "integer", range = c(1L, 4L), inclusive = c(TRUE, TRUE),
    label = c(degree = "Degree")))
  search <- list(method = "bayes", initial = 2L, iter = 2L, param_info = params,
    resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
    control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 2L))
  config <- list(strategy = "two_stage", model = search)
  first <- generics::fit(calibrated_workflow(base, "truth", tuning = config, seed = 7), data)
  second <- generics::fit(calibrated_workflow(base, "truth", tuning = config, seed = 7), data)
  choices <- venn_abers_tuning(first)$models$Inductive1
  testthat::expect_true(choices$chosen %in% choices$candidates$id)
  testthat::expect_true(all(choices$candidates$status %in% c("eligible", "excluded")))
  testthat::expect_true(choices$stop_reason %in% c("budget", "no_improve", "time_limit", "repeated_candidate"))
  testthat::expect_identical(choices$selected, venn_abers_tuning(second)$models$Inductive1$selected)
  testthat::expect_true(length(choices$steps) >= 1L)
  search$control <- tune::control_bayes(seed = 4, save_pred = TRUE, uncertain = 1)
  testthat::expect_error(calibrated_workflow(base, "truth",
    tuning = list(strategy = "two_stage", model = search)), "uncertain = Inf")
})

testthat::test_that("sized pair and search fits replay on outer workers and restore without learning", {
  for (pkg in c("tune", "yardstick", "future", "recipes", "dials")) testthat::skip_if_not_installed(pkg)
  testthat::skip_if_not(worker_replay_available())
  # Outputs plus every stored plan with U and the resolved counts; for searches
  # also the selections, roles and costs.
  snapshot <- function(fit, new, schema = NULL) {
    state <- workflows::extract_fit_engine(fit)
    tuning <- if (!is.null(state$tuning)) venn_abers_tuning(fit)
    selection <- function(x) if (!is.null(x)) x[c("chosen", "selected", "candidates", "grid", "resamples", "seed",
      "budget", "cost", "scores", "selection_rows", "final_rows")]
    list(outputs = lapply(c(0L, 1L, 7L), function(n) {
      batch <- new[seq_len(n), , drop = FALSE]
      list(prediction = predict(fit, batch, type = "prob"), details = venn_abers_details(fit, batch))
    }), plans = if (is.null(state$pair_models)) state$folds else lapply(state$pair_models, `[[`, "folds"),
      effective = state$parameters$effective[c("training_size", "calibration_size", "unassigned_size")],
      models = lapply(tuning$models, selection), calibration = selection(tuning$calibration),
      roles = tuning$fit_roles, cost = tuning$cost)
  }
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  folds <- function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)
  data <- withr::with_seed(91, data.frame(x = stats::rnorm(240), truth = factor(rep(c("a", "b", "c"), 80))))
  data$x <- data$x + as.integer(data$truth)
  cases <- list(
    pair = list(data = data, args = list(classification = "ovo", model_training = "pairwise",
      training_size = 30, calibration_size = 20, prediction_model = "full")),
    search = list(data = droplevels(data[data$truth != "c", ]), args = list(training_size = .5, calibration_size = .25,
      loss = tune::tune("va_loss"), tuning = list(strategy = "joint", selection_prop = .5,
        model = list(method = "grid", grid = data.frame(degree = 1:2), resamples = folds),
        calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")), resamples = folds)))))
  for (case in names(cases)) {
    d <- cases[[case]]$data
    base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, d) |>
      recipes::step_poly(x, degree = if (case == "search") tune::tune("degree") else 1L)) |>
      workflows::add_model(parsnip::logistic_reg(engine = "glm"))
    model <- do.call(calibrated_workflow, c(list(workflow = base, outcome = "truth", seed = 31),
      cases[[case]]$args), quote = TRUE)
    outer <- withr::with_seed(92, rsample::vfold_cv(d, v = 2))
    saveRDS(list(snapshot = snapshot, replay = list(model = model, outer = outer, new = d[1:7, ], schema = NULL)), saved)
    run_guarded_replay(saved, "yardstick::metric_set(yardstick::mn_log_loss)", c("vennabers", "recipes"), info = case,
      checks = "for (x in restored) stopifnot(all(unlist(s$snapshot(x$workflow, r$new)$effective$unassigned_size) > 0L))")
  }
})

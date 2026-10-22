testthat::test_that("tuned postprocessing examples evaluate natively and predict from RDS in fresh R without learning", {
  for (pkg in c("tailor", "probably", "recipes", "tune", "yardstick", "nnet")) testthat::skip_if_not_installed(pkg)
  example <- new.env()
  sys.source(system.file("examples", "postprocessing-grid.R", package = "vennabers"), example)
  testthat::expect_named(example$cases, c("binary_sigmoid_ensemble", "pairwise_tailor_joint",
    "ovr_multinomial_options", "regression_isotonic_joint", "regression_linear_ensemble"))
  snapshot <- function(fit, new) lapply(c(0L, 1L, 7L), function(n) {
    batch <- new[seq_len(n), , drop = FALSE]
    regression <- identical(workflows::extract_fit_engine(fit)$parameters$requested$mode, "regression")
    list(prediction = predict(fit, batch, type = if (regression) "numeric" else "prob"),
      details = venn_abers_details(fit, batch), tuning = venn_abers_tuning(fit)[c("models", "calibration", "joint", "cost")])
  })
  records <- lapply(example$cases, function(r) {
    r$expected <- snapshot(r$fit, r$new)
    r
  })
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  saveRDS(records, saved)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "library(vennabers); loadNamespace('workflows')",
    "guard <- function(...) stop('Unexpected learning, planning or search')",
    "registerS3method('fit', 'workflow', guard, envir=asNamespace('generics'))",
    "registerS3method('fit', 'tailor', guard, envir=asNamespace('generics'))",
    "for (f in c('optim', 'glm.fit', 'isoreg')) trace(f, where=asNamespace('stats'), tracer=quote(guard()), print=FALSE)",
    "trace('prep.recipe', where=asNamespace('recipes'), tracer=quote(guard()), print=FALSE)",
    "for (f in c('initial_split', 'vfold_cv')) trace(f, where=asNamespace('rsample'), tracer=quote(guard()), print=FALSE)",
    "for (f in c('tune_grid', 'tune_bayes')) trace(f, where=asNamespace('tune'), tracer=quote(guard()), print=FALSE)",
    paste0("snapshot <- ", paste(deparse(snapshot), collapse = "\n")),
    paste0("records <- readRDS(", deparse(saved), ")"),
    "for (r in records) stopifnot(identical(r$expected, snapshot(r$fit, r$new)))"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

testthat::test_that("outer parallel workers reproduce sequential tuned postprocessing", {
  for (pkg in c("tailor", "probably", "recipes", "tune", "yardstick", "future")) testthat::skip_if_not_installed(pkg)
  script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(script, log)), add = TRUE)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "Sys.setenv(OMP_NUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1', MKL_NUM_THREADS = '1')",
    "d <- withr::with_seed(8, { x <- stats::rnorm(600, sd = 1.5); data.frame(x = x,",
    "  truth = factor(ifelse(stats::runif(600) < stats::plogis(x - .4 * x^2 + .5), 'yes', 'no'))) })",
    "base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, d) |>",
    "  recipes::step_poly(x, degree = tune::tune('degree'))) |> workflows::add_model(parsnip::logistic_reg(engine = 'glm'))",
    "stratified <- function(x) rsample::vfold_cv(x, v = 2, strata = truth, pool = 0)",
    "model <- vennabers::calibrated_workflow(base, 'truth', method = 'cross', v = 2, prediction_model = 'full', seed = 3,",
    "  tuning = list(strategy = 'joint', model = list(method = 'grid', grid = data.frame(degree = 1:2), resamples = stratified)),",
    "  postprocessing = list(method = 'sigmoid', sigmoid_input = 'probability', strategy = 'holdout', prop = .3))",
    "outer <- withr::with_seed(31, rsample::vfold_cv(d, v = 3))",
    "control <- tune::control_resamples(save_pred = TRUE, extract = vennabers::venn_abers_tuning, parallel_over = 'resamples', pkgs = c('vennabers', 'stats'))",
    "run <- function() withr::with_seed(303, tune::fit_resamples(model, outer, metrics = yardstick::metric_set(yardstick::mn_log_loss), control = control))",
    "serial <- run()",
    "future::plan(future::multisession, workers = 2)",
    "parallel <- run()",
    "future::plan(future::sequential)",
    "stopifnot(nrow(tune::collect_notes(serial)) == 0L)",
    "stopifnot(nrow(tune::collect_notes(parallel)) == 0L)",
    "stopifnot(identical(tune::collect_predictions(serial), tune::collect_predictions(parallel)))",
    "choices <- function(r) lapply(tune::collect_extracts(r)$.extracts, function(t) lapply(t$models, function(p) lapply(p, `[[`, 'chosen')))",
    "stopifnot(length(unlist(choices(parallel))) == 9L)",
    "stopifnot(identical(choices(serial), choices(parallel)))"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

testthat::test_that("postprocessed Grid ensembles preserve complete searches and costs on outer workers", {
  for (pkg in c("tune", "yardstick", "future", "withr", "recipes", "rpart", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  source("helper-helpers.R")
  # Nested joint ensembles need B large enough for a nonconstant tree and
  # full-rank Tailor P; a small constant-tree fixture correctly fails instead.
  data <- withr::with_seed(81, {
    x <- rnorm(1920)
    data.frame(x = x, truth = 3 * x + rnorm(1920, sd = .5))
  })
  base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, data) |>
    recipes::step_normalize(x) |> recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::decision_tree(tree_depth = tune::tune("depth"), min_n = 10) |>
      parsnip::set_mode("regression") |> parsnip::set_engine("rpart"))
  outer <- withr::with_seed(82, rsample::vfold_cv(data, v = 2, repeats = 2))
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  for (method in c("isotonic_numeric", "tailor_linear")) for (strategy in c("joint", "two_stage")) {
    model <- calibrated_workflow(base, "truth", method = "cross", v = 2, prediction_model = "full", seed = 31,
      tuning = list(strategy = strategy, model = list(method = "grid", grid = data.frame(degree = 1:2, depth = 1:2),
        resamples = function(d) rsample::vfold_cv(d, v = 2))),
      postprocessing = list(method = method, strategy = "ensemble", v = 2))
    saveRDS(list(snapshot = postprocessing_replay_snapshot, replay = list(model = model, outer = outer, new = data[1:7, ], schema = NULL)), saved)
    run_guarded_replay(saved, "yardstick::metric_set(yardstick::rmse)", c("vennabers", "rpart", "probably"),
      info = paste(method, strategy))
  }
})

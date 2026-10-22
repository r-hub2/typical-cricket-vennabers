testthat::test_that("Bayes postprocessing examples evaluate natively and predict from RDS in fresh R without learning", {
  for (pkg in c("tailor", "probably", "recipes", "tune", "dials", "yardstick", "nnet")) testthat::skip_if_not_installed(pkg)
  example <- new.env()
  sys.source(system.file("examples", "postprocessing-bayes.R", package = "vennabers"), example)
  testthat::expect_named(example$cases, c("binary_sigmoid_joint", "pairwise_tailor_two_stage",
    "ovr_multinomial_options", "regression_linear_joint"))
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

testthat::test_that("outer parallel workers reproduce sequential Bayes postprocessing with failed candidates", {
  for (pkg in c("tailor", "probably", "recipes", "tune", "dials", "yardstick", "future")) testthat::skip_if_not_installed(pkg)
  script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(script, log)), add = TRUE)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "Sys.setenv(OMP_NUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1', MKL_NUM_THREADS = '1')",
    "d <- withr::with_seed(8, { x <- stats::rnorm(900, sd = 1.5); z <- sample(1:4, 900, TRUE); data.frame(x = x, z = z,",
    "  truth = factor(ifelse(stats::runif(900) < stats::plogis(x - .3 * z + .5), 'yes', 'no'))) })",
    "base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x + z, d) |>",
    "  recipes::step_poly(z, degree = tune::tune('degree'))) |> workflows::add_model(parsnip::logistic_reg(engine = 'glm'))",
    "stratified <- function(x) rsample::vfold_cv(x, v = 2, strata = truth, pool = 0)",
    "params <- dials::parameters(dials::new_quant_param(type = 'integer', range = c(1L, 6L), inclusive = c(TRUE, TRUE), label = c(degree = 'Degree')))",
    "search <- list(method = 'bayes', initial = 3L, iter = 3L, param_info = params, resamples = stratified,",
    "  control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 3L))",
    "model <- vennabers::calibrated_workflow(base, 'truth', method = 'cross', v = 2, prediction_model = 'full', seed = 3,",
    "  tuning = list(strategy = 'joint', model = search),",
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
    "history <- function(r) lapply(tune::collect_extracts(r)$.extracts, function(t) lapply(t$models, function(p)",
    "  lapply(p, function(s) list(s$chosen, s$candidates[c('id', 'degree', 'status')], s$stop_reason, s$attempts))))",
    "stopifnot(length(unlist(lapply(history(parallel), function(t) lapply(t, names)))) == 9L)",
    "stopifnot(identical(history(serial), history(parallel)))",
    "statuses <- unlist(lapply(history(serial), function(t) lapply(t, function(p) lapply(p, function(s) s[[2L]]$status))))",
    "stopifnot(any(statuses == 'excluded'))"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

testthat::test_that("postprocessed Bayes ensembles replay search, failure and cost histories on outer workers", {
  for (pkg in c("tune", "dials", "yardstick", "future", "withr", "recipes", "rpart", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  source("helper-helpers.R")
  # z has four distinct values, so step_poly(z) fails for degree 4 in every B.
  # The space-filling initial points 1/1, 2/2, 3/1 and 4/2 (degree/depth)
  # leave the three complete points the GP needs; proposal 4/1 fails as well.
  data <- withr::with_seed(81, {
    x <- rnorm(1920)
    z <- rep(1:4, 480)
    data.frame(x = x, z = z, truth = 3 * x + z + rnorm(1920, sd = .5))
  })
  base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x + z, data) |>
    recipes::step_normalize(x) |> recipes::step_poly(z, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::decision_tree(tree_depth = tune::tune("depth"), min_n = 10) |>
      parsnip::set_mode("regression") |> parsnip::set_engine("rpart"))
  integer_param <- function(id, upper) dials::new_quant_param(type = "integer", range = c(1L, upper),
    inclusive = c(TRUE, TRUE), label = stats::setNames(id, id))
  search <- list(method = "bayes", initial = 4L, iter = 2L,
    param_info = dials::parameters(integer_param("degree", 4L), integer_param("depth", 2L)),
    resamples = function(d) rsample::vfold_cv(d, v = 2),
    control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 2L))
  outer <- withr::with_seed(82, rsample::vfold_cv(data, v = 2, repeats = 2))
  # Every restored local search keeps its failed degree-4 candidates, with
  # reasons and notes, outside the GP input.
  failed <- c("for (x in restored) for (path in venn_abers_tuning(x$workflow)$models) for (search in path) {",
    "  bad <- search$candidates$degree == 4L",
    "  stopifnot(any(bad), all(search$candidates$status[bad] == 'excluded'), all(!is.na(search$candidates$reason[bad])),",
    "    !any(tune::collect_metrics(search$native)$degree == 4L), grepl('unique points', paste(search$notes$note, collapse = ' ')))",
    "}")
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  for (method in c("isotonic_numeric", "tailor_linear")) for (strategy in c("joint", "two_stage")) {
    model <- calibrated_workflow(base, "truth", method = "cross", v = 2, prediction_model = "full", seed = 31,
      tuning = list(strategy = strategy, model = search),
      postprocessing = list(method = method, strategy = "ensemble", v = 2))
    saveRDS(list(snapshot = postprocessing_replay_snapshot,
      replay = list(model = model, outer = outer, new = data[1:7, ], schema = NULL)), saved)
    run_guarded_replay(saved, "yardstick::metric_set(yardstick::rmse)", c("vennabers", "rpart", "probably"),
      info = paste(method, strategy), checks = failed)
  }
})

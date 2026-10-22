# Replays the saved `s$replay` on two fresh workers in a fresh process, then
# restores every worker fit from RDS under guards against learning, planning
# or search; `s$snapshot` must match the sequential reference throughout.
# Optional `checks` are script lines run last on the `restored` worker fits.
run_guarded_replay <- function(saved, metrics, packages, info = "", checks = character()) {
  fitted <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(fitted, script, log)), add = TRUE)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"), "library(vennabers)",
    paste0("s <- readRDS(", deparse(saved), ")"),
    worker_replay(metrics, packages, character()),
    "checker <- new.env(); sys.source(system.file('examples', 'check-encoded-evaluation.R', package = 'vennabers'), checker)",
    "checker$check_outer_identities(parallel, r$outer, metrics = 1L)",
    paste0("saveRDS(found, ", deparse(fitted), ")"),
    paste0("restored <- readRDS(", deparse(fitted), ")"),
    "guard <- function(...) stop('Unexpected learning, planning or search')",
    "registerS3method('fit', 'workflow', guard, envir = asNamespace('generics'))",
    "registerS3method('fit', 'tailor', guard, envir = asNamespace('generics'))",
    "for (f in c('venn_abers', 'venn_abers_regression')) trace(f, where = asNamespace('vennabers'), tracer = quote(guard()), print = FALSE)",
    "trace('isoreg', where = asNamespace('stats'), tracer = quote(guard()), print = FALSE)",
    "trace('prep.recipe', where = asNamespace('recipes'), tracer = quote(guard()), print = FALSE)",
    "for (f in c('initial_split', 'vfold_cv')) trace(f, where = asNamespace('rsample'), tracer = quote(guard()), print = FALSE)",
    "for (f in c('tune_grid', 'tune_bayes')) trace(f, where = asNamespace('tune'), tracer = quote(guard()), print = FALSE)",
    "for (i in seq_along(restored)) stopifnot(identical(s$snapshot(restored[[i]]$workflow, r$new), s$snapshot(expected[[i]], r$new)))",
    checks
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(info, paste(readLines(log), collapse = "\n")))
}

# Complete stored state of a tuned postprocessed fit for run_guarded_replay():
# outputs for 0/1/7 rows, every local Grid or Bayes search record with its
# native predictions, B/P/C plans with seeds and RNG kind, costs and roles.
postprocessing_replay_snapshot <- function(fit, new, schema = NULL) {
  state <- workflows::extract_fit_engine(fit)
  tuning <- venn_abers_tuning(fit)
  selection <- function(x) c(x[c("chosen", "selected", "candidates", "grid", "resamples", "seed",
    "search_seed", "budget", "cost", "scores", "report_metrics", "decisions",
    "bayes_seed", "steps", "notes", "attempts", "stop_reason", "no_improve")],
    list(predictions = tune::collect_predictions(x$native)))
  list(outputs = lapply(c(0L, 1L, 7L), function(n) {
    batch <- new[seq_len(n), , drop = FALSE]
    list(prediction = predict(fit, batch), details = venn_abers_details(fit, batch))
  }), models = lapply(tuning$models, function(path) lapply(path, selection)),
    plans = state$postprocessing, cost = tuning$cost, roles = tuning$fit_roles, sources = tuning$sources)
}

# Shared checks for the independently scheduled search/domain replay tests.
check_tuning_replay <- function(search, domain) {
  for (pkg in c("tune", "yardstick", "future", "recipes", "rpart")) testthat::skip_if_not_installed(pkg)
  source("helper-helpers.R")
  snapshot <- function(fit, new, schema = NULL) {
    tuning <- venn_abers_tuning(fit)
    selection <- function(x) {
      if (is.null(x)) return(NULL)
      x[c("chosen", "selected", "candidates", "grid", "resamples", "seed", "budget", "cost",
        "bayes_seed", "steps", "notes", "attempts", "stop_reason", "no_improve", "scores", "report_metrics")]
    }
    models <- if (identical(tuning$model_training, "pairwise")) lapply(tuning$models, function(x) lapply(x, selection)) else
      lapply(tuning$models, selection)
    list(outputs = lapply(c(0L, 1L, 7L), function(n) {
      batch <- new[seq_len(n), , drop = FALSE]
      list(prediction = predict(fit, batch, type = if (is.numeric(new$truth)) "numeric" else "prob"),
        details = venn_abers_details(fit, batch))
    }), models = models, calibration = selection(tuning$calibration), joint = selection(tuning$joint),
      cost = tuning$cost, roles = tuning$fit_roles, sources = tuning$sources)
  }
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  for (strategy in c("joint", "two_stage")) {
    data <- withr::with_seed(81, data.frame(x = stats::rnorm(720),
      truth = if (domain == "regression") stats::rnorm(720) else
        factor(rep(if (domain == "binary") c("a", "b") else c("a", "b", "c"), length.out = 720))))
    base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, data) |>
      recipes::step_poly(x, degree = tune::tune("degree"))) |>
      workflows::add_model(parsnip::decision_tree(tree_depth = tune::tune("depth"), min_n = 10) |>
        parsnip::set_mode(if (domain == "regression") "regression" else "classification") |>
        parsnip::set_engine("rpart"))
    folds <- function(d) rsample::vfold_cv(d, v = 2)
    grid <- data.frame(degree = 1:2, depth = 1:2)
    # The Cartesian pair route is already covered with two model candidates;
    # one here keeps worker replay small while still tuning calibration.
    if (search == "grid" && domain == "pairwise" && strategy == "joint") grid <- grid[1L, , drop = FALSE]
    tuning <- list(strategy = strategy, selection_prop = .5,
      model = list(method = "grid", grid = grid, resamples = folds),
      calibration = list(method = "grid", resamples = folds,
        grid = if (domain == "regression") data.frame(moderation = c(.2, .3)) else
          data.frame(va_loss = c("log", "brier"))))
    if (search == "grid" && domain == "pairwise" && strategy == "joint") tuning$model$resamples <- NULL
    if (search == "bayes") {
      params <- dials::parameters(
        dials::new_quant_param(type = "integer", range = c(1L, 2L), inclusive = c(TRUE, TRUE),
          label = c(degree = "Degree")),
        dials::new_quant_param(type = "integer", range = c(1L, 3L), inclusive = c(TRUE, TRUE),
          label = c(depth = "Depth")))
      options <- if (domain == "regression") dials::new_quant_param(type = "double",
        range = c(.2, .3), inclusive = c(TRUE, TRUE), label = c(moderation = "Moderation")) else
        dials::new_qual_param(type = "character", values = c("log", "brier"), label = c(va_loss = "Loss"))
      tuning$model <- list(method = "bayes", initial = 3L, iter = 2L, param_info = params,
        resamples = folds, control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 2L))
      tuning$calibration <- list(method = "bayes", initial = 2L, iter = 1L,
        param_info = dials::parameters(options), resamples = folds,
        control = tune::control_bayes(seed = 5, save_pred = TRUE, no_improve = 2L))
    }
    args <- list(workflow = base, outcome = "truth", method = "cross", v = 2,
      prediction_model = "full", seed = 31, tuning = tuning)
    if (domain == "regression") args$epsilon <- tune::tune("moderation") else args$loss <- tune::tune("va_loss")
    if (domain %in% c("ovo", "ovr", "pairwise")) args$classification <- if (domain == "pairwise") "ovo" else domain
    if (domain == "pairwise") args$model_training <- "pairwise"
    model <- do.call(calibrated_workflow, args, quote = TRUE)
    outer <- withr::with_seed(82, rsample::vfold_cv(data, v = 2, repeats = 2))
    saveRDS(list(snapshot = snapshot, replay = list(model = model, outer = outer, new = data[1:7, ], schema = NULL)), saved)
    metrics <- if (domain == "regression") "yardstick::metric_set(yardstick::rmse)" else
      "yardstick::metric_set(yardstick::mn_log_loss)"
    run_guarded_replay(saved, metrics, c("vennabers", "rpart"), info = paste(search, domain, strategy))
  }
}

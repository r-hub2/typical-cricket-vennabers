testthat::test_that("schema-multiclass-resampling", {
library(vennabers)
if (all(vapply(c("tune", "yardstick", "rsample", "workflows", "recipes", "nnet", "future"),
               requireNamespace, logical(1), quietly = TRUE))) local({
  example <- new.env()
  sys.source(system.file("examples", "encoded-multiclass-regression.R", package = "vennabers"), example)
  check <- example$check_encoded_evaluation
  error <- function(expr) {
    e <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(e, "error"))
  }
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); output <- tempfile(fileext = ".log")
  on.exit(unlink(c(saved, script, output)), add = TRUE)
  for (case in example$examples) {
    schema <- case$schema; data <- case$data; outcome <- schema$encoded_outcome
    base <- case$base; outer <- case$outer
    regression <- schema$mode == "regression"
    method <- workflows::extract_fit_engine(case$fit)$parameters$requested$method
    scoring <- if (regression) yardstick::metric_set(yardstick::rmse, yardstick::mae) else
      yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
    control <- tune::control_resamples(save_pred = TRUE, extract = identity, allow_par = FALSE)
    run <- function(model, metrics = scoring) withr::with_seed(43,
      tune::fit_resamples(model, outer, metrics = metrics, control = control))
    automatic <- case$evaluation
    checked <- check(automatic, outer, schema)
    states <- tune::collect_extracts(automatic)$.extracts
    prepared <- lapply(seq_along(states), function(i) {
      state <- workflows::extract_fit_engine(states[[i]])
      local_data <- rsample::analysis(outer$splits[[i]])
      held <- rsample::assessment(outer$splits[[i]])
      direct <- venn_abers_decode(venn_abers_details(states[[i]], held), schema)
      rows <- checked$decoded$identity$id == outer$id[i] & checked$decoded$identity$id2 == outer$id2[i]
      stopifnot(identical(state$schema, schema), identical(names(state$models), names(state$folds)))
      if (regression) {
        stopifnot(identical(direct$prediction, checked$decoded$prediction[rows]))
      } else {
        stopifnot(identical(unname(direct$probability), unname(checked$decoded$probability[rows, ])))
        stopifnot(abs(brier_score(checked$decoded$truth[rows], checked$decoded$probability[rows, ]) / 2 -
          checked$metrics$.estimate[checked$metrics$id == outer$id[i] & checked$metrics$id2 == outer$id2[i] &
                                      checked$metrics$.metric == "brier_class"]) <= 1e-10)
      }
      splits <- lapply(names(state$folds), function(id) {
        f <- state$folds[[id]]
        recipe <- workflows::extract_recipe(state$models[[id]])
        stopifnot(abs(recipe$steps[[1]]$means[[1]] - mean(local_data[[1]][f$analysis])) <= 1e-10)
        rsample::make_splits(f[c("analysis", "assessment")], local_data)
      })
      if (method == "cross") rsample::manual_rset(splits, names(state$folds)) else splits[[1]]
    })
    resolver <- local_calibration_splits(prepared, names(data))
    seen <- list()
    callback <- function(local) {
      stopifnot(setequal(names(local), names(data)))
      if (!regression) stopifnot(identical(levels(local[[outcome]]), unname(schema$classes)))
      seen[[length(seen) + 1L]] <<- local[[1]]
      resolver(local)
    }
    planned <- calibrated_workflow(base, outcome, schema = schema, method = method, splits = callback, seed = 42)
    explicit <- run(planned)
    check(explicit, outer, schema)
    stopifnot(identical(seen, lapply(outer$splits, function(s) rsample::analysis(s)[[1]])),
              identical(tune::collect_predictions(automatic, summarize = FALSE), tune::collect_predictions(explicit, summarize = FALSE)),
              identical(tune::collect_metrics(automatic, summarize = FALSE), tune::collect_metrics(explicit, summarize = FALSE)))
    global <- if (method == "cross") rsample::vfold_cv(data, v = 2) else rsample::initial_split(data)
    bad <- suppressWarnings(run(calibrated_workflow(base, outcome, schema = schema, method = method, splits = global)))
    error(check(bad, outer, schema))
    stopifnot(nrow(tune::collect_notes(bad)) == nrow(outer))
    fail_metric <- if (regression) yardstick::new_numeric_metric(function(data, ...) stop("metric failed"), direction = "minimize") else
      yardstick::new_prob_metric(function(data, ...) stop("metric failed"), direction = "minimize")
    bad <- suppressWarnings(run(case$model, if (regression) yardstick::metric_set(yardstick::rmse, fail_metric) else
                                 yardstick::metric_set(yardstick::brier_class, fail_metric)))
    error(check(bad, outer, schema))
    stopifnot(nrow(tune::collect_notes(bad)) > 0L)
    first <- rsample::analysis(outer$splits[[1]])[[1]]
    partial <- function(local) {
      if (identical(local[[1]], first)) stop("local failure")
      resolver(local)
    }
    bad <- suppressWarnings(run(calibrated_workflow(base, outcome, schema = schema, method = method, splits = partial)))
    error(check(bad, outer, schema))
    stopifnot(nrow(tune::collect_notes(bad)) > 0L, nrow(tune::collect_predictions(bad)) < nrow(checked$decoded$identity))
    # Failed inner models cannot be replaced by surviving fold predictions.
    broken <- if (regression) parsnip::linear_reg(engine = "lm") |> parsnip::set_engine("lm", qr = NA) else
      parsnip::multinom_reg(engine = "nnet") |> parsnip::set_engine("nnet", MaxNWts = 1)
    bad <- suppressWarnings(run(calibrated_workflow(workflows::update_model(base, broken), outcome,
                                                   schema = schema, method = method, seed = 42)))
    error(check(bad, outer, schema))
    stopifnot(nrow(tune::collect_notes(bad)) == nrow(outer))
    for (change in c("missing", "duplicate", "config", "metric", "truth")) {
      bad <- automatic
      if (change == "missing") bad$.predictions[[1]] <- bad$.predictions[[1]][-1, ]
      if (change == "duplicate") bad$.predictions[[1]] <- rbind(bad$.predictions[[1]], bad$.predictions[[1]][1, ])
      if (change == "config") bad$.predictions[[1]]$.config <- "unexpected"
      if (change == "metric") bad$.metrics[[1]] <- bad$.metrics[[1]][-1, ]
      if (change == "truth") bad$.predictions[[1]][[outcome]][1] <- NA
      error(check(bad, outer, schema))
    }
    # Full native and decoded roles survive fresh-process restore and two workers.
    snapshot <- function(fit, data, schema) lapply(c(0L, 1L, 3L), function(n) {
      new <- data[seq_len(n), setdiff(names(data), schema$encoded_outcome), drop = FALSE]
      native <- list(prediction = predict(fit, new, type = if (schema$mode == "regression") "numeric" else "prob"),
                     details = venn_abers_details(fit, new))
      if (schema$mode != "regression") native$class <- predict(fit, new, type = "class")
      list(native = native, decoded = lapply(native, venn_abers_decode, schema = schema))
    })
    saveRDS(list(case = case, snapshot = snapshot, expected = snapshot(case$fit, data, schema), checked = checked), saved)
    writeLines(c(
      paste0('.libPaths(', paste(deparse(.libPaths()), collapse = ''), ')'),
      'library(vennabers); loadNamespace("workflows"); loadNamespace("stats"); loadNamespace("nnet")',
      paste0('s <- readRDS(', deparse(saved), '); c <- s$case'),
      'source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"))',
      'stopifnot(identical(s$snapshot(c$fit, c$data, c$schema), s$expected),',
      '  identical(check_encoded_evaluation(c$evaluation, c$outer, c$schema), s$checked))',
      'Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")',
      'metrics <- if (c$schema$mode == "regression") yardstick::metric_set(yardstick::rmse, yardstick::mae) else',
      '  yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)',
      'control <- tune::control_resamples(save_pred = TRUE, extract = identity, parallel_over = "resamples",',
      '  pkgs = c("vennabers", "stats", "nnet"))',
      'future::plan(future::multisession, workers = 2)',
      'one <- withr::with_seed(43, tune::fit_resamples(c$model, c$outer, metrics = metrics, control = control))',
      'future::plan(future::sequential)',
      'stopifnot(identical(tune::collect_predictions(one), tune::collect_predictions(c$evaluation)),',
      '  identical(tune::collect_metrics(one), tune::collect_metrics(c$evaluation)),',
      '  identical(check_encoded_evaluation(one, c$outer, c$schema), s$checked))',
      'a <- tune::collect_extracts(one)$.extracts; b <- tune::collect_extracts(c$evaluation)$.extracts',
      'for (i in seq_along(a)) {',
      '  stopifnot(identical(workflows::extract_fit_engine(a[[i]])$schema, c$schema),',
      '    identical(s$snapshot(a[[i]], c$data, c$schema), s$snapshot(b[[i]], c$data, c$schema)))',
      '}'
    ), script)
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = output, stderr = output)
    if (status != 0L) stop(paste(readLines(output), collapse = "\n"))
  }
  # Preserve metric conventions even at exact zero/one probability extremes.
  for (k in 2:4) {
    original <- data.frame(y = factor(seq_len(k), labels = c("class", ".pred_yes", "a b", "a.b")[seq_len(k)]))
    schema <- venn_abers_schema(original, "y", classification = "ovo")
    p <- diag(k)[, rev(seq_len(k)), drop = FALSE]
    native <- as.data.frame(p); names(native) <- paste0(".pred_", unname(schema$classes))
    decoded <- venn_abers_decode(native, schema)$probability
    estimate <- if (k == 2) decoded[, 1] else decoded
    stopifnot(is.infinite(log_loss(original$y, decoded)),
              is.finite(yardstick::mn_log_loss_vec(original$y, estimate)),
              abs(yardstick::brier_class_vec(original$y, estimate) - brier_score(original$y, decoded) / 2) <= 1e-10,
              identical(unname(decoded), p))
  }
})

  testthat::expect_true(TRUE)
})

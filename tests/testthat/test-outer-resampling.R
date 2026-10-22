testthat::test_that("outer-resampling", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")
if (all(vapply(c("tune", "yardstick", "rsample", "workflows", "recipes"),
               requireNamespace, logical(1), quietly = TRUE))) {
  # Run the installed example, including its public-collector completeness audit.
  example <- new.env()
  sys.source(system.file("examples", "calibrated-resampling.R", package = "vennabers"), example)
  check <- example$check_evaluation
  error <- function(expr, pattern = "") {
    result <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
  }
  audit <- new.env()
  audit$fits <- 0L
  audit_fit <- function(x, y, fail_id = NULL) {
    audit$fits <- audit$fits + 1L
    if (!is.null(fail_id) && any(x$id %in% fail_id)) stop("instrumented inner fit failure")
    list(x = x, y = y, levels = levels(y), draw = runif(1))
  }
  audit_prob <- function(object, new_data) {
    p <- plogis(new_data$x)
    stats::setNames(data.frame(p, 1 - p), paste0(".pred_", object$levels))
  }
  register_reference_engine("va_audit", "classification", "audit_fit", "audit_prob", "prob", "audit")
  spec <- parsnip::new_model_spec("va_audit", args = list(), eng_args = NULL,
                                 mode = "classification", method = NULL, engine = "audit")
  d <- data.frame(x = sin(1:120), id = 1:120, truth = factor(rep(c("yes", "no", "no", "yes"), 30)))
  rec <- recipes::recipe(truth ~ x + id, d) |> recipes::step_normalize(x)
  base <- workflows::workflow() |> workflows::add_recipe(rec) |> workflows::add_model(spec)
  scoring <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
  control <- function(event = "first", extract = workflows::extract_fit_engine) {
    tune::control_resamples(save_pred = TRUE, extract = extract, event_level = event,
                            parallel_over = "resamples", allow_par = FALSE, pkgs = "vennabers")
  }
  run <- function(model, outer, event = "first", metrics = scoring, extract = workflows::extract_fit_engine) {
    withr::with_seed(201, tune::fit_resamples(model, outer, metrics = metrics, control = control(event, extract)))
  }
  for (method in c("inductive", "cross")) {
    ordinary <- withr::with_seed(200, rsample::vfold_cv(d, v = 5))
    repeated <- withr::with_seed(200, rsample::vfold_cv(d, v = 3, repeats = 2))
    holdout <- rsample::manual_rset(list(rsample::make_splits(list(analysis = 1:80, assessment = 81:120), d)), "Test")
    for (outer in list(ordinary, repeated, holdout)) {
      for (event in c("first", "second")) for (loss in c("log", "brier")) {
        before <- audit$fits
        model <- calibrated_workflow(base, "truth", method = method, seed = 42, loss = loss, event_level = event)
        result <- run(model, outer, event)
        checked <- check(result, outer, event)
        states <- tune::collect_extracts(result)
        stopifnot(audit$fits == before + nrow(outer) * if (method == "cross") 5L else 1L)
        for (i in seq_len(nrow(outer))) {
          local_data <- rsample::analysis(outer$splits[[i]])
          held <- rsample::assessment(outer$splits[[i]])
          state <- states$.extracts[[i]]
          calibration <- truth <- fold <- list()
          new_p <- matrix(NA_real_, nrow(held), length(state$folds), dimnames = list(NULL, names(state$folds)))
          for (id in names(state$folds)) {
            inner <- state$folds[[id]]
            engine <- workflows::extract_fit_engine(state$models[[id]])
            stopifnot(identical(engine$x$id, local_data$id[inner$analysis]),
                      identical(as.vector(engine$y), as.vector(local_data$truth[inner$analysis])),
                      !any(engine$x$id %in% held$id),
                      !any(local_data$id[inner$assessment] %in% held$id))
            learned <- workflows::extract_recipe(state$models[[id]])
            stopifnot(abs(learned$steps[[1]]$means[["x"]] - mean(local_data$x[inner$analysis])) <= 1e-10)
            calibration[[id]] <- predict(state$models[[id]], local_data[inner$assessment, ], type = "prob")[[paste0(".pred_", state$event_class)]]
            truth[[id]] <- local_data$truth[inner$assessment] == state$event_class
            fold[[id]] <- rep(id, length(inner$assessment))
            new_p[, id] <- predict(state$models[[id]], held, type = "prob")[[paste0(".pred_", state$event_class)]]
          }
          if (method == "cross") {
            core <- cross_venn_abers(unlist(truth), unlist(calibration), unlist(fold))
            expected <- predict(core, new_p, loss = loss)
          } else {
            core <- venn_abers(truth[[1L]], calibration[[1L]])
            expected <- predict(core, new_p[, 1L], loss = loss)
          }
          ids <- setdiff(names(outer), "splits")
          rows <- Reduce(`&`, lapply(ids, function(id) checked$predictions[[id]] == outer[[id]][i]))
          actual <- checked$predictions[rows, ]
          stopifnot(max(abs(actual[[paste0(".pred_", state$event_class)]] - expected$probability)) <= 1e-10)
          direct_splits <- lapply(state$folds, function(f) rsample::make_splits(f[c("analysis", "assessment")], local_data))
          direct_plan <- if (method == "cross") rsample::manual_rset(direct_splits, names(state$folds)) else direct_splits[[1L]]
          direct <- generics::fit(calibrated_workflow(base, "truth", method = method, splits = direct_plan, loss = loss, event_level = event), local_data)
          details <- venn_abers_details(direct, held)
          stopifnot(max(abs(details$prediction[[paste0(".pred_", state$event_class)]] - expected$probability)) <= 1e-10)

        }
      }
    }
    # Automatic plans and identical prepared plans give the same complete results.
    prepared <- lapply(repeated$splits, function(split) {
      withr::with_seed(42,
        if (method == "cross") rsample::vfold_cv(rsample::analysis(split), v = 5, strata = truth, pool = 0) else
        rsample::initial_split(rsample::analysis(split), .8, strata = truth, pool = 0), .rng_kind = "L'Ecuyer-CMRG")
    })
    resolver <- local_calibration_splits(prepared, names(d))
    planned <- calibrated_workflow(base, "truth", method = method, splits = resolver, seed = 42, loss = "brier")
    explicit <- run(planned, repeated)
    automatic <- run(calibrated_workflow(base, "truth", method = method, seed = 42, loss = "brier"), repeated)
    stopifnot(identical(tune::collect_predictions(explicit), tune::collect_predictions(automatic)),
              identical(tune::collect_metrics(explicit), tune::collect_metrics(automatic)))
    # A callback records only the local analysis rows supplied by the engine.
    seen <- list()
    callback <- function(data) { seen[[length(seen) + 1L]] <<- data$id; resolver(data) }
    run(calibrated_workflow(base, "truth", method = method, splits = callback, seed = 42), repeated)
    stopifnot(identical(seen, lapply(repeated$splits, function(s) rsample::analysis(s)$id)))

    # Changing only outer truth leaves all fitted state and predictions unchanged.
    changed <- d
    changed$truth[81:120] <- factor(rep("yes", 40), levels = levels(d$truth))
    changed_outer <- rsample::manual_rset(list(rsample::make_splits(list(analysis = 1:80, assessment = 81:120), changed)), "Test")
    model <- calibrated_workflow(base, "truth", method = method, seed = 42, loss = "brier")
    before <- run(model, holdout)
    after <- run(model, changed_outer)
    check(before, holdout)
    check(after, changed_outer)
    s1 <- tune::collect_extracts(before)$.extracts[[1]]
    s2 <- tune::collect_extracts(after)$.extracts[[1]]
    assert_conditions(identical(s1[setdiff(names(s1), "models")], s2[setdiff(names(s2), "models")]),
              identical(lapply(s1$models, workflows::extract_fit_engine), lapply(s2$models, workflows::extract_fit_engine)),
              identical(lapply(s1$models, function(w) workflows::extract_recipe(w)$steps),
                        lapply(s2$models, function(w) workflows::extract_recipe(w)$steps)),
              identical(tune::collect_predictions(before)[c(".pred_no", ".pred_yes")],
                        tune::collect_predictions(after)[c(".pred_no", ".pred_yes")]),
              !identical(tune::collect_metrics(before)$mean, tune::collect_metrics(after)$mean))
    # Full workflow and compact diagnostic extraction both work through collectors.
    full <- run(model, holdout, extract = identity)
    full_fit <- tune::collect_extracts(full)$.extracts[[1]]
    stopifnot(inherits(full_fit, "workflow"), identical(predict(full_fit, d, type = "prob"), venn_abers_details(full_fit, d)$prediction))
    compact <- run(model, holdout, extract = function(w) workflows::extract_fit_engine(w)$folds)
    stopifnot(identical(tune::collect_extracts(compact)$.extracts[[1]], s1$folds))

    # An invalid global binding or reserved name fails before any base fit.
    before_count <- audit$fits
    global <- if (method == "cross") rsample::vfold_cv(d, v = 5) else rsample::initial_split(d)
    bad <- suppressWarnings(run(calibrated_workflow(base, "truth", method = method, splits = global), ordinary))
    error(check(bad, ordinary))
    stopifnot(audit$fits == before_count, nrow(tune::collect_notes(bad)) == nrow(ordinary))
    bad_names <- d
    levels(bad_names$truth) <- c("class", "yes")
    bad_outer <- rsample::manual_rset(list(rsample::make_splits(list(analysis = 1:80, assessment = 81:120), bad_names)), "Names")
    error(run(model, bad_outer), "reserved")
    stopifnot(audit$fits == before_count)

    # One inner failure, a metric-only error, silent NA removal, and total failure
    # must never pass the completeness audit just because a collector returns data.
    fail_base <- workflows::update_model(base, parsnip::set_engine(spec, "audit", fail_id = 1L))
    partial_outer <- rsample::manual_rset(list(holdout$splits[[1]],
      rsample::make_splits(list(analysis = 41:120, assessment = 1:40), d)), c("Fails", "Works"))
    fixed <- function(data) {
      split <- rsample::make_splits(list(analysis = 1:60, assessment = 61:80), data)
      if (method == "inductive") return(split)
      rsample::manual_rset(list(split, rsample::make_splits(list(analysis = 61:80, assessment = 1:60), data)), c("First", "Second"))
    }
    fail_model <- calibrated_workflow(fail_base, "truth", method = method, splits = fixed, seed = 42)
    fail <- suppressWarnings(run(fail_model, partial_outer))
    error(check(fail, partial_outer))
    stopifnot(nrow(tune::collect_notes(fail)) == 1L,
              nrow(tune::collect_predictions(fail)) == 40L)
    metric_error <- function(data, ...) stop("instrumented metric failure")
    metric_error <- yardstick::new_prob_metric(metric_error, direction = "minimize")
    bad_metrics <- yardstick::metric_set(yardstick::brier_class, metric_error)
    bad <- suppressWarnings(run(model, ordinary, metrics = bad_metrics))
    error(check(bad, ordinary))
    stopifnot(nrow(tune::collect_notes(bad)) > 0L)
    missing_truth <- d
    missing_truth$truth[81] <- NA
    missing_outer <- rsample::manual_rset(list(rsample::make_splits(list(analysis = 1:80, assessment = 81:120), missing_truth)), "Test")
    error(run(model, missing_outer), "nonmissing")
    brier_class <- function(data, ...) {
      data$truth[1] <- NA
      yardstick::brier_class(data, ...)
    }
    brier_class <- yardstick::new_prob_metric(brier_class, direction = "minimize")
    removed <- run(model, holdout, metrics = yardstick::metric_set(brier_class, yardstick::mn_log_loss))
    stopifnot(all(is.finite(tune::collect_metrics(removed)$mean)))
    error(check(removed, holdout))
    total_base <- workflows::update_model(base, parsnip::set_engine(spec, "audit", fail_id = 1:120))
    bad <- suppressWarnings(run(calibrated_workflow(total_base, "truth", method = method, seed = 42), ordinary))
    error(check(bad, ordinary))
    stopifnot(nrow(tune::collect_notes(bad)) == nrow(ordinary))
    # Deliberately corrupt successful collector tables to test missing/duplicate keys.
    for (change in c("missing", "duplicate", "config", "metric")) {
      bad <- before
      if (change == "missing") bad$.predictions[[1]] <- bad$.predictions[[1]][-1, ]
      if (change == "duplicate") bad$.predictions[[1]] <- rbind(bad$.predictions[[1]], bad$.predictions[[1]][1, ])
      if (change == "config") bad$.predictions[[1]]$.config <- "unexpected"
      if (change == "metric") bad$.metrics[[1]] <- bad$.metrics[[1]][-1, ]
      error(check(bad, holdout))
    }
    overlap <- rsample::manual_rset(list(rsample::make_splits(list(analysis = 1:80, assessment = 80:120), d)), "Overlap")
    error(example$check_outer_splits(overlap))

  }

  # Hand calculations for both event conventions and 0/1 extremes. Yardstick
  # clips only when scoring; the supplied probabilities stay exactly 0 and 1.
  truth <- factor(c("yes", "no", "yes", "no"), levels = c("yes", "no"))
  p <- c(0, 1, .25, .75)
  for (event in c("first", "second")) {
    event_p <- if (event == "first") p else 1 - p
    stopifnot(abs(yardstick::brier_class_vec(truth, event_p, event_level = event) - .78125) < 1e-10,
              abs(yardstick::mn_log_loss_vec(truth, event_p, event_level = event) -
                    (-log(.Machine$double.eps) - log(.25)) / 2) < 1e-10,
              identical(p, c(0, 1, .25, .75)))
  }
}

  testthat::expect_true(TRUE)
})

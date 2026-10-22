testthat::test_that("reference-splits", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")
if (all(vapply(c("parsnip", "workflows", "rsample", "recipes", "nnet", "withr", "tune", "yardstick"),
               requireNamespace, logical(1), quietly = TRUE))) {
  example <- new.env()
  sys.source(system.file("examples", "reference-splits.R", package = "vennabers"), example)
  ordered <- example$ordered_holdout
  error <- function(expr, pattern) {
    e <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
  }
  partition <- function(data, held = list(1:20, 21:40, 41:60)) {
    rsample::manual_rset(lapply(held, function(h) rsample::make_splits(
      list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data)),
      paste0("Block", seq_along(held)))
  }
  # Independent literal memberships, with descending values and odd sizes.
  small <- data.frame(x = 11:1)
  for (q in c(.25, .5)) {
    s <- ordered(small, q)
    expected <- if (q == .25) list(1:8, 9:11) else list(1:5, 6:11)
    stopifnot(identical(as.integer(s, data = "analysis"), expected[[1L]]),
              identical(as.integer(s, data = "assessment"), expected[[2L]]),
              identical(rsample::analysis(s)$x, small$x[expected[[1L]]]))
  }
  for (q in list(0, 1, NA_real_, Inf, .99, c(.1, .2))) error(ordered(small, q), "")
  error(ordered(small[1, , drop = FALSE]), "")

  fits <- 0L
  reference_fit <- function(x, y, fail = FALSE) {
    fits <<- fits + 1L
    if (fail) stop("instrumented untrainable base model")
    list(x = x, y = y)
  }
  reference_prob <- function(object, new_data) {
    k <- nlevels(object$y)
    p <- exp(base::outer(new_data$x, seq_len(k) / k))
    p <- sweep(p, 2L, 1 + tabulate(object$y, nbins = k), "*")
    stats::setNames(as.data.frame(p / rowSums(p)), paste0(".pred_", levels(object$y)))
  }
  reference_numeric <- function(object, new_data) mean(object$y) + new_data$x
  parsnip::set_new_model("split_reference")
  for (mode in c("classification", "regression")) {
    register_reference_engine("split_reference", mode, "reference_fit",
      if (mode == "regression") "reference_numeric" else "reference_prob",
      if (mode == "regression") "numeric" else "prob", new_model = FALSE)
  }
  # Public fit/predict/extract seams: observe actual recipe and engine inputs,
  # then compare the core on exactly the outputs of the retained fold models.
  for (kind in c("binary", "ovo", "ovr", "regression")) {
    regression <- kind == "regression"
    mode <- if (regression) "regression" else "classification"
    data <- data.frame(x = c(100, NA, sin(3:72)), id = 72:1)
    data$truth <- if (regression) 10 + cos(1:72) else
      factor(rep(if (kind == "binary") c("b", "a") else c("b", "a", "c"), length.out = 72))
    spec <- parsnip::new_model_spec("split_reference", args = list(), eng_args = NULL,
                                   mode = mode, method = NULL, engine = "reference")
    base <- workflows::workflow() |>
      workflows::add_recipe(recipes::recipe(truth ~ x + id, data) |>
        recipes::step_impute_median(x) |> recipes::step_normalize(x)) |>
      workflows::add_model(spec)
    options <- if (regression) list(m = 1L) else list(classification = kind, precision = 2)
    make_model <- function(method, splits = NULL, ...) do.call(calibrated_workflow,
      c(list(workflow = base, outcome = "truth", method = method, splits = splits, seed = 42), options, list(...)))
    # Local order deliberately differs from the original/global positions.
    development <- data[c(21:60, 1:20), ]
    new <- data[61:72, ]
    for (method in c("inductive", "cross")) {
      resolve <- if (method == "inductive") ordered else partition
      plan <- resolve(development)
      prepared <- local_calibration_splits(list(plan, plan), names(data))
      fit <- generics::fit(make_model(method, prepared), development)
      state <- workflows::extract_fit_engine(fit)
      before <- fits
      assert_conditions(identical(names(state$models), names(state$folds)),
                is.null(state$parameters$effective$calibration_prop),
                is.null(state$parameters$effective$v), is.null(state$parameters$effective$strata),
                state$parameters$explicit[["splits"]])
      calibration <- truth <- inputs <- ids <- list()
      probability <- function(w, d) {
        p <- predict(w, d, type = if (regression) "numeric" else "prob")
        if (regression) return(p$.pred)
        if (kind == "binary") return(p[[paste0(".pred_", state$event_class)]])
        p <- as.matrix(p)
        colnames(p) <- sub("^\\.pred_", "", colnames(p))
        p
      }
      expected_folds <- if (method == "inductive") list(plan) else plan$splits
      for (i in seq_along(state$folds)) {
        id <- names(state$folds)[i]
        f <- state$folds[[i]]
        assert_conditions(identical(f$analysis, as.integer(expected_folds[[i]], data = "analysis")),
                  identical(f$assessment, as.integer(expected_folds[[i]], data = "assessment")),
                  setequal(f$analysis, setdiff(1:60, f$assessment)),
                  f$training_size == length(f$analysis), f$calibration_size == length(f$assessment))
        engine <- workflows::extract_fit_engine(state$models[[id]])
        recipe <- workflows::extract_recipe(state$models[[id]])
        x <- development$x[f$analysis]
        median <- median(x, na.rm = TRUE)
        x[is.na(x)] <- median
        stopifnot(identical(engine$x$id, development$id[f$analysis]),
                  identical(as.vector(engine$y), as.vector(development$truth[f$analysis])),
                  abs(recipe$steps[[1L]]$medians[["x"]] - median) < 1e-10,
                  abs(recipe$steps[[2L]]$means[["x"]] - mean(x)) < 1e-10)
        calibration[[id]] <- probability(state$models[[id]], development[f$assessment, ])
        truth[[id]] <- if (kind == "binary") development$truth[f$assessment] == state$event_class else
          development$truth[f$assessment]
        inputs[[id]] <- probability(state$models[[id]], new)
        ids[[id]] <- rep(id, length(f$assessment))
      }
      if (method == "cross") stopifnot(identical(sort(unlist(lapply(state$folds, `[[`, "assessment"), use.names = FALSE)), 1:60))
      core_args <- list(truth = do.call(c, truth))
      core_args[[if (regression) "prediction" else "probability"]] <-
        if (regression || kind == "binary") unlist(calibration) else do.call(rbind, calibration)
      if (regression) core_args$m <- 1L else {
        core_args$precision <- 2
        if (kind != "binary") core_args$method <- kind
      }
      if (method == "cross") core_args$fold <- unlist(ids)
      core_fun <- if (regression) {
        if (method == "cross") cross_venn_abers_regression else venn_abers_regression
      } else if (method == "cross") cross_venn_abers else venn_abers
      core <- do.call(core_fun, core_args)
      input <- if (method == "inductive") inputs[[1L]] else if (regression || kind == "binary") do.call(cbind, inputs) else inputs
      for (loss in if (regression) "unused" else c("log", "brier")) {
        actual <- if (regression) predict(fit, new) else predict(fit, new, type = "prob", loss = loss)
        expected <- if (regression) predict(core, input)$prediction else {
          p <- predict(core, input, loss = loss)
          if (kind == "binary") p$probability else p
        }
        actual <- if (regression) actual$.pred else if (kind == "binary") actual[[paste0(".pred_", state$event_class)]] else as.matrix(actual)
        tolerance <- if (regression) 1e-10 + 1e-10 * abs(expected) else 1e-10
        stopifnot(all(abs(actual - expected) <= tolerance))
      }
      stopifnot(fits == before) # Prediction uses stored models, never learns.

      # Change only one held fold's truth: its model/recipe must be unchanged;
      # calibration must see the change. Cross sibling models may legitimately change.
      changed <- development
      h <- state$folds[[1L]]$assessment
      changed$truth[h] <- if (regression) changed$truth[h] + 100 else
        factor(rep(levels(data$truth)[1L], length(h)), levels = levels(data$truth))
      if (!regression) changed$truth[tail(h, nlevels(data$truth) - 1L)] <- levels(data$truth)[-1L]
      changed_fit <- generics::fit(make_model(method, resolve), changed)
      changed_state <- workflows::extract_fit_engine(changed_fit)
      stopifnot(identical(workflows::extract_fit_engine(state$models[[1L]]),
                          workflows::extract_fit_engine(changed_state$models[[1L]])),
                identical(workflows::extract_recipe(state$models[[1L]])$steps,
                          workflows::extract_recipe(changed_state$models[[1L]])$steps),
                !identical(state$calibrator, changed_state$calibrator))
      error(generics::fit(make_model(method, prepared), changed), "No prepared plan")
      # Full transfer keeps the matched held-out values and calibrators; only its
      # extra Full1 fit learns the complete local domain, changed held truth included.
      before <- fits
      full <- lapply(list(development, changed), function(d) workflows::extract_fit_engine(
        generics::fit(make_model(method, resolve, prediction_model = "full"), d)))
      engines <- lapply(full, function(s) workflows::extract_fit_engine(s$models$Full1))
      assert_conditions(fits == before + 2L * (length(state$folds) + 1L),
                identical(full[[1L]]$calibration, state$calibration), identical(full[[1L]]$calibrator, state$calibrator),
                identical(full[[2L]]$calibration, changed_state$calibration),
                identical(full[[2L]]$calibrator, changed_state$calibrator),
                identical(names(full[[2L]]$models), "Full1"), all(full[[2L]]$sources == "Full1"),
                identical(engines[[1L]]$x, engines[[2L]]$x), identical(engines[[2L]]$x$id, development$id),
                identical(as.vector(engines[[2L]]$y), as.vector(changed$truth)),
                !identical(engines[[1L]]$y, engines[[2L]]$y))

      # Complete outer evaluation sees the same locally ordered plan. Changes
      # solely in outer assessment truth cannot affect learned state/predictions.
      outer <- function(d) rsample::manual_rset(list(rsample::make_splits(
        list(analysis = c(21:60, 1:20), assessment = 61:72), d)), "Held")
      metrics <- if (regression) yardstick::metric_set(yardstick::rmse) else
        yardstick::metric_set(yardstick::brier_class)
      run <- function(d) tune::fit_resamples(make_model(method, prepared), outer(d), metrics = metrics,
        control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, allow_par = FALSE))
      first <- run(data)
      changed <- data
      changed$truth[61:72] <- if (regression) data$truth[61:72] + 50 else
        factor(rep(levels(data$truth)[1L], 12), levels = levels(data$truth))
      second <- run(changed)
      stopifnot(nrow(tune::collect_notes(first)) == 0L, nrow(tune::collect_notes(second)) == 0L)
      a <- tune::collect_extracts(first)$.extracts[[1L]]
      b <- tune::collect_extracts(second)$.extracts[[1L]]
      assert_conditions(identical(a$folds, state$folds), identical(a$calibrator, state$calibrator),
                identical(a$folds, b$folds), identical(a$calibrator, b$calibrator),
                identical(lapply(a$models, workflows::extract_fit_engine), lapply(b$models, workflows::extract_fit_engine)))
      p <- tune::collect_predictions(first)
      q <- tune::collect_predictions(second)
      columns <- names(p)[startsWith(names(p), ".pred")]
      stopifnot(identical(p[columns], q[columns]), identical(sort(p$.row), 61:72),
                !identical(p$truth, q$truth))

      before <- fits
      if (method == "cross") {
        bad <- plan
        # A forward-style analysis prefix is not the exact third-fold complement.
        bad$splits[[3L]]$in_id <- 1:10
        error(generics::fit(make_model(method, bad), development), "Block3.*partition")
        for (held in list(list(1:20, 20:40, 41:60), list(1:20, 22:40, 41:60),
                          list(1:20, 21:40, 21:40))) {
          error(generics::fit(make_model(method, partition(development, held)), development), "exactly once")
        }
      } else {
        bad <- plan
        bad$in_id <- integer()
        error(generics::fit(make_model(method, bad), development), "Inductive1.*partition")
      }
      if (regression) {
        too_small <- if (method == "inductive") ordered(development, .02) else
          partition(development, list(1:2, 3:30, 31:60))
        error(generics::fit(make_model(method, too_small), development), "calibration size 2")
      } else {
        missing <- development
        missing$truth[] <- levels(data$truth)[1L]
        error(generics::fit(make_model(method, resolve), missing), "all classes")
        if (kind != "binary") {
          absent <- development
          absent$truth[state$folds[[1L]]$assessment] <- levels(data$truth)[1L]
          if (method == "inductive" || kind == "ovo") {
            error(generics::fit(make_model(method, resolve), absent),
                  if (method == "inductive") "class" else "Block1")
          }
        }
      }
      stopifnot(fits == before) # Invalid plans/classes/sizes fail before learning.
      failed_base <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail = TRUE))
      failed_model <- do.call(calibrated_workflow, c(list(workflow = failed_base,
        outcome = "truth", method = method, splits = plan), options))
      error(generics::fit(failed_model, development), paste0(names(state$folds)[1L], " base workflow training"))
      stopifnot(fits == before + 1L) # No retry, reshuffle or ensemble reduction.

      # Mixed automatic plans remain reproducible for ordinary/stratified paths.
      for (strata in if (regression) list(FALSE) else list(NULL, FALSE)) {
        one <- workflows::extract_fit_engine(generics::fit(make_model(method, strata = strata), development))
        two <- workflows::extract_fit_engine(generics::fit(make_model(method, strata = strata), development))
        stopifnot(identical(one$folds, two$folds), identical(one$calibrator, two$calibrator))
      }
    }
    if (kind %in% c("ovo", "ovr")) {
      strat_data <- data[1:11, ]
      strat_data$truth <- example$small$truth
      strat_plan <- partition(strat_data, example$stratified_held)
      strat_fit <- generics::fit(make_model("cross", strat_plan), strat_data)
      f <- workflows::extract_fit_engine(strat_fit)$folds
      stopifnot(identical(unname(lapply(f, `[[`, "assessment")), example$stratified_held),
                identical(unname(vapply(f, `[[`, 0L, "calibration_size")), c(4L, 4L, 3L)))
    }
  }

  # Restore real workflows AND prepared resolvers in a fresh R process for all
  # eight paths; current-process success alone cannot prove namespace restoration.
  payload <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  output <- tempfile(fileext = ".log")
  saveRDS(example$reference_results, payload)
  writeLines(c(paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    'library(vennabers)', 'library(workflows)',
    paste0('saved <- readRDS(', deparse(payload), ')'),
    'for (entry in saved) {',
    '  type <- if (is.factor(entry$development$truth)) "prob" else "numeric"',
    '  stopifnot(identical(predict(entry$fit, entry$assessment, type = type), entry$prediction))',
    '  refit <- generics::fit(entry$model, entry$development)',
    '  stopifnot(identical(workflows::extract_fit_engine(refit)$folds, entry$folds),',
    '            identical(predict(refit, entry$assessment, type = type), entry$prediction))',
    '}'), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = output, stderr = output)
  if (status != 0L) stop(paste(readLines(output), collapse = "\n"))
  unlink(c(payload, script, output))
}

  testthat::expect_true(TRUE)
})

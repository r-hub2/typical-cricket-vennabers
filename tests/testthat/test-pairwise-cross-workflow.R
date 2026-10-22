testthat::test_that("pairwise-cross-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")

packages <- c("parsnip", "workflows", "recipes", "rsample", "withr")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  data <- data.frame(x = sin(1:90), row_id = 1:90,
    truth = factor(rep(c("z", "a", "b"), c(24, 30, 36)), levels = c("z", "a", "b")))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  fit <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", method = "cross", v = 5, seed = 42), data)
  state <- workflows::extract_fit_engine(fit)
  stopifnot(sum(vapply(state$pair_models, function(p) length(p$models), 1L)) == 15L,
    nrow(venn_abers_details(fit, data[1:7, ])$details) == 105L)
}

if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  close <- function(actual, expected) stopifnot(identical(dim(actual), dim(expected)),
    all(is.finite(actual)), all(abs(actual - expected) <= 1e-10))
  error <- function(expr, text = "") {
    e <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(e, "error"), grepl(text, conditionMessage(e), fixed = TRUE))
  }
  list2env(pairwise_reference("va_pair_cross_reference"), environment())
  # Ordinary non-shuffled folds: one-sided assessment sets are intentional.
  partition <- function(data, v = 3L) {
    held <- split(seq_len(nrow(data)), cut(seq_len(nrow(data)), v, labels = FALSE))
    rsample::manual_rset(lapply(held, function(h) rsample::make_splits(
      list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data)),
      paste0("Fold", seq_len(v)))
  }
  for (k in 2:4) {
    classes <- c("z", "a b", "a.b", "4")[seq_len(k)]
    n <- sum(18L + 6L * seq_len(k))
    data <- data.frame(x = c(NA, sin(2:(n - 6L)), 1000 + 1:6), row_id = seq_len(n),
      truth = factor(rep(classes, 18L + 6L * seq_len(k)), levels = classes))
    base <- make_base(data)
    locals <- pair_data(data)
    plans <- setNames(lapply(seq_along(locals), function(i) partition(locals[[i]], 3L + i %% 3L)), names(locals))
    cost <- sum(vapply(plans, nrow, 1L))
    new <- data[c(1:7, n), ]
    for (precision in list(NULL, 0, 2, -1)) {
      before <- length(audit$fits)
      calls <- length(audit$predictions)
      model <- calibrated_workflow(base, "truth", classification = "ovo", method = "cross",
        model_training = "pairwise", splits = plans[rev(names(plans))], seed = 51, precision = precision)
      fit <- generics::fit(model, data)
      state <- workflows::extract_fit_engine(fit)
      assert_conditions(length(audit$fits) == before + cost, length(audit$predictions) == calls + cost,
        identical(state$class_levels, classes),
        identical(state$rng_kind, RNGkind()), identical(state$parameters$effective$precision, precision))
      for (j in seq_along(locals)) {
        component <- state$pair_models[[j]]
        stopifnot(identical(names(component$models), plans[[j]]$id),
          identical(component$calibrator$folds, plans[[j]]$id),
          identical(data$row_id[component$rows], locals[[j]]$row_id))
        for (h in seq_len(nrow(plans[[j]]))) {
          fold <- component$folds[[h]]
          engine <- workflows::extract_fit_engine(component$models[[h]])
          training <- rsample::analysis(plans[[j]]$splits[[h]])
          held <- rsample::assessment(plans[[j]]$splits[[h]])
          assert_conditions(identical(engine$x$row_id, training$row_id), identical(engine$y, training$truth),
            identical(engine$levels, component$class_levels), !any(engine$x$row_id %in% held$row_id),
            identical(fold$analysis, as.integer(plans[[j]]$splits[[h]], data = "analysis")),
            identical(fold$assessment, as.integer(plans[[j]]$splits[[h]], data = "assessment")),
            fold$training_size == nrow(training), fold$calibration_size == nrow(held))
          x <- training$x
          x[is.na(x)] <- median(x, na.rm = TRUE)
          close(engine$x$x, (x - mean(x)) / sd(x))
        }
      }
      stopifnot(any(vapply(state$pair_models, function(p)
        any(vapply(p$calibration, function(f) length(unique(f$truth)) == 1L, logical(1))), logical(1))))
      for (loss in c("log", "brier")) {
        calls <- length(audit$predictions)
        before <- length(audit$fits)
        details <- venn_abers_details(fit, new, loss = loss)
        assert_conditions(length(audit$fits) == before, length(audit$predictions) == calls + cost,
          all(details$details$model_training == "pairwise"), all(details$details$prediction_model == "matched"),
          identical(details$details$model_source, details$details$pair_id),
          all(details$details$event_class == details$details$class_2))
        q <- matrix(0, nrow(new), length(locals))
        for (j in seq_along(locals)) {
          local <- locals[[j]]
          event <- levels(local$truth)[2L]
          column <- paste0(".pred_", event)
          lower <- upper <- point <- matrix(0, nrow(new), nrow(plans[[j]]))
          for (h in seq_len(nrow(plans[[j]]))) {
            split <- plans[[j]]$splits[[h]]
            binary <- generics::fit(make_base(local), rsample::analysis(split))
            held <- rsample::assessment(split)
            pcal <- predict(binary, held, type = "prob")[[column]]
            cal <- venn_abers(held$truth == event, pcal, precision = precision)
            expected <- predict(cal, predict(binary, new, type = "prob")[[column]], loss = loss)
            part <- details$details[details$details$pair_id == names(locals)[j] &
              details$details$fold == plans[[j]]$id[h], ]
            close(as.matrix(part[c("lower", "upper", "probability")]), as.matrix(expected))
            close(state$pair_models[[j]]$calibration[[h]]$probability, pcal)
            lower[, h] <- expected$lower
            upper[, h] <- expected$upper
            point[, h] <- expected$probability
          }
          # Independent geometric-products reference, without the Cross merger.
          q[, j] <- if (loss == "brier") rowMeans(point) else {
            a <- apply(upper, 1, prod)^(1 / ncol(upper))
            b <- apply(1 - lower, 1, prod)^(1 / ncol(lower))
            a / (a + b)
          }
        }
        close(as.matrix(details$prediction), couple(q, classes))
        for (n in c(0L, 1L, nrow(new))) {
          input <- new[seq_len(n), , drop = FALSE]
          before <- length(audit$fits)
          calls <- length(audit$predictions)
          result <- predict(fit, input, type = "prob", loss = loss)
          cl <- predict(fit, input, type = "class", loss = loss)
          d <- venn_abers_details(fit, input, loss = loss)
          assert_conditions(identical(result, details$prediction[seq_len(n), ]), identical(result, d$prediction),
            identical(names(d$details), names(details$details)), nrow(d$details) == n * cost,
            identical(cl$.pred_class, factor(classes[max.col(as.matrix(result), ties.method = "first")], levels = classes)),
            length(audit$fits) == before, length(audit$predictions) == calls + if (n) 3L * cost else 0L)
        }
      }
    }
    for (strata in list(NULL, FALSE)) {
      model <- calibrated_workflow(base, "truth", classification = "ovo", method = "cross",
        model_training = "pairwise", v = 5, strata = strata, seed = 42)
      set.seed(23)
      rng <- .Random.seed
      before <- length(audit$fits)
      automatic <- generics::fit(model, data)
      stopifnot(identical(.Random.seed, rng), length(audit$fits) == before + 5L * length(locals))
      a <- workflows::extract_fit_engine(automatic)$pair_models
      b <- workflows::extract_fit_engine(generics::fit(model, data))$pair_models
      for (id in names(a)) {
        stopifnot(identical(a[[id]]$folds, b[[id]]$folds), identical(a[[id]]$calibration, b[[id]]$calibration),
          identical(a[[id]]$calibrator, b[[id]]$calibrator))
        for (h in names(a[[id]]$models)) stopifnot(identical(workflows::extract_fit_engine(a[[id]]$models[[h]]),
          workflows::extract_fit_engine(b[[id]]$models[[h]])))
      }
      prepared <- lapply(locals, function(local) withr::with_seed(42,
        if (is.null(strata)) rsample::vfold_cv(local, v = 5, strata = truth, pool = 0) else rsample::vfold_cv(local, v = 5)))
      for (j in seq_along(prepared)) for (h in seq_len(5L)) stopifnot(
        identical(a[[j]]$folds[[h]]$analysis, as.integer(prepared[[j]]$splits[[h]], data = "analysis")))
      resolver <- local_calibration_splits(prepared, names(data))
      explicit <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo", method = "cross",
        model_training = "pairwise", splits = resolver, seed = 42), data)
      stopifnot(identical(venn_abers_details(automatic, new), venn_abers_details(explicit, new)))
    }
  # Reject malformed later pairs before *any* training.
  before <- length(audit$fits)
  run <- function(plans) generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", method = "cross", splits = plans), data)
  malformed <- if (length(plans) > 1L) list(plans[-1L], c(plans, plans[1L]),
    setNames(plans, rev(names(plans))), setNames(plans, rep("pair_1", length(plans))), unname(plans)) else list()
  for (bad in malformed) error(run(bad))
  for (change in c("values", "type", "levels", "order", "missing", "overlap", "duplicate", "fraction",
                   "empty", "outside", "na", "fold_id", "missing_id", "repeated", "forward", "coverage")) {
    bad <- plans
    j <- length(bad)
    f <- bad[[j]]$splits[[1L]]
    if (change == "values") f$data$x[1L] <- 456
    if (change == "type") f$data$row_id <- as.double(f$data$row_id)
    if (change == "levels") f$data$truth <- factor(f$data$truth, levels = rev(levels(f$data$truth)))
    if (change == "order") f$data <- f$data[rev(seq_len(nrow(f$data))), ]
    if (change %in% c("missing", "forward")) f$in_id <- f$in_id[-1L]
    if (change == "overlap") f$in_id <- c(f$in_id, f$out_id[1L])
    if (change == "duplicate") f$in_id <- c(f$in_id, f$in_id[1L])
    if (change == "fraction") f$in_id[1L] <- 1.5
    if (change == "empty") f$out_id <- integer()
    if (change == "outside") f$in_id[1L] <- nrow(f$data) + 1L
    if (change == "na") f$in_id[1L] <- NA_integer_
    bad[[j]]$splits[[1L]] <- f
    if (change == "fold_id") bad[[j]]$id[1L] <- bad[[j]]$id[2L]
    if (change == "missing_id") bad[[j]]$id[1L] <- NA_character_
    if (change == "repeated") bad[[j]] <- rsample::vfold_cv(locals[[j]], v = 3, repeats = 2)
    if (change == "coverage") bad[[j]]$splits[[1L]] <- bad[[j]]$splits[[2L]]
    error(run(bad), paste0("pair_", length(plans)))
  }
  error(generics::fit(calibrated_workflow(base, "truth", classification = "ovo", method = "cross",
    model_training = "pairwise", v = nrow(data) + 1L), data), "not reduced")
  error(calibrated_workflow(base, "truth", classification = "ovo", method = "cross",
    model_training = "pairwise", splits = plans, v = 3), "Explicit")
  stopifnot(length(audit$fits) == before)
  conflict <- prepared[[1L]]
  conflict$splits[[1L]]$in_id <- rev(conflict$splits[[1L]]$in_id)
  error(local_calibration_splits(c(prepared, list(conflict)), names(data)), "Conflicting")
  if (length(plans) > 1L) for (change in c("pair", "classes", "pair_model", "fold_names")) {
    bad <- fit
    s <- bad$fit$fit$fit
    if (change == "pair") s$pair_models[c(1, 2)] <- s$pair_models[c(2, 1)]
    if (change == "classes") s$pair_models[[1L]]$class_levels <- s$pair_models[[2L]]$class_levels
    if (change == "pair_model") s$pair_models[[1L]]$models[[1L]] <- s$pair_models[[2L]]$models[[1L]]
    if (change == "fold_names") names(s$pair_models[[1L]]$models)[1L] <- "wrong"
    bad$fit$fit$fit <- s
    error(predict(bad, new, type = "prob"), "mapping")
  }
  for (bad in list(list(fail_size = length(plans[[1L]]$splits[[1L]]$in_id)), list(bad = TRUE))) {
    broken <- workflows::update_model(base, do.call(parsnip::set_engine, c(list(spec, "reference"), bad)))
    error(generics::fit(calibrated_workflow(broken, "truth", classification = "ovo", method = "cross",
      model_training = "pairwise", splits = plans), data), "pair_")
  }
  # An untrainable complement fails; folds are never silently removed.
  bad <- plans
  local <- locals[[1L]]
  h <- which(local$truth == levels(local$truth)[1L])
  a <- setdiff(seq_len(nrow(local)), h)
  bad[[1L]] <- rsample::manual_rset(list(rsample::make_splits(list(analysis = a, assessment = h), local),
    rsample::make_splits(list(analysis = h, assessment = a), local)), c("A", "B"))
  error(run(bad), "pair_1")
  tied <- data.frame(x = rep(1, 24), row_id = rep(1L, 24), truth = factor(rep(c("z", "a"), 12), levels = c("z", "a")))
  tied_base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |> workflows::add_model(spec)
  tied_fit <- generics::fit(calibrated_workflow(tied_base, "truth", classification = "ovo", method = "cross",
    model_training = "pairwise", splits = partition), tied)
  for (loss in c("log", "brier")) stopifnot(all(predict(tied_fit, tied, type = "prob", loss = loss) == .5),
    all(predict(tied_fit, tied, type = "class", loss = loss)$.pred_class == "z"))
  }
}

if (exists("pair_fit") && all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
  for (k in 2:4) {
    classes <- c("z", "a b", "a.b", "4")[seq_len(k)]
    n <- sum(18L + 6L * seq_len(k))
    data <- data.frame(x = c(NA, sin(2:(n - 6L)), 1000 + 1:6), row_id = seq_len(n),
      truth = factor(rep(classes, 18L + 6L * seq_len(k)), levels = classes))
    base <- make_base(data)
    locals <- pair_data(data)
    plans <- setNames(lapply(seq_along(locals), function(i) partition(locals[[i]], 3L + i %% 3L)), names(locals))
    pair_count <- length(locals)
    fold_count <- nrow(partition(locals[[1L]]))
  checker <- new.env()
  sys.source(system.file("examples", "calibrated-resampling.R", package = "vennabers"), checker)
  check <- checker$check_evaluation
  scoring <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
  control <- tune::control_resamples(save_pred = TRUE, extract = identity, allow_par = FALSE)
  # Interleave classes so every outer domain retains all four classes.
  data <- data[order(ave(data$row_id, data$truth, FUN = seq_along)), ]
  outer <- withr::with_seed(201, rsample::vfold_cv(data, v = 2, repeats = 2, strata = truth))
  h <- seq.int(5L, nrow(data), 5L)
  holdout <- rsample::make_splits(list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data)
  for (loss in c("log", "brier")) {
    model <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
      method = "cross", splits = partition, seed = 42, loss = loss)
    before <- length(audit$fits)
    result <- tune::fit_resamples(model, outer, metrics = scoring, control = control)
    stopifnot(length(audit$fits) == before + nrow(outer) * fold_count * pair_count)
    checked <- check(result, outer, classification = "ovo")
    fits <- tune::collect_extracts(result)$.extracts
    prepared <- list()
    for (i in seq_along(fits)) {
      local <- rsample::analysis(outer$splits[[i]])
      held <- rsample::assessment(outer$splits[[i]])
      locals <- pair_data(local)
      state <- workflows::extract_fit_engine(fits[[i]])
      rows <- checked$predictions$id == outer$id[i] & checked$predictions$id2 == outer$id2[i]
      direct <- predict(fits[[i]], held, type = "prob")
      close(as.matrix(checked$predictions[rows, names(direct)]), as.matrix(direct))
      for (id in names(locals)) {
        s <- state$pair_models[[id]]
        for (fold_id in names(s$folds)) {
          engine <- workflows::extract_fit_engine(s$models[[fold_id]])
          stopifnot(identical(engine$x$row_id, locals[[id]]$row_id[s$folds[[fold_id]]$analysis]),
            !any(engine$x$row_id %in% held$row_id), identical(local$row_id[s$rows], locals[[id]]$row_id),
            !any(local$row_id[s$rows][s$folds[[fold_id]]$assessment] %in% held$row_id))
        }
        prepared[[length(prepared) + 1L]] <- rsample::manual_rset(lapply(s$folds, function(f)
          rsample::make_splits(f[c("analysis", "assessment")], locals[[id]])), names(s$folds))
      }
    }
    explicit <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
      method = "cross", splits = local_calibration_splits(prepared, names(data)), seed = 42, loss = loss)
    replay <- tune::fit_resamples(explicit, outer, metrics = scoring, control = control)
    check(replay, outer, classification = "ovo")
    stopifnot(identical(tune::collect_predictions(result), tune::collect_predictions(replay)),
      identical(tune::collect_metrics(result), tune::collect_metrics(replay)))
    # Assessment truth mutations cannot change any learned pair/fold state.
    for (resampling in c(FALSE, TRUE)) {
      split <- if (resampling) outer$splits[[1L]] else holdout
      changed <- data
      changed$truth[as.integer(split, data = "assessment")] <- classes[1L]
      altered <- rsample::make_splits(list(analysis = as.integer(split, data = "analysis"),
        assessment = as.integer(split, data = "assessment")), changed)
      evaluate <- function(split) if (resampling) tune::fit_resamples(model,
        rsample::manual_rset(list(split), "Outer1"), metrics = scoring, control = control) else
        tune::last_fit(model, split, metrics = scoring)
      before <- length(audit$fits)
      final <- evaluate(split)
      mutated <- evaluate(altered)
      stopifnot(length(audit$fits) == before + 2L * fold_count * pair_count, nrow(tune::collect_notes(final)) == 0L,
        nrow(tune::collect_notes(mutated)) == 0L)
      p <- tune::collect_predictions(final)
      m <- tune::collect_metrics(final)
      estimate <- if (resampling) "mean" else ".estimate"
      assert_conditions(identical(p$.row, as.integer(split, data = "assessment")), !anyDuplicated(p$.row),
        all(p$.config == "pre0_mod0_post0"), identical(p$truth, data$truth[p$.row]),
        nrow(m) == 2L, estimate %in% names(m), all(is.finite(m[[estimate]])))
      one <- if (resampling) tune::collect_extracts(final)$.extracts[[1L]] else tune::extract_workflow(final)
      two <- if (resampling) tune::collect_extracts(mutated)$.extracts[[1L]] else tune::extract_workflow(mutated)
      a <- workflows::extract_fit_engine(one)$pair_models
      b <- workflows::extract_fit_engine(two)$pair_models
      for (id in names(a)) {
        stopifnot(identical(a[[id]]$calibration, b[[id]]$calibration),
          identical(a[[id]]$calibrator, b[[id]]$calibrator), identical(a[[id]]$folds, b[[id]]$folds))
        for (fold in names(a[[id]]$models)) stopifnot(
          identical(workflows::extract_fit_engine(a[[id]]$models[[fold]]),
                    workflows::extract_fit_engine(b[[id]]$models[[fold]])))
      }
      stopifnot(identical(venn_abers_details(one, data), venn_abers_details(two, data)))
    }
    if (loss == "log") {
      for (change in c("missing", "duplicate", "config", "metric", "repeat", "truth")) {
        bad <- result
        if (change == "missing") bad$.predictions[[1L]] <- bad$.predictions[[1L]][-1L, ]
        if (change == "duplicate") bad$.predictions[[1L]] <- rbind(bad$.predictions[[1L]], bad$.predictions[[1L]][1L, ])
        if (change == "config") bad$.predictions[[1L]]$.config <- "wrong"
        if (change == "metric") bad$.metrics[[1L]] <- bad$.metrics[[1L]][-1L, ]
        if (change == "repeat") bad$id2[1L] <- "wrong"
        if (change == "truth") bad$.predictions[[1L]]$truth[1L] <- NA
        error(check(bad, outer, classification = "ovo"))
      }
      # One failed inner fold invalidates its complete outer result.
      first_size <- nrow(rsample::analysis(partition(pair_data(rsample::analysis(outer$splits[[1L]]))[[1L]])$splits[[1L]]))
      broken <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail_size = first_size))
      bad_model <- calibrated_workflow(broken, "truth", classification = "ovo", method = "cross",
        model_training = "pairwise", splits = partition)
      bad <- suppressWarnings(tune::fit_resamples(bad_model, outer, metrics = scoring, control = control))
      stopifnot(nrow(tune::collect_notes(bad)) > 0L)
      error(check(bad, outer, classification = "ovo"))
      bad <- suppressWarnings(tune::last_fit(bad_model, outer$splits[[1L]], metrics = scoring))
      stopifnot(nrow(tune::collect_notes(bad)) > 0L, is.null(bad[[".workflow"]]))
      # Development-wide pair lists cannot bind to smaller outer domains.
      bad_model <- calibrated_workflow(base, "truth", classification = "ovo", method = "cross",
        model_training = "pairwise", splits = plans)
      bad <- suppressWarnings(tune::fit_resamples(bad_model, outer, metrics = scoring, control = control))
      stopifnot(nrow(tune::collect_notes(bad)) == nrow(outer))
      error(check(bad, outer, classification = "ovo"))
    }
  }
  }
}

  testthat::expect_true(TRUE)
})

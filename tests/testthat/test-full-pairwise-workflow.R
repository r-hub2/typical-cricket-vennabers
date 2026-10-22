testthat::test_that("full-pairwise-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")
packages <- c("parsnip", "workflows", "recipes", "rsample", "withr")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  data <- data.frame(x = sin(1:90), truth = factor(rep(c("z", "a", "b"), 30)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (method in c("inductive", "cross")) {
    fits <- lapply(c("matched", "full"), function(mode) generics::fit(calibrated_workflow(
      base, "truth", classification = "ovo", model_training = "pairwise",
      method = method, prediction_model = mode, seed = 42), data))
    states <- lapply(fits, workflows::extract_fit_engine)
    for (id in names(states[[1L]]$pair_models)) {
      source <- states[[1L]]$pair_models[[id]]
      full <- states[[2L]]$pair_models[[id]]
      stopifnot(identical(source$calibration, full$calibration),
        identical(source$calibrator, full$calibrator), identical(names(full$models), "Full1"))
    }
    stopifnot(all(venn_abers_details(fits[[2L]], data)$details$prediction_model == "full"))
  }
}

if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  close <- function(a, b) stopifnot(identical(dim(a), dim(b)), all(is.finite(a)), all(abs(a - b) <= 1e-10))
  error <- function(expr, text = "") {
    e <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(e, "error"), grepl(text, conditionMessage(e), fixed = TRUE))
  }
  audit <- new.env()
  audit$fits <- list()
  audit$predict <- 0L
  full_pair_fit <- function(x, y, fail_size = integer(), bad_size = integer(), random = FALSE) {
    audit$fits[[length(audit$fits) + 1L]] <- list(x = x, y = y)
    if (nrow(x) %in% fail_size) stop("planned training failure")
    list(x = x, y = y, levels = levels(y), bad_size = bad_size, draw = runif(1), random = random)
  }
  full_pair_prob <- function(object, new_data) {
    stopifnot(nrow(new_data) > 0L)
    audit$predict <- audit$predict + 1L
    p <- plogis(new_data$p + nrow(object$x) / 10 - 3 + if (object$random) object$draw else 0)
    if (nrow(object$x) %in% object$bad_size) p[] <- NA_real_
    # Reverse columns to verify named event selection, not positional selection.
    setNames(data.frame(p, 1 - p), paste0(".pred_", rev(object$levels)))
  }
  register_reference_engine("va_full_pair", "classification", "full_pair_fit", "full_pair_prob", "prob")
  spec <- parsnip::new_model_spec("va_full_pair", args = list(), eng_args = NULL,
    mode = "classification", method = NULL, engine = "reference")
  make_base <- function(data, model = spec) workflows::workflow() |> workflows::add_recipe(
    recipes::recipe(truth ~ ., data) |> recipes::step_impute_median(x) |> recipes::step_normalize(x)) |>
    workflows::add_model(model)
  pair_data <- function(data) setNames(lapply(combn(levels(data$truth), 2, simplify = FALSE), function(pair) {
    d <- data[data$truth %in% pair, , drop = FALSE]
    d$truth <- factor(d$truth, levels = pair, ordered = is.ordered(data$truth))
    d
  }), paste0("pair_", seq_len(choose(nlevels(data$truth), 2))))
  plan <- function(data, method, v = 3L) {
    held <- if (method == "inductive") list(seq.int(nrow(data) - 7L, nrow(data))) else
      split(seq_len(nrow(data)), cut(seq_len(nrow(data)), v, labels = FALSE))
    splits <- lapply(held, function(h) rsample::make_splits(
      list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data))
    if (method == "inductive") splits[[1L]] else rsample::manual_rset(splits, paste0("Fold", seq_along(splits)))
  }
  couple <- function(q, classes) {
    pairs <- combn(classes, 2, simplify = FALSE)
    p <- vapply(classes, function(cl) {
      relevant <- which(vapply(pairs, function(pair) cl %in% pair, logical(1)))
      odds <- lapply(relevant, function(j) {
        probability <- if (cl == pairs[[j]][2L]) q[, j] else 1 - q[, j]
        (1 - probability) / probability
      })
      1 / (1 + Reduce(`+`, odds))
    }, numeric(nrow(q)))
    p <- matrix(p, nrow(q), length(classes))
    p / rowSums(p)
  }
  for (k in 2:4) for (method in c("inductive", "cross")) {
    classes <- c("z", "a b", "a.b", "4")[seq_len(k)]
    n <- 24L * k
    data <- data.frame(p = seq(-8, 0, length.out = n), row_id = seq_len(n),
      x = c(NA, sin(2:(n - 8L)), 10000 + 1:8),
      truth = factor(rep(classes, 24), levels = classes, ordered = k == 2L))
    base <- make_base(data)
    locals <- pair_data(data)
    plans <- setNames(lapply(seq_along(locals), function(j) plan(locals[[j]], method, 2L + j %% 4L)), names(locals))
    cost <- if (method == "inductive") length(plans) else sum(vapply(plans, nrow, 1L))
    new <- data[c(1:7, n), ]
    new$p <- seq(-3, 1, length.out = nrow(new))
    for (precision in list(NULL, 0, 2, -1)) {
      fits <- states <- list()
      for (mode in c("matched", "full")) {
        before <- length(audit$fits)
        set.seed(77)
        rng <- .Random.seed
        fits[[mode]] <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
          model_training = "pairwise", prediction_model = mode, method = method,
          splits = plans[rev(names(plans))], precision = precision, seed = 51), data)
        states[[mode]] <- workflows::extract_fit_engine(fits[[mode]])
        stopifnot(identical(rng, .Random.seed), length(audit$fits) == before + cost + (mode == "full") * length(locals))
      }
      state <- states$full
      stopifnot(state$parameters$effective$prediction_model == "full",
        state$parameters$explicit[["prediction_model"]], identical(state$rng_kind, RNGkind()),
        all(c("R", "vennabers", "recipes", "workflows") %in% names(state$versions)))
      for (id in names(locals)) {
        a <- states$matched$pair_models[[id]]
        b <- state$pair_models[[id]]
        assert_conditions(identical(a$calibrator, b$calibrator), identical(a$calibration, b$calibration),
          identical(a$folds, b$folds), identical(a$rows, b$rows), identical(names(b$models), "Full1"))
        engine <- workflows::extract_fit_engine(b$models$Full1)
        stopifnot(identical(engine$x$row_id, locals[[id]]$row_id), identical(engine$y, locals[[id]]$truth),
          identical(levels(engine$y), b$class_levels), identical(data$row_id[b$rows], locals[[id]]$row_id))
        x <- locals[[id]]$x
        x[is.na(x)] <- median(x, na.rm = TRUE)
        close(workflows::extract_recipe(b$models$Full1)$steps[[2L]]$means[["x"]], mean(x))
        for (fold in names(a$folds)) {
          f <- a$folds[[fold]]
          source <- workflows::extract_fit_engine(a$models[[fold]])
          stopifnot(identical(source$x$row_id, locals[[id]]$row_id[f$analysis]),
            !any(source$x$row_id %in% locals[[id]]$row_id[f$assessment]),
            identical(b$calibration[[fold]]$row, f$assessment),
            identical(b$calibration[[fold]]$truth, locals[[id]]$truth[f$assessment] == b$event_class))
          close(b$calibration[[fold]]$probability, plogis(locals[[id]]$p[f$assessment] + length(f$analysis) / 10 - 3))
          x <- locals[[id]]$x[f$analysis]
          x[is.na(x)] <- median(x, na.rm = TRUE)
          close(workflows::extract_recipe(a$models[[fold]])$steps[[2L]]$means[["x"]], mean(x))
        }
        if (method == "inductive") stopifnot(abs(workflows::extract_recipe(b$models$Full1)$steps[[2L]]$means[["x"]] -
          workflows::extract_recipe(a$models$Inductive1)$steps[[2L]]$means[["x"]]) > 100)
      }
      for (loss in c("log", "brier")) {
        before <- length(audit$fits)
        calls <- audit$predict
        d <- venn_abers_details(fits$full, new, loss = loss)
        assert_conditions(audit$predict == calls + length(locals), all(d$details$prediction_model == "full"),
          all(d$details$model_training == "pairwise"), all(d$details$event_class == d$details$class_2),
          identical(d$details$model_source, paste0(d$details$pair_id, "/Full1")))
        q <- matrix(0, nrow(new), length(locals))
        for (j in seq_along(locals)) {
          component <- state$pair_models[[j]]
          p <- plogis(new$p + nrow(locals[[j]]) / 10 - 3)
          lower <- upper <- point <- matrix(0, nrow(new), length(component$folds))
          for (h in seq_along(component$folds)) {
            cal <- component$calibration[[h]]
            expected <- predict(venn_abers(cal$truth, cal$probability, precision = precision), p, loss = loss)
            part <- d$details[d$details$pair_id == names(locals)[j] & d$details$fold == names(component$folds)[h], ]
            close(as.matrix(part[c("lower", "upper", "probability")]), as.matrix(expected))
            lower[, h] <- expected$lower
            upper[, h] <- expected$upper
            point[, h] <- expected$probability
          }
          q[, j] <- if (loss == "brier") rowMeans(point) else {
            a <- apply(upper, 1, prod)^(1 / ncol(upper))
            b <- apply(1 - lower, 1, prod)^(1 / ncol(lower))
            a / (a + b)
          }
        }
        close(as.matrix(d$prediction), couple(q, classes))
        for (size in c(0L, 1L, nrow(new))) {
          x <- new[seq_len(size), , drop = FALSE]
          calls <- audit$predict
          prob <- predict(fits$full, x, type = "prob", loss = loss)
          cls <- predict(fits$full, x, type = "class", loss = loss)
          small <- venn_abers_details(fits$full, x, loss = loss)
          assert_conditions(audit$predict == calls + if (size) 3L * length(locals) else 0L,
            identical(small$prediction, prob), identical(prob, d$prediction[seq_len(size), ]),
            identical(names(small$details), names(d$details)), nrow(small$details) == size * cost,
            identical(levels(cls$.pred_class), classes),
            identical(as.character(cls$.pred_class), classes[max.col(as.matrix(prob), ties.method = "first")]))
        }
        singles <- do.call(rbind, lapply(seq_len(nrow(new)), function(i) predict(fits$full, new[i, ], type = "prob", loss = loss)))
        stopifnot(identical(singles, d$prediction), length(audit$fits) == before)
        if (is.null(precision)) stopifnot(max(abs(as.matrix(d$prediction) -
          as.matrix(predict(fits$matched, new, type = "prob", loss = loss)))) > 1e-5)
      }
    }
    # Foreign-class rows never reach a pair. Held pair truth reaches its calibrator
    # and Full1, but not the retained source predictions used as calibration values.
    fit_full <- function(d) workflows::extract_fit_engine(generics::fit(calibrated_workflow(base, "truth",
      classification = "ovo", model_training = "pairwise", prediction_model = "full", method = method,
      splits = function(local) plan(local, method), seed = 51), d))
    kept <- fit_full(data)
    if (k > 2L) {
      foreign <- data
      rows <- data$truth == classes[k]
      foreign$x[rows] <- foreign$x[rows] + 50
      foreign$p[rows] <- -foreign$p[rows]
      other <- fit_full(foreign)
      for (j in seq_along(locals)) {
        a <- kept$pair_models[[j]]
        b <- other$pair_models[[j]]
        same <- !classes[k] %in% a$class_levels
        assert_conditions(identical(a$rows, b$rows), identical(a$folds, b$folds),
          identical(a$calibration, b$calibration) == same, identical(a$calibrator, b$calibrator) == same,
          identical(workflows::extract_fit_engine(a$models$Full1), workflows::extract_fit_engine(b$models$Full1)) == same)
      }
    }
    a <- kept$pair_models$pair_1
    rows <- a$rows[a$folds[[1L]]$assessment]
    moved <- data
    moved$truth[rows] <- rev(moved$truth[rows])
    b <- fit_full(moved)$pair_models$pair_1
    source <- workflows::extract_fit_engine(a$models$Full1)
    engine <- workflows::extract_fit_engine(b$models$Full1)
    assert_conditions(!identical(moved$truth, data$truth), identical(a$rows, b$rows), identical(a$folds, b$folds),
      identical(a$calibration[[1L]]$probability, b$calibration[[1L]]$probability),
      !identical(a$calibration[[1L]]$truth, b$calibration[[1L]]$truth), !identical(a$calibrator, b$calibrator),
      identical(engine$x, source$x), identical(as.character(engine$y), as.character(moved$truth[a$rows])),
      !identical(engine$y, source$y))
    # Automatic plans and source predictions also stay unchanged for random engines.
    random_base <- make_base(data, parsnip::set_engine(spec, "reference", random = TRUE))
    for (strata in list(NULL, FALSE)) {
      before <- length(audit$fits)
      auto <- lapply(c("matched", "full"), function(mode) generics::fit(calibrated_workflow(
        random_base, "truth", classification = "ovo", model_training = "pairwise", method = method,
        prediction_model = mode, strata = strata, seed = 42), data))
      stopifnot(length(audit$fits) == before + (if (method == "cross") 11L else 3L) * length(locals))
      a <- workflows::extract_fit_engine(auto[[1L]])$pair_models
      b <- workflows::extract_fit_engine(auto[[2L]])$pair_models
      for (id in names(a)) stopifnot(identical(a[[id]]$folds, b[[id]]$folds),
        identical(a[[id]]$calibration, b[[id]]$calibration), identical(a[[id]]$calibrator, b[[id]]$calibrator))
      # Rebuild the exact prepared local plans and replay them in full mode.
      prepared <- lapply(names(b), function(id) {
        fs <- lapply(b[[id]]$folds, function(f) rsample::make_splits(f[c("analysis", "assessment")], locals[[id]]))
        if (method == "inductive") fs[[1L]] else rsample::manual_rset(fs, names(b[[id]]$folds))
      })
      replay <- generics::fit(calibrated_workflow(random_base, "truth", classification = "ovo",
        model_training = "pairwise", method = method, prediction_model = "full", seed = 42,
        splits = local_calibration_splits(prepared, names(data))), data)
      stopifnot(identical(venn_abers_details(auto[[2L]], new), venn_abers_details(replay, new)))
    }
    error(predict(fits$full, new, prediction_model = "matched"), "prediction_model")
    error(venn_abers_details(fits$full, new, prediction_model = "matched"), "prediction_model")
    stopifnot(any(grepl("full-model transfer", capture.output(print(state)), fixed = TRUE)),
      any(grepl("Full1", capture.output(print(state)), fixed = TRUE)))
    for (size in c(nrow(locals[[1L]]), length(state$pair_models[[1L]]$folds[[1L]]$analysis))) {
      for (option in c("fail_size", "bad_size")) {
        bad_spec <- do.call(parsnip::set_engine, c(list(spec, "reference"), setNames(list(size), option)))
        error(generics::fit(calibrated_workflow(make_base(data, bad_spec), "truth", classification = "ovo",
          model_training = "pairwise", method = method, prediction_model = "full", splits = plans), data), "pair_1")
      }
    }
    before <- length(audit$fits)
    bad <- plans
    bad[[length(bad)]] <- plans[[1L]]
    if (k > 2L) error(generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
      model_training = "pairwise", method = method, prediction_model = "full", splits = bad), data), "data binding")
    error(generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
      model_training = "pairwise", method = method, prediction_model = "full", splits = plans), data[n:1, ]), "data binding")
    stopifnot(length(audit$fits) == before)
    if (k > 2L) {
      # A foreign pair's model is rejected through its fitted two-class mold.
      bad <- fits$full
      s <- bad$fit$fit$fit
      s$pair_models[[1L]]$models$Full1 <- states$matched$pair_models[[2L]]$models[[1L]]
      bad$fit$fit$fit <- s
      error(predict(bad, new, type = "prob"), "mapping")
    }
  }
  # Identical inputs to the pinned Python full-transfer component fixtures.
  reference <- read.csv(file.path("extdata", "full-model-reference.csv"))
  data <- data.frame(p = seq(-3, 3, length.out = 40), x = sin(1:40), row_id = 1:40,
    truth = factor(rep(c("event", "other", "other", "event"), 10), levels = c("other", "event")))
  for (method in c("inductive", "cross")) for (precision in list(NULL, 0, 2, -1)) {
    held <- if (method == "inductive") list(33:40) else list(1:8, 9:20, 21:40)
    fs <- lapply(held, function(h) rsample::make_splits(list(analysis = setdiff(1:40, h), assessment = h), data))
    reference_plan <- if (method == "inductive") fs[[1L]] else rsample::manual_rset(fs, c("small", "middle", "large"))
    fit <- generics::fit(calibrated_workflow(make_base(data), "truth", classification = "ovo",
      model_training = "pairwise", prediction_model = "full", method = method, precision = precision,
      splits = list(pair_1 = reference_plan)), data)
    for (loss in c("log", "brier")) {
      expected <- reference[reference$method == method & reference$loss == loss &
        reference$precision == if (is.null(precision)) "none" else as.character(precision), ]
      new <- data[1:7, ]
      new$p <- expected$x[expected$fold == 1L]
      d <- venn_abers_details(fit, new, loss = loss)
      close(as.matrix(d$details[c("lower", "upper", "probability")]), as.matrix(expected[c("lower", "upper", "probability")]))
      close(d$prediction$.pred_event, expected$combined[expected$fold == 1L])
    }
  }
}

if (exists("full_pair_fit") && all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
  checker <- new.env()
  sys.source(system.file("examples", "calibrated-resampling.R", package = "vennabers"), checker)
  check <- checker$check_evaluation
  scoring <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
  control <- tune::control_resamples(save_pred = TRUE, extract = identity, allow_par = FALSE)
  data <- data.frame(p = sin(1:120), x = cos(1:120), row_id = 1:120,
    truth = factor(rep(c("z", "a b", "a.b"), 40), levels = c("z", "a b", "a.b")))
  base <- make_base(data)
  outer <- withr::with_seed(201, rsample::vfold_cv(data, v = 2, repeats = 2, strata = truth))
  holdout <- rsample::make_splits(list(analysis = 1:90, assessment = 91:120), data)
  for (method in c("inductive", "cross")) for (loss in c("log", "brier")) {
    model <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
      prediction_model = "full", method = method, seed = 42, loss = loss)
    cost <- if (method == "cross") 18L else 6L
    before <- length(audit$fits)
    result <- tune::fit_resamples(model, outer, metrics = scoring, control = control)
    stopifnot(length(audit$fits) == before + nrow(outer) * cost)
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
        engine <- workflows::extract_fit_engine(s$models$Full1)
        stopifnot(identical(engine$x$row_id, locals[[id]]$row_id), !any(engine$x$row_id %in% held$row_id),
          identical(local$row_id[s$rows], locals[[id]]$row_id), identical(names(s$models), "Full1"))
        fs <- lapply(s$folds, function(f) rsample::make_splits(f[c("analysis", "assessment")], locals[[id]]))
        prepared[[length(prepared) + 1L]] <- if (method == "cross") rsample::manual_rset(fs, names(s$folds)) else fs[[1L]]
      }
    }
    replay_model <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
      prediction_model = "full", method = method, seed = 42, loss = loss,
      splits = local_calibration_splits(prepared, names(data)))
    replay <- tune::fit_resamples(replay_model, outer, metrics = scoring, control = control)
    check(replay, outer, classification = "ovo")
    stopifnot(identical(tune::collect_predictions(result), tune::collect_predictions(replay)),
      identical(tune::collect_metrics(result), tune::collect_metrics(replay)))
    for (resampling in c(FALSE, TRUE)) {
      split <- if (resampling) outer$splits[[1L]] else holdout
      changed <- data
      changed$truth[as.integer(split, data = "assessment")] <- levels(data$truth)[1L]
      altered <- rsample::make_splits(list(analysis = as.integer(split, data = "analysis"),
        assessment = as.integer(split, data = "assessment")), changed)
      evaluate <- function(split) if (resampling) tune::fit_resamples(model,
        rsample::manual_rset(list(split), "Outer1"), metrics = scoring, control = control) else
        tune::last_fit(model, split, metrics = scoring)
      before <- length(audit$fits)
      final <- evaluate(split)
      mutated <- evaluate(altered)
      stopifnot(length(audit$fits) == before + 2L * cost, nrow(tune::collect_notes(final)) == 0L,
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
          identical(a[[id]]$calibrator, b[[id]]$calibrator), identical(a[[id]]$folds, b[[id]]$folds),
          identical(workflows::extract_fit_engine(a[[id]]$models$Full1), workflows::extract_fit_engine(b[[id]]$models$Full1)))
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
      # Full-model failure invalidates an entire outer fit, without a source fallback.
      size <- nrow(pair_data(rsample::analysis(outer$splits[[1L]]))[[1L]])
      broken <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail_size = size))
      bad_model <- calibrated_workflow(broken, "truth", classification = "ovo", method = method,
        model_training = "pairwise", prediction_model = "full", seed = 42)
      bad <- suppressWarnings(tune::fit_resamples(bad_model, outer, metrics = scoring, control = control))
      stopifnot(nrow(tune::collect_notes(bad)) > 0L)
      error(check(bad, outer, classification = "ovo"))
      bad <- suppressWarnings(tune::last_fit(bad_model, outer$splits[[1L]], metrics = scoring))
      stopifnot(nrow(tune::collect_notes(bad)) > 0L, is.null(bad[[".workflow"]]))
    }
  }
}

  testthat::expect_true(TRUE)
})

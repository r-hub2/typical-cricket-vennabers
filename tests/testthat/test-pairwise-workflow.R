testthat::test_that("pairwise-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")

packages <- c("parsnip", "workflows", "recipes", "rsample", "withr", "tune")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  close <- function(actual, expected) stopifnot(identical(dim(actual), dim(expected)),
    all(is.finite(actual)), all(abs(actual - expected) <= 1e-10))
  error <- function(expr, text = "") {
    e <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(e, "error"), grepl(text, conditionMessage(e), fixed = TRUE))
  }
  data <- data.frame(x = sin(1:120), z = cos(1:120),
    truth = factor(rep(c("z", "a", "b"), 40), levels = c("z", "a", "b")))
  base <- workflows::workflow() |> workflows::add_recipe(
    recipes::recipe(truth ~ ., data) |> recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  ordered <- function(data) {
    n <- nrow(data)
    held <- seq.int(n - ceiling(n / 4) + 1L, n)
    rsample::make_splits(list(analysis = setdiff(seq_len(n), held), assessment = held), data)
  }
  pairs <- combn(levels(data$truth), 2, simplify = FALSE)
  new <- data[1:7, ]
  for (mode in c("matched", "full")) {
    fit <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
      model_training = "pairwise", prediction_model = mode, splits = ordered, seed = 42), data)
    d <- venn_abers_details(fit, new)
    if (mode == "full") {
      stopifnot(all(vapply(workflows::extract_fit_engine(fit)$pair_models,
        function(x) identical(names(x$models), "Full1"), logical(1))),
        all(d$details$prediction_model == mode))
      next
    }
  q <- matrix(0, nrow(new), length(pairs))
  for (j in seq_along(pairs)) {
    pair <- pairs[[j]]
    local <- data[data$truth %in% pair, ]
    local$truth <- factor(local$truth, levels = pair)
    plan <- ordered(local)
    reference_base <- workflows::workflow() |> workflows::add_recipe(
      recipes::recipe(truth ~ ., local) |> recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(parsnip::logistic_reg(engine = "glm"))
    binary <- generics::fit(reference_base, rsample::analysis(plan))
    held <- rsample::assessment(plan)
    column <- paste0(".pred_", pair[2L])
    cal <- venn_abers(held$truth == pair[2L], predict(binary, held, type = "prob")[[column]])
    expected <- predict(cal, predict(binary, new, type = "prob")[[column]])
    part <- d$details[d$details$pair_id == paste0("pair_", j), ]
    close(as.matrix(part[c("lower", "upper", "probability")]), as.matrix(expected))
    q[, j] <- expected$probability
  }
  # Independent odds formulation of reciprocal OvO coupling.
  p <- vapply(levels(data$truth), function(cl) {
    relevant <- which(vapply(pairs, function(pair) cl %in% pair, logical(1)))
    odds <- vapply(relevant, function(j) {
      prob <- if (cl == pairs[[j]][2L]) q[, j] else 1 - q[, j]
      (1 - prob) / prob
    }, numeric(nrow(new)))
    1 / (1 + rowSums(odds))
  }, numeric(nrow(new)))
  close(as.matrix(d$prediction), p / rowSums(p))
  }
}

if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  list2env(pairwise_reference(), environment())
  for (k in 2:4) {
    classes <- c("z", "a b", "a.b", "4")[seq_len(k)]
    data <- data.frame(x = c(NA, sin(2:90), 1000 + 1:30), row_id = 1:120,
      truth = factor(rep(classes, length.out = 120), levels = classes))
    base <- make_base(data)
    locals <- pair_data(data)
    plans <- lapply(locals, ordered)
    new <- data[c(1:7, 31, 60), ]
    for (precision in list(NULL, 0, 2, -1)) {
      before <- length(audit$fits)
      calls <- length(audit$predictions)
      fit <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
        model_training = "pairwise", splits = plans[rev(names(plans))], seed = 51, precision = precision), data)
      state <- workflows::extract_fit_engine(fit)
      assert_conditions(length(audit$fits) == before + length(locals),
        length(audit$predictions) == calls + length(locals),
        identical(state$class_levels, classes),
        identical(names(state$pair_models), names(locals)), is.null(state$event_class))
      for (id in names(locals)) {
        s <- state$pair_models[[id]]
        engine <- workflows::extract_fit_engine(s$models$Inductive1)
        training <- rsample::analysis(plans[[id]])
        held <- rsample::assessment(plans[[id]])
        assert_conditions(identical(engine$x$row_id, training$row_id), identical(engine$y, training$truth),
          nlevels(engine$y) == 2L, identical(s$event_class, tail(levels(training$truth), 1)),
          !any(engine$x$row_id %in% held$row_id), identical(data$row_id[s$rows], locals[[id]]$row_id),
          identical(s$folds$Inductive1$assessment, as.integer(plans[[id]], data = "assessment")))
        x <- training$x
        x[is.na(x)] <- median(x, na.rm = TRUE)
        close(engine$x$x, (x - mean(x)) / sd(x))
      }
      for (loss in c("log", "brier")) {
        before <- length(audit$fits)
        calls <- length(audit$predictions)
        details <- venn_abers_details(fit, new, loss = loss)
        assert_conditions(length(audit$fits) == before, length(audit$predictions) == calls + length(locals),
          all(details$details$model_training == "pairwise"), all(details$details$prediction_model == "matched"),
          identical(details$details$model_source, details$details$pair_id),
          all(details$details$fold == "Inductive1"), all(details$details$event_class == details$details$class_2))
        q <- matrix(0, nrow(new), length(locals))
        for (j in seq_along(locals)) {
          local <- locals[[j]]
          binary <- generics::fit(make_base(local), rsample::analysis(plans[[j]]))
          held <- rsample::assessment(plans[[j]])
          event <- levels(local$truth)[2L]
          column <- paste0(".pred_", event)
          pcal <- predict(binary, held, type = "prob")[[column]]
          cal <- venn_abers(held$truth == event, pcal, precision = precision)
          expected <- predict(cal, predict(binary, new, type = "prob")[[column]], loss = loss)
          part <- details$details[details$details$pair_id == names(locals)[j], ]
          close(as.matrix(part[c("lower", "upper", "probability")]), as.matrix(expected))
          close(state$pair_models[[j]]$calibration$Inductive1$probability, pcal)
          q[, j] <- expected$probability
        }
        close(as.matrix(details$prediction), couple(q, classes))
        for (n in c(0L, 1L, nrow(new))) {
          input <- new[seq_len(n), , drop = FALSE]
          calls <- length(audit$predictions)
          before <- length(audit$fits)
          result <- predict(fit, input, type = "prob", loss = loss)
          cl <- predict(fit, input, type = "class", loss = loss)
          d <- venn_abers_details(fit, input, loss = loss)
          assert_conditions(identical(result, details$prediction[seq_len(n), ]), identical(result, d$prediction),
            identical(names(d$details), names(details$details)), nrow(d$details) == n * length(locals),
            identical(cl$.pred_class, factor(classes[max.col(as.matrix(result), ties.method = "first")], levels = classes)),
            length(audit$fits) == before, length(audit$predictions) == calls + if (n) 3L * length(locals) else 0L)
        }
      }
    }
    automatic_model <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise", seed = 42)
    set.seed(23)
    rng <- .Random.seed
    automatic <- generics::fit(automatic_model, data)
    stopifnot(identical(.Random.seed, rng))
    again <- generics::fit(automatic_model, data)
    a <- workflows::extract_fit_engine(automatic)$pair_models
    b <- workflows::extract_fit_engine(again)$pair_models
    for (id in names(a)) stopifnot(identical(a[[id]]$folds, b[[id]]$folds),
      identical(a[[id]]$calibration, b[[id]]$calibration), identical(a[[id]]$calibrator, b[[id]]$calibrator),
      identical(workflows::extract_fit_engine(a[[id]]$models$Inductive1),
                workflows::extract_fit_engine(b[[id]]$models$Inductive1)))
    state <- workflows::extract_fit_engine(automatic)
    prepared <- lapply(names(locals), function(id) {
      plan <- withr::with_seed(42, rsample::initial_split(locals[[id]], .8, strata = truth, pool = 0))
      stopifnot(identical(state$pair_models[[id]]$folds$Inductive1$analysis, as.integer(plan, data = "analysis")))
      plan
    })
    resolver <- local_calibration_splits(prepared, names(data))
    explicit <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
      model_training = "pairwise", splits = resolver, seed = 42), data)
    stopifnot(identical(venn_abers_details(automatic, new), venn_abers_details(explicit, new)))
    # A real shared model sees every class; pair models do not substitute its outputs.
    shared <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo", splits = ordered), data)
    stopifnot(nlevels(workflows::extract_fit_engine(workflows::extract_fit_engine(shared)$models$Inductive1)$y) == k)
    if (k > 2L) stopifnot(any(abs(as.matrix(predict(shared, new, type = "prob")) -
      as.matrix(predict(fit, new, type = "prob"))) > 1e-6))
  # Every malformed pair plan fails before the first model fit, including later pairs.
  before <- length(audit$fits)
  malformed <- if (length(plans) > 1L) list(plans[-1L], c(plans, plans[1L]),
    setNames(plans, rev(names(plans))), setNames(plans, rep("pair_1", length(plans))), unname(plans)) else list()
  for (bad in malformed) {
    error(generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
      model_training = "pairwise", splits = bad), data))
  }
  for (change in c("values", "type", "levels", "order", "missing", "overlap", "duplicate", "fraction", "empty", "cal_class", "train_class")) {
    bad <- plans
    j <- length(bad)
    if (change == "values") bad[[j]]$data$x[1L] <- 456
    if (change == "type") bad[[j]]$data$row_id <- as.double(bad[[j]]$data$row_id)
    if (change == "levels") bad[[j]]$data$truth <- factor(bad[[j]]$data$truth, levels = rev(levels(bad[[j]]$data$truth)))
    if (change == "order") bad[[j]]$data <- bad[[j]]$data[rev(seq_len(nrow(bad[[j]]$data))), ]
    if (change == "missing") bad[[j]]$in_id <- bad[[j]]$in_id[-1L]
    if (change == "overlap") bad[[j]]$in_id <- c(bad[[j]]$in_id, bad[[j]]$out_id[1L])
    if (change == "duplicate") bad[[j]]$in_id <- c(bad[[j]]$in_id, bad[[j]]$in_id[1L])
    if (change == "fraction") bad[[j]]$in_id[1L] <- 1.5
    if (change == "empty") bad[[j]]$in_id <- integer()
    if (change %in% c("cal_class", "train_class")) {
      h <- which(locals[[j]]$truth == levels(locals[[j]]$truth)[1L])
      a <- setdiff(seq_len(nrow(locals[[j]])), h)
      bad[[j]] <- rsample::make_splits(if (change == "cal_class") list(analysis = a, assessment = h) else
        list(analysis = h, assessment = a), locals[[j]])
    }
    error(generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
      model_training = "pairwise", splits = bad), data),
      paste0("pair_", length(plans), " Inductive1"))
  }
  for (args in c(if (k > 2L) list(list(classification = "binary")) else list(), list(list(classification = "ovr")),
    list(classification = "ovo", event_level = "second"), list(classification = "ovo", m = NULL),
    list(classification = "ovo", splits = ordered(data)), list(classification = "ovo", splits = plans, strata = FALSE))) {
    error(do.call(calibrated_workflow, c(list(base, "truth", model_training = "pairwise"), args)))
  }
  regression <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  error(calibrated_workflow(regression, "truth", classification = "ovo", model_training = "pairwise"))
  error(calibrated_workflow(state$pair_models[[1L]]$models$Inductive1, "truth", classification = "ovo", model_training = "pairwise"))
  placeholder <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(penalty = tune::tune(), engine = "glm"))
  error(calibrated_workflow(placeholder, "truth", classification = "ovo", model_training = "pairwise"), "placeholders")
  stopifnot(length(audit$fits) == before)
  conflict <- prepared[[1L]]
  conflict$in_id <- rev(conflict$in_id)
  error(local_calibration_splits(c(prepared, list(conflict)), names(data)), "Conflicting")
  if (length(plans) > 1L) for (change in c("pair", "classes", "model")) {
    bad <- fit
    s <- bad$fit$fit$fit
    if (change == "pair") s$pair_models[c(1, 2)] <- s$pair_models[c(2, 1)]
    if (change == "classes") s$pair_models[[1L]]$class_levels <- s$pair_models[[2L]]$class_levels
    if (change == "model") s$pair_models[[1L]]$models <- s$pair_models[[2L]]$models
    bad$fit$fit$fit <- s
    error(predict(bad, new, type = "prob"), "mapping")
  }
  for (bad in list(list(fail_size = length(plans[[1L]]$in_id)), list(bad = TRUE))) {
    broken <- workflows::update_model(base, do.call(parsnip::set_engine, c(list(spec, "reference"), bad)))
    error(generics::fit(calibrated_workflow(broken, "truth", classification = "ovo",
      model_training = "pairwise", splits = plans), data), "pair_1")
  }
  # A two-class global tie picks the first level, although the pair event is second.
  tied <- data.frame(x = rep(1, 20), row_id = rep(1L, 20), truth = factor(rep(c("z", "a"), 10), levels = c("z", "a")))
  tied_base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |> workflows::add_model(spec)
  tied_fit <- generics::fit(calibrated_workflow(tied_base, "truth", classification = "ovo",
    model_training = "pairwise", splits = function(data) rsample::make_splits(list(analysis = 1:12, assessment = 13:20), data)), tied)
  for (loss in c("log", "brier")) stopifnot(all(predict(tied_fit, tied, type = "prob", loss = loss) == .5),
    all(predict(tied_fit, tied, type = "class", loss = loss)$.pred_class == "z"))
  }
}

# Minimal and unequal calibration sets on an imbalanced global domain.
if (exists("pair_fit")) local({
  data <- data.frame(x = sin(1:24),
    truth = factor(rep(c("z", "a", "b"), c(6, 8, 10)), levels = c("z", "a", "b")))
  base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, data)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (minimum in c(TRUE, FALSE)) {
    plan <- function(local) {
      held <- unlist(lapply(levels(local$truth), function(cl) head(which(local$truth == cl),
        if (minimum) 1L else match(cl, levels(data$truth)))))
      rsample::make_splits(list(analysis = setdiff(seq_len(nrow(local)), held), assessment = held), local)
    }
    for (precision in list(NULL, 0, 2, -1)) {
      fit <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
        model_training = "pairwise", splits = plan, seed = 10, precision = precision), data)
      state <- workflows::extract_fit_engine(fit)
      sizes <- vapply(state$pair_models, function(s) s$folds$Inductive1$calibration_size, 1L)
      stopifnot(identical(unname(sizes), if (minimum) c(2L, 2L, 2L) else c(3L, 4L, 5L)))
      for (loss in c("log", "brier")) {
        actual <- venn_abers_details(fit, data, loss = loss)
        q <- matrix(0, nrow(data), length(state$pair_models))
        for (j in seq_along(state$pair_models)) {
          component <- state$pair_models[[j]]
          local <- data[component$rows, ]
          local$truth <- factor(local$truth, levels = component$class_levels)
          split <- plan(local)
          reference <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
            workflows::add_model(parsnip::logistic_reg(engine = "glm"))
          reference <- generics::fit(reference, rsample::analysis(split))
          held <- rsample::assessment(split)
          column <- paste0(".pred_", component$event_class)
          cal <- venn_abers(held$truth == component$event_class,
            predict(reference, held, type = "prob")[[column]], precision = precision)
          expected <- predict(cal, predict(reference, data, type = "prob")[[column]], loss = loss)
          close(as.matrix(actual$details[actual$details$pair_id == component$pair_id,
            c("lower", "upper", "probability")]), as.matrix(expected))
          q[, j] <- expected$probability
        }
        close(as.matrix(actual$prediction), couple(q, levels(data$truth)))
      }
    }
  }
})

if (exists("pair_fit") && all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
  for (k in 2:4) {
    classes <- c("z", "a b", "a.b", "4")[seq_len(k)]
    data <- data.frame(x = c(NA, sin(2:90), 1000 + 1:30), row_id = 1:120,
      truth = factor(rep(classes, length.out = 120), levels = classes))
    base <- make_base(data)
    pair_count <- length(pair_data(data))
    plans <- lapply(pair_data(data), ordered)
  checker <- new.env()
  sys.source(system.file("examples", "calibrated-resampling.R", package = "vennabers"), checker)
  check <- checker$check_evaluation
  scoring <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
  control <- tune::control_resamples(save_pred = TRUE, extract = identity, allow_par = FALSE)
  outer <- withr::with_seed(201, rsample::vfold_cv(data, v = 2, repeats = 2, strata = truth))
  assessment <- tail(seq_len(nrow(data)), ceiling(nrow(data) / 5L))
  holdout <- rsample::make_splits(list(analysis = setdiff(seq_len(nrow(data)), assessment), assessment = assessment), data)
  for (loss in c("log", "brier")) {
    model <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise", seed = 42, loss = loss)
    before <- length(audit$fits)
    result <- tune::fit_resamples(model, outer, metrics = scoring, control = control)
    stopifnot(length(audit$fits) == before + nrow(outer) * choose(nlevels(data$truth), 2))
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
        engine <- workflows::extract_fit_engine(s$models$Inductive1)
        stopifnot(identical(engine$x$row_id, locals[[id]]$row_id[s$folds$Inductive1$analysis]),
          !any(engine$x$row_id %in% held$row_id), identical(local$row_id[s$rows], locals[[id]]$row_id))
        prepared[[length(prepared) + 1L]] <- rsample::make_splits(
          s$folds$Inductive1[c("analysis", "assessment")], locals[[id]])
      }
    }
    resolver <- local_calibration_splits(prepared, names(data))
    explicit <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
      splits = resolver, seed = 42, loss = loss)
    replay <- tune::fit_resamples(explicit, outer, metrics = scoring, control = control)
    check(replay, outer, classification = "ovo")
    stopifnot(identical(tune::collect_predictions(result), tune::collect_predictions(replay)),
      identical(tune::collect_metrics(result), tune::collect_metrics(replay)))
    model <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
      splits = ordered, seed = 42, loss = loss)
    ordered_result <- tune::fit_resamples(model, outer, metrics = scoring, control = control)
    check(ordered_result, outer, classification = "ovo")
    before <- length(audit$fits)
    final <- tune::last_fit(model, holdout, metrics = scoring)
    p <- tune::collect_predictions(final)
    m <- tune::collect_metrics(final)
    assert_conditions(length(audit$fits) == before + pair_count, nrow(tune::collect_notes(final)) == 0L,
      identical(p$.row, assessment), !anyDuplicated(p$.row), all(p$.config == "pre0_mod0_post0"),
      identical(p$truth, data$truth[assessment]), nrow(m) == 2L,
      setequal(m$.metric, c("brier_class", "mn_log_loss")), all(is.finite(m$.estimate)))
    changed <- data
    changed$truth[assessment] <- classes[1L]
    mutated <- tune::last_fit(model, rsample::make_splits(
      list(analysis = setdiff(seq_len(nrow(data)), assessment), assessment = assessment), changed), metrics = scoring)
    one <- tune::extract_workflow(final)
    two <- tune::extract_workflow(mutated)
    a <- workflows::extract_fit_engine(one)$pair_models
    b <- workflows::extract_fit_engine(two)$pair_models
    for (id in names(a)) stopifnot(identical(a[[id]]$calibration, b[[id]]$calibration),
      identical(a[[id]]$calibrator, b[[id]]$calibrator), identical(a[[id]]$folds, b[[id]]$folds),
      identical(workflows::extract_fit_engine(a[[id]]$models$Inductive1),
                workflows::extract_fit_engine(b[[id]]$models$Inductive1)))
    stopifnot(identical(venn_abers_details(one, data), venn_abers_details(two, data)))
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
      partial <- rsample::manual_rset(list(holdout,
        rsample::make_splits(list(analysis = 1:80, assessment = 81:120), data)), c("First", "Second"))
      sizes <- vapply(partial$splits, function(split) nrow(rsample::analysis(
        ordered(pair_data(rsample::analysis(split))[[1L]]))), 1L)
      for (sizes in list(sizes[1L], sizes)) {
        broken <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail_size = sizes))
        bad_model <- calibrated_workflow(broken, "truth", classification = "ovo", model_training = "pairwise", splits = ordered)
        bad <- suppressWarnings(tune::fit_resamples(bad_model, partial, metrics = scoring, control = control))
        stopifnot(nrow(tune::collect_notes(bad)) == length(sizes))
        error(check(bad, partial, classification = "ovo"))
      }
      fail_metric <- yardstick::new_prob_metric(function(data, ...) stop("metric failed"), direction = "minimize")
      bad <- suppressWarnings(tune::fit_resamples(model, outer,
        metrics = yardstick::metric_set(yardstick::brier_class, fail_metric), control = control))
      stopifnot(nrow(tune::collect_notes(bad)) > 0L)
      error(check(bad, outer, classification = "ovo"))
      # Development-wide pair plans cannot leak into smaller outer domains.
      bad_model <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise", splits = plans)
      bad <- suppressWarnings(tune::fit_resamples(bad_model, outer, metrics = scoring, control = control))
      stopifnot(nrow(tune::collect_notes(bad)) == nrow(outer))
      error(check(bad, outer, classification = "ovo"))
    }
  }
  }
}

  testthat::expect_true(TRUE)
})

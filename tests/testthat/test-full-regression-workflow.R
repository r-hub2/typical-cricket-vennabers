testthat::test_that("full-regression-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")

packages <- c("parsnip", "workflows", "rsample", "recipes", "withr")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  # Fixed before comparison: one outcome unit; scale s with any unit change.
  close <- function(actual, expected, s = 1) assert_conditions(
    identical(dim(actual), dim(expected)), length(actual) == length(expected),
    all(is.finite(actual)), all(is.finite(expected)),
    all(abs(actual - expected) <= 1e-10 * s + 1e-10 * abs(expected)))
  data <- data.frame(x = sin(1:80), truth = 10 + cos(1:80) + 3 * sin(1:80))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  for (method in c("inductive", "cross")) {
    matched <- generics::fit(calibrated_workflow(base, "truth", method = method, seed = 42), data)
    full <- generics::fit(calibrated_workflow(base, "truth", method = method,
      prediction_model = "full", seed = 42), data)
    source <- workflows::extract_fit_engine(matched)
    state <- workflows::extract_fit_engine(full)
    assert_conditions(identical(state$calibrator, source$calibrator),
      identical(state$calibration, source$calibration), identical(state$folds, source$folds),
      identical(names(state$models), "Full1"), state$parameters$effective$prediction_model == "full")
    input <- predict(state$models$Full1, data)$.pred
    p <- if (method == "cross") matrix(input, nrow(data), length(state$folds),
      dimnames = list(NULL, names(state$folds))) else input
    expected <- predict(source$calibrator, p)
    actual <- venn_abers_details(full, data)
    close(actual$prediction$.pred, expected$prediction)
    stopifnot(all(actual$details$prediction_model == "full"), all(actual$details$model_source == "Full1"))
  }
}

if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  error <- function(expr, text = "") {
    e <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(e, "error"), grepl(text, conditionMessage(e), fixed = TRUE))
  }
  audit <- new.env()
  audit$fits <- audit$predictions <- list()
  transfer_fit <- function(x, y, fail_size = integer(), bad_size = integer(), raw = FALSE) {
    audit$fits[[length(audit$fits) + 1L]] <- list(x = x, y = y)
    if (nrow(x) %in% fail_size) stop("untrainable regression model")
    list(x = x, y = y, offset = mean(y), scale = sd(y), bad_size = bad_size, raw = raw)
  }
  transfer_predict <- function(object, new_data) {
    stopifnot(nrow(new_data) > 0L)
    audit$predictions[[length(audit$predictions) + 1L]] <- list(train = object$x$id, rows = new_data$id)
    p <- if (object$raw) new_data$x else object$offset + object$scale * new_data$x
    if (nrow(object$x) %in% object$bad_size) p[1L] <- Inf
    p
  }
  register_reference_engine("va_regression_transfer", "regression", "transfer_fit", "transfer_predict", "numeric")
  spec <- parsnip::new_model_spec("va_regression_transfer", args = list(), eng_args = NULL,
    mode = "regression", method = NULL, engine = "reference")
  make_base <- function(data, raw = FALSE) {
    pre <- recipes::recipe(truth ~ x + id, data)
    if (!raw) pre <- pre |> recipes::step_impute_median(x) |> recipes::step_normalize(x)
    workflows::workflow() |> workflows::add_recipe(pre) |>
      workflows::add_model(parsnip::set_engine(spec, "reference", raw = raw))
  }
  make_plan <- function(data, method, held = list(1:3, 4:11, 12:30)) {
    if (method == "inductive") held <- held[2L]
    splits <- lapply(held, function(h) rsample::make_splits(
      list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data))
    if (method == "inductive") splits[[1L]] else rsample::manual_rset(splits, paste0("part", seq_along(splits)))
  }
  # Independent max-min isotonic characterization, including the hypothetical
  # observation and equal-score pooling. No production calibrator or hull scan.
  isotonic <- function(y, p, new, outer) {
    y <- c(y, outer); p <- c(p, new)
    knots <- sort(unique(p)); at <- match(new, knots)
    max(vapply(seq_len(at), function(left) min(vapply(seq.int(at, length(knots)), function(right)
      mean(y[p >= knots[left] & p <= knots[right]]), 0)), 0))
  }
  reference <- function(cal, new, m) {
    ranked <- sort(cal$truth); k <- length(ranked)
    a <- ranked[m]; b <- ranked[k - m + 1L]
    lo <- pmin(pmax(cal$truth, a), ranked[k - m])
    hi <- pmin(pmax(cal$truth, ranked[m + 1L]), b)
    lower <- vapply(new, function(p) isotonic(lo, cal$prediction, p, a), 0)
    upper <- vapply(new, function(p) isotonic(hi, cal$prediction, p, b), 0)
    weight <- lower - a + b - upper
    point <- (lower * (lower - a) + upper * (b - upper)) / weight
    point[weight == 0] <- ((lower + upper) / 2)[weight == 0]
    data.frame(prediction = point, lower = lower, upper = upper)
  }
  original <- data.frame(x = c(100, NA, sin(3:30)), id = 1:30, truth = 2 + cos(1:30) + sin(1:30))
  for (scale in c(1e-12, 1, 1e120)) for (shift in c(0, 100)) {
    data <- original; data$truth <- scale * (data$truth + shift)
    base <- make_base(data)
    new <- data[c(1, 2, 8, 15, 30), ]
    for (method in c("inductive", "cross")) for (option in list(list(), list(m = 1), list(epsilon = .5), list(m = 2))) {
      held <- if (identical(option, list(m = 2))) list(1:5, 6:12, 13:30) else list(1:3, 4:11, 12:30)
      plan <- make_plan(data, method, held)
      fits <- states <- list()
      for (mode in c("matched", "full")) {
        before <- length(audit$fits); calls <- length(audit$predictions)
        set.seed(901); rng <- .Random.seed
        fits[[mode]] <- generics::fit(do.call(calibrated_workflow, c(list(base, "truth", method = method,
          prediction_model = mode, splits = plan, seed = 42), option)), data)
        s <- states[[mode]] <- workflows::extract_fit_engine(fits[[mode]])
        assert_conditions(length(audit$fits) == before + length(s$folds) + (mode == "full"),
          length(audit$predictions) == calls + length(s$folds) + (mode == "full"), identical(rng, .Random.seed),
          s$parameters$effective$prediction_model == mode, is.null(s$class_levels), is.null(s$loss))
        for (i in seq_along(s$folds)) {
          f <- s$folds[[i]]; learned <- audit$fits[[before + i]]
          x <- data$x[f$analysis]; x[is.na(x)] <- median(x, na.rm = TRUE)
          close(learned$x$x, (x - mean(x)) / sd(x))
          assert_conditions(identical(learned$x$id, data$id[f$analysis]),
            identical(as.vector(learned$y), data$truth[f$analysis]),
            identical(audit$predictions[[calls + i]], list(train = data$id[f$analysis], rows = data$id[f$assessment])),
            identical(s$calibration[[i]]$row, f$assessment), identical(s$calibration[[i]]$truth, data$truth[f$assessment]))
        }
      }
      s <- states$full
      stopifnot(identical(s$calibrator, states$matched$calibrator), identical(s$calibration, states$matched$calibration),
        identical(s$folds, states$matched$folds), identical(names(s$models), "Full1"))
      engine <- workflows::extract_fit_engine(s$models$Full1)
      recipe <- workflows::extract_recipe(s$models$Full1)
      x <- data$x; x[is.na(x)] <- median(x, na.rm = TRUE)
      stopifnot(identical(engine$x$id, data$id), identical(as.vector(engine$y), data$truth))
      close(recipe$steps[[1L]]$medians[["x"]], median(x))
      close(recipe$steps[[2L]]$means[["x"]], mean(x))
      close(engine$x$x, (x - mean(x)) / sd(x))
      input <- predict(s$models$Full1, new)$.pred
      before <- length(audit$fits); calls <- length(audit$predictions)
      actual <- venn_abers_details(fits$full, new)
      stopifnot(length(audit$predictions) == calls + 1L,
        all(actual$details$prediction_model == "full"), all(actual$details$model_source == "Full1"))
      expected <- lapply(names(s$folds), function(id) {
        f <- s$folds[[id]]; cal <- s$calibration[[id]]
        close(cal$prediction, predict(states$matched$models[[id]], data[f$assessment, ])$.pred, scale)
        m <- if (length(option$epsilon)) floor(option$epsilon * (nrow(cal) + 1) / 2) else if (length(option$m)) option$m else 1
        stopifnot(f$moderation$m == m, s$parameters$effective$m[[id]] == m,
          f$moderation$epsilon == 2 * m / (nrow(cal) + 1), identical(f$moderation$requested, option))
        ref <- reference(cal, input, m)
        close(as.matrix(actual$details[actual$details$fold == id, c("prediction", "lower", "upper")]), as.matrix(ref), scale)
        ref
      })
      mean_points <- Reduce(`+`, lapply(expected, `[[`, "prediction")) / length(expected)
      close(actual$prediction$.pred, mean_points, scale)
      stopifnot(identical(actual$prediction, predict(fits$full, new)),
        identical(actual$prediction, predict(fits$full, new, type = "numeric")))
      if (method == "cross") {
        close(as.matrix(actual$aggregate[c("lower", "upper")]),
          as.matrix(Reduce(`+`, expected)[c("lower", "upper")] / length(expected)), scale)
        stopifnot(all(actual$aggregate$prediction_model == "full"),
          identical(actual$aggregate$model_source, rep(list("Full1"), nrow(new))))
      }
      if (scale == 1 && shift == 0 && !length(option)) {
        stopifnot(any(abs(actual$prediction$.pred - predict(fits$matched, new)$.pred) > 1e-5))
        if (method == "cross") {
          weights <- vapply(s$folds, `[[`, 0, "calibration_size")
          wrong <- as.vector(do.call(cbind, lapply(expected, `[[`, "prediction")) %*% weights / sum(weights))
          pooled <- do.call(rbind, s$calibration)
          stopifnot(any(abs(mean_points - wrong) > 1e-5),
            any(abs(mean_points - reference(pooled, input, 1)$prediction) > 1e-5))
          # Re-merging averaged endpoints, even with averaged a/b, is wrong.
          bounds <- Reduce(`+`, expected) / length(expected)
          a <- mean(vapply(s$calibrator$calibrators, `[[`, 0, "a"))
          b <- mean(vapply(s$calibrator$calibrators, `[[`, 0, "b"))
          wrong <- with(bounds, (lower * (lower - a) + upper * (b - upper)) / (lower - a + b - upper))
          stopifnot(any(abs(mean_points - wrong) > 1e-5))
        }
      }
      for (n in 0:1) {
        calls <- length(audit$predictions)
        small <- venn_abers_details(fits$full, new[seq_len(n), ])
        stopifnot(identical(small$prediction, predict(fits$full, new[seq_len(n), ])),
          length(audit$predictions) == calls + 2L * n,
          nrow(small$details) == n * length(s$folds),
          identical(vapply(small$details, typeof, ""), vapply(actual$details, typeof, "")))
        if (method == "cross") stopifnot(nrow(small$aggregate) == n,
          identical(vapply(small$aggregate, typeof, ""), vapply(actual$aggregate, typeof, "")))
      }
      singles <- do.call(rbind, lapply(seq_len(nrow(new)), function(i) predict(fits$full, new[i, ])))
      stopifnot(identical(actual$prediction, singles), length(audit$fits) == before)
    }
  }
}

if (exists("transfer_fit")) {
  # Offline Python 1.5.4 components: feed identical held-out and new values
  # through a raw numeric engine. Invalid direct Python cells are excluded;
  # the independent repaired oracle and the exact zero-weight cases remain.
  read_fixture <- function(name, text_columns) {
    d <- read.csv(file.path("extdata", name), colClasses = "character")
    for (column in setdiff(names(d), text_columns)) d[[column]] <- as.numeric(d[[column]])
    d
  }
  calibration <- read_fixture("regression-calibration.csv", "case_id")
  expected <- read_fixture("regression-predictions.csv", "case_id")
  options <- read_fixture("regression-options.csv", c("case_id", "parameter"))
  for (i in which(!startsWith(options$case_id, "random_") | startsWith(options$case_id, "random_00_"))) {
    option <- options[i, ]
    cal <- calibration[calibration$case_id == option$case_id, ]
    wanted <- expected[expected$case_id == option$case_id, ]
    data <- data.frame(x = rep(cal$prediction, 2), id = seq_len(2 * nrow(cal)), truth = rep(cal$truth, 2))
    new <- data.frame(x = wanted$input, id = seq_len(nrow(wanted)))
    held <- list(seq_len(nrow(cal)), nrow(cal) + seq_len(nrow(cal)))
    for (method in c("inductive", "cross")) {
      model <- do.call(calibrated_workflow, c(list(make_base(data, raw = TRUE), "truth", method = method,
        prediction_model = "full", splits = make_plan(data, method, held)), setNames(list(option$requested), option$parameter)))
      fit <- generics::fit(model, data)
      actual <- venn_abers_details(fit, new)
      for (column in c("prediction", "lower", "upper")) {
        value <- if (column == "prediction") actual$prediction$.pred else if (method == "cross")
          actual$aggregate[[column]] else actual$details[[column]]
        close(value, wanted[[column]], option$scale)
        valid <- wanted$direct_valid == 1
        close(value[valid], wanted[[paste0("direct_", column)]][valid], option$scale)
      }
    }
  }
  data <- original; base <- make_base(data)
  before <- length(audit$fits)
  for (option in list(list(classification = "binary"), list(event_level = "first"), list(precision = NULL),
    list(loss = "log"), list(model_training = "shared"), list(m = 1, epsilon = .5), list(m = TRUE),
    list(m = 1.5), list(epsilon = 0), list(epsilon = 1), list(epsilon = NA_real_))) {
    error(do.call(calibrated_workflow, c(list(base, "truth", prediction_model = "full"), option)))
  }
  for (bad in list(NULL, "f", "FULL", NA_character_, c("full", "matched"), TRUE))
    error(calibrated_workflow(base, "truth", prediction_model = bad), "prediction_model")
  error(calibrated_workflow(base, "truth", prediction_m = "full"), "prediction_m")
  for (method in c("inductive", "cross")) {
    plan <- make_plan(data, method)
    error(generics::fit(calibrated_workflow(base, "truth", method = method,
      prediction_model = "full", splits = plan), data[30:1, ]), "data binding")
    for (option in list(list(m = 9), list(epsilon = .01))) {
      error(generics::fit(do.call(calibrated_workflow, c(list(base, "truth", method = method,
        prediction_model = "full", splits = plan), option)), data), "calibration size")
    }
    invalid <- make_plan(data, method, list(1:28, 29:30))
    error(generics::fit(calibrated_workflow(base, "truth", method = method,
      prediction_model = "full", splits = invalid), data), "calibration size 2")
  }
  stopifnot(length(audit$fits) == before)
  for (method in c("inductive", "cross")) {
    for (size in c(30L, if (method == "cross") 27L else 22L)) {
      for (arg in c("fail_size", "bad_size")) {
        engine <- do.call(parsnip::set_engine, c(list(spec, "reference"), setNames(list(size), arg)))
        broken <- workflows::update_model(base, engine)
        error(generics::fit(calibrated_workflow(broken, "truth", method = method, prediction_model = "full",
          splits = make_plan(data, method)), data), if (size == 30L) "Full1" else if (method == "cross") "part1" else "Inductive1")
      }
    }
    for (truth in list(rep(7, 30), rep(c(0, 1), 15))) {
      d <- data; d$truth <- truth
      fit <- generics::fit(calibrated_workflow(make_base(d), "truth", method = method,
        prediction_model = "full", seed = 42), d)
      stopifnot(is.double(predict(fit, d)$.pred), is.null(workflows::extract_fit_engine(fit)$class_levels))
      if (all(truth == 7)) stopifnot(all(predict(fit, d)$.pred == 7))
    }
    zero <- data.frame(x = rep(as.double(1:4), 2), id = 1:8, truth = rep(c(0, 0, 2, 2), 2))
    fit <- generics::fit(calibrated_workflow(make_base(zero, raw = TRUE), "truth", method = method,
      prediction_model = "full", splits = make_plan(zero, method, list(1:4, 5:8))), zero)
    d <- venn_abers_details(fit, data.frame(x = 2.5, id = 1L))
    stopifnot(all(d$details$lower == 0), all(d$details$upper == 2), all(d$details$prediction == 1))
    auto <- generics::fit(calibrated_workflow(base, "truth", method = method, prediction_model = "full", seed = 42), data)
    plan <- withr::with_seed(42, if (method == "cross") rsample::vfold_cv(data, v = 5) else rsample::initial_split(data, .8))
    explicit <- generics::fit(calibrated_workflow(base, "truth", method = method, prediction_model = "full",
      splits = local_calibration_splits(list(plan), names(data)), seed = 42), data)
    stopifnot(identical(venn_abers_details(auto, data), venn_abers_details(explicit, data)),
      any(grepl("full-model transfer", capture.output(print(workflows::extract_fit_engine(auto))), fixed = TRUE)))
    for (type in c("prob", "class", "conf_int", "pred_int")) error(predict(auto, data, type = type), "type")
    error(predict(auto, data, prediction_model = "matched"), "prediction_model")
    error(venn_abers_details(auto, data, prediction_model = "matched"), "prediction_model")
    error(predict(auto, data, loss = "log"), "Regression")
    error(venn_abers_details(auto, data, loss = NULL), "Regression")
    error(calibrate_workflow(workflows::extract_fit_engine(auto)$models$Full1, data, "truth",
      prediction_model = "full"), "prediction_model")
    # Change predictors exclusively in one calibration fold. Its source recipe
    # and fitted engine are unchanged; the full recipe learns the changed data.
    altered <- data; altered$x[4:11] <- 1e6
    before <- length(audit$fits)
    one <- generics::fit(calibrated_workflow(base, "truth", method = method,
      prediction_model = "full", splits = make_plan(data, method)), data)
    two <- generics::fit(calibrated_workflow(base, "truth", method = method,
      prediction_model = "full", splits = make_plan(altered, method)), altered)
    source_index <- if (method == "cross") 2L else 1L
    count <- if (method == "cross") 4L else 2L
    stopifnot(identical(audit$fits[[before + source_index]], audit$fits[[before + count + source_index]]))
    a <- workflows::extract_fit_engine(one); b <- workflows::extract_fit_engine(two)
    stopifnot(!identical(workflows::extract_recipe(a$models$Full1)$steps,
      workflows::extract_recipe(b$models$Full1)$steps))
  }
}

if (exists("transfer_fit") && all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
  example <- new.env()
  sys.source(system.file("examples", "calibrated-resampling.R", package = "vennabers"), example)
  check <- function(result, outer) example$check_evaluation(result, outer, classification = "regression")
  data <- data.frame(x = c(100, NA, sin(3:60)), id = 1:60, truth = 2 + cos(1:60) + sin(1:60))
  base <- make_base(data)
  outer <- withr::with_seed(201, rsample::vfold_cv(data, v = 2, repeats = 2))
  holdout <- rsample::make_splits(list(analysis = 1:48, assessment = 49:60), data)
  scoring <- yardstick::metric_set(yardstick::rmse, yardstick::mae)
  control <- tune::control_resamples(save_pred = TRUE, extract = identity, allow_par = FALSE)
  for (method in c("inductive", "cross")) for (option in list(list(m = 2), list(epsilon = .5))) {
    options <- c(list(base, "truth", method = method, prediction_model = "full", seed = 42), option)
    model <- do.call(calibrated_workflow, options)
    before <- length(audit$fits)
    result <- tune::fit_resamples(model, outer, metrics = scoring, control = control)
    checked <- check(result, outer)
    fits <- tune::collect_extracts(result)$.extracts
    stopifnot(length(audit$fits) == before + nrow(outer) * if (method == "cross") 6L else 2L)
    offset <- before
    prepared <- lapply(seq_along(fits), function(i) {
      state <- workflows::extract_fit_engine(fits[[i]])
      local <- rsample::analysis(outer$splits[[i]])
      held <- rsample::assessment(outer$splits[[i]])
      for (f in c(state$folds, list(list(analysis = seq_len(nrow(local)))))) {
        offset <<- offset + 1L
        engine <- audit$fits[[offset]]
        x <- local$x[f$analysis]; x[is.na(x)] <- median(x, na.rm = TRUE)
        stopifnot(identical(engine$x$id, local$id[f$analysis]), !any(engine$x$id %in% held$id),
          identical(as.vector(engine$y), local$truth[f$analysis]))
        close(engine$x$x, (x - mean(x)) / sd(x))
      }
      rows <- checked$predictions$id == outer$id[i] & checked$predictions$id2 == outer$id2[i]
      p <- checked$predictions[rows, ]
      close(p$.pred, predict(fits[[i]], held)$.pred)
      input <- predict(state$models$Full1, held)$.pred
      expected <- lapply(names(state$folds), function(id)
        reference(state$calibration[[id]], input, state$parameters$effective$m[[id]])$prediction)
      close(p$.pred, Reduce(`+`, expected) / length(expected))
      parts <- lapply(state$folds, function(f) rsample::make_splits(f[c("analysis", "assessment")], local))
      if (method == "cross") rsample::manual_rset(parts, names(state$folds)) else parts[[1L]]
    })
    explicit <- do.call(calibrated_workflow, c(options, list(splits = local_calibration_splits(prepared, names(data)))))
    replay <- tune::fit_resamples(explicit, outer, metrics = scoring, control = control)
    check(replay, outer)
    stopifnot(identical(tune::collect_predictions(result), tune::collect_predictions(replay)),
      identical(tune::collect_metrics(result), tune::collect_metrics(replay)))
    ordered <- function(local) make_plan(local, method, if (method == "inductive")
      list(NULL, seq.int(nrow(local) - 7L, nrow(local))) else
        split(seq_len(nrow(local)), cut(seq_len(nrow(local)), 3L, labels = FALSE)))
    ordered_model <- do.call(calibrated_workflow, c(options, list(splits = ordered)))
    ordered_result <- tune::fit_resamples(ordered_model, outer, metrics = scoring, control = control)
    check(ordered_result, outer)
    # last_fit covers automatic, explicit local bindings and ordered plans.
    for (plan in list(NULL, local_calibration_splits(list(ordered(rsample::analysis(holdout))), names(data)), ordered)) {
      model <- do.call(calibrated_workflow, c(options, list(splits = plan)))
      before <- length(audit$fits)
      final <- tune::last_fit(model, holdout, metrics = scoring)
      stopifnot(length(audit$fits) == before + if (method == "cross") (if (is.null(plan)) 6L else 4L) else 2L)
      p <- tune::collect_predictions(final); metrics <- tune::collect_metrics(final)
      assert_conditions(nrow(tune::collect_notes(final)) == 0L, identical(p$.row, 49:60), !anyDuplicated(p$.row),
        all(p$.config == "pre0_mod0_post0"), identical(p$truth, data$truth[49:60]), nrow(metrics) == 2L,
        setequal(metrics$.metric, c("rmse", "mae")), all(is.finite(metrics$.estimate)),
        all(metrics$.config == "pre0_mod0_post0"))
      changed <- data; changed$truth[49:60] <- changed$truth[49:60] + 1000
      altered <- rsample::make_splits(list(analysis = 1:48, assessment = 49:60), changed)
      mutated <- tune::last_fit(model, altered, metrics = scoring)
      one <- tune::extract_workflow(final); two <- tune::extract_workflow(mutated)
      a <- workflows::extract_fit_engine(one); b <- workflows::extract_fit_engine(two)
      assert_conditions(identical(a$calibration, b$calibration), identical(a$calibrator, b$calibrator),
        identical(a$folds, b$folds), identical(a$parameters, b$parameters),
        identical(workflows::extract_fit_engine(a$models$Full1), workflows::extract_fit_engine(b$models$Full1)),
        identical(workflows::extract_fit_engine(a$models$Full1)$x$id, data$id[1:48]),
        identical(workflows::extract_recipe(a$models$Full1)$steps, workflows::extract_recipe(b$models$Full1)$steps),
        identical(venn_abers_details(one, data), venn_abers_details(two, data)))
    }
    heldout <- rsample::manual_rset(list(holdout), "Test")
    altered_outer <- rsample::manual_rset(list(altered), "Test")
    first <- tune::fit_resamples(ordered_model, heldout, metrics = scoring, control = control)
    second <- tune::fit_resamples(ordered_model, altered_outer, metrics = scoring, control = control)
    check(first, heldout); check(second, altered_outer)
    a <- workflows::extract_fit_engine(tune::collect_extracts(first)$.extracts[[1L]])
    b <- workflows::extract_fit_engine(tune::collect_extracts(second)$.extracts[[1L]])
    assert_conditions(identical(tune::collect_predictions(first)$.pred, tune::collect_predictions(second)$.pred),
      identical(a$calibration, b$calibration), identical(a$calibrator, b$calibrator),
      identical(a$folds, b$folds), identical(a$parameters, b$parameters),
      identical(workflows::extract_fit_engine(a$models$Full1), workflows::extract_fit_engine(b$models$Full1)),
      identical(workflows::extract_recipe(a$models$Full1)$steps, workflows::extract_recipe(b$models$Full1)$steps))
    if (!is.null(option$m)) {
      for (change in c("missing", "duplicate", "config", "metric", "repeat", "truth")) {
        bad <- result
        if (change == "missing") bad$.predictions[[1L]] <- bad$.predictions[[1L]][-1L, ]
        if (change == "duplicate") bad$.predictions[[1L]] <- rbind(bad$.predictions[[1L]], bad$.predictions[[1L]][1L, ])
        if (change == "config") bad$.predictions[[1L]]$.config <- "wrong"
        if (change == "metric") bad$.metrics[[1L]] <- bad$.metrics[[1L]][-1L, ]
        if (change == "repeat") bad$id2[1L] <- "wrong"
        if (change == "truth") bad$.predictions[[1L]]$truth[1L] <- NA_real_
        error(check(bad, outer))
      }
      fail_metric <- yardstick::new_numeric_metric(function(data, ...) stop("metric failed"), direction = "minimize")
      bad <- suppressWarnings(tune::fit_resamples(ordered_model, outer,
        metrics = yardstick::metric_set(yardstick::rmse, fail_metric), control = control))
      stopifnot(nrow(tune::collect_notes(bad)) > 0L)
      error(check(bad, outer))
      partial <- rsample::manual_rset(list(holdout,
        rsample::make_splits(list(analysis = 1:36, assessment = 37:60), data)), c("First", "Second"))
      for (sizes in list(48L, c(36L, 48L))) {
        broken <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail_size = sizes))
        bad_model <- calibrated_workflow(broken, "truth", method = method, prediction_model = "full", splits = ordered)
        bad <- suppressWarnings(tune::fit_resamples(bad_model, partial, metrics = scoring, control = control))
        stopifnot(nrow(tune::collect_notes(bad)) == length(sizes))
        error(check(bad, partial))
      }
    }
  }
}

  testthat::expect_true(TRUE)
})

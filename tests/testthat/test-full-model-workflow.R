testthat::test_that("full-model-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")

packages <- c("parsnip", "workflows", "rsample", "recipes", "withr")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  data <- data.frame(x = sin(1:60), truth = factor(rep(c("yes", "no", "no"), 20)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (method in c("inductive", "cross")) {
    matched <- generics::fit(calibrated_workflow(base, "truth", method = method, seed = 42), data)
    full <- generics::fit(calibrated_workflow(base, "truth", method = method, seed = 42,
                                             prediction_model = "full"), data)
    source <- workflows::extract_fit_engine(matched)
    state <- workflows::extract_fit_engine(full)
    stopifnot(identical(state$calibrator, source$calibrator),
              identical(state$folds, source$folds), identical(names(state$models), "Full1"),
              state$parameters$effective$prediction_model == "full")
    full_base <- generics::fit(base, data)
    p <- predict(full_base, data, type = "prob")[[paste0(".pred_", state$event_class)]]
    input <- if (method == "cross") {
      matrix(p, nrow(data), length(state$folds), dimnames = list(NULL, names(state$folds)))
    } else p
    expected <- predict(source$calibrator, input)$probability
    stopifnot(max(abs(predict(full, data, type = "prob")[[paste0(".pred_", state$event_class)]] - expected)) <= 1e-10)
  }
}

if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  error <- function(expr, text) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"), grepl(text, conditionMessage(result), fixed = TRUE))
  }
  counts <- new.env()
  counts$fit <- counts$predict <- 0L
  counts$training <- list()
  transfer_fit <- function(x, y, fail_size = 0L, tie = FALSE, bad_size = 0L) {
    counts$fit <- counts$fit + 1L
    if (nrow(x) %in% fail_size) stop("instrumented training failure")
    counts$training[[length(counts$training) + 1L]] <- x
    list(x = x, y = y, levels = levels(y), tie = tie, bad_size = bad_size, draw = runif(1))
  }
  transfer_prob <- function(object, new_data) {
    counts$predict <- counts$predict + 1L
    stopifnot(nrow(new_data) > 0L)
    # Training-domain size deliberately changes the new model's probabilities.
    p <- if (object$tie) rep(.5, nrow(new_data)) else plogis(new_data$p + nrow(object$x) / 10 - 3)
    if (nrow(object$x) %in% object$bad_size) p[] <- NA_real_
    stats::setNames(data.frame(1 - p, p), paste0(".pred_", rev(object$levels)))
  }
  register_reference_engine("va_transfer", "classification", "transfer_fit", "transfer_prob", "prob")
  spec <- parsnip::new_model_spec("va_transfer", args = list(), eng_args = NULL,
                                 mode = "classification", method = NULL, engine = "reference")
  data <- data.frame(p = seq(-3, 3, length.out = 40), id = 1:40,
                     x = c(NA, sin(2:32), 10000 + 1:8),
                     truth = factor(rep(c("a b", "a.b", "a.b", "a b"), 10), levels = c("a b", "a.b")))
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ ., data) |>
      recipes::step_impute_median(x) |> recipes::step_normalize(x)) |>
    workflows::add_model(spec)
  make_plan <- function(data, method) {
    held <- if (method == "inductive") list(33:40) else list(1:8, 9:20, 21:40)
    splits <- lapply(held, function(h) rsample::make_splits(
      list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data))
    if (method == "inductive") splits[[1L]] else rsample::manual_rset(splits, c("small", "middle", "large"))
  }
  for (method in c("inductive", "cross")) for (event in c("first", "second")) {
    for (precision in list(NULL, 0, 2, -1)) {
      plan <- make_plan(data, method)
      fits <- states <- list()
      k <- if (method == "inductive") 1L else 3L
      for (mode in c("matched", "full")) {
        before <- counts$fit
        set.seed(77)
        rng <- .Random.seed
        fits[[mode]] <- generics::fit(calibrated_workflow(base, "truth", method = method,
          splits = plan, prediction_model = mode, event_level = event, precision = precision,
          loss = "brier", seed = 13), data)
        state <- states[[mode]] <- workflows::extract_fit_engine(fits[[mode]])
        assert_conditions(counts$fit == before + k + (mode == "full"), identical(rng, .Random.seed),
                  length(state$models) == if (mode == "full") 1L else k,
                  state$parameters$requested$prediction_model == mode,
                  state$parameters$explicit[["prediction_model"]],
                  identical(state$parameters$requested$precision, precision),
                  identical(names(state$calibration), names(state$folds)))
      }
      source <- states$matched
      state <- states$full
      assert_conditions(identical(source$calibrator, state$calibrator), identical(source$calibration, state$calibration),
                identical(source$folds, state$folds), identical(state$rng_kind, RNGkind()),
                all(c("R", "vennabers", "workflows", "recipes") %in% names(state$versions)))
      full_engine <- workflows::extract_fit_engine(state$models$Full1)
      stopifnot(identical(full_engine$x$id, data$id), identical(as.vector(full_engine$y), as.vector(data$truth)))
      for (id in names(source$folds)) {
        f <- source$folds[[id]]
        learned <- workflows::extract_recipe(source$models[[id]])
        x <- data$x[f$analysis]
        x[is.na(x)] <- median(x, na.rm = TRUE)
        stopifnot(abs(learned$steps[[2L]]$means[["x"]] - mean(x)) < 1e-10,
                  identical(state$calibration[[id]]$row, f$assessment),
                  identical(state$calibration[[id]]$truth, data$truth[f$assessment] == state$event_class),
                  identical(state$calibration[[id]]$probability,
                    predict(source$models[[id]], data[f$assessment, ], type = "prob")[[paste0(".pred_", state$event_class)]]))
      }
      x <- data$x
      x[is.na(x)] <- median(x, na.rm = TRUE)
      learned <- workflows::extract_recipe(state$models$Full1)
      stopifnot(abs(learned$steps[[2L]]$means[["x"]] - mean(x)) < 1e-10)
      if (method == "inductive") stopifnot(abs(learned$steps[[2L]]$means[["x"]] -
        workflows::extract_recipe(source$models$Inductive1)$steps[[2L]]$means[["x"]]) > 1000)
      new <- data[c(1, 8, 15, 22, 29, 33, 35, 38, 40), ]
      p <- predict(state$models$Full1, new, type = "prob")[[paste0(".pred_", state$event_class)]]
      for (loss in c("log", "brier")) {
        input <- if (method == "cross") matrix(p, nrow(new), k, dimnames = list(NULL, names(state$folds))) else p
        expected <- predict(source$calibrator, input, loss = loss)$probability
        before_fit <- counts$fit
        before <- counts$predict
        details <- venn_abers_details(fits$full, new, loss = loss)
        assert_conditions(counts$predict == before + 1L,
                  all(details$details$prediction_model == "full"), all(details$details$model_source == "Full1"),
                  identical(unique(details$details$fold), names(state$folds)),
                  max(abs(details$prediction[[paste0(".pred_", state$event_class)]] - expected)) <= 1e-10)
        before <- counts$predict
        prob <- predict(fits$full, new, type = "prob", loss = loss)
        cls <- predict(fits$full, new, type = "class", loss = loss)$.pred_class
        stopifnot(counts$predict == before + 2L, identical(details$prediction, prob),
                  identical(names(prob), paste0(".pred_", levels(data$truth))),
                  identical(as.character(cls), ifelse(expected >= .5, state$event_class, setdiff(state$class_levels, state$event_class))))
        lower <- matrix(details$details$lower, nrow(new))
        upper <- matrix(details$details$upper, nrow(new))
        points <- matrix(details$details$probability, nrow(new))
        reconstructed <- if (loss == "log") plogis(rowMeans(log(upper)) - rowMeans(log1p(-lower))) else rowMeans(points)
        stopifnot(max(abs(expected - reconstructed)) <= 1e-10)
        for (n in 0:1) {
          before <- counts$predict
          small <- venn_abers_details(fits$full, new[seq_len(n), ], loss = loss)
          small_p <- predict(fits$full, new[seq_len(n), ], type = "prob", loss = loss)
          small_c <- predict(fits$full, new[seq_len(n), ], type = "class", loss = loss)
          assert_conditions(counts$predict == before + 3L * n, identical(names(small$details), names(details$details)),
                    identical(small$prediction, small_p), nrow(small_p) == n, nrow(small$details) == k * n,
                    is.integer(small$details$row), is.character(small$details$model_source),
                    identical(levels(small_c$.pred_class), levels(data$truth)))
        }
        singles <- do.call(rbind, lapply(seq_len(nrow(new)), function(i) predict(fits$full, new[i, ], type = "prob", loss = loss)))
        stopifnot(identical(prob, singles), counts$fit == before_fit,
                  workflows::extract_fit_engine(fits$full)$loss == "brier")
        if (is.null(precision)) {
          stopifnot(max(abs(as.matrix(prob) - as.matrix(predict(fits$matched, new, type = "prob", loss = loss)))) > 1e-5)
          if (method == "cross") {
            pooled <- do.call(rbind, state$calibration)
            wrong <- predict(venn_abers(pooled$truth, pooled$probability), p, loss = loss)$probability
            stopifnot(max(abs(expected - wrong)) > 1e-5)
          }
        }
      }
      stopifnot(any(grepl("full-model transfer", capture.output(print(fits$full)), fixed = TRUE)),
                any(grepl("Full1", capture.output(print(state)), fixed = TRUE)))
    }
  }
  before <- counts$fit
  for (bad in list(NULL, "f", "FULL", NA_character_, c("full", "matched"), TRUE)) {
    error(calibrated_workflow(base, "truth", prediction_model = bad), "prediction_model")
  }
  error(calibrated_workflow(base, "truth", prediction_m = "full"), "prediction_m")
  reg <- workflows::workflow() |> workflows::add_formula(truth ~ p) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  stopifnot(inherits(calibrated_workflow(reg, "truth", prediction_model = "full"), "workflow"))
  error(predict(fits$full, data, prediction_model = "matched"), "prediction_model")
  error(venn_abers_details(fits$full, data, prediction_model = "matched"), "prediction_model")
  error(calibrate_workflow(state$models$Full1, data, "truth", prediction_model = "full"), "prediction_model")
  stopifnot(counts$fit == before)
  for (method in c("inductive", "cross")) {
    plan <- make_plan(data, method)
    for (fail_size in c(40L, if (method == "inductive") 32L else 20L)) {
      bad <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail_size = fail_size))
      error(generics::fit(calibrated_workflow(bad, "truth", method = method, splits = plan,
        prediction_model = "full"), data), "base workflow training")
    }
    bad <- workflows::update_model(base, parsnip::set_engine(spec, "reference", bad_size = 40L))
    error(generics::fit(calibrated_workflow(bad, "truth", method = method, splits = plan,
      prediction_model = "full"), data), "Full1 base workflow must return finite probabilities")
    tied <- workflows::update_model(base, parsnip::set_engine(spec, "reference", tie = TRUE))
    for (event in c("first", "second")) {
      fit <- generics::fit(calibrated_workflow(tied, "truth", method = method, seed = 42,
        prediction_model = "full", event_level = event), data)
      for (loss in c("log", "brier")) stopifnot(all(predict(fit, data, type = "prob", loss = loss) == .5),
        all(predict(fit, data, type = "class", loss = loss)$.pred_class == levels(data$truth)[if (event == "first") 1L else 2L]))
    }
    # Same prepared plan binding and fold validation as matching models, before learning.
    before <- counts$fit
    error(generics::fit(calibrated_workflow(base, "truth", method = method, splits = plan,
      prediction_model = "full"), data[40:1, ]), "data binding")
    stopifnot(counts$fit == before)
  }

  reference <- read.csv(file.path("extdata", "full-model-reference.csv"))
  stopifnot(nrow(reference) == 224L)
  maximum_error <- 0
  for (method in c("inductive", "cross")) for (precision in list(NULL, 0, 2, -1)) {
    fitted <- generics::fit(calibrated_workflow(base, "truth", method = method, prediction_model = "full",
      splits = make_plan(data, method), precision = precision), data)
    for (loss in c("log", "brier")) {
      expected <- reference[reference$method == method & reference$loss == loss &
        reference$precision == if (is.null(precision)) "none" else as.character(precision), ]
      new <- data[1:7, ]
      new$p <- expected$x[expected$fold == 1L]
      actual <- venn_abers_details(fitted, new, loss = loss)
      stopifnot(identical(actual$details$row, expected$row))
      maximum_error <- max(maximum_error,
        abs(as.matrix(actual$details[c("lower", "upper", "probability")]) - as.matrix(expected[c("lower", "upper", "probability")])),
        abs(actual$prediction[[".pred_a b"]] - expected$combined[expected$fold == 1L]))
    }
  }
  stopifnot(maximum_error <= 1e-10)
  cat("Full-transfer Python component maximum absolute error:", maximum_error, "\n")

  if (all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
    # Reuse the shipped completeness audit, including repeat/row/configuration keys.
    example <- new.env()
    sys.source(system.file("examples", "calibrated-resampling.R", package = "vennabers"), example)
    check <- example$check_evaluation
    scoring <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
    outer <- withr::with_seed(201, rsample::vfold_cv(data, v = 2, repeats = 2))
    holdout <- rsample::make_splits(list(analysis = 1:32, assessment = 33:40), data)
    local_plan <- function(data, method) {
      h <- seq.int(nrow(data) - 3L, nrow(data))
      split <- rsample::make_splits(list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data)
      if (method == "inductive") split else rsample::manual_rset(list(split,
        rsample::make_splits(list(analysis = h, assessment = setdiff(seq_len(nrow(data)), h)), data)), c("A", "B"))
    }
    for (method in c("inductive", "cross")) for (event in c("first", "second")) for (loss in c("log", "brier")) {
      prepared <- lapply(outer$splits, function(s) local_plan(rsample::analysis(s), method))
      resolver <- local_calibration_splits(prepared, names(data))
      model <- calibrated_workflow(base, "truth", method = method, splits = resolver,
        prediction_model = "full", seed = 42, event_level = event, loss = loss)
      control <- tune::control_resamples(save_pred = TRUE, extract = identity, event_level = event, allow_par = FALSE)
      before <- counts$fit
      result <- tune::fit_resamples(model, outer, metrics = scoring, control = control)
      checked <- check(result, outer, event)
      stopifnot(counts$fit == before + nrow(outer) * if (method == "cross") 3L else 2L)
      fits <- tune::collect_extracts(result)$.extracts
      for (i in seq_along(fits)) {
        state <- workflows::extract_fit_engine(fits[[i]])
        local <- rsample::analysis(outer$splits[[i]])
        held <- rsample::assessment(outer$splits[[i]])
        engine <- workflows::extract_fit_engine(state$models$Full1)
        stopifnot(identical(engine$x$id, local$id), !any(engine$x$id %in% held$id),
                  identical(as.vector(engine$y), as.vector(local$truth)), length(state$models) == 1L)
        rows <- checked$predictions$id == outer$id[i] & checked$predictions$id2 == outer$id2[i]
        direct <- predict(fits[[i]], held, type = "prob")
        stopifnot(identical(as.matrix(checked$predictions[rows, names(direct)]), as.matrix(direct)))
      }
      # Assessment targets are inaccessible to either inner model domain.
      changed <- data
      changed$truth[33:40] <- rev(changed$truth[33:40])
      changed$truth[33L] <- if (changed$truth[33L] == "a b") "a.b" else "a b"
      altered <- rsample::make_splits(list(analysis = 1:32, assessment = 33:40), changed)
      model <- calibrated_workflow(base, "truth", method = method,
        splits = function(data) local_plan(data, method), prediction_model = "full", seed = 42,
        event_level = event, loss = loss)
      before <- counts$fit
      final <- tune::last_fit(model, holdout, metrics = scoring,
        control = tune::control_last_fit(event_level = event))
      stopifnot(counts$fit == before + if (method == "cross") 3L else 2L)
      final_p <- tune::collect_predictions(final)
      final_m <- tune::collect_metrics(final)
      assert_conditions(nrow(tune::collect_notes(final)) == 0L,
                identical(final_p$.row, 33:40), !anyDuplicated(final_p$.row),
                all(final_p$.config == "pre0_mod0_post0"),
                identical(final_p$truth, data$truth[33:40]),
                setequal(final_m$.metric, c("brier_class", "mn_log_loss")),
                nrow(final_m) == 2L, all(is.finite(final_m$.estimate)),
                all(final_m$.config == "pre0_mod0_post0"))
      mutated <- tune::last_fit(model, altered, metrics = scoring,
        control = tune::control_last_fit(event_level = event))
      state <- workflows::extract_fit_engine(tune::extract_workflow(final))
      other <- workflows::extract_fit_engine(tune::extract_workflow(mutated))
      assert_conditions(identical(state$calibration, other$calibration), identical(state$calibrator, other$calibrator),
                identical(state$folds, other$folds),
                identical(workflows::extract_fit_engine(state$models$Full1), workflows::extract_fit_engine(other$models$Full1)),
                identical(workflows::extract_fit_engine(state$models$Full1)$x$id, data$id[1:32]),
                identical(predict(tune::extract_workflow(final), data, type = "prob"),
                          predict(tune::extract_workflow(mutated), data, type = "prob")))
      if (event == "first" && loss == "log") {
        for (change in c("missing", "duplicate", "config", "metric", "repeat")) {
          bad <- result
          if (change == "missing") bad$.predictions[[1L]] <- bad$.predictions[[1L]][-1L, ]
          if (change == "duplicate") bad$.predictions[[1L]] <- rbind(bad$.predictions[[1L]], bad$.predictions[[1L]][1L, ])
          if (change == "config") bad$.predictions[[1L]]$.config <- "wrong"
          if (change == "metric") bad$.metrics[[1L]] <- bad$.metrics[[1L]][-1L, ]
          if (change == "repeat") bad$id2[1L] <- "wrong"
          error(check(bad, outer, event), "")
        }
        partial_outer <- rsample::manual_rset(list(holdout,
          rsample::make_splits(list(analysis = 1:20, assessment = 21:40), data)), c("First", "Second"))
        # Fail only the extra full fit, with every source model still trainable.
        for (sizes in list(32L, c(20L, 32L))) {
          bad_base <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail_size = sizes))
          bad_model <- calibrated_workflow(bad_base, "truth", method = method,
            splits = function(data) local_plan(data, method), prediction_model = "full")
          bad <- suppressWarnings(tune::fit_resamples(bad_model, partial_outer, metrics = scoring, control = control))
          stopifnot(nrow(tune::collect_notes(bad)) == length(sizes))
          error(check(bad, partial_outer), "")
          if (length(sizes) == 1L) stopifnot(nrow(tune::collect_predictions(bad)) == 20L)
        }
      }
    }
  }
}

  testthat::expect_true(TRUE)
})

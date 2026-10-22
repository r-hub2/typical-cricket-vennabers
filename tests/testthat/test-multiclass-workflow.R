testthat::test_that("multiclass-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")

if (all(vapply(c("parsnip", "workflows", "rsample", "nnet"),
               requireNamespace, logical(1), quietly = TRUE))) {
  # Both explicit multiclass contracts use the same fitted model and holdout.
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  split <- withr::with_seed(32, rsample::initial_split(data, .8, strata = truth, pool = 0))
  for (method in c("ovo", "ovr")) for (mode in c("matched", "full")) {
    fit <- generics::fit(calibrated_workflow(base, "truth", classification = method,
                                            prediction_model = mode, splits = split, seed = 32), data)
    state <- workflows::extract_fit_engine(fit)
    if (mode == "full") {
      stopifnot(identical(names(state$models), "Full1"),
        all(venn_abers_details(fit, data)$details$prediction_model == mode))
      next
    }
    probabilities <- function(new_data) {
      p <- as.matrix(predict(state$models$Inductive1, new_data, type = "prob"))
      colnames(p) <- sub("^\\.pred_", "", colnames(p))
      p
    }
    held <- rsample::assessment(split)
    manual <- venn_abers(held$truth, probabilities(held), method = method)
    stopifnot(max(abs(as.matrix(predict(fit, data, type = "prob")) -
                       predict(manual, probabilities(data)))) <= 1e-10)
  }
}

if (all(vapply(c("parsnip", "workflows", "rsample", "recipes"),
               requireNamespace, logical(1), quietly = TRUE))) {
  error <- function(expr, text) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"), grepl(text, conditionMessage(result), fixed = TRUE))
    conditionMessage(result)
  }
  close <- function(actual, expected) {
    stopifnot(identical(dim(actual), dim(expected)), all(is.finite(actual)),
              all(abs(actual - expected) <= 1e-10))
  }
  # Public reference engine: known outputs, original names, and observable inputs.
  audit <- new.env()
  audit$fits <- audit$predictions <- 0L
  multi_fit <- function(x, y, output = "reverse", fail = FALSE) {
    audit$fits <- audit$fits + 1L
    if (fail || length(unique(y)) < 2L) stop("untrainable reference model")
    list(x = x, y = y, levels = levels(y), output = output, draw = runif(1))
  }
  multi_prob <- function(object, new_data) {
    audit$predictions <- audit$predictions + 1L
    stopifnot(nrow(new_data) > 0L)
    p <- as.data.frame(new_data[paste0("p", seq_along(object$levels))])
    names(p) <- paste0(".pred_", object$levels)
    if (object$output == "missing") return(p[-1L])
    if (object$output == "duplicate") names(p)[2L] <- names(p)[1L]
    if (object$output == "unknown") names(p)[1L] <- ".pred_unknown"
    if (object$output == "nonfinite") p[1L, 1L] <- Inf
    if (object$output == "nonnumeric") p[[1L]] <- as.character(p[[1L]])
    if (object$output == "range") p[1L, 1L] <- -0.1
    if (object$output == "sum") p <- p * .9
    if (object$output == "tolerance") p <- p * (1 + 5e-9)
    p[rev(seq_along(p))]
  }
  register_reference_engine("va_multi", "classification", "multi_fit", "multi_prob", "prob")
  spec <- parsnip::new_model_spec("va_multi", args = list(), eng_args = NULL,
                                  mode = "classification", method = NULL, engine = "reference")
  probabilities <- function(model, data) {
    p <- as.matrix(predict(model, data, type = "prob"))
    colnames(p) <- sub("^\\.pred_", "", colnames(p))
    p
  }
  for (k in 2:4) {
    classes <- c("z", "a b", "a.b", "A")[seq_len(k)]
    p <- matrix(abs(sin(seq_len(120 * k))), 120, k)
    p[seq(1, 120, 7), 1:2] <- 0
    if (k == 2) p[seq(1, 120, 7), 2L] <- 1
    p <- p / rowSums(p)
    data <- as.data.frame(p)
    names(data) <- paste0("p", seq_len(k))
    data$id <- 1:120
    data$x <- sin(1:120)
    data$truth <- factor(rep(classes, length.out = 120), levels = classes)
    base <- workflows::workflow() |> workflows::add_variables(outcomes = truth, predictors = tidyselect::everything()) |>
      workflows::add_model(spec)
    new <- data[1:7, ]
    for (method in c("ovo", "ovr")) for (precision in list(NULL, 0, 2, -1)) {
      before <- audit$fits
      set.seed(700)
      rng <- .Random.seed
      fit <- generics::fit(calibrated_workflow(base, "truth", classification = method,
                                               seed = 32, precision = precision, loss = "brier"), data)
      stopifnot(identical(rng, .Random.seed), audit$fits == before + 1L)
      state <- workflows::extract_fit_engine(fit)
      split <- state$folds$Inductive1
      assert_conditions(identical(state$class_levels, classes), is.null(state$event_class),
                is.null(state$parameters$effective$event_level), state$calibrator$method == method,
                length(state$models) == 1L,
                identical(workflows::extract_fit_engine(state$models$Inductive1)$x$id, data$id[split$analysis]))
      held <- data[split$assessment, ]
      manual <- venn_abers(held$truth, probabilities(state$models$Inductive1, held),
                           precision = precision, method = method)
      input <- probabilities(state$models$Inductive1, new)
      for (loss in c("log", "brier")) {
        expected <- predict(manual, input, type = "details", loss = loss)
        details <- venn_abers_details(fit, new, loss = loss)
        stopifnot(identical(details$prediction, predict(fit, new, type = "prob", loss = loss)),
                  identical(names(details$prediction), paste0(".pred_", classes)),
                  identical(predict(fit, new, type = "class", loss = loss)$.pred_class,
                            predict(manual, input, type = "class", loss = loss)))
        close(as.matrix(details$prediction), expected$prediction)
        d <- details$details
        part <- if (method == "ovo") expected$pairs else expected$classes
        close(as.matrix(d[c("lower", "upper", "probability")]),
              as.matrix(part[c("lower", "upper", "probability")]))
        stopifnot(all(d$fold == "Inductive1"), all(d$method == method),
                  identical(d$event_class, part$event))
        reconstructed <- matrix(if (method == "ovo") 2 - k else 0, nrow(new), k,
                                dimnames = list(NULL, classes))
        for (row in seq_len(nrow(d))) {
          i <- d$row[row]
          q <- d$probability[row]
          if (method == "ovo") {
            stopifnot(d$class_1[row] == part$non_event[row], d$class_2[row] == part$event[row])
            reconstructed[i, d$class_1[row]] <- reconstructed[i, d$class_1[row]] + 1 / (1 - q)
            reconstructed[i, d$class_2[row]] <- reconstructed[i, d$class_2[row]] + 1 / q
          } else {
            stopifnot(d$class[row] == d$event_class[row])
            reconstructed[i, d$class[row]] <- q
          }
        }
        if (method == "ovo") reconstructed <- 1 / reconstructed
        close(reconstructed / rowSums(reconstructed), expected$prediction)
        for (n in 0:1) {
          small <- new[seq_len(n), , drop = FALSE]
          calls <- audit$predictions
          small_details <- venn_abers_details(fit, small, loss = loss)
          prob <- predict(fit, small, type = "prob", loss = loss)
          cls <- predict(fit, small, type = "class", loss = loss)
          if (!n) stopifnot(audit$predictions == calls)
          assert_conditions(identical(small_details$prediction, prob), nrow(prob) == n,
                    nrow(small_details$details) == n * if (method == "ovo") choose(k, 2) else k,
                    identical(levels(cls$.pred_class), classes),
                    identical(names(small_details$details), names(d)),
                    identical(vapply(small_details$details, typeof, ""), vapply(d, typeof, "")))
        }
      }
      before <- audit$fits
      batch <- predict(fit, new, type = "prob")
      singles <- do.call(rbind, lapply(seq_len(nrow(new)), function(i) predict(fit, new[i, ], type = "prob")))
      stopifnot(identical(batch, singles), audit$fits == before,
                identical(batch, predict(fit, new, type = "prob", loss = "brier")),
                workflows::extract_fit_engine(fit)$loss == "brier")
      plan <- rsample::make_splits(split[c("analysis", "assessment")], data)
      resolver <- local_calibration_splits(list(plan, plan), names(data))
      explicit <- generics::fit(calibrated_workflow(base, "truth", classification = method,
        splits = resolver, seed = 32, precision = precision, loss = "brier"), data)
      stopifnot(identical(batch, predict(explicit, new, type = "prob")))
    }
  }

  for (method in c("ovo", "ovr")) {
    # Calibration classes are checked before model learning; no split repair.
    before <- audit$fits
    for (name in c("class", ".pred_yes")) {
      bad <- data
      levels(bad$truth)[1L] <- name
      error(generics::fit(calibrated_workflow(base, "truth", classification = method), bad), name)
      rec <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ ., bad)) |>
        workflows::add_model(spec)
      error(calibrated_workflow(rec, "truth", classification = method), name)
    }
    error(calibrated_workflow(base, "truth", classification = method, event_level = "first"), "event_level")
    held <- which(data$truth == classes[1L])
    plan <- rsample::make_splits(list(analysis = setdiff(1:120, held), assessment = held), data)
    message <- error(generics::fit(calibrated_workflow(base, "truth", classification = method, splits = plan), data), "missing")
    p <- as.matrix(data[held, paste0("p", seq_along(classes))])
    colnames(p) <- classes
    stopifnot(identical(message, error(venn_abers(data$truth[held], p, method = method), "missing")))
    absent <- data[data$truth != classes[1L], ]
    error(generics::fit(calibrated_workflow(base, "truth", classification = method), absent), "all classes")
    foreign <- data
    foreign$x <- -foreign$x
    error(generics::fit(calibrated_workflow(base, "truth", classification = method, splits = plan), foreign), "binding")
    stopifnot(audit$fits == before)
    one_class <- which(data$truth == classes[1L])[1:4]
    plan <- rsample::make_splits(list(analysis = one_class, assessment = setdiff(1:120, one_class)), data)
    error(generics::fit(calibrated_workflow(base, "truth", classification = method, splits = plan), data),
          "Inductive1 base workflow training")
    for (output in c("missing", "duplicate", "unknown", "nonfinite", "nonnumeric", "range", "sum")) {
      bad <- workflows::update_model(base, parsnip::set_engine(spec, "reference", output = output))
      error(generics::fit(calibrated_workflow(bad, "truth", classification = method, seed = 32), data), "probabilit")
    }
    # Within-tolerance row sums are normalized before calibration and prediction.
    near <- data
    near[1:4] <- (near[1:4] + .1) / 1.4
    noisy_base <- workflows::update_model(base, parsnip::set_engine(spec, "reference", output = "tolerance"))
    noisy <- generics::fit(calibrated_workflow(noisy_base, "truth", classification = method, seed = 32), near)
    s <- workflows::extract_fit_engine(noisy)
    held <- near[s$folds$Inductive1$assessment, ]
    core <- venn_abers(held$truth, probabilities(s$models$Inductive1, held), method = method)
    close(as.matrix(predict(noisy, near, type = "prob")), predict(core, probabilities(s$models$Inductive1, near)))
    # Exactly balanced two-class ties favor the first level for either loss.
    tied <- data.frame(p1 = rep(.5, 40), p2 = rep(.5, 40),
                        truth = factor(rep(c("a b", "a.b"), 20), levels = c("a b", "a.b")))
    tied_base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |> workflows::add_model(spec)
    fit <- generics::fit(calibrated_workflow(tied_base, "truth", classification = method, seed = 32), tied)
    for (loss in c("log", "brier")) {
      stopifnot(all(predict(fit, tied, type = "prob", loss = loss) == .5),
                all(predict(fit, tied, type = "class", loss = loss)$.pred_class == "a b"))
    }
  }
}

if (exists("multi_fit") && all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
  example <- new.env()
  sys.source(system.file("examples", "calibrated-resampling.R", package = "vennabers"), example)
  check <- example$check_evaluation
  scoring <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
  control <- tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine,
                                     allow_par = FALSE, parallel_over = "resamples", pkgs = "vennabers")
  run <- function(model, outer) withr::with_seed(401, tune::fit_resamples(model, outer, metrics = scoring, control = control))
  for (k in 2:4) {
    classes <- c("z", "a b", "a.b", "A")[seq_len(k)]
    p <- matrix(abs(sin(seq_len(120 * k))), 120, k)
    data <- as.data.frame(p / rowSums(p))
    names(data) <- paste0("p", seq_len(k))
    data$id <- 1:120
    data$x <- sin(1:120)
    data$truth <- factor(rep(classes, length.out = 120), levels = classes)
    rec <- recipes::recipe(truth ~ ., data) |> recipes::step_normalize(x)
    base <- workflows::workflow() |> workflows::add_recipe(rec) |> workflows::add_model(spec)
    outer <- withr::with_seed(400, rsample::vfold_cv(data, v = 3, repeats = 2, strata = truth))
    for (method in c("ovo", "ovr")) for (loss in c("log", "brier")) {
      model <- calibrated_workflow(base, "truth", classification = method, seed = 32, loss = loss)
      before <- audit$fits
      result <- run(model, outer)
      checked <- check(result, outer, classification = method)
      states <- tune::collect_extracts(result)$.extracts
      stopifnot(audit$fits == before + nrow(outer), length(states) == nrow(outer))
      prepared <- vector("list", nrow(outer))
      for (i in seq_len(nrow(outer))) {
        local_data <- rsample::analysis(outer$splits[[i]])
        assessment <- rsample::assessment(outer$splits[[i]])
        s <- states[[i]]
        f <- s$folds$Inductive1
        engine <- workflows::extract_fit_engine(s$models$Inductive1)
        learned <- workflows::extract_recipe(s$models$Inductive1)
        assert_conditions(identical(engine$x$id, local_data$id[f$analysis]),
                  identical(as.vector(engine$y), as.vector(local_data$truth[f$analysis])),
                  !any(engine$x$id %in% assessment$id),
                  !any(local_data$id[f$assessment] %in% assessment$id), s$loss == loss,
                  abs(learned$steps[[1L]]$means[["x"]] - mean(local_data$x[f$analysis])) <= 1e-10)
        held <- local_data[f$assessment, ]
        core <- venn_abers(held$truth, probabilities(s$models$Inductive1, held), method = method)
        expected <- predict(core, probabilities(s$models$Inductive1, assessment), loss = loss)
        ids <- setdiff(names(outer), "splits")
        rows <- Reduce(`&`, lapply(ids, function(id) checked$predictions[[id]] == outer[[id]][i]))
        actual <- checked$predictions[rows, ]
        close(as.matrix(actual[paste0(".pred_", classes)]), expected)
        stopifnot(abs(yardstick::brier_class_vec(assessment$truth,
                    if (k == 2) expected[, 1L] else expected) - brier_score(assessment$truth, expected) / 2) <= 1e-10)
        prepared[[i]] <- rsample::make_splits(f[c("analysis", "assessment")], local_data)
      }
      resolver <- local_calibration_splits(prepared, names(data))
      seen <- list()
      callback <- function(data) { seen[[length(seen) + 1L]] <<- data$id; resolver(data) }
      planned <- run(calibrated_workflow(base, "truth", classification = method,
                     splits = callback, seed = 32, loss = loss), outer)
      check(planned, outer, classification = method)
      stopifnot(identical(seen, lapply(outer$splits, function(s) rsample::analysis(s)$id)),
                identical(tune::collect_predictions(result), tune::collect_predictions(planned)),
                identical(tune::collect_metrics(result), tune::collect_metrics(planned)))
    }
  }
  for (method in c("ovo", "ovr")) {
    # Modify only the outer assessment labels: fit, preprocessing and predictions
    # must be identical while scores change. Large calibration x cannot leak.
    make_outer <- function(data) rsample::manual_rset(list(rsample::make_splits(
      list(analysis = 1:80, assessment = 81:120), data)), "Test")
    changed <- data
    changed$truth[81:120] <- classes[1L]
    model <- calibrated_workflow(base, "truth", classification = method, seed = 32, loss = "brier")
    before <- run(model, make_outer(data))
    after <- run(model, make_outer(changed))
    check(before, make_outer(data), classification = method)
    check(after, make_outer(changed), classification = method)
    s1 <- tune::collect_extracts(before)$.extracts[[1L]]
    s2 <- tune::collect_extracts(after)$.extracts[[1L]]
    assert_conditions(identical(s1[setdiff(names(s1), "models")], s2[setdiff(names(s2), "models")]),
              identical(workflows::extract_fit_engine(s1$models$Inductive1), workflows::extract_fit_engine(s2$models$Inductive1)),
              identical(workflows::extract_recipe(s1$models$Inductive1)$steps, workflows::extract_recipe(s2$models$Inductive1)$steps),
              identical(tune::collect_predictions(before)[paste0(".pred_", classes)],
                        tune::collect_predictions(after)[paste0(".pred_", classes)]),
              !identical(tune::collect_metrics(before)$mean, tune::collect_metrics(after)$mean))
    extreme <- data
    plan <- withr::with_seed(32, rsample::initial_split(extreme, .8, strata = truth, pool = 0))
    train <- as.integer(plan, data = "analysis")
    held <- as.integer(plan, data = "assessment")
    extreme$x[held] <- 1e10
    fitted <- generics::fit(model, extreme)
    learned <- workflows::extract_recipe(workflows::extract_fit_engine(fitted)$models$Inductive1)
    stopifnot(abs(learned$steps[[1L]]$means[["x"]] - mean(extreme$x[train])) <= 1e-10)
    failed_base <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail = TRUE))
    failed <- suppressWarnings(run(calibrated_workflow(failed_base, "truth", classification = method, seed = 32), make_outer(data)))
    stopifnot(nrow(tune::collect_notes(failed)) == 1L)
    error(check(failed, make_outer(data), classification = method), "")
  }
  # Metric conventions affect scoring only, including the two-column matrix.
  for (k in 2:4) {
    truth <- factor(seq_len(k))
    p <- diag(k)[, rev(seq_len(k)), drop = FALSE]
    colnames(p) <- levels(truth)
    original <- p
    stopifnot(is.infinite(log_loss(truth, p)),
              is.finite(yardstick::mn_log_loss_vec(truth, if (k == 2) p[, 1L] else p)),
              abs(yardstick::brier_class_vec(truth, if (k == 2) p[, 1L] else p) - brier_score(truth, p) / 2) <= 1e-10,
              identical(p, original))
  }
}

  testthat::expect_true(TRUE)
})

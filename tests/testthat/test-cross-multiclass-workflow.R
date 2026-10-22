testthat::test_that("cross-multiclass-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")

if (all(vapply(c("parsnip", "workflows", "rsample", "nnet"),
               requireNamespace, logical(1), quietly = TRUE))) {
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  for (method in c("ovo", "ovr")) {
    fit <- generics::fit(calibrated_workflow(base, "truth", method = "cross",
                                            classification = method, seed = 33), data)
    state <- workflows::extract_fit_engine(fit)
    probability <- function(model, data) {
      p <- as.matrix(predict(model, data, type = "prob"))
      colnames(p) <- sub("^\\.pred_", "", colnames(p))
      p
    }
    held <- lapply(state$folds, function(f) data[f$assessment, ])
    p <- Map(probability, state$models, held)
    truth <- do.call(c, lapply(held, `[[`, "truth"))
    ids <- rep(names(held), vapply(held, nrow, 0L))
    core <- cross_venn_abers(truth, do.call(rbind, p), ids, method = method)
    inputs <- lapply(state$models, probability, data = data)
    for (loss in c("log", "brier")) {
      details <- venn_abers_details(fit, data, loss = loss)
      expected <- predict(core, inputs, loss = loss)
      stopifnot(length(state$models) == 5L,
                all(abs(as.matrix(details$prediction) - expected) <= 1e-10),
                identical(details$prediction, predict(fit, data, type = "prob", loss = loss)),
                setequal(details$details$fold, names(state$models)))
    }
  }
}

if (all(vapply(c("parsnip", "workflows", "rsample", "recipes"),
               requireNamespace, logical(1), quietly = TRUE))) {
  close <- function(actual, expected) {
    stopifnot(identical(dim(actual), dim(expected)), all(is.finite(actual)),
              all(abs(actual - expected) <= 1e-10))
  }
  error <- function(expr, text) {
    e <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(e, "error"), grepl(text, conditionMessage(e), fixed = TRUE))
    conditionMessage(e)
  }
  audit <- new.env()
  audit$fits <- 0L
  audit$predictions <- list()
  cross_multi_fit <- function(x, y, output = "reverse", fail = FALSE) {
    audit$fits <- audit$fits + 1L
    if (fail || length(unique(y)) < 2L) stop("untrainable reference model")
    list(x = x, y = y, levels = levels(y), output = output, draw = runif(1))
  }
  cross_multi_prob <- function(object, new_data) {
    stopifnot(nrow(new_data) > 0L)
    audit$predictions[[length(audit$predictions) + 1L]] <- list(train = object$x$id, rows = new_data$id)
    p <- as.matrix(new_data[paste0("p", seq_along(object$levels))])
    # Different stored models yield different scores, preserving exact zero pairs.
    p <- p ^ (0.5 + mean(object$x$id) / 30)
    p <- as.data.frame(p / rowSums(p))
    names(p) <- paste0(".pred_", object$levels)
    if (object$output == "missing") return(p[-1L])
    if (object$output == "tolerance") p <- p * (1 + 5e-9)
    p[rev(seq_along(p))]
  }
  register_reference_engine("va_cross_multi", "classification", "cross_multi_fit", "cross_multi_prob", "prob")
  spec <- parsnip::new_model_spec("va_cross_multi", args = list(), eng_args = NULL,
                                  mode = "classification", method = NULL, engine = "reference")
  make_plan <- function(data, held = list(1:7, 8:24, 25:60)) {
    rsample::manual_rset(lapply(held, function(h) rsample::make_splits(
      list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data)),
      paste0("part", seq_along(held)))
  }
  probability <- function(model, data) {
    p <- as.matrix(predict(model, data, type = "prob"))
    colnames(p) <- sub("^\\.pred_", "", colnames(p))
    p
  }
  manual_inputs <- function(state, data, new) {
    held <- lapply(state$folds, function(f) data[f$assessment, , drop = FALSE])
    list(truth = do.call(c, lapply(held, `[[`, "truth")),
         p = do.call(rbind, Map(probability, state$models[names(held)], held)),
         fold = rep(names(held), vapply(held, nrow, 0L)),
         new = lapply(state$models, probability, data = new))
  }
  for (k in 2:4) {
    classes <- c("z", "a b", "a.b", "A")[seq_len(k)]
    p <- matrix(abs(sin(seq_len(60 * k))), 60, k)
    p[seq(1, 60, 7), 1:2] <- 0
    if (k == 2) p[seq(1, 60, 7), 2L] <- 1
    data <- as.data.frame(p / rowSums(p))
    names(data) <- paste0("p", seq_len(k))
    data$id <- 1:60
    data$x <- c(1000, NA, sin(3:60))
    data$truth <- factor(rep(classes, length.out = 60), levels = classes)
    rec <- recipes::recipe(truth ~ ., data) |> recipes::step_impute_median(x) |> recipes::step_normalize(x)
    base <- workflows::workflow() |> workflows::add_recipe(rec) |> workflows::add_model(spec)
    plan <- make_plan(data)
    new <- data[1:7, ]
    for (method in c("ovo", "ovr")) for (precision in list(NULL, 0, 2, -1)) {
      before <- audit$fits
      audit$predictions <- list()
      set.seed(700)
      rng <- .Random.seed
      fit <- generics::fit(calibrated_workflow(base, "truth", method = "cross", classification = method,
        splits = plan, seed = 33, precision = precision, loss = "brier"), data)
      stopifnot(identical(rng, .Random.seed), audit$fits == before + 3L, length(audit$predictions) == 3L)
      state <- workflows::extract_fit_engine(fit)
      assert_conditions(identical(names(state$models), plan$id), identical(state$class_levels, classes), is.null(state$event_class),
          is.null(state$parameters$effective$event_level), identical(state$parameters$requested$precision,
              precision))
      for (i in seq_along(state$folds)) {
        f <- state$folds[[i]]
        model <- workflows::extract_fit_engine(state$models[[i]])
        learned <- workflows::extract_recipe(state$models[[i]])
        x <- data$x[f$analysis]
        x[is.na(x)] <- median(x, na.rm = TRUE)
        assert_conditions(identical(model$x$id, data$id[f$analysis]), identical(as.vector(model$y), as.vector(data$truth[f$analysis])),
            identical(audit$predictions[[i]], list(train = data$id[f$analysis], rows = data$id[f$assessment])),
            abs(learned$steps[[1]]$medians[["x"]] - median(x)) <= 1e-10, abs(learned$steps[[2]]$means[["x"]] -
                mean(x)) <= 1e-10)
      }
      inputs <- manual_inputs(state, data, new)
      core <- cross_venn_abers(inputs$truth, inputs$p, inputs$fold, precision = precision, method = method)
      for (loss in c("log", "brier")) {
        expected <- predict(core, inputs$new, loss = loss, type = "details")
        actual <- venn_abers_details(fit, new, loss = loss)
        close(as.matrix(actual$prediction), expected$prediction)
        stopifnot(identical(actual$prediction, predict(fit, new, type = "prob", loss = loss)),
                  identical(names(actual$prediction), paste0(".pred_", classes)),
                  identical(predict(fit, new, type = "class", loss = loss)$.pred_class,
                            predict(core, inputs$new, type = "class", loss = loss)))
        d <- actual$details
        part <- if (method == "ovo") expected$pairs else expected$classes
        close(as.matrix(d[c("lower", "upper", "probability")]), as.matrix(part[c("lower", "upper", "probability")]))
        stopifnot(identical(d$row, part$row), identical(d$fold, part$fold),
                  identical(d$event_class, part$event), all(d$method == method))
        problems <- if (method == "ovo") combn(classes, 2, simplify = FALSE) else as.list(classes)
        reconstructed <- matrix(if (method == "ovo") 2 - k else 0, nrow(new), k, dimnames = list(NULL, classes))
        for (problem in problems) {
          rows <- if (method == "ovo") d$class_1 == problem[1] & d$class_2 == problem[2] else d$class == problem
          parts <- d[rows, ]
          q <- vapply(seq_len(nrow(new)), function(i) {
            f <- parts[parts$row == i, ]
            stopifnot(setequal(f$fold, plan$id), nrow(f) == 3L)
            if (loss == "log") plogis(mean(log(f$upper)) - mean(log1p(-f$lower))) else mean(f$probability)
          }, 0.0)
          event <- tail(problem, 1)
          if (method == "ovo") {
            reconstructed[, problem[1]] <- reconstructed[, problem[1]] + 1 / (1 - q)
            reconstructed[, event] <- reconstructed[, event] + 1 / q
          } else reconstructed[, event] <- q
        }
        if (method == "ovo") reconstructed <- 1 / reconstructed
        close(reconstructed / rowSums(reconstructed), expected$prediction)
        if (is.null(precision)) {
          swapped <- rev(inputs$new)
          names(swapped) <- names(inputs$new)
          stopifnot(max(abs(predict(core, swapped, loss = loss) - expected$prediction)) > 1e-5)
          if (k > 2L) {
            coupled <- lapply(plan$id, function(id) {
              h <- inputs$fold == id
              one <- venn_abers(inputs$truth[h], inputs$p[h, ], method = method)
              predict(one, inputs$new[[id]], loss = loss)
            })
            stopifnot(max(abs(Reduce(`+`, coupled) / 3 - expected$prediction)) > 1e-5)
          }
        }
        for (n in 0:1) {
          calls <- length(audit$predictions)
          small <- venn_abers_details(fit, new[seq_len(n), ], loss = loss)
          cls <- predict(fit, new[seq_len(n), ], type = "class", loss = loss)
          assert_conditions(nrow(small$prediction) == n, nrow(cls) == n, nrow(small$details) == n * 3 * length(problems),
              identical(levels(cls$.pred_class), classes), identical(names(small$details), names(d)), identical(vapply(small$details,
                  typeof, ""), vapply(d, typeof, "")))
          if (!n) stopifnot(length(audit$predictions) == calls)
        }
      }
      before <- audit$fits
      batch <- predict(fit, new, type = "prob")
      singles <- do.call(rbind, lapply(seq_len(nrow(new)), function(i) predict(fit, new[i, ], type = "prob")))
      stopifnot(identical(batch, singles), audit$fits == before,
                identical(batch, predict(fit, new, type = "prob", loss = "brier")),
                workflows::extract_fit_engine(fit)$loss == "brier")
    }
  }

  # A late empty pair must fail before even the first base fit.
  before <- audit$fits
  empty <- make_plan(data, list(2:30, 31:60, 1L))
  message <- error(generics::fit(calibrated_workflow(base, "truth", method = "cross", classification = "ovo",
    splits = empty), data), "pair [a b, a.b] in fold [part3]")
  held <- lapply(empty$splits, rsample::assessment)
  truth <- do.call(c, lapply(held, `[[`, "truth"))
  p <- as.matrix(do.call(rbind, held)[paste0("p", seq_along(classes))])
  colnames(p) <- classes
  ids <- rep(empty$id, vapply(held, nrow, 0L))
  stopifnot(identical(message, error(cross_venn_abers(truth, p, ids), "Empty calibration")))
  stopifnot(audit$fits == before)
}

if (exists("cross_multi_fit")) {
  for (method in c("ovo", "ovr")) {
    before <- audit$fits
    absent <- data[data$truth != classes[1L], ]
    error(generics::fit(calibrated_workflow(base, "truth", method = "cross", classification = method), absent), "all classes")
    for (name in c("class", ".pred_yes")) {
      bad <- data
      levels(bad$truth)[1L] <- name
      error(generics::fit(calibrated_workflow(base, "truth", method = "cross", classification = method), bad), name)
    }
    stopifnot(audit$fits == before)
    broken <- workflows::update_model(base, parsnip::set_engine(spec, "reference", output = "missing"))
    error(generics::fit(calibrated_workflow(broken, "truth", method = "cross", classification = method,
      splits = plan), data), "part1 base workflow must return")
    one_class <- which(data$truth == classes[1L])
    untrainable <- make_plan(data, list(setdiff(1:60, one_class), one_class))
    # OvR permits singleton/one-class calibration; base trainability is separate.
    if (method == "ovr") error(generics::fit(calibrated_workflow(base, "truth", method = "cross",
      classification = method, splits = untrainable), data), "part1 base workflow training")
    auto <- generics::fit(calibrated_workflow(base, "truth", method = "cross", classification = method,
      seed = 33), data)
    s <- workflows::extract_fit_engine(auto)
    prepared <- rsample::manual_rset(lapply(s$folds, function(f) {
      rsample::make_splits(f[c("analysis", "assessment")], data)
    }), names(s$folds))
    resolver <- local_calibration_splits(list(prepared, prepared), names(data))
    explicit <- generics::fit(calibrated_workflow(base, "truth", method = "cross", classification = method,
      splits = resolver, seed = 33), data)
    stopifnot(identical(venn_abers_details(auto, new), venn_abers_details(explicit, new)))
    # Legal incomplete folds: OvO needs at least one member of every pair;
    # OvR can even calibrate a singleton in a four-class domain.
    small <- if (method == "ovo") data[data$truth != classes[4L], ] else data
    small$truth <- droplevels(small$truth)
    small <- small[c(paste0("p", seq_len(nlevels(small$truth))), "id", "x", "truth")]
    small[seq_len(nlevels(small$truth))] <- small[seq_len(nlevels(small$truth))] + .1
    small[seq_len(nlevels(small$truth))] <- small[seq_len(nlevels(small$truth))] / rowSums(small[seq_len(nlevels(small$truth))])
    small_base <- workflows::workflow() |> workflows::add_variables(outcomes = truth, predictors = tidyselect::everything()) |>
      workflows::add_model(spec)
    held <- if (method == "ovo") 1:2 else 1L
    incomplete <- make_plan(small, list(held, setdiff(seq_len(nrow(small)), held)))
    # Singleton complements need an engine that accepts one-class training.
    if (method == "ovr") incomplete <- make_plan(small, list(1L, 2:30, 31:nrow(small)))
    fitted <- generics::fit(calibrated_workflow(small_base, "truth", method = "cross", classification = method,
      splits = incomplete), small)
    s <- workflows::extract_fit_engine(fitted)
    input <- manual_inputs(s, small, small)
    manual <- cross_venn_abers(input$truth, input$p, input$fold, method = method)
    close(as.matrix(predict(fitted, small, type = "prob")), predict(manual, input$new))
    # Cross ties use the first multiclass level, including a two-class domain.
    tied <- data.frame(p1 = rep(.5, 40), p2 = rep(.5, 40), id = 1:40,
      truth = factor(rep(c("a b", "a.b"), 20), levels = c("a b", "a.b")))
    tied_base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |> workflows::add_model(spec)
    fitted <- generics::fit(calibrated_workflow(tied_base, "truth", method = "cross", classification = method, seed = 33), tied)
    for (loss in c("log", "brier")) {
      stopifnot(all(predict(fitted, tied, type = "prob", loss = loss) == .5),
                all(predict(fitted, tied, type = "class", loss = loss)$.pred_class == "a b"))
    }
    single <- generics::fit(calibrated_workflow(tied_base, "truth", method = "cross", classification = method,
      splits = make_plan(tied, list(1L, 2:20, 21:40))), tied)
    stopifnot(nrow(predict(single, tied, type = "prob")) == 40L)
  }
}

if (exists("cross_multi_fit") && all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
  example <- new.env()
  sys.source(system.file("examples", "calibrated-resampling.R", package = "vennabers"), example)
  check <- example$check_evaluation
  control <- tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine,
    allow_par = FALSE, parallel_over = "resamples", pkgs = "vennabers")
  run <- function(model, outer) withr::with_seed(331, tune::fit_resamples(model, outer,
    metrics = yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss), control = control))
  for (method in c("ovo", "ovr")) for (loss in c("log", "brier")) {
    model <- calibrated_workflow(base, "truth", method = "cross", classification = method, seed = 33, loss = loss)
    outer <- withr::with_seed(330, rsample::vfold_cv(data, v = 5, strata = truth))
    before <- audit$fits
    result <- run(model, outer)
    checked <- check(result, outer, classification = method)
    states <- tune::collect_extracts(result)$.extracts
    stopifnot(audit$fits == before + 25L, length(states) == 5L)
    prepared <- vector("list", nrow(outer))
    for (i in seq_len(nrow(outer))) {
      local_data <- rsample::analysis(outer$splits[[i]])
      assessment <- rsample::assessment(outer$splits[[i]])
      s <- states[[i]]
      for (id in names(s$folds)) {
        f <- s$folds[[id]]
        engine <- workflows::extract_fit_engine(s$models[[id]])
        learned <- workflows::extract_recipe(s$models[[id]])
        x <- local_data$x[f$analysis]
        x[is.na(x)] <- median(x, na.rm = TRUE)
        assert_conditions(identical(engine$x$id, local_data$id[f$analysis]), identical(as.vector(engine$y), as.vector(local_data$truth[f$analysis])),
            !any(engine$x$id %in% assessment$id), !any(local_data$id[f$assessment] %in% assessment$id), abs(learned$steps[[2]]$means[["x"]] -
                mean(x)) <= 1e-10)
      }
      inputs <- manual_inputs(s, local_data, assessment)
      core <- cross_venn_abers(inputs$truth, inputs$p, inputs$fold, method = method)
      expected <- predict(core, inputs$new, loss = loss)
      actual <- checked$predictions[checked$predictions$id == outer$id[i], ]
      close(as.matrix(actual[paste0(".pred_", classes)]), expected)
      prepared[[i]] <- rsample::manual_rset(lapply(s$folds, function(f) {
        rsample::make_splits(f[c("analysis", "assessment")], local_data)
      }), names(s$folds))
    }
    resolver <- local_calibration_splits(prepared, names(data))
    explicit <- calibrated_workflow(base, "truth", method = "cross", classification = method,
      splits = resolver, seed = 33, loss = loss)
    # Repeat identical local data with the same prepared plan, preserving id/id2.
    repeated <- rsample::manual_rset(rep(outer$splits, 2),
      tibble::tibble(id = rep(c("Repeat1", "Repeat2"), each = 5), id2 = rep(outer$id, 2)))
    planned <- run(explicit, repeated)
    checked_repeat <- check(planned, repeated, classification = method)
    first <- checked_repeat$predictions[checked_repeat$predictions$id == "Repeat1", ]
    second <- checked_repeat$predictions[checked_repeat$predictions$id == "Repeat2", ]
    stopifnot(identical(first[setdiff(names(first), "id")], second[setdiff(names(second), "id")]),
              identical(checked$predictions[paste0(".pred_", classes)], first[paste0(".pred_", classes)]))
  }
  for (method in c("ovo", "ovr")) {
    make_outer <- function(data) rsample::manual_rset(list(rsample::make_splits(
      list(analysis = 1:40, assessment = 41:60), data)), "Test")
    changed <- data
    changed$truth[41:60] <- classes[1L]
    model <- calibrated_workflow(base, "truth", method = "cross", classification = method, seed = 33, loss = "brier")
    before <- run(model, make_outer(data))
    after <- run(model, make_outer(changed))
    check(before, make_outer(data), classification = method)
    check(after, make_outer(changed), classification = method)
    s1 <- tune::collect_extracts(before)$.extracts[[1L]]
    s2 <- tune::collect_extracts(after)$.extracts[[1L]]
    stopifnot(identical(s1[setdiff(names(s1), "models")], s2[setdiff(names(s2), "models")]))
    for (id in names(s1$models)) {
      stopifnot(identical(workflows::extract_fit_engine(s1$models[[id]]), workflows::extract_fit_engine(s2$models[[id]])),
                identical(workflows::extract_recipe(s1$models[[id]])$steps, workflows::extract_recipe(s2$models[[id]])$steps))
    }
    stopifnot(identical(tune::collect_predictions(before)[paste0(".pred_", classes)],
                        tune::collect_predictions(after)[paste0(".pred_", classes)]),
              !identical(tune::collect_metrics(before)$mean, tune::collect_metrics(after)$mean))
    failed_base <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail = TRUE))
    failed <- suppressWarnings(run(calibrated_workflow(failed_base, "truth", method = "cross",
      classification = method, seed = 33), make_outer(data)))
    stopifnot(nrow(tune::collect_notes(failed)) == 1L)
    error(check(failed, make_outer(data), classification = method), "")
  }
}

  testthat::expect_true(TRUE)
})

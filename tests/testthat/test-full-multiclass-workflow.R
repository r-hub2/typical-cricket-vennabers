testthat::test_that("full-multiclass-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")

packages <- c("parsnip", "workflows", "rsample", "recipes", "nnet")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |>
    workflows::add_model(parsnip::multinom_reg(penalty = .1) |>
      parsnip::set_engine("nnet", trace = FALSE))
  for (method in c("inductive", "cross")) for (classification in c("ovo", "ovr")) {
    matched <- generics::fit(calibrated_workflow(base, "truth", method = method,
      classification = classification, seed = 42), data)
    full <- generics::fit(calibrated_workflow(base, "truth", method = method,
      classification = classification, prediction_model = "full", seed = 42), data)
    source <- workflows::extract_fit_engine(matched)
    state <- workflows::extract_fit_engine(full)
    assert_conditions(identical(state$calibrator, source$calibrator),
              identical(state$folds, source$folds), identical(names(state$models), "Full1"),
              state$parameters$effective$prediction_model == "full",
              state$parameters$effective$model_training == "shared")
    p <- as.matrix(predict(state$models$Full1, data, type = "prob"))
    colnames(p) <- state$class_levels
    input <- if (method == "cross") setNames(rep(list(p), length(state$folds)), names(state$folds)) else p
    expected <- predict(source$calibrator, input)
    stopifnot(max(abs(as.matrix(predict(full, data, type = "prob")) - expected)) <= 1e-10)
  }
}

if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  close <- function(actual, expected) stopifnot(identical(dim(actual), dim(expected)),
    all(is.finite(actual)), all(abs(actual - expected) <= 1e-10))
  error <- function(expr, text = "") {
    e <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(e, "error"), grepl(text, conditionMessage(e), fixed = TRUE))
  }
  audit <- new.env()
  audit$fits <- list()
  audit$predictions <- list()
  shared_fit <- function(x, y, fail_size = integer(), bad_size = integer(), output = "reverse") {
    audit$fits[[length(audit$fits) + 1L]] <- list(x = x, y = y)
    if (nrow(x) %in% fail_size || length(unique(y)) < 2L) stop("untrainable shared model")
    list(x = x, y = y, levels = levels(y), bad_size = bad_size, output = output, draw = runif(1))
  }
  shared_prob <- function(object, new_data) {
    stopifnot(nrow(new_data) > 0L)
    audit$predictions[[length(audit$predictions) + 1L]] <- list(train = object$x$id, rows = new_data$id)
    p <- as.matrix(new_data[paste0("p", seq_along(object$levels))])
    # Training size changes full/source outputs; zero pairs remain exactly zero.
    p <- p ^ (nrow(object$x) / 10)
    p <- as.data.frame(p / rowSums(p))
    names(p) <- paste0(".pred_", object$levels)
    if (nrow(object$x) %in% object$bad_size) {
      if (object$output == "missing") return(p[-1L])
      if (object$output == "duplicate") names(p)[2L] <- names(p)[1L]
      if (object$output == "unknown") names(p)[1L] <- ".pred_unknown"
      if (object$output == "nonfinite") p[1L, 1L] <- Inf
      if (object$output == "range") p[1L, 1L] <- -0.1
      if (object$output == "sum") p <- p * (1 - 2e-8)
      if (object$output == "tolerance") p <- p * (1 - 5e-9)
      if (object$output == "rows") p <- p[-1L, , drop = FALSE]
    }
    p[rev(seq_along(p))]
  }
  register_reference_engine("va_shared_transfer", "classification", "shared_fit", "shared_prob", "prob")
  spec <- parsnip::new_model_spec("va_shared_transfer", args = list(), eng_args = NULL,
    mode = "classification", method = NULL, engine = "reference")
  make_plan <- function(data, method, held = NULL) {
    if (is.null(held)) held <- if (method == "inductive") list(49:60) else list(1:8, 9:24, 25:60)
    splits <- lapply(held, function(h) rsample::make_splits(
      list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data))
    if (method == "inductive") splits[[1L]] else rsample::manual_rset(splits, paste0("part", seq_along(held)))
  }
  probability <- function(model, new) {
    p <- predict(model, new, type = "prob")
    classes <- workflows::extract_fit_engine(model)$levels
    p <- as.matrix(p[paste0(".pred_", classes)])
    colnames(p) <- classes
    p
  }
  # Binary reference components plus independent odds coupling. No multiclass
  # core is used to compute either the expected endpoints or final probabilities.
  reference <- function(state, input, classification, precision, loss, details) {
    classes <- state$class_levels
    problems <- if (classification == "ovo") combn(classes, 2, simplify = FALSE) else as.list(classes)
    normalize <- function(p) p / Reduce(`+`, as.data.frame(p))
    score <- function(p, problem) {
      p <- normalize(p)
      if (length(problem) == 1L) return(p[, problem])
      total <- p[, problem[1L]] + p[, problem[2L]]
      ifelse(total == 0, .5, p[, problem[2L]] / total)
    }
    q <- matrix(0, nrow(input), length(problems))
    fold_points <- vector("list", length(problems))
    for (j in seq_along(problems)) {
      problem <- problems[[j]]
      event <- tail(problem, 1L)
      parts <- details[details$event_class == event, ]
      if (classification == "ovo") parts <- parts[parts$class_1 == problem[1L], ]
      endpoints <- lapply(names(state$folds), function(id) {
        cal <- state$calibration[[id]]
        rows <- if (classification == "ovo") cal$truth %in% problem else rep(TRUE, length(cal$truth))
        binary <- venn_abers(cal$truth[rows] == event, score(cal$probability, problem)[rows], precision)
        expected <- predict(binary, score(input, problem), loss = loss)
        close(as.matrix(parts[parts$fold == id, c("lower", "upper", "probability")]), as.matrix(expected))
        expected
      })
      lower <- vapply(endpoints, `[[`, numeric(nrow(input)), "lower")
      upper <- vapply(endpoints, `[[`, numeric(nrow(input)), "upper")
      points <- vapply(endpoints, `[[`, numeric(nrow(input)), "probability")
      q[, j] <- if (loss == "log") plogis(rowMeans(log(upper)) - rowMeans(log1p(-lower))) else rowMeans(points)
      fold_points[[j]] <- points
    }
    couple <- function(q) {
      p <- vapply(classes, function(cl) {
        if (classification == "ovr") return(q[, match(cl, classes)])
        relevant <- which(vapply(problems, function(pair) cl %in% pair, logical(1)))
        odds <- vapply(relevant, function(j) {
          prob <- if (cl == problems[[j]][1L]) 1 - q[, j] else q[, j]
          (1 - prob) / prob
        }, numeric(nrow(input)))
        1 / (1 + rowSums(odds))
      }, numeric(nrow(input)))
      p / rowSums(p)
    }
    list(prediction = couple(q), coupled_first = Reduce(`+`, lapply(seq_along(state$folds), function(i) {
      couple(vapply(fold_points, function(p) p[, i], numeric(nrow(input))))
    })) / length(state$folds))
  }
  maximum_error <- 0
  for (k in 2:4) {
    classes <- c("z", "a b", "a.b", "4")[seq_len(k)]
    p <- matrix(abs(sin(seq_len(60 * k))), 60, k)
    p[seq(1, 60, 7), 1:2] <- 0
    if (k == 2L) p[seq(1, 60, 7), 2L] <- 1
    data <- as.data.frame(p / rowSums(p))
    names(data) <- paste0("p", seq_len(k))
    data$id <- 1:60
    data$x <- c(NA, sin(2:48), 10000 + 1:12)
    data$truth <- factor(rep(classes, length.out = 60), levels = classes)
    base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ ., data) |>
      recipes::step_impute_median(x) |> recipes::step_normalize(x)) |> workflows::add_model(spec)
    new <- data[c(1:7, 17, 31, 60), ]
    for (method in c("inductive", "cross")) for (classification in c("ovo", "ovr")) {
      for (precision in list(NULL, 0, 2, -1)) {
        plan <- make_plan(data, method)
        fits <- states <- list()
        folds <- if (method == "inductive") 1L else 3L
        for (mode in c("matched", "full")) {
          before <- length(audit$fits)
          calls <- length(audit$predictions)
          set.seed(700)
          rng <- .Random.seed
          fits[[mode]] <- generics::fit(calibrated_workflow(base, "truth", method = method,
            classification = classification, model_training = "shared", prediction_model = mode,
            splits = plan, seed = 33, precision = precision, loss = "brier"), data)
          s <- states[[mode]] <- workflows::extract_fit_engine(fits[[mode]])
          assert_conditions(length(audit$fits) == before + folds + (mode == "full"), identical(.Random.seed, rng),
            length(audit$predictions) == calls + folds + (mode == "full"),
            s$parameters$effective$model_training == "shared", s$parameters$effective$prediction_model == mode,
            s$parameters$explicit[["model_training"]], is.null(s$event_class),
            identical(names(s$calibration), names(s$folds)), is.null(s$parameters$effective$event_level))
          for (i in seq_along(s$folds)) {
            f <- s$folds[[i]]
            learned <- audit$fits[[before + i]]
            x <- data$x[f$analysis]
            x[is.na(x)] <- median(x, na.rm = TRUE)
            assert_conditions(identical(learned$x$id, data$id[f$analysis]),
              identical(as.vector(learned$y), as.vector(data$truth[f$analysis])),
              identical(audit$predictions[[calls + i]], list(train = data$id[f$analysis], rows = data$id[f$assessment])),
              identical(s$calibration[[i]]$row, f$assessment), identical(s$calibration[[i]]$truth, data$truth[f$assessment]))
            close(learned$x$x, (x - mean(x)) / sd(x))
          }
        }
        s <- states$full
        assert_conditions(identical(s$calibrator, states$matched$calibrator), identical(s$calibration, states$matched$calibration),
          identical(s$folds, states$matched$folds), identical(names(s$models), "Full1"),
          identical(workflows::extract_fit_engine(s$models$Full1)$x$id, data$id))
        learned <- workflows::extract_recipe(s$models$Full1)
        x <- data$x
        x[is.na(x)] <- median(x, na.rm = TRUE)
        close(unname(learned$steps[[2L]]$means[["x"]]), mean(x))
        for (id in names(s$folds)) close(s$calibration[[id]]$probability,
          probability(states$matched$models[[id]], data[s$folds[[id]]$assessment, ]))
        input <- probability(s$models$Full1, new)
        for (loss in c("log", "brier")) {
          before <- length(audit$fits)
          calls <- length(audit$predictions)
          actual <- venn_abers_details(fits$full, new, loss = loss)
          d <- actual$details
          assert_conditions(length(audit$predictions) == calls + 1L, all(d$model_training == "shared"),
            all(d$prediction_model == "full"), all(d$model_source == "Full1"),
            identical(unique(d$fold), names(s$folds)), identical(names(actual$prediction), paste0(".pred_", classes)))
          expected <- reference(s, input, classification, precision, loss, d)
          close(as.matrix(actual$prediction), expected$prediction)
          maximum_error <- max(maximum_error, abs(as.matrix(actual$prediction) - expected$prediction))
          stopifnot(identical(actual$prediction, predict(fits$full, new, type = "prob", loss = loss)),
            identical(as.character(predict(fits$full, new, type = "class", loss = loss)$.pred_class),
              classes[max.col(expected$prediction, ties.method = "first")]))
          for (n in 0:1) {
            calls <- length(audit$predictions)
            small <- venn_abers_details(fits$full, new[seq_len(n), ], loss = loss)
            prob <- predict(fits$full, new[seq_len(n), ], type = "prob", loss = loss)
            cls <- predict(fits$full, new[seq_len(n), ], type = "class", loss = loss)
            assert_conditions(length(audit$predictions) == calls + 3L * n, identical(small$prediction, prob),
              nrow(small$details) == n * folds * if (classification == "ovo") choose(k, 2) else k,
              identical(levels(cls$.pred_class), classes), identical(names(small$details), names(d)),
              identical(vapply(small$details, typeof, ""), vapply(d, typeof, "")))
          }
          singles <- do.call(rbind, lapply(seq_len(nrow(new)), function(i) predict(fits$full, new[i, ], type = "prob", loss = loss)))
          stopifnot(identical(actual$prediction, singles), length(audit$fits) == before)
          if (is.null(precision)) {
            stopifnot(max(abs(as.matrix(actual$prediction) - as.matrix(predict(fits$matched, new, type = "prob", loss = loss)))) > 1e-5)
            if (method == "cross" && k > 2L) {
              stopifnot(max(abs(expected$prediction - expected$coupled_first)) > 1e-5)
              pooled <- venn_abers(do.call(c, lapply(s$calibration, `[[`, "truth")),
                do.call(rbind, lapply(s$calibration, `[[`, "probability")), method = classification)
              stopifnot(max(abs(expected$prediction - predict(pooled, input, loss = loss))) > 1e-5)
            }
          }
        }
      }
    }
  }
  cat("Shared full-model binary-component maximum absolute error:", maximum_error, "\n")
}

if (exists("shared_fit")) {
  before <- length(audit$fits)
  for (bad in list(NULL, "s", "SHARED", NA_character_, c("shared", "shared"), TRUE)) {
    for (classification in c("ovo", "ovr")) error(calibrated_workflow(base, "truth",
      classification = classification, model_training = bad), "model_training")
  }
  error(calibrated_workflow(base, "truth", classification = "ovo", model_train = "shared"), "model_train")
  binary <- workflows::workflow() |> workflows::add_formula(truth ~ x) |> workflows::add_model(spec)
  error(calibrated_workflow(binary, "truth", model_training = "shared"), "model_training")
  reg <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  error(calibrated_workflow(reg, "truth", model_training = "shared"), "model_training")
  error(predict(fits$full, new, model_training = "shared"), "model_training")
  error(predict(fits$full, new, prediction_model = "matched"), "prediction_model")
  error(venn_abers_details(fits$full, new, model_training = "shared"), "model_training")
  error(calibrate_workflow(s$models$Full1, data, "truth", classification = "ovo",
    model_training = "shared"), "model_training")
  stopifnot(length(audit$fits) == before)
  for (method in c("inductive", "cross")) for (classification in c("ovo", "ovr")) {
    plan <- make_plan(data, method)
    before <- length(audit$fits)
    error(generics::fit(calibrated_workflow(base, "truth", method = method, classification = classification,
      prediction_model = "full", splits = plan), data[60:1, ]), "data binding")
    absent <- data[data$truth != classes[1L], ]
    error(generics::fit(calibrated_workflow(base, "truth", method = method, classification = classification,
      prediction_model = "full"), absent), "all classes")
    invalid <- make_plan(data, method, if (method == "inductive") list(1L) else list(1L, 2:30, 31:60))
    if (method == "inductive" || classification == "ovo") error(generics::fit(calibrated_workflow(base, "truth",
      method = method, classification = classification, prediction_model = "full", splits = invalid), data), "calibration")
    stopifnot(length(audit$fits) == before)
    for (size in c(60L, if (method == "inductive") 48L else 24L)) {
      broken <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail_size = size))
      error(generics::fit(calibrated_workflow(broken, "truth", method = method, classification = classification,
        prediction_model = "full", splits = plan), data), "base workflow training")
      for (output in c("missing", "duplicate", "unknown", "nonfinite", "range", "sum", "rows")) {
        broken <- workflows::update_model(base, parsnip::set_engine(spec, "reference", bad_size = size, output = output))
        error(generics::fit(calibrated_workflow(broken, "truth", method = method, classification = classification,
          prediction_model = "full", splits = plan), data), "probabilit")
      }
    }
    tolerant <- workflows::update_model(base, parsnip::set_engine(spec, "reference", bad_size = 60L, output = "tolerance"))
    fit <- generics::fit(calibrated_workflow(tolerant, "truth", method = method, classification = classification,
      prediction_model = "full", splits = plan), data)
    state <- workflows::extract_fit_engine(fit)
    d <- venn_abers_details(fit, new)
    close(as.matrix(d$prediction), reference(state, probability(state$models$Full1, new), classification,
      NULL, "log", d$details)$prediction)
    # Automatic split and an explicitly replayed plan agree, including diagnostics.
    automatic <- generics::fit(calibrated_workflow(base, "truth", method = method, classification = classification,
      prediction_model = "full", seed = 33), data)
    state <- workflows::extract_fit_engine(automatic)
    prepared <- make_plan(data, method, lapply(state$folds, `[[`, "assessment"))
    if (method == "cross") prepared$id <- names(state$folds)
    resolver <- local_calibration_splits(list(prepared), names(data))
    explicit <- generics::fit(calibrated_workflow(base, "truth", method = method, classification = classification,
      model_training = "shared", prediction_model = "full", splits = resolver, seed = 33), data)
    stopifnot(identical(venn_abers_details(automatic, new), venn_abers_details(explicit, new)),
      !state$parameters$explicit[["model_training"]],
      any(grepl("Model training: shared", capture.output(print(state)), fixed = TRUE)),
      any(grepl("full-model transfer", capture.output(print(state)), fixed = TRUE)))
    tied <- data.frame(p1 = rep(.5, 60), p2 = rep(.5, 60), id = 1:60,
      truth = factor(rep(c("z", "a"), 30), levels = c("z", "a")))
    tied_base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |> workflows::add_model(spec)
    fit <- generics::fit(calibrated_workflow(tied_base, "truth", method = method, classification = classification,
      prediction_model = "full", splits = make_plan(tied, method)), tied)
    for (loss in c("log", "brier")) stopifnot(all(predict(fit, tied, type = "prob", loss = loss) == .5),
      all(predict(fit, tied, type = "class", loss = loss)$.pred_class == "z"))
  }
}

if (exists("shared_fit") && all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
  example <- new.env()
  sys.source(system.file("examples", "calibrated-resampling.R", package = "vennabers"), example)
  check <- example$check_evaluation
  scoring <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
  control <- tune::control_resamples(save_pred = TRUE, extract = identity, allow_par = FALSE)
  outer <- withr::with_seed(201, rsample::vfold_cv(data, v = 2, repeats = 2, strata = truth))
  holdout <- rsample::make_splits(list(analysis = 1:48, assessment = 49:60), data)
  for (method in c("inductive", "cross")) for (classification in c("ovo", "ovr")) {
    for (loss in c("log", "brier")) {
      options <- list(base, "truth", method = method, classification = classification,
        prediction_model = "full", seed = 33, loss = loss)
      automatic <- do.call(calibrated_workflow, options)
      result <- tune::fit_resamples(automatic, outer, metrics = scoring, control = control)
      checked <- check(result, outer, classification = classification)
      fits <- tune::collect_extracts(result)$.extracts
      prepared <- lapply(seq_along(fits), function(i) {
        state <- workflows::extract_fit_engine(fits[[i]])
        local <- rsample::analysis(outer$splits[[i]])
        held <- rsample::assessment(outer$splits[[i]])
        engine <- workflows::extract_fit_engine(state$models$Full1)
        stopifnot(identical(engine$x$id, local$id), !any(engine$x$id %in% held$id),
          identical(as.vector(engine$y), as.vector(local$truth)), identical(names(state$models), "Full1"))
        rows <- checked$predictions$id == outer$id[i] & checked$predictions$id2 == outer$id2[i]
        direct <- predict(fits[[i]], held, type = "prob")
        close(as.matrix(checked$predictions[rows, names(direct)]), as.matrix(direct))
        p <- make_plan(local, method, lapply(state$folds, `[[`, "assessment"))
        if (method == "cross") p$id <- names(state$folds)
        p
      })
      model <- do.call(calibrated_workflow, c(options, list(splits = local_calibration_splits(prepared, names(data)))))
      before <- length(audit$fits)
      calls <- length(audit$predictions)
      explicit <- tune::fit_resamples(model, outer, metrics = scoring, control = control)
      check(explicit, outer, classification = classification)
      stopifnot(length(audit$fits) == before + nrow(outer) * if (method == "cross") 6L else 2L,
        identical(tune::collect_predictions(result), tune::collect_predictions(explicit)),
        identical(tune::collect_metrics(result), tune::collect_metrics(explicit)))
      fitted <- tune::collect_extracts(explicit)$.extracts
      # Every source fit is confined to its local analysis rows, with its own
      # recipe; the extra full fit uses all and only the outer analysis rows.
      offset <- before
      for (i in seq_along(fitted)) {
        state <- workflows::extract_fit_engine(fitted[[i]])
        local <- rsample::analysis(outer$splits[[i]])
        for (f in c(state$folds, list(list(analysis = seq_len(nrow(local)))))) {
          offset <- offset + 1L
          engine <- audit$fits[[offset]]
          x <- local$x[f$analysis]
          x[is.na(x)] <- median(x, na.rm = TRUE)
          stopifnot(identical(engine$x$id, local$id[f$analysis]),
            identical(as.vector(engine$y), as.vector(local$truth[f$analysis])))
          close(engine$x$x, (x - mean(x)) / sd(x))
        }
      }
      # Ordered plans are resolved on the current local domain, also for last_fit.
      ordered <- function(local) make_plan(local, method, if (method == "inductive")
        list(seq.int(nrow(local) - 7L, nrow(local))) else split(seq_len(nrow(local)),
          cut(seq_len(nrow(local)), 3L, labels = FALSE)))
      model <- do.call(calibrated_workflow, c(options, list(splits = ordered)))
      ordered_result <- tune::fit_resamples(model, outer, metrics = scoring, control = control)
      check(ordered_result, outer, classification = classification)
      before <- length(audit$fits)
      final <- tune::last_fit(model, holdout, metrics = scoring)
      stopifnot(length(audit$fits) == before + if (method == "cross") 4L else 2L)
      final_p <- tune::collect_predictions(final)
      final_m <- tune::collect_metrics(final)
      assert_conditions(nrow(tune::collect_notes(final)) == 0L, identical(final_p$.row, 49:60),
        !anyDuplicated(final_p$.row), all(final_p$.config == "pre0_mod0_post0"),
        identical(final_p$truth, data$truth[49:60]), nrow(final_m) == 2L,
        setequal(final_m$.metric, c("brier_class", "mn_log_loss")), all(is.finite(final_m$.estimate)))
      changed <- data
      changed$truth[49:60] <- classes[1L]
      altered <- rsample::make_splits(list(analysis = 1:48, assessment = 49:60), changed)
      mutated <- tune::last_fit(model, altered, metrics = scoring)
      one <- tune::extract_workflow(final)
      two <- tune::extract_workflow(mutated)
      s1 <- workflows::extract_fit_engine(one)
      s2 <- workflows::extract_fit_engine(two)
      assert_conditions(identical(s1$calibration, s2$calibration), identical(s1$calibrator, s2$calibrator),
        identical(s1$folds, s2$folds),
        identical(workflows::extract_fit_engine(s1$models$Full1), workflows::extract_fit_engine(s2$models$Full1)),
        identical(workflows::extract_fit_engine(s1$models$Full1)$x$id, data$id[1:48]),
        identical(venn_abers_details(one, data), venn_abers_details(two, data)))
      if (loss == "log") {
        for (change in c("missing", "duplicate", "config", "metric", "repeat", "truth")) {
          bad <- result
          if (change == "missing") bad$.predictions[[1L]] <- bad$.predictions[[1L]][-1L, ]
          if (change == "duplicate") bad$.predictions[[1L]] <- rbind(bad$.predictions[[1L]], bad$.predictions[[1L]][1L, ])
          if (change == "config") bad$.predictions[[1L]]$.config <- "wrong"
          if (change == "metric") bad$.metrics[[1L]] <- bad$.metrics[[1L]][-1L, ]
          if (change == "repeat") bad$id2[1L] <- "wrong"
          if (change == "truth") bad$.predictions[[1L]]$truth[1L] <- NA
          error(check(bad, outer, classification = classification))
        }
        fail_metric <- yardstick::new_prob_metric(function(data, ...) stop("metric failed"), direction = "minimize")
        bad <- suppressWarnings(tune::fit_resamples(model, outer,
          metrics = yardstick::metric_set(yardstick::brier_class, fail_metric), control = control))
        stopifnot(nrow(tune::collect_notes(bad)) > 0L)
        error(check(bad, outer, classification = classification))
        partial <- rsample::manual_rset(list(holdout,
          rsample::make_splits(list(analysis = 1:36, assessment = 37:60), data)), c("First", "Second"))
        for (sizes in list(48L, c(36L, 48L))) {
          broken <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail_size = sizes))
          bad_model <- calibrated_workflow(broken, "truth", method = method, classification = classification,
            prediction_model = "full", splits = ordered)
          bad <- suppressWarnings(tune::fit_resamples(bad_model, partial, metrics = scoring, control = control))
          stopifnot(nrow(tune::collect_notes(bad)) == length(sizes))
          error(check(bad, partial, classification = classification))
        }
      }
    }
  }
}

  testthat::expect_true(TRUE)
})

testthat::test_that("frozen-multiclass-regression", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

packages <- c("parsnip", "workflows", "hardhat", "rsample", "recipes", "nnet",
              "rlang", "tibble", "withr", "generics", "tidyselect")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  error <- function(expr, text) {
    e <- tryCatch(force(expr), error = identity)
    if (!inherits(e, "error") || !grepl(text, conditionMessage(e), fixed = TRUE)) {
      stop("Expected ", text, ": ", if (inherits(e, "error")) conditionMessage(e) else "no error")
    }
  }
  close <- function(actual, expected, s = NULL) {
    stopifnot(identical(dim(actual), dim(expected)), length(actual) == length(expected),
              all(is.finite(actual)), all(is.finite(expected)))
    tolerance <- if (is.null(s)) 1e-10 else 1e-10 * s + 1e-10 * abs(expected)
    stopifnot(all(abs(actual - expected) <= tolerance))
  }
  matrix_prob <- function(workflow, data) {
    p <- as.matrix(predict(workflow, data, type = "prob"))
    colnames(p) <- sub("^\\.pred_", "", colnames(p))
    p
  }
  set.seed(904)
  data <- data.frame(x = rnorm(360), z = rnorm(360))
  data$truth <- factor(rep(c("third", "first", "second"), 120),
                       levels = c("third", "first", "second"))
  data$x[c(2, 155, 255, 355)] <- NA_real_
  post_available <- all(vapply(c("tailor", "probably", "mgcv"), requireNamespace, logical(1), quietly = TRUE))
  cases <- list()
  for (k in 2:3) {
    d <- data
    if (k == 2L) d$truth <- factor(rep(c("third", "first"), 180), levels = c("third", "first"))
    base <- workflows::workflow() |>
      workflows::add_recipe(recipes::recipe(truth ~ ., d[1:120, ]) |>
        recipes::step_impute_median(x) |>
        recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(if (k == 2L) parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
    frozen <- generics::fit(base, d[1:120, ])
    for (post in if (post_available) c(FALSE, TRUE) else FALSE) {
      model <- if (!post) frozen else generics::fit(workflows::add_tailor(base, tailor::tailor() |>
        tailor::adjust_probability_calibration(method = if (k == 2L) "logistic" else "multinomial", smooth = (k == 3L))),
        d[1:120, ], data_calibration = d[121:200, ])
      calibration <- d[201:320, ]
      new <- d[321:360, c("x", "z")]
      if (post) stopifnot(max(abs(matrix_prob(model, new) - matrix_prob(frozen, new))) > 1e-4)
      cases[[length(cases) + 1L]] <- list(model = model, calibration = calibration, new = new)
    }
  }
  regression_data <- data
  regression_data$x[is.na(regression_data$x)] <- 0
  regression_data$truth <- 2 + regression_data$x + cos(1:360)
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ ., regression_data[1:120, ]) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  regressors <- list(generics::fit(base, regression_data[1:120, ]))
  if (post_available) {
    regressors[[2L]] <- generics::fit(workflows::add_tailor(base, tailor::tailor() |>
      tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE)),
      regression_data[1:120, ], data_calibration = regression_data[121:200, ])
    stopifnot(max(abs(predict(regressors[[1L]], regression_data)$.pred -
                       predict(regressors[[2L]], regression_data)$.pred)) > 1e-4)
  }

  # Guard learning only after the complete original procedures have been fitted.
  guards <- list(c("stats", "glm.fit"), c("stats", "lm.fit"), c("nnet", "multinom"),
                 c("workflows", "fit.workflow"), c("recipes", "prep.recipe"))
  if (post_available) guards <- c(guards, list(c("tailor", "fit.tailor"), c("mgcv", "gam")))
  for (g in guards) trace(g[2L], where = asNamespace(g[1L]), tracer = quote(stop("Unexpected learning")), print = FALSE)
  on.exit(for (g in guards) untrace(g[2L], where = asNamespace(g[1L])), add = TRUE)

  for (case in cases) {
    model <- case$model
    calibration <- case$calibration
    new <- case$new
    classes <- levels(calibration$truth)
    final_cal <- matrix_prob(model, calibration[c("x", "z")])
    final_new <- matrix_prob(model, new)
    # Input factor order may differ; only the stored workflow defines class order.
    calibration$truth <- factor(calibration$truth, levels = rev(classes))
    for (method in c("ovo", "ovr")) for (precision in list(NULL, 0, 2, -1)) {
      for (default_loss in c("log", "brier")) {
        fit <- calibrate_workflow(model, calibration, "truth", classification = method,
                                  precision = precision, loss = default_loss)
        manual <- venn_abers(factor(calibration$truth, levels = classes), final_cal,
                             method = method, precision = precision)
        trace("venn_abers", where = asNamespace("vennabers"),
              tracer = quote(stop("Unexpected calibration learning")), print = FALSE)
        for (loss in c("log", "brier")) {
          expected <- predict(manual, final_new, type = "details", loss = loss)
          actual <- venn_abers_details(fit, new, loss = loss)
          close(as.matrix(actual$prediction), expected$prediction)
          parts <- if (method == "ovo") expected$pairs else expected$classes
          close(as.matrix(actual$details[c("lower", "upper", "probability")]),
                as.matrix(parts[c("lower", "upper", "probability")]))
          cls <- predict(fit, new, type = "class", loss = loss)$.pred_class
          assert_conditions(identical(actual$prediction, predict(fit, new, loss = loss)),
                    identical(names(actual$prediction), paste0(".pred_", classes)),
                    identical(cls, predict(manual, final_new, type = "class", loss = loss)),
                    identical(actual$details$row, parts$row),
                    identical(actual$details$event_class, parts$event),
                    all(actual$details$method == method), all(actual$details$fold == "Inductive1"))
          if (method == "ovo") stopifnot(identical(actual$details$class_1, parts$non_event),
                                         identical(actual$details$class_2, parts$event)) else
            stopifnot(identical(actual$details$class, parts$event))
          for (n in 0:1) {
            small <- new[seq_len(n), , drop = FALSE]
            detail <- venn_abers_details(fit, small, loss = loss)
            assert_conditions(nrow(detail$prediction) == n, nrow(detail$details) == n * length(manual$calibrators),
                      identical(names(detail$details), names(actual$details)),
                      identical(detail$prediction, predict(fit, small, loss = loss)),
                      all(vapply(detail$prediction, is.double, logical(1))),
                      identical(levels(predict(fit, small, type = "class")$.pred_class), classes))
          }
        }
        assert_conditions(identical(fit$models[[1L]], model), identical(fit$calibrator, manual),
                  identical(fit$parameters$requested$precision, precision),
                  identical(fit$parameters$effective$precision, precision), fit$loss == default_loss,
                  is.null(fit$event_class), is.null(fit$parameters$effective$event_level),
                  identical(predict(fit, new), do.call(rbind, lapply(seq_len(nrow(new)), function(i) predict(fit, new[i, ])))))
        untrace("venn_abers", where = asNamespace("vennabers"))
      }
      changed <- calibration
      changed$truth <- rev(changed$truth)
      other <- calibrate_workflow(model, changed, "truth", classification = method)
      stopifnot(identical(other$models[[1L]], model), !identical(other$calibrator,
        calibrate_workflow(model, calibration, "truth", classification = method)$calibrator),
        identical(predict(other, new), predict(other, transform(new, truth = NA))))
    }
    for (method in c("ovo", "ovr")) {
      error(calibrate_workflow(model, calibration[calibration$truth != classes[1L], ], "truth", classification = method), "Every class")
      error(calibrate_workflow(model, calibration, "truth", classification = method, event_level = "first"), "event_level")
    }
    if (length(classes) == 3L) error(calibrate_workflow(model, calibration, "truth"), "exactly two")
  }

  # Unit scale s is fixed before each comparison, including shifted units.
  calibration <- regression_data[201:320, ]
  new <- regression_data[321:360, c("x", "z")]
  for (model in regressors) for (scale in c(1e-9, 1, 1e9)) for (shift in c(0, 100)) {
    cal <- calibration
    cal$truth <- scale * (cal$truth + shift)
    for (options in list(list(), list(m = NULL, epsilon = NULL), list(m = 1), list(m = 59), list(epsilon = .1))) {
      fit <- do.call(calibrate_workflow, c(list(model, cal, "truth"), options))
      manual <- do.call(venn_abers_regression, c(list(cal$truth, predict(model, cal[c("x", "z")])$.pred),
                                               options[!vapply(options, is.null, logical(1))]))
      trace("venn_abers_regression", where = asNamespace("vennabers"),
            tracer = quote(stop("Unexpected calibration learning")), print = FALSE)
      expected <- predict(manual, predict(model, new)$.pred)
      actual <- venn_abers_details(fit, new)
      close(actual$prediction$.pred, expected$prediction, scale)
      close(as.matrix(actual$details[c("prediction", "lower", "upper")]), as.matrix(expected), scale)
      assert_conditions(identical(actual$prediction, predict(fit, new)),
                identical(actual$prediction, predict(fit, new, type = "numeric")),
                identical(names(actual$prediction), ".pred"), is.double(actual$prediction$.pred),
                identical(actual$details$row, seq_len(nrow(new))), all(actual$details$fold == "Inductive1"),
                identical(fit$calibrator, manual), identical(fit$models[[1L]], model),
                identical(fit$parameters$requested$m, options$m),
                identical(fit$parameters$requested$epsilon, options$epsilon),
                identical(fit$parameters$effective$m, manual$m), identical(fit$parameters$effective$epsilon, manual$epsilon),
                is.null(fit$class_levels), is.null(fit$event_class), is.null(fit$loss),
                is.null(fit$parameters$effective$classification), is.null(fit$parameters$effective$precision),
                identical(predict(fit, new), predict(fit, transform(new, truth = NA))),
                identical(predict(fit, new), do.call(rbind, lapply(seq_len(nrow(new)), function(i) predict(fit, new[i, ])))))
      for (n in 0:1) {
        small <- new[seq_len(n), , drop = FALSE]
        detail <- venn_abers_details(fit, small)
        stopifnot(identical(detail$prediction, predict(fit, small)), nrow(detail$details) == n,
                  identical(names(detail$details), names(actual$details)), is.double(detail$prediction$.pred))
      }
      untrace("venn_abers_regression", where = asNamespace("vennabers"))
    }
  }
  model <- regressors[[1L]]
  fit <- calibrate_workflow(model, calibration, "truth")
  changed <- calibration
  changed$truth <- rev(changed$truth)
  other <- calibrate_workflow(model, changed, "truth")
  stopifnot(!identical(other$calibrator, fit$calibrator), identical(other$models, fit$models))
  for (name in c("classification", "event_level", "loss", "precision")) {
    value <- list(classification = "binary", event_level = "first", loss = "log", precision = NULL)[name]
    error(do.call(calibrate_workflow, c(list(model, calibration, "truth"), value)), "Regression does not accept")
  }
  for (value in list(NULL, "log", "brier")) {
    error(predict(fit, new, loss = value), "Regression does not accept")
    error(venn_abers_details(fit, new, loss = value), "Regression does not accept")
  }
  for (type in c("prob", "class", "conf_int", "pred_int")) error(predict(fit, new, type = type), "type")
  for (m in list(0, -1, 1.5, 60, Inf, NA_real_, TRUE, "1", numeric(), c(1, 2))) {
    error(calibrate_workflow(model, calibration, "truth", m = m), "m")
  }
  for (epsilon in list(0, 1, -1, 1e-10, Inf, NA_real_, TRUE, "1", numeric(), c(.1, .2))) {
    error(calibrate_workflow(model, calibration, "truth", epsilon = epsilon), "`")
  }
  error(calibrate_workflow(model, calibration, "truth", m = 1, epsilon = .1), "only one")
  for (k in 1:2) error(calibrate_workflow(model, calibration[seq_len(k), ], "truth"), "at least three")
  for (truth in list(rep(7, 3), c(0, 1, 0), c(-2, 0, 3))) {
    cal <- calibration[1:3, ]
    cal$truth <- truth
    f <- calibrate_workflow(model, cal, "truth")
    expected <- predict(venn_abers_regression(truth, predict(model, cal)$.pred), predict(model, new)$.pred)
    close(as.matrix(venn_abers_details(f, new)$details[c("prediction", "lower", "upper")]), as.matrix(expected), 1)
    if (all(truth == 7)) stopifnot(all(predict(f, new)$.pred == 7))
  }
  for (truth in list(factor(rep("a", 120)), rep(NA_real_, 120), rep(Inf, 120), as.character(calibration$truth))) {
    bad <- calibration
    bad$truth <- truth
    error(calibrate_workflow(model, bad, "truth"), "regression truth")
  }

  # Inject outputs at the original workflow's public predict boundary. This also
  # verifies outcome removal, exact row counts, and no call for empty inputs.
  calls <- 0L
  output <- "reverse"
  predict_reference <- function(object, new_data, type, ...) {
    calls <<- calls + 1L
    stopifnot(!"truth" %in% names(new_data), nrow(new_data) > 0L)
    if (output == "error") stop("Stored prediction failed")
    p <- NextMethod()
    if (output == "rows") return(p[-1L, , drop = FALSE])
    if (output == "missing") return(p[-1L])
    if (output == "extra") p$unknown <- 0
    if (output == "unknown") names(p)[1L] <- ".pred_unknown"
    if (output == "duplicate") p <- cbind(p, p[1L])
    if (output == "nan") p[[1L]][1L] <- NaN
    if (output == "inf") p[[1L]][1L] <- Inf
    if (output == "string") p[[1L]] <- as.character(p[[1L]])
    if (type == "numeric") {
      if (output == "scores") p$.pred <- new_data$x
      return(p)
    }
    if (output == "negative") p[[1L]][1L] <- -.1
    if (output == "range") p[[1L]][1L] <- 2
    if (output == "sum") p[[1L]][1L] <- .5
    if (output == "tolerance") p <- p * (1 + 5e-9)
    if (output == "zero_pair") {
      p[,] <- 0
      p[[1L]] <- rep(1, nrow(p))
    }
    if (output == "tie") p[,] <- 1 / ncol(p)
    p[rev(seq_len(ncol(p)))]
  }
  registerS3method("predict", "va_frozen_multi_reference", predict_reference, envir = asNamespace("stats"))
  for (case in cases[c(1L, if (post_available) 3L else 2L)]) for (method in c("ovo", "ovr")) {
    reference <- case$model
    class(reference) <- c("va_frozen_multi_reference", class(reference))
    cal <- case$calibration
    new <- case$new
    output <- "reverse"
    before <- calls
    fit <- calibrate_workflow(reference, cal, "truth", classification = method)
    stopifnot(calls == before + 1L)
    before <- calls
    for (type in c("prob", "class")) predict(fit, new[0, ], type = type)
    venn_abers_details(fit, new[0, ])
    stopifnot(calls == before)
    for (output in c("rows", "missing", "extra", "unknown", "duplicate", "nan", "inf", "string",
                     "negative", "range", "sum", "error")) {
      error(calibrate_workflow(reference, cal, "truth", classification = method), "Frozen calibration")
      error(predict(fit, new), "base")
      error(venn_abers_details(fit, new), "base")
    }
    for (output in c("zero_pair", "tie", "tolerance")) for (precision in list(NULL, 0, 2, -1)) {
      # One observation per class exercises minimal nonempty subproblems.
      small <- cal[match(levels(cal$truth), cal$truth), ]
      fit <- calibrate_workflow(reference, small, "truth", classification = method, precision = precision)
      manual <- venn_abers(small$truth, matrix_prob(reference, small[c("x", "z")]), method = method, precision = precision)
      for (loss in c("log", "brier")) {
        actual <- venn_abers_details(fit, new, loss = loss)
        expected <- predict(manual, matrix_prob(reference, new), loss = loss, type = "details")
        close(as.matrix(actual$prediction), expected$prediction)
        parts <- if (method == "ovo") expected$pairs else expected$classes
        close(as.matrix(actual$details[c("lower", "upper", "probability")]),
              as.matrix(parts[c("lower", "upper", "probability")]))
        if (output == "tie" && is.null(precision)) stopifnot(all(predict(fit, new, type = "class", loss = loss)$.pred_class == levels(cal$truth)[1L]))
      }
    }
  }
  reference <- regressors[[1L]]
  class(reference) <- c("va_frozen_multi_reference", class(reference))
  output <- "reverse"
  fit <- calibrate_workflow(reference, calibration, "truth")
  before <- calls
  predict(fit, new[0, ])
  venn_abers_details(fit, new[0, ])
  stopifnot(calls == before)
  for (output in c("rows", "missing", "extra", "unknown", "duplicate", "nan", "inf", "string", "error")) {
    error(calibrate_workflow(reference, calibration, "truth"), "Frozen calibration")
    error(predict(fit, new), "base")
    error(venn_abers_details(fit, new), "base")
  }
  output <- "scores"
  # Nonconstant truth with both merger weights zero, plus scaled/tied scores.
  cal <- data.frame(x = 1:3, z = 0, truth = c(-1e16, -3, -3))
  fit <- calibrate_workflow(reference, cal, "truth")
  actual <- venn_abers_details(fit, data.frame(x = 0, z = 0))
  close(as.matrix(actual$details[c("prediction", "lower", "upper")]),
        matrix(c(-1e16 / 2 - 3 / 2, -1e16, -3), 1), 1)
  for (s in c(1e-9, 1, 1e9)) for (shift in c(0, 100)) {
    cal <- data.frame(x = rep(s * shift, 5), z = 0, truth = s * (c(0, 0, 0, 4, 10) + shift))
    fit <- calibrate_workflow(reference, cal, "truth")
    actual <- venn_abers_details(fit, cal[1L, c("x", "z")])
    close(as.matrix(actual$details[c("prediction", "lower", "upper")]),
          matrix(s * (c(116 / 33, 4 / 3, 4) + shift), 1), s)
  }
})

  testthat::expect_true(TRUE)
})

testthat::test_that("fit diagnoses of frozen and trained workflows show their VA calibrators", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

packages <- c("parsnip", "workflows", "hardhat", "rsample", "recipes", "nnet",
              "rlang", "tibble", "withr", "generics", "tidyselect")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  error <- function(expr, text) {
    e <- tryCatch(force(expr), error = identity)
    if (!inherits(e, "error") || !grepl(text, conditionMessage(e), fixed = TRUE)) {
      stop("Expected ", text, ": ", if (inherits(e, "error")) conditionMessage(e) else "no error")
    }
  }
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  layout <- graphics::par("mfrow")
  # Every new panel records the grid it is drawn in.
  grids <- list()
  setHook("plot.new", function() grids[[length(grids) + 1L]] <<- graphics::par("mfrow"))
  on.exit({ setHook("plot.new", NULL, "replace"); grDevices::dev.off(); unlink(path) })
  panels <- function(expr) {
    grids <<- list()
    force(expr)
    grids
  }
  values <- c("score", "lower", "upper", "probability", "width")

  set.seed(911)
  d <- data.frame(x = rnorm(180), z = rnorm(180))
  d$truth <- factor(rep(c("third", "first", "second"), 60), levels = c("third", "first", "second"))
  d$outcome <- 2 * d$x + d$z + rnorm(180)
  two <- d[d$truth != "second", ]
  two$truth <- droplevels(two$truth)
  logistic <- workflows::workflow() |> workflows::add_formula(truth ~ x + z) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  multinomial <- workflows::workflow() |> workflows::add_formula(truth ~ x + z) |>
    workflows::add_model(parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
  linear <- workflows::workflow() |> workflows::add_formula(outcome ~ x + z) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))

  # Frozen fits plot their stored calibrator directly, with the stored loss
  # as the default; trained fits plot through extract_fit_engine().
  frozen <- list(
    binary = calibrate_workflow(generics::fit(logistic, two[1:60, ]), two[61:120, ], "truth", loss = "brier"),
    ovr = calibrate_workflow(generics::fit(multinomial, d[1:90, ]), d[91:180, ], "truth", classification = "ovr"),
    ovo = calibrate_workflow(generics::fit(multinomial, d[1:90, ]), d[91:180, ], "truth", classification = "ovo",
                             loss = "brier"))
  trained <- list(
    binary = generics::fit(calibrated_workflow(logistic, "truth", loss = "brier", seed = 7), two),
    cross = generics::fit(calibrated_workflow(logistic, "truth", method = "cross", v = 3, seed = 7), two),
    ovr = generics::fit(calibrated_workflow(multinomial, "truth", classification = "ovr", seed = 7), d),
    ovo = generics::fit(calibrated_workflow(multinomial, "truth", classification = "ovo", method = "cross",
                                            v = 3, loss = "brier", seed = 7), d))
  fits <- c(frozen, lapply(trained, workflows::extract_fit_engine))
  for (k in seq_along(fits)) {
    fit <- fits[[k]]
    stored <- fit$loss
    other <- setdiff(c("log", "brier"), stored)
    assert_conditions(inherits(fit, "vennabers_fit"), identical(k <= length(frozen), inherits(fit, "vennabers_frozen")),
                      !identical(plot(fit), plot(fit, loss = other)))
    drawn <- withVisible(plot(fit, c("width", "map")))
    assert_conditions(!drawn$visible, identical(drawn$value, plot(fit$calibrator, "map", loss = stored)),
                      identical(plot(fit, loss = other), plot(fit$calibrator, loss = other)),
                      identical(graphics::par("mfrow"), layout))
  }

  # Regression delegates without a loss; passing one is an error.
  regression <- list(
    frozen = calibrate_workflow(generics::fit(linear, d[1:90, ]), d[91:180, ], "outcome"),
    trained = workflows::extract_fit_engine(generics::fit(
      calibrated_workflow(linear, "outcome", method = "cross", v = 3, seed = 7), d)))
  for (fit in regression) {
    assert_conditions(identical(plot(fit, c("map", "width")), plot(fit$calibrator)),
                      identical(plot(fit, xlim = c(-1, 1)), plot(fit$calibrator, xlim = c(-1, 1))))
    for (bad in list("log", "brier", NULL)) error(plot(fit, loss = bad), "Regression does not accept `loss`.")
  }

  # Postprocessing models are source models before the VA calibrator: the
  # map is that of the stored calibrator.
  sigmoid <- workflows::extract_fit_engine(generics::fit(calibrated_workflow(logistic, "truth", seed = 7,
    postprocessing = list(method = "sigmoid", sigmoid_input = "probability", strategy = "holdout", prop = .3)), two))
  assert_conditions(!is.null(sigmoid$postprocessing), identical(plot(sigmoid), plot(sigmoid$calibrator)))

  # Pairwise fits draw one panel per pair component; `pair_id` leads the
  # identity columns, each pair's rows are the map of its own calibrator.
  for (method in c("inductive", "cross")) {
    cross <- method == "cross"
    options <- list(logistic, "truth", classification = "ovo", model_training = "pairwise",
                    method = method, seed = 7)
    if (cross) options$v <- 3
    fit <- workflows::extract_fit_engine(generics::fit(do.call(calibrated_workflow, options), d))
    for (loss in c("log", "brier")) {
      drawn <- panels(maps <- plot(fit, loss = loss))
      assert_conditions(identical(names(maps), c("pair_id", if (cross) "fold", "class_1", "class_2", values)),
                        identical(unique(maps$pair_id), names(fit$pair_models)),
                        identical(rownames(maps), as.character(seq_len(nrow(maps)))),
                        length(drawn) == 3L, identical(drawn[[1L]], c(2L, 2L)),
                        identical(graphics::par("mfrow"), layout))
      for (part in fit$pair_models) {
        own <- maps[maps$pair_id == part$pair_id, ]
        rownames(own) <- NULL
        assert_conditions(all(own$class_1 == part$class_levels[1L]), all(own$class_2 == part$class_levels[2L]),
                          identical(own[setdiff(names(own), c("pair_id", "class_1", "class_2"))],
                                    plot(part$calibrator, loss = loss)))
      }
    }
    assert_conditions(identical(plot(fit), plot(fit, loss = fit$loss)))
    drawn <- panels(both <- plot(fit, c("width", "map")))
    assert_conditions(identical(both, plot(fit)), length(drawn) == 7L,
                      all(vapply(drawn, identical, TRUE, c(4L, 2L))),
                      identical(graphics::par("mfrow"), layout))
    error(plot(fit, ylim = "bad"), "ylim")
    assert_conditions(identical(graphics::par("mfrow"), layout))
    for (bad in list("hinge", NA_character_)) error(plot(fit, loss = bad), "`loss` must be one of: log, brier.")
    error(plot(fit, which = "fit"), "`which` must be one or more distinct values of: map, width.")
    broken <- fit
    broken$pair_models[[2L]]$calibrator <- NULL
    error(plot(broken), "`x` must be a fitted calibrated_workflow() or calibrate_workflow().")
    broken$pair_models <- NULL
    error(plot(broken), "`x` must be a fitted calibrated_workflow() or calibrate_workflow().")
  }

  # A fit without a VA calibrator stops clearly.
  for (fit in c(fits[c(1L, length(frozen) + 2L)], regression)) {
    broken <- fit
    broken$calibrator <- NULL
    error(plot(broken), "`x` must be a fitted calibrated_workflow() or calibrate_workflow().")
    broken$calibrator <- unclass(fit$calibrator)
    error(plot(broken), "`x` must be a fitted calibrated_workflow() or calibrate_workflow().")
  }
  error(plot(structure(list(), class = "vennabers_fit")),
        "`x` must be a fitted calibrated_workflow() or calibrate_workflow().")
})

testthat::expect_true(TRUE)
})

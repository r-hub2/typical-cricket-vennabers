testthat::test_that("frozen-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

packages <- c("parsnip", "workflows", "hardhat", "rsample", "recipes",
              "rlang", "tibble", "withr", "generics", "tidyselect")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  set.seed(903)
  data <- data.frame(x = rnorm(240), z = rnorm(240))
  data$truth <- factor(ifelse(runif(240) < plogis(data$x), "yes", "no"),
                       levels = c("yes", "no"))
  data$x[c(2, 105, 205)] <- NA_real_
  training <- data[1:100, ]
  calibration <- data[101:200, ]
  new <- data[201:240, c("x", "z")]
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ ., training) |>
                           recipes::step_impute_median(x) |>
                           recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  frozen <- generics::fit(base, training)
  original <- frozen
  fitted <- calibrate_workflow(frozen, calibration, "truth")
  manual <- venn_abers(calibration$truth == "yes",
                       predict(frozen, calibration[c("x", "z")], type = "prob")$.pred_yes)
  expected <- predict(manual, predict(frozen, new, type = "prob")$.pred_yes)
  actual <- venn_abers_details(fitted, new)
  assert_conditions(max(abs(as.matrix(actual$details[c("lower", "upper", "probability")]) -
                     as.matrix(expected))) <= 1e-10,
            identical(actual$prediction, predict(fitted, new, type = "prob")),
            identical(original, fitted$models[[1L]]),
            !inherits(fitted, "workflow"), !inherits(fitted, "model_spec"))

  if (all(vapply(c("tailor", "probably"), requireNamespace, logical(1), quietly = TRUE))) {
    post_base <- workflows::add_tailor(base, tailor::tailor() |>
      tailor::adjust_probability_calibration(method = "logistic", smooth = FALSE))
    post_frozen <- generics::fit(post_base, training, data_calibration = data[101:140, ])
    calibration <- data[141:200, ]
    workflows_to_test <- list(frozen, post_frozen)
    stopifnot(max(abs(as.matrix(predict(frozen, new, type = "prob")) -
                       as.matrix(predict(post_frozen, new, type = "prob")))) > .01)
  } else workflows_to_test <- list(frozen)

  regression <- generics::fit(workflows::workflow() |>
    workflows::add_formula(x ~ z) |> workflows::add_model(parsnip::linear_reg(engine = "lm")),
    training[!is.na(training$x), ])

  # Activate after the original fits. All subsequent calibration/prediction must
  # retain these learned components; predicting must not even fit another VA.
  guards <- list(c("stats", "glm.fit"), c("workflows", "fit.workflow"),
                 c("recipes", "prep.recipe"))
  if (length(workflows_to_test) == 2L) guards <- c(guards, list(c("tailor", "fit.tailor")))
  for (g in guards) trace(g[2L], where = asNamespace(g[1L]),
                          tracer = quote(stop("Unexpected learning")), print = FALSE)
  on.exit(for (g in guards) untrace(g[2L], where = asNamespace(g[1L])), add = TRUE)
  expect_error <- function(expr, text) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    if (!inherits(error, "error") || !grepl(text, conditionMessage(error), fixed = TRUE)) {
      stop("Expected ", text, ": ", if (is.null(error)) "no error" else conditionMessage(error))
    }
  }
  for (workflow in workflows_to_test) {
    final_calibration <- predict(workflow, calibration[c("x", "z")], type = "prob")
    final_new <- predict(workflow, new, type = "prob")
    for (event in c("first", "second")) for (precision in list(NULL, 1L, 3L)) {
      event_class <- levels(training$truth)[if (event == "first") 1L else 2L]
      column <- paste0(".pred_", event_class)
      for (default_loss in c("log", "brier")) {
        calibrated <- calibrate_workflow(workflow, calibration, "truth",
          event_level = event, precision = precision, loss = default_loss)
        manual <- venn_abers(calibration$truth == event_class, as.double(final_calibration[[column]]), precision = precision)
        trace("venn_abers", where = asNamespace("vennabers"),
              tracer = quote(stop("Unexpected calibration learning")), print = FALSE)
        for (loss in c("log", "brier")) {
          expected <- predict(manual, final_new[[column]], loss = loss)
          actual <- venn_abers_details(calibrated, new, loss = loss)
          prob <- predict(calibrated, new, type = "prob", loss = loss)
          cls <- predict(calibrated, new, type = "class", loss = loss)$.pred_class
          assert_conditions(max(abs(as.matrix(actual$details[c("lower", "upper", "probability")]) -
                             as.matrix(expected))) <= 1e-10,
                    identical(prob, actual$prediction), identical(names(prob), c(".pred_yes", ".pred_no")),
                    identical(levels(cls), levels(training$truth)),
                    identical(as.character(cls), ifelse(expected$probability >= .5, event_class,
                      setdiff(levels(training$truth), event_class))),
                    all(actual$details$event_class == event_class),
                    identical(actual$details$row, seq_len(nrow(new))),
                    all(actual$details$fold == "Inductive1"))
          for (n in 0:1) {
            small <- new[seq_len(n), , drop = FALSE]
            detail <- venn_abers_details(calibrated, small, loss = loss)
            assert_conditions(identical(detail$prediction, predict(calibrated, small, loss = loss)),
                      identical(names(detail$details), names(actual$details)),
                      nrow(detail$details) == n,
                      all(vapply(detail$prediction, is.double, logical(1))),
                      identical(levels(predict(calibrated, small, type = "class", loss = loss)$.pred_class),
                                levels(training$truth)))
          }
        }
        assert_conditions(identical(calibrated$loss, default_loss),
                  identical(calibrated$models[[1L]], workflow), identical(calibrated$calibrator, manual),
                  identical(calibrated$parameters$requested$precision, precision),
                  identical(calibrated$parameters$effective$event_level, event),
                  identical(predict(calibrated, new),
                    do.call(rbind, lapply(seq_len(nrow(new)), function(i) predict(calibrated, new[i, ])))))
        untrace("venn_abers", where = asNamespace("vennabers"))
      }
    }
    changed <- calibration
    changed$x <- 1e6
    changed$truth <- factor(rev(changed$truth), levels = rev(levels(training$truth)))
    changed_fit <- calibrate_workflow(workflow, changed, "truth")
    stopifnot(identical(changed_fit$models[[1L]], workflow),
              identical(predict(workflow, new, type = "prob"), final_new),
              !identical(changed_fit$calibrator, calibrate_workflow(workflow, calibration, "truth")$calibrator),
              identical(predict(changed_fit, new), predict(changed_fit, transform(new, truth = NA))))
    one_sided <- calibration
    one_sided$truth[] <- "yes"
    stopifnot(nrow(predict(calibrate_workflow(workflow, one_sided, "truth"), new)) == nrow(new))
  }

  for (bad in list(base, workflows::extract_fit_parsnip(frozen), workflows::extract_fit_engine(frozen), NULL)) {
    expect_error(calibrate_workflow(bad, calibration, "truth"), "fully fitted workflow")
  }
  expect_error(calibrate_workflow(regression, calibration, "truth"), "outcome")
  expect_error(calibrate_workflow(frozen, calibration, "wrong"), "outcome")
  for (outcome in list(c("truth", "x"), "", NA_character_, 1)) {
    expect_error(calibrate_workflow(frozen, calibration, outcome), "outcome")
  }
  for (option in c("workflow", "data", "out", "classification_method", "method", "v", "splits", "seed",
                   "calibration_prop", "strata", "predictors", "m", "epsilon", "full_model", "train", "fit")) {
    expect_error(do.call(calibrate_workflow, c(list(frozen, calibration, "truth"), setNames(list(NULL), option))), option)
  }
  for (name in c("classification", "event_level", "loss")) for (value in list(NULL, NA, "bad", c("first", "second"))) {
    expect_error(do.call(calibrate_workflow, c(list(frozen, calibration, "truth"), setNames(list(value), name))), name)
  }
  for (precision in list(TRUE, .5, Inf, NA_real_, numeric(), "2", matrix(2))) {
    expect_error(calibrate_workflow(frozen, calibration, "truth", precision = precision), "precision")
  }
  expect_error(calibrate_workflow(frozen, calibration[0, ], "truth"), "nonempty")
  for (bad in list(new, calibration["truth"], as.matrix(calibration))) {
    expect_error(calibrate_workflow(frozen, bad, "truth"), "Frozen")
  }
  for (truth in list(as.character(calibration$truth), as.integer(calibration$truth),
                    factor(rep("yes", nrow(calibration))),
                    factor(calibration$truth, levels = c("yes", "no", "unknown")),
                    factor(rep(NA_character_, nrow(calibration)), levels = levels(training$truth)))) {
    bad <- calibration
    bad$truth <- truth
    expect_error(calibrate_workflow(frozen, bad, "truth"), "nonmissing factor with trained classes")
  }
  bad <- calibration
  names(bad)[2L] <- "x"
  expect_error(calibrate_workflow(frozen, bad, "truth"), "x")
  bad <- frozen
  bad$pre$mold$blueprint$recipe$steps[[1L]]$trained <- FALSE
  expect_error(calibrate_workflow(bad, calibration, "truth"), "recipe")
  bad <- frozen
  bad$fit$fit$lvl <- rev(bad$fit$fit$lvl)
  expect_error(calibrate_workflow(bad, calibration, "truth"), "trained outcome order")
  if (length(workflows_to_test) == 2L) {
    bad <- post_frozen
    bad$post$fit <- NULL
    expect_error(calibrate_workflow(bad, calibration, "truth"), "not been fit")
    bad <- post_frozen
    bad$post$fit$adjustments[[1L]]$trained <- FALSE
    expect_error(calibrate_workflow(bad, calibration, "truth"), "postprocessor")
  }

  # Instrument the public workflow boundary, including malformed engine outputs.
  # Calling an extracted engine would bypass both this guard and these outputs.
  calls <- 0L
  output <- "reverse"
  predict_reference <- function(object, new_data, type, ...) {
    calls <<- calls + 1L
    stopifnot(!"truth" %in% names(new_data), nrow(new_data) > 0L, type == "prob")
    if (output == "error") stop("Original predictor x failed in stored preprocessing")
    p <- NextMethod()
    if (output == "missing") return(p[1L])
    if (output == "additional") p$.pred_unknown <- 0
    if (output == "unknown") names(p)[1L] <- ".pred_unknown"
    if (output == "duplicate") names(p)[2L] <- names(p)[1L]
    if (output == "nonfinite") p[[1L]][1L] <- NA_real_
    if (output == "infinite") p[[1L]][1L] <- Inf
    if (output == "range") p[[1L]][1L] <- 2
    if (output == "negative") p[[1L]][1L] <- -.1
    if (output == "sum") p[[1L]] <- p[[2L]] <- rep(.2, nrow(p))
    if (output == "string") p[[1L]] <- as.character(p[[1L]])
    if (output == "rows") return(p[-1L, ])
    if (output == "class") return(tibble::tibble(.pred_class = factor(rep("yes", nrow(p)))))
    if (output == "tie") p[[1L]] <- p[[2L]] <- rep(.5, nrow(p))
    p[c(2L, 1L, setdiff(seq_len(ncol(p)), 1:2))]
  }
  registerS3method("predict", "va_frozen_reference", predict_reference, envir = asNamespace("stats"))
  reference <- frozen
  class(reference) <- c("va_frozen_reference", class(reference))
  reference_fit <- calibrate_workflow(reference, calibration, "truth")
  stopifnot(calls == 1L)
  expected_fit <- calibrate_workflow(frozen, calibration, "truth")
  stopifnot(identical(predict(reference_fit, new), predict(expected_fit, new)))
  before <- calls
  for (type in c("prob", "class")) predict(reference_fit, new[0, ], type = type)
  venn_abers_details(reference_fit, new[0, ])
  stopifnot(calls == before)
  for (output in c("missing", "additional", "unknown", "duplicate", "nonfinite", "infinite", "range",
                   "negative", "sum", "string", "rows", "class", "error")) {
    expect_error(calibrate_workflow(reference, calibration, "truth"), "Frozen calibration")
    expect_error(predict(reference_fit, new), "base")
    expect_error(venn_abers_details(reference_fit, new), "base")
  }
  output <- "tie"
  for (event in c("first", "second")) {
    balanced <- calibration[1:2, ]
    balanced$truth <- factor(c("yes", "no"), levels = levels(training$truth))
    tie <- calibrate_workflow(reference, balanced, "truth", event_level = event)
    for (loss in c("log", "brier")) {
      stopifnot(all(predict(tie, new, loss = loss)[[1L]] == .5),
                all(predict(tie, new, type = "class", loss = loss)$.pred_class ==
                      levels(training$truth)[if (event == "first") 1L else 2L]))
    }
  }
  output <- "reverse"
  for (value in list(NULL, "l", "Brier", NA_character_, c("log", "brier"))) {
    expect_error(predict(reference_fit, new, loss = value), "loss")
    expect_error(venn_abers_details(reference_fit, new, loss = value), "loss")
  }
  for (type in c("p", "numeric", "raw")) expect_error(predict(reference_fit, new, type = type), "type")
  expect_error(predict(reference_fit, new, los = "brier"), "los")
  expect_error(venn_abers_details(reference_fit, new, los = "brier"), "los")
  for (bad in list(new["z"], transform(new, x = as.character(x)))) {
    expect_error(predict(reference_fit, bad), "x")
    expect_error(venn_abers_details(reference_fit, bad), "x")
  }
  expect_error(generics::fit(reference_fit, training), "fit")
})

  testthat::expect_true(TRUE)
})

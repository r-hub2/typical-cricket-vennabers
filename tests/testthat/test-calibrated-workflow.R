testthat::test_that("calibrated-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")

integration_packages <- c("parsnip", "workflows", "hardhat", "rsample",
                          "recipes", "rlang", "tibble", "withr", "generics")
if (all(vapply(integration_packages, requireNamespace, logical(1), quietly = TRUE))) local({
  data <- data.frame(x = sin(seq_len(100)),
                     truth = factor(rep(c("yes", "no", "no", "yes"), 25)))
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  model <- calibrated_workflow(base, "truth", seed = 42)
  fitted <- generics::fit(model, data)
  state <- workflows::extract_fit_engine(fitted)
  split <- withr::with_seed(42, rsample::initial_split(data, prop = .8,
                                                      strata = truth, pool = 0))
  train <- as.integer(split, data = "analysis")
  held <- as.integer(split, data = "assessment")
  manual_model <- generics::fit(base, data[train, ])
  p <- predict(manual_model, data[held, ], type = "prob")[[paste0(".pred_", state$event_class)]]
  manual <- venn_abers(data$truth[held] == state$event_class, p)
  new_p <- predict(manual_model, data, type = "prob")[[paste0(".pred_", state$event_class)]]
  expected <- predict(manual, new_p)
  actual <- venn_abers_details(fitted, data)
  assert_conditions(max(abs(actual$details$probability - expected$probability)) <= 1e-10, max(abs(actual$details$lower -
      expected$lower)) <= 1e-10, max(abs(actual$details$upper - expected$upper)) <= 1e-10, identical(actual$prediction,
      predict(fitted, data, type = "prob")), identical(state$folds$Inductive1$analysis, train), identical(state$folds$Inductive1$assessment,
      held))
})

if (all(vapply(integration_packages, requireNamespace, logical(1), quietly = TRUE))) {
  expect_error <- function(expr, text) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    if (!inherits(error, "error") || !grepl(text, conditionMessage(error), fixed = TRUE)) {
      stop("Expected ", text, ": ", if (is.null(error)) "no error" else conditionMessage(error))
    }
  }
  # A deterministic public parsnip reference engine preserves arbitrary names,
  # accepts a one-class training subset, and exposes its actual training input.
  counts <- new.env()
  counts$fit <- counts$predict <- 0L
  reference_fit <- function(x, y, output = "reverse") {
    counts$fit <- counts$fit + 1L
    structure(list(x = x, y = y, levels = levels(y), output = output), class = "va_reference")
  }
  reference_prob <- function(object, new_data) {
    counts$predict <- counts$predict + 1L
    stopifnot(nrow(new_data) > 0L)
    p <- if ("p" %in% names(new_data)) new_data$p else rep(.5, nrow(new_data))
    result <- setNames(data.frame(p, 1 - p), paste0(".pred_", object$levels))
    if (object$output == "missing") return(result[1L])
    if (object$output == "duplicate") names(result)[2L] <- names(result)[1L]
    if (object$output == "unknown") names(result)[1L] <- ".pred_unknown"
    if (object$output == "nonfinite") result[[1L]][1L] <- NA_real_
    if (object$output == "nonnumeric") result[[1L]] <- as.character(result[[1L]])
    if (object$output == "range") result[[1L]][1L] <- 2
    if (object$output == "sum") result[[1L]] <- result[[2L]] <- rep(.2, nrow(new_data))
    result[2:1]
  }
  register_reference_engine("va_reference", "classification", "reference_fit", "reference_prob", "prob")
  reference <- parsnip::new_model_spec("va_reference", args = list(), eng_args = NULL,
                                      mode = "classification", method = NULL, engine = "reference")
  data <- data.frame(p = seq(.013, .987, length.out = 40),
                     id = seq_len(40), z = factor(rep(c("a", "b"), 20)),
                     truth = factor(rep(c("a b", "a.b", "a.b", "a b"), 10),
                                    levels = c("a b", "a.b")))
  base <- workflows::workflow() |>
    workflows::add_variables(outcomes = truth, predictors = c(p, id, z),
                            blueprint = hardhat::default_xy_blueprint(allow_novel_levels = TRUE)) |>
    workflows::add_model(reference)
  for (event_level in c("first", "second")) for (precision in list(NULL, 0, 2, -1)) {
    fitted <- generics::fit(calibrated_workflow(base, "truth", seed = 73,
                                               event_level = event_level, precision = precision), data)
    state <- workflows::extract_fit_engine(fitted)
    held <- state$folds$Inductive1$assessment
    train <- state$folds$Inductive1$analysis
    assert_conditions(identical(tibble::as_tibble(workflows::extract_fit_engine(state$models$Inductive1)$x),
        tibble::as_tibble(data[train, c("p", "id", "z")])), identical(state$columns, names(data)), identical(state$class_levels,
        levels(data$truth)), identical(state$parameters$requested$precision, precision), state$parameters$explicit[["precision"]],
        !state$parameters$explicit[["loss"]], length(state$models) == 1L)
    event <- state$event_class
    calibration <- predict(state$models$Inductive1, data[held, ], type = "prob")
    manual <- venn_abers(data$truth[held] == event, calibration[[paste0(".pred_", event)]], precision)
    input <- predict(state$models$Inductive1, data, type = "prob")[[paste0(".pred_", event)]]
    for (loss in c("log", "brier")) {
      expected <- predict(manual, input, loss = loss)
      details <- venn_abers_details(fitted, data, loss = loss)
      stopifnot(max(abs(as.matrix(details$details[c("lower", "upper", "probability")]) -
                            as.matrix(expected))) <= 1e-10,
                identical(details$prediction, predict(fitted, data, type = "prob", loss = loss)),
                identical(names(details$prediction), paste0(".pred_", levels(data$truth))),
                identical(as.character(predict(fitted, data, type = "class", loss = loss)$.pred_class),
                          ifelse(expected$probability >= .5, event, setdiff(levels(data$truth), event))))
      for (n in 0:1) {
        new <- data[seq_len(n), , drop = FALSE]
        calls <- counts$predict
        small <- venn_abers_details(fitted, new, loss = loss)
        prob <- predict(fitted, new, type = "prob", loss = loss)
        cls <- predict(fitted, new, type = "class", loss = loss)
        if (!n) stopifnot(counts$predict == calls)
        assert_conditions(identical(small$prediction, prob), nrow(prob) == n, identical(levels(cls$.pred_class),
            levels(data$truth)), identical(names(small$details), c("row", "fold", "event_class", "lower", "upper",
            "probability", "prediction_model", "model_source")), all(vapply(prob, is.double, logical(1))), is.integer(small$details$row),
            is.character(small$details$fold), is.character(small$details$event_class))
      }
    }
    before <- counts$fit
    batch <- predict(fitted, data, type = "prob")
    singles <- do.call(rbind, lapply(seq_len(nrow(data)), function(i) predict(fitted, data[i, ], type = "prob")))
    stopifnot(isTRUE(all.equal(batch, singles)), counts$fit == before,
              identical(state$loss, "log"),
              identical(workflows::extract_fit_engine(fitted)$loss, "log"))
  }

  # Local RNG restoration; NULL consumes the caller's stream, with no retries.
  set.seed(88)
  rng <- .Random.seed
  fitted <- generics::fit(calibrated_workflow(base, "truth", seed = 5), data)
  stopifnot(identical(rng, .Random.seed))
  for (strata in list(NULL, FALSE)) {
    set.seed(3)
    expected <- if (is.null(strata)) rsample::initial_split(data, .8, strata = truth, pool = 0) else rsample::initial_split(data, .8)
    rng_after <- .Random.seed
    set.seed(3)
    actual <- generics::fit(calibrated_workflow(base, "truth", strata = strata), data)
    positions <- workflows::extract_fit_engine(actual)$folds$Inductive1
    stopifnot(identical(positions$analysis, as.integer(expected, data = "analysis")),
              identical(positions$assessment, as.integer(expected, data = "assessment")),
              identical(rng_after, .Random.seed))
  }
  new <- data[1:2, ]
  new$z <- factor(c("novel", "a"))
  stopifnot(identical(predict(fitted, new, type = "prob"), venn_abers_details(fitted, new)$prediction))
  for (new in list(data[c("id", "z")], transform(data, p = as.character(p)))) {
    expect_error(predict(fitted, new, type = "prob"), "p")
    expect_error(venn_abers_details(fitted, new), "p")
  }
  for (bad in list(NULL, "l", "Brier", NA_character_, c("log", "brier"))) {
    expect_error(predict(fitted, data, type = "prob", loss = bad), "loss")
    expect_error(venn_abers_details(fitted, data, loss = bad), "loss")
  }
  expect_error(predict(fitted, data, type = "prob", los = "brier"), "los")
  expect_error(venn_abers_details(fitted, data, los = "brier"), "los")
  expect_error(predict(fitted, data, type = "p"), "type")

  # All rejected constructor/fit inputs must precede inner training.
  before <- counts$fit
  for (name in c(".pred", ".pred_class", ".pred_yes", ".row", ".config", "id", "id1", "id2", "id12")) {
    expect_error(calibrated_workflow(base, name), name)
  }
  for (name in c("class", ".pred_yes")) {
    bad <- data
    levels(bad$truth)[1L] <- name
    expect_error(generics::fit(calibrated_workflow(base, "truth"), bad), name)
  }
  for (outcome in list(c("truth", "p"), "", NA_character_)) expect_error(calibrated_workflow(base, outcome), "outcome")
  expect_error(calibrated_workflow(base, "truth", predictors = c("p", "truth")), "truth")
  expect_error(calibrated_workflow(base, "truth", predictors = c("p", "p")), "p")
  for (option in c("v", "m", "epsilon", "calibration_p", "clas", "unknown")) {
    expect_error(do.call(calibrated_workflow, c(list(base, "truth"), setNames(list(NULL), option))), option)
  }
  expect_error(calibrated_workflow(work = base, outcome = "truth"), "work")
  expect_error(calibrated_workflow(base, out = "truth"), "out")
  for (method in c("cro", "ind", "inductive ")) expect_error(calibrated_workflow(base, "truth", method = method), "method")
  for (classification in c("ov", "OVR", "bin")) expect_error(calibrated_workflow(base, "truth", classification = classification), "classification")
  for (prop in list(0, 1, -.1, Inf, NA_real_, TRUE, numeric(), c(.1, .2))) {
    expect_error(calibrated_workflow(base, "truth", calibration_prop = prop), "calibration_prop")
  }
  for (precision in list(TRUE, .5, Inf, NA_real_, numeric(), "2", matrix(2))) {
    expect_error(calibrated_workflow(base, "truth", precision = precision), "precision")
  }
  for (seed in list(TRUE, 1.5, Inf, NA_real_, 2^32, numeric())) {
    expect_error(calibrated_workflow(base, "truth", seed = seed), "seed")
  }
  expect_error(calibrated_workflow(base, "truth", strata = "truth"), "strata")
  for (truth in list(as.character(data$truth), as.integer(data$truth),
                    factor(rep("a", 40)), factor(rep(c("a", "b", "c", "a"), 10)),
                    factor(rep("a", 40), levels = c("a", "b")))) {
    bad <- data
    bad$truth <- truth
    expect_error(generics::fit(calibrated_workflow(base, "truth"), bad), "")
  }
  bad <- data
  names(bad)[2L] <- "p"
  expect_error(generics::fit(calibrated_workflow(base, "truth"), bad), "p")
  stopifnot(counts$fit == before)

  trained <- generics::fit(base, data)
  expect_error(calibrated_workflow(trained, "truth"), "fitted")
  reg <- workflows::workflow() |> workflows::add_formula(truth ~ p) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  stopifnot(workflows::extract_spec_parsnip(calibrated_workflow(reg, "truth"))$mode == "regression")
  expect_error(generics::fit(calibrated_workflow(reg, "truth"), data), "regression truth")
  recipe <- recipes::recipe(truth ~ p + id, data)
  prepped <- recipes::prep(recipe)
  # workflows itself rejects prepped recipes, before the constructor can see one.
  expect_error(workflows::add_recipe(workflows::workflow(), prepped), "trained")
  tuned <- workflows::update_model(base, parsnip::set_engine(reference, "reference", output = quote(tune())))
  expect_error(calibrated_workflow(tuned, "truth"), "tune")
  if (requireNamespace("tailor", quietly = TRUE)) {
    post <- workflows::add_tailor(base, tailor::tailor() |> tailor::adjust_probability_threshold(.5))
    expect_error(calibrated_workflow(post, "truth"), "postprocessor")
  }
  for (output in c("missing", "duplicate", "unknown", "nonfinite", "nonnumeric", "range", "sum")) {
    bad <- workflows::update_model(base, parsnip::set_engine(reference, "reference", output = output))
    expect_error(generics::fit(calibrated_workflow(bad, "truth", seed = 42), data), "probabilit")
  }

  # Balanced tied calibration gives exactly .5: both losses favor the event.
  tied <- data
  tied$p <- .5
  for (level in c("first", "second")) {
    fit <- generics::fit(calibrated_workflow(base, "truth", seed = 42, event_level = level, loss = "brier"), tied)
    event <- levels(tied$truth)[if (level == "first") 1L else 2L]
    for (loss in c("log", "brier")) {
      stopifnot(all(predict(fit, tied, type = "prob", loss = loss) == .5),
                all(predict(fit, tied, type = "class", loss = loss)$.pred_class == event),
                identical(venn_abers_details(fit, tied, loss = loss)$prediction, predict(fit, tied, type = "prob", loss = loss)))
    }
    stopifnot(identical(predict(fit, tied, type = "prob"), predict(fit, tied, type = "prob", loss = "brier")),
              workflows::extract_fit_engine(fit)$loss == "brier")
  }

  # Singleton calibration is valid even when the other class is only in training.
  tiny <- data[1:2, ]
  fit <- generics::fit(calibrated_workflow(base, "truth", calibration_prop = .5, strata = FALSE, seed = 4), tiny)
  stopifnot(workflows::extract_fit_engine(fit)$calibrator$observations == 1L)
  before <- counts$fit
  expect_error(generics::fit(calibrated_workflow(base, "truth", calibration_prop = .99, seed = 4), tiny), "nonempty")
  stopifnot(counts$fit == before)

  # Held-out extremes never affect imputation or scaling learned by the recipe.
  split <- withr::with_seed(42, rsample::initial_split(data, .8, strata = truth, pool = 0))
  held <- as.integer(split, data = "assessment")
  train <- as.integer(split, data = "analysis")
  data$id[train[1L]] <- NA_integer_
  data$id[held] <- 1000000L
  rec <- recipes::recipe(truth ~ p + id, data) |>
    recipes::step_impute_median(recipes::all_numeric_predictors()) |>
    recipes::step_normalize(id)
  recipe_base <- workflows::workflow() |> workflows::add_recipe(rec) |> workflows::add_model(reference)
  before <- counts$fit
  fit <- generics::fit(calibrated_workflow(recipe_base, "truth", seed = 42), data)
  learned <- workflows::extract_recipe(workflows::extract_fit_engine(fit)$models$Inductive1)
  expected <- data$id[train]
  expected[is.na(expected)] <- median(expected, na.rm = TRUE)
  stopifnot(counts$fit == before + 1L,
            abs(learned$steps[[1L]]$medians[["id"]] - median(expected)) < 1e-10,
            abs(learned$steps[[2L]]$means[["id"]] - mean(expected)) < 1e-10,
            abs(learned$steps[[2L]]$sds[["id"]] - sd(expected)) < 1e-10)
}

  testthat::expect_true(TRUE)
})

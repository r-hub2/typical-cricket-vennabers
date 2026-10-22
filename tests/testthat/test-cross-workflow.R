testthat::test_that("cross-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")
if (all(vapply(c("parsnip", "workflows", "rsample", "recipes", "withr"),
               requireNamespace, logical(1), quietly = TRUE))) {
  data <- data.frame(x = sin(1:80), truth = factor(rep(c("yes", "no", "no", "yes"), 20)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  fitted <- generics::fit(calibrated_workflow(base, "truth", method = "cross", seed = 42), data)
  state <- workflows::extract_fit_engine(fitted)
  plan <- withr::with_seed(42, rsample::vfold_cv(data, v = 5, strata = truth, pool = 0))
  stopifnot(identical(names(state$models), plan$id), identical(names(state$folds), plan$id))
  p <- truth <- fold <- list()
  new_p <- matrix(NA_real_, nrow(data), nrow(plan), dimnames = list(NULL, plan$id))
  for (i in seq_len(nrow(plan))) {
    id <- plan$id[i]
    held <- as.integer(plan$splits[[i]], data = "assessment")
    train <- as.integer(plan$splits[[i]], data = "analysis")
    stopifnot(identical(state$folds[[id]]$analysis, train), identical(state$folds[[id]]$assessment, held))
    model <- state$models[[id]]
    p[[i]] <- predict(model, data[held, ], type = "prob")[[paste0(".pred_", state$event_class)]]
    truth[[i]] <- data$truth[held] == state$event_class
    fold[[i]] <- rep(id, length(held))
    new_p[, id] <- predict(model, data, type = "prob")[[paste0(".pred_", state$event_class)]]
  }
  core <- cross_venn_abers(unlist(truth), unlist(p), unlist(fold))
  for (loss in c("log", "brier")) {
    expected <- predict(core, new_p, loss = loss, type = "details")
    actual <- venn_abers_details(fitted, data, loss = loss)
    stopifnot(max(abs(actual$prediction[[paste0(".pred_", state$event_class)]] - expected$prediction$probability)) <= 1e-10,
              identical(actual$prediction, predict(fitted, data, type = "prob", loss = loss)),
              identical(as.data.frame(actual$details[names(expected$folds)]), expected$folds))
  }
  # Normalized prepared rsets preserve fold identity, order and local data.
  manual_plan <- rsample::manual_rset(lapply(state$folds, function(f) {
    rsample::make_splits(f[c("analysis", "assessment")], data)
  }), names(state$folds))
  resolve <- local_calibration_splits(list(plan, manual_plan), names(data))
  explicit <- generics::fit(calibrated_workflow(base, "truth", method = "cross", splits = resolve), data)
  stopifnot(identical(venn_abers_details(fitted, data), venn_abers_details(explicit, data)))
}

if (all(vapply(c("parsnip", "workflows", "rsample", "recipes", "withr"),
               requireNamespace, logical(1), quietly = TRUE))) {
  error <- function(expr, text = "") {
    result <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(result, "error"), grepl(text, conditionMessage(result), fixed = TRUE))
  }
  counts <- new.env()
  counts$fit <- counts$predict <- 0L
  cross_fit <- function(x, y, output = "reverse", fail = FALSE) {
    counts$fit <- counts$fit + 1L
    if (fail) stop("instrumented training failure")
    structure(list(x = x, y = y, levels = levels(y), output = output,
                   offset = mean(as.integer(y)) - 1.5, draw = runif(1)), class = "cross_reference")
  }
  cross_prob <- function(object, new_data) {
    counts$predict <- counts$predict + 1L
    stopifnot(nrow(new_data) > 0L)
    p <- if (object$output == "tie") rep(.5, nrow(new_data)) else plogis(new_data$x + object$offset)
    result <- setNames(data.frame(p, 1 - p), paste0(".pred_", object$levels))
    if (object$output == "missing") return(result[1L])
    if (object$output == "unknown") names(result)[1L] <- ".pred_unknown"
    if (object$output == "duplicate") names(result)[1L] <- names(result)[2L]
    if (object$output == "nonfinite") result[[1L]][1L] <- NA_real_
    result[2:1]
  }
  register_reference_engine("cross_reference", "classification", "cross_fit", "cross_prob", "prob")
  spec <- parsnip::new_model_spec("cross_reference", args = list(), eng_args = NULL,
                                 mode = "classification", method = NULL, engine = "reference")
  data <- data.frame(x = c(1000, NA, sin(3:40)), id = 1:40,
                     truth = factor(c(rep("a b", 9), rep(c("a.b", "a.b", "a b"), 10), "a.b"),
                                    levels = c("a b", "a.b")))
  recipe <- recipes::recipe(truth ~ x + id, data) |>
    recipes::step_impute_median(x) |> recipes::step_normalize(x)
  base <- workflows::workflow() |> workflows::add_recipe(recipe) |> workflows::add_model(spec)
  make_plan <- function(data, held = list(1:7, 8:20, 21:40), ids = c("small", "middle", "large")) {
    rsample::manual_rset(lapply(held, function(h) {
      rsample::make_splits(list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data)
    }), ids)
  }
  plan <- make_plan(data)
  model <- function(...) calibrated_workflow(base, "truth", method = "cross", ...)
  for (level in c("first", "second")) for (precision in list(NULL, 0, 2, -1)) {
    before <- counts$fit
    fitted <- generics::fit(model(splits = plan, precision = precision, event_level = level, loss = "brier"), data)
    state <- workflows::extract_fit_engine(fitted)
    assert_conditions(counts$fit == before + 3L, identical(names(state$models), plan$id),
              identical(names(state$calibrator$calibrators), plan$id),
              is.null(state$parameters$effective$v), is.null(state$parameters$effective$calibration_prop),
              identical(state$parameters$requested$precision, precision))
    calibration <- truth <- fold <- list()
    new_p <- matrix(NA_real_, nrow(data), 3, dimnames = list(NULL, plan$id))
    for (i in seq_len(nrow(plan))) {
      id <- plan$id[i]
      inner <- state$folds[[id]]
      train <- inner$analysis
      held <- inner$assessment
      engine <- workflows::extract_fit_engine(state$models[[id]])
      learned <- workflows::extract_recipe(state$models[[id]])
      x <- data$x[train]
      x[is.na(x)] <- median(x, na.rm = TRUE)
      assert_conditions(identical(engine$x$id, data$id[train]), identical(as.vector(engine$y), as.vector(data$truth[train])),
                abs(learned$steps[[1]]$medians[["x"]] - median(x)) <= 1e-10,
                abs(learned$steps[[2]]$means[["x"]] - mean(x)) <= 1e-10,
                abs(learned$steps[[2]]$sds[["x"]] - sd(x)) <= 1e-10,
                inner$training_size == length(train), inner$calibration_size == length(held))
      calibration[[i]] <- predict(state$models[[id]], data[held, ], type = "prob")[[paste0(".pred_", state$event_class)]]
      truth[[i]] <- data$truth[held] == state$event_class
      fold[[i]] <- rep(id, length(held))
      new_p[, id] <- predict(state$models[[id]], data, type = "prob")[[paste0(".pred_", state$event_class)]]
    }
    core <- cross_venn_abers(unlist(truth), unlist(calibration), unlist(fold), precision = precision)
    for (loss in c("log", "brier")) {
      expected <- predict(core, new_p, loss = loss, type = "details")
      actual <- venn_abers_details(fitted, data, loss = loss)
      p <- actual$prediction[[paste0(".pred_", state$event_class)]]
      assert_conditions(max(abs(p - expected$prediction$probability)) <= 1e-10,
                identical(as.data.frame(actual$details[names(expected$folds)]), expected$folds),
                identical(actual$prediction, predict(fitted, data, type = "prob", loss = loss)),
                identical(names(actual$prediction), paste0(".pred_", levels(data$truth))),
                identical(as.character(predict(fitted, data, type = "class", loss = loss)$.pred_class),
                          ifelse(p >= .5, state$event_class, setdiff(state$class_levels, state$event_class))))
      lower <- matrix(actual$details$lower, nrow(data))
      upper <- matrix(actual$details$upper, nrow(data))
      points <- matrix(actual$details$probability, nrow(data))
      reconstructed <- if (loss == "log") plogis(rowMeans(log(upper)) - rowMeans(log1p(-lower))) else rowMeans(points)
      stopifnot(max(abs(p - reconstructed)) <= 1e-10)
      if (is.null(precision)) {
        swapped <- new_p[, 3:1]
        colnames(swapped) <- colnames(new_p)
        stopifnot(max(abs(p - predict(core, swapped, loss = loss)$probability)) > 1e-5,
                  max(abs(p - as.vector(points %*% c(7, 13, 20) / 40))) > 1e-5)
      }
      before_predict <- counts$predict
      for (n in 0:1) {
        small <- data[seq_len(n), ]
        details <- venn_abers_details(fitted, small, loss = loss)
        prob <- predict(fitted, small, type = "prob", loss = loss)
        cls <- predict(fitted, small, type = "class", loss = loss)
        assert_conditions(nrow(prob) == n, nrow(cls) == n, nrow(details$details) == 3L * n,
                  identical(details$prediction, prob), identical(levels(cls$.pred_class), levels(data$truth)),
                  identical(names(details$details), c("row", "fold", "event_class", "lower", "upper", "probability", "prediction_model", "model_source")),
                  is.integer(details$details$row), is.character(details$details$fold),
                  is.character(details$details$event_class), all(vapply(prob, is.double, logical(1))))
        if (!n) stopifnot(counts$predict == before_predict)
      }
    }
    before <- counts$fit
    singles <- do.call(rbind, lapply(seq_len(nrow(data)), function(i) predict(fitted, data[i, ], type = "prob")))
    stopifnot(isTRUE(all.equal(singles, predict(fitted, data, type = "prob"))), counts$fit == before,
              identical(predict(fitted, data, type = "prob"), predict(fitted, data, type = "prob", loss = "brier")),
              workflows::extract_fit_engine(fitted)$loss == "brier")
  }

  # Validate the entire partition before any inner model learns, even a bad last fold.
  before <- counts$fit
  for (v in list(NULL, TRUE, 1, 1.5, Inf, NA_real_, c(2, 3), 41)) error(generics::fit(model(v = v), data), "v")
  error(model(calibration_prop = .2), "calibration_prop")
  error(model(splits = plan, v = 5), "cannot be combined")
  error(model(splits = plan, strata = NULL), "cannot be combined")
  error(model(splits = plan$splits[[1]]), "rset")
  for (ids in list(c("a", "a", "b"), c("a", "", "b"), c("a", NA, "b"))) {
    bad <- plan
    bad$id <- ids
    error(generics::fit(model(splits = bad), data), "fold ID")
  }
  for (indices in list(c(1, 1), c(0, 1), c(1, 41), c(1, NA), c(1, Inf), c(1, 2.5), numeric())) {
    for (field in c("in_id", "out_id")) {
      bad <- plan
      bad$splits[[3]][[field]] <- indices
      error(generics::fit(model(splits = bad), data), "partition")
    }
  }
  for (held in list(list(1:7, 7:20, 21:40), list(1:7, 9:20, 21:40))) {
    error(generics::fit(model(splits = make_plan(data, held)), data), "exactly once")
  }
  bad <- plan
  bad$splits[[3]]$in_id <- 1:19
  error(generics::fit(model(splits = bad), data), "partition")
  bad <- plan
  bad$splits[[3]] <- bad$splits[[2]]
  error(generics::fit(model(splits = bad), data), "exactly once")
  bad <- plan
  class(bad) <- c("group_vfold_cv", class(bad))
  error(generics::fit(model(splits = bad), data), "simple Cross")
  bad <- plan
  class(bad$splits[[3]]) <- c("boot_split", class(bad$splits[[3]]))
  error(generics::fit(model(splits = bad), data), "simple Cross")
  error(generics::fit(model(splits = rsample::vfold_cv(data, v = 2, repeats = 2)), data), "simple Cross")
  error(generics::fit(model(splits = rsample::bootstraps(data, times = 2)), data), "simple Cross")
  error(generics::fit(model(splits = function(data) 1:10), data), "simple Cross")
  for (change in c("values", "type", "levels", "order", "target")) {
    foreign <- data
    if (change == "values") foreign$x[1] <- 99
    if (change == "type") foreign$id <- as.double(foreign$id)
    if (change == "levels") levels(foreign$truth) <- rev(levels(foreign$truth))
    if (change == "order") foreign <- foreign[40:1, ]
    if (change == "target") foreign$truth[1] <- "a.b"
    error(generics::fit(model(splits = plan), foreign), "data binding")
  }
  stopifnot(counts$fit == before)

  # Prepared assignments check every entry and include fold order and IDs in identity.
  resolver <- local_calibration_splits(list(plan, plan), names(data))
  same <- tibble::as_tibble(data[c("truth", "id", "x")])
  same$unused <- 1:40
  same_fit <- generics::fit(calibrated_workflow(base, "truth", method = "cross", predictors = c("x", "id"), splits = resolver), same)
  original_fit <- generics::fit(model(splits = plan), data)
  stopifnot(identical(venn_abers_details(same_fit, data), venn_abers_details(original_fit, data)))
  conflict <- plan
  conflict$id <- rev(conflict$id)
  error(local_calibration_splits(list(plan, conflict), names(data)), "entries 1 and 2")
  conflict <- rsample::manual_rset(rev(plan$splits), rev(plan$id))
  error(local_calibration_splits(list(plan, conflict), names(data)), "Conflicting")
  error(local_calibration_splits(list(plan, make_plan(data, list(1:8, 9:20, 21:40))), names(data)), "Conflicting")
  error(local_calibration_splits(list(plan, plan$splits[[1]]), names(data)), "Conflicting")
  error(resolver(transform(data, id = id + 1L)), "No prepared plan")
  partial <- local_calibration_splits(list(plan), c("x", "id"))
  error(generics::fit(model(splits = partial), data), "canonical selected")

  # One RNG scope covers resolver and all fold fits; errors restore it as well.
  set.seed(66)
  saved_rng <- .Random.seed
  draws <- function(fit) vapply(workflows::extract_fit_engine(fit)$models,
                               function(w) workflows::extract_fit_engine(w)$draw, numeric(1))
  one <- generics::fit(model(seed = 42), data)
  two <- generics::fit(model(seed = 42), data)
  stopifnot(identical(saved_rng, .Random.seed), identical(draws(one), draws(two)), !anyDuplicated(draws(one)))
  for (strata in list(NULL, FALSE)) {
    set.seed(93)
    expected <- if (is.null(strata)) rsample::vfold_cv(data, v = 5, strata = truth, pool = 0) else rsample::vfold_cv(data, v = 5)
    expected_draws <- runif(5)
    expected_rng <- .Random.seed
    set.seed(93)
    fit <- generics::fit(model(strata = strata), data)
    stopifnot(identical(unname(draws(fit)), expected_draws), identical(.Random.seed, expected_rng),
              identical(lapply(workflows::extract_fit_engine(fit)$folds, `[[`, "assessment"),
                        setNames(lapply(expected$splits, as.integer, data = "assessment"), expected$id)))
  }
  calls <- 0L
  broken <- function(data) { calls <<- calls + 1L; runif(1); stop("resolver failure") }
  saved_rng <- .Random.seed
  error(generics::fit(model(splits = broken, seed = 42), data), "resolver failure")
  stopifnot(calls == 1L, identical(saved_rng, .Random.seed))

  # Invalid engine probabilities and a single engine failure cannot shrink the ensemble.
  for (output in c("missing", "unknown", "duplicate", "nonfinite")) {
    bad_base <- workflows::update_model(base, parsnip::set_engine(spec, "reference", output = output))
    error(generics::fit(calibrated_workflow(bad_base, "truth", method = "cross", splits = plan), data), "probabilit")
  }
  bad_base <- workflows::update_model(base, parsnip::set_engine(spec, "reference", fail = TRUE))
  error(generics::fit(calibrated_workflow(bad_base, "truth", method = "cross", splits = plan), data), "small base workflow training")

  # One-class/singleton calibrators are legal with a capable base engine.
  plain <- workflows::workflow() |> workflows::add_formula(truth ~ x) |> workflows::add_model(spec)
  tiny <- data[c(1, 10), ]
  tiny_plan <- make_plan(tiny, list(1L, 2L), c("one", "two"))
  tiny_fit <- generics::fit(calibrated_workflow(plain, "truth", method = "cross", splits = tiny_plan), tiny)
  stopifnot(all(vapply(workflows::extract_fit_engine(tiny_fit)$calibrator$calibrators, `[[`, integer(1), "observations") == 1L))
  tied <- data
  tied$truth <- factor(rep(c("a b", "a.b"), 20), levels = levels(data$truth))
  tied_base <- workflows::update_model(base, parsnip::set_engine(spec, "reference", output = "tie"))
  for (level in c("first", "second")) {
    fit <- generics::fit(calibrated_workflow(tied_base, "truth", method = "cross", seed = 42, event_level = level), tied)
    event <- levels(tied$truth)[if (level == "first") 1L else 2L]
    for (loss in c("log", "brier")) {
      stopifnot(all(predict(fit, tied, type = "prob", loss = loss) == .5),
                all(predict(fit, tied, type = "class", loss = loss)$.pred_class == event))
    }
  }
}

  testthat::expect_true(TRUE)
})

testthat::test_that("an oversized fold count names the rows it splits", {
  skip_if_no_integration("tune", "yardstick")
  d <- withr::with_seed(3, data.frame(x = stats::rnorm(120), truth = factor(rep(c("yes", "no", "maybe"), 40))))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  binary <- d[d$truth != "maybe", ]
  binary$truth <- droplevels(binary$truth)
  binary_base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  testthat::expect_error(generics::fit(calibrated_workflow(binary_base, "truth", method = "cross", v = 81), binary),
    "^`v` cannot exceed the number of local rows; the fold count is not reduced\\.$")
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", classification = "ovo", method = "cross",
    model_training = "pairwise", v = 81), d),
    "^pair_1 Cross: `v` cannot exceed the number of local pair rows; the fold count is not reduced\\.$")
  # After calibration selection, the shared split runs on the final-fit F.
  ten <- function(data) rsample::make_splits(list(analysis = 1:10, assessment = 11:nrow(data)), data)
  testthat::expect_error(generics::fit(calibrated_workflow(binary_base, "truth", method = "cross", v = 11,
    loss = tune::tune("va_loss"), tuning = list(strategy = "joint", selection_split = ten,
      calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")),
        resamples = function(data) rsample::vfold_cv(data, v = 2)))), binary),
    "^`v` cannot exceed the number of final-fit F rows; the fold count is not reduced\\.$")
})

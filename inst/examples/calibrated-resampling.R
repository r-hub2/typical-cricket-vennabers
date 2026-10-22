# Run against an installed vennabers package. Helpers here are an executable
source(system.file("examples", "assertions.R", package = "vennabers"))
# audit of these fixed-configuration examples, not an evaluation API.
check_outer_splits <- function(outer) {
  ids <- setdiff(names(outer), "splits")
  stopifnot(!anyDuplicated(outer[ids]))
  for (split in outer$splits) {
    a <- as.integer(split, data = "analysis")
    b <- as.integer(split, data = "assessment")
    assert_conditions(length(a) > 0L, length(b) > 0L, !anyDuplicated(a), !anyDuplicated(b),
              !anyNA(c(a, b)), all(c(a, b) == floor(c(a, b))),
              all(c(a, b) >= 1L & c(a, b) <= nrow(split$data)), !length(intersect(a, b)))
  }
  invisible(NULL)
}

check_evaluation <- function(result, outer, event_level = "first", classification = "binary") {
  check_outer_splits(outer)
  regression <- classification == "regression"
  # A collector can itself fail when all fits fail. That is an incomplete run.
  notes <- tune::collect_notes(result)
  predictions <- tune::collect_predictions(result, summarize = FALSE)
  metrics <- tune::collect_metrics(result, summarize = FALSE)
  stopifnot(nrow(notes) == 0L)
  ids <- setdiff(names(outer), "splits")
  config <- "pre0_mod0_post0" # tune 2.1: the one fixed fit_resamples configuration
  expected <- do.call(rbind, lapply(seq_len(nrow(outer)), function(i) {
    rows <- as.integer(outer$splits[[i]], data = "assessment")
    cbind(as.data.frame(outer[rep(i, length(rows)), ids, drop = FALSE]),
          .row = rows, .config = config)
  }))
  expected_metrics <- do.call(rbind, lapply(seq_len(nrow(outer)), function(i) {
    cbind(as.data.frame(outer[rep(i, 2L), ids, drop = FALSE]),
          .metric = if (regression) c("rmse", "mae") else c("brier_class", "mn_log_loss"), .config = config)
  }))
  same_keys <- function(actual, expected) {
    actual <- as.data.frame(actual[names(expected)])
    expected <- as.data.frame(expected)
    stopifnot(!anyNA(actual), !anyDuplicated(actual), nrow(actual) == nrow(expected))
    actual <- actual[do.call(order, actual), , drop = FALSE]
    expected <- expected[do.call(order, expected), , drop = FALSE]
    rownames(actual) <- rownames(expected) <- NULL
    stopifnot(identical(actual, expected))
  }
  same_keys(predictions, expected)
  same_keys(metrics, expected_metrics)
  estimator <- if (regression) "standard" else if (nlevels(outer$splits[[1L]]$data$truth) == 2L) "binary" else "multiclass"
  stopifnot(all(is.finite(metrics$.estimate)), all(metrics$.estimator == estimator))
  for (i in seq_len(nrow(outer))) {
    matches <- function(table) Reduce(`&`, lapply(ids, function(id) table[[id]] == outer[[id]][i]))
    p <- predictions[matches(predictions), ]
    truth <- outer$splits[[i]]$data$truth[p$.row]
    stopifnot(identical(p$truth, truth), !anyNA(truth))
    if (regression) {
      stopifnot(is.double(p$.pred), all(is.finite(p$.pred)), all(is.finite(truth)))
      manual <- c(rmse = sqrt(mean((truth - p$.pred)^2)), mae = mean(abs(truth - p$.pred)))
      m <- metrics[matches(metrics), ]
      stopifnot(all(abs(m$.estimate - manual[m$.metric]) <= 1e-10 + 1e-10 * abs(manual[m$.metric])))
      next
    }
    event <- levels(truth)[if (event_level == "first") 1L else 2L]
    probability <- p[[paste0(".pred_", event)]]
    all_p <- as.matrix(p[paste0(".pred_", levels(truth))])
    stopifnot(all(is.finite(all_p)), all(all_p >= 0 & all_p <= 1),
              all(abs(rowSums(all_p) - 1) <= 1e-10))
    # yardstick binary vector Brier; only the metric clips for log loss.
    clipped <- pmin(pmax(probability, .Machine$double.eps), 1 - .Machine$double.eps)
    manual <- c(brier_class = mean(((truth == event) - probability)^2),
                mn_log_loss = -mean(ifelse(truth == event, log(clipped), log1p(-clipped))))
    if (classification != "binary") {
      colnames(all_p) <- levels(truth)
      true_p <- all_p[cbind(seq_along(truth), as.integer(truth))]
      manual <- c(brier_class = vennabers::brier_score(truth, all_p) / 2,
                  mn_log_loss = -mean(log(pmin(pmax(true_p, .Machine$double.eps), 1 - .Machine$double.eps))))
    }
    m <- metrics[matches(metrics), ]
    stopifnot(all(abs(m$.estimate - manual[m$.metric]) <= 1e-10))
  }
  invisible(list(predictions = predictions, metrics = metrics, notes = notes))
}

if (all(vapply(c("vennabers", "tune", "yardstick", "rsample", "workflows"),
               requireNamespace, logical(1), quietly = TRUE))) {
  development <- data.frame(x = sin(1:120), truth = factor(rep(c("yes", "no", "no", "yes"), 30)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  # Direct explicit calibration uses a bound local rsplit.
  direct_split <- withr::with_seed(101, rsample::initial_split(development, .8))
  direct <- generics::fit(vennabers::calibrated_workflow(base, "truth", splits = direct_split), development)
  stopifnot(nrow(predict(direct, development, type = "prob")) == nrow(development))

  # Seed the outer plan separately. Prepared inner assignments see only analysis.
  outer <- withr::with_seed(102, rsample::vfold_cv(development, v = 3, repeats = 2))
  check_outer_splits(outer)
  prepared <- withr::with_seed(103, lapply(outer$splits, function(split) {
    rsample::initial_split(rsample::analysis(split), .8)
  }))
  resolver <- vennabers::local_calibration_splits(prepared, c("x", "truth"))
  model <- vennabers::calibrated_workflow(base, "truth", splits = resolver,
                                         seed = 104, event_level = "first", loss = "brier")
  # For compact diagnostics use extract = function(w) {
  #   s <- workflows::extract_fit_engine(w)
  #   list(folds = s$folds, parameters = s$parameters, rng_kind = s$rng_kind)
  # }. extract = identity retains every entire fitted outer workflow instead.
  control <- tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine,
                                     event_level = "first", parallel_over = "resamples",
                                     pkgs = c("vennabers", "stats"))
  result <- withr::with_seed(105, tune::fit_resamples(model, outer,
    metrics = yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss), control = control))
  checked <- check_evaluation(result, outer, event_level = "first")
  states <- tune::collect_extracts(result)
  stopifnot(nrow(states) == nrow(outer), all(vapply(states$.extracts, inherits, logical(1), "vennabers_fit")))
}

if (all(vapply(c("vennabers", "tune", "yardstick", "rsample", "workflows"),
               requireNamespace, logical(1), quietly = TRUE))) {
  # Parameters were fixed independently; a final test set remains untouched.
  # Cross keeps each learned recipe/model paired with its held-out calibrator.
  cross <- vennabers::calibrated_workflow(base, "truth", method = "cross", seed = 106)
  fitted_cross <- generics::fit(cross, development)
  new_data <- data.frame(x = c(-.75, .25, .8))
  cross_probability <- predict(fitted_cross, new_data, type = "prob")
  cross_class <- predict(fitted_cross, new_data, type = "class")
  cross_details <- vennabers::venn_abers_details(fitted_cross, new_data, loss = "brier")
  stopifnot(identical(cross_details$prediction, predict(fitted_cross, new_data, type = "prob", loss = "brier")),
            nrow(cross_details$details) == 5L * nrow(new_data))

  # Calibration is not independent evaluation. Refit the full procedure inside
  # every outer analysis set and score only its untouched assessment observations.
  cross_outer <- withr::with_seed(107, rsample::vfold_cv(development, v = 5))
  cross_result <- withr::with_seed(108, tune::fit_resamples(cross, cross_outer,
    metrics = yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss), control = control))
  cross_checked <- check_evaluation(cross_result, cross_outer)
  cross_states <- tune::collect_extracts(cross_result)$.extracts
  stopifnot(length(cross_states) == 5L,
            sum(vapply(cross_states, function(s) length(s$models), integer(1))) == 25L)
  # These extracts retain 25 full workflows, recipes and calibrators. There is no
  # hidden full-data or final fit. Fitting a final Cross workflow costs five more.
  # For prepared inner Cross plans use local_calibration_splits() on a list of
  # vfold_cv(analysis(outer_split), v = 5) rsets, then pass that resolver as splits.
}

if (all(vapply(c("vennabers", "tune", "yardstick", "rsample", "workflows", "nnet"),
               requireNamespace, logical(1), quietly = TRUE))) {
  # One shared multinomial model supplies every OvO pair or OvR class score.
  # Keep final test observations outside this development/evaluation example.
  multiclass_data <- iris
  names(multiclass_data)[5L] <- "truth"
  multiclass_base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ ., multiclass_data) |>
                            recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  multiclass_outer <- withr::with_seed(201, rsample::vfold_cv(multiclass_data, v = 3, repeats = 2, strata = truth))
  for (method in c("inductive", "cross")) for (classification in c("ovo", "ovr")) {
    model <- vennabers::calibrated_workflow(multiclass_base, "truth",
      classification = classification, method = method, seed = 202, precision = 2, loss = "brier")
    fitted <- generics::fit(model, multiclass_data)
    new_data <- multiclass_data[1:3, 1:4]
    probability <- predict(fitted, new_data, type = "prob")
    predicted_class <- predict(fitted, new_data, type = "class")
    details <- vennabers::venn_abers_details(fitted, new_data)
    stopifnot(identical(probability, details$prediction), nrow(predicted_class) == 3L)
    # details$details$probability contains binary subproblem points. OvO couples the
    # pair points; OvR normalizes its class-against-rest points to sum to one.
    # Cross first merges each binary subproblem across matching folds (Log/Brier),
    # then couples or normalizes. Never average completed multiclass fold tables.
    state <- workflows::extract_fit_engine(fitted)
    direct_plan <- if (method == "cross") rsample::manual_rset(lapply(state$folds, function(f) {
      rsample::make_splits(f[c("analysis", "assessment")], multiclass_data)
    }), names(state$folds)) else
      rsample::make_splits(state$folds$Inductive1[c("analysis", "assessment")], multiclass_data)
    explicit <- generics::fit(vennabers::calibrated_workflow(multiclass_base, "truth",
      classification = classification, method = method, splits = direct_plan, seed = 202, precision = 2, loss = "brier"), multiclass_data)
    stopifnot(identical(probability, predict(explicit, new_data, type = "prob")))
    result <- withr::with_seed(203, tune::fit_resamples(model, multiclass_outer,
      metrics = yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss),
      control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine,
        parallel_over = "resamples", pkgs = c("vennabers", "nnet"))))
    checked <- check_evaluation(result, multiclass_outer, classification = classification)
    # The audit checks every repeat/resample/configuration/row and both metrics.
    # Yardstick Brier is half the manual matrix score, even with two classes.
    stopifnot(nrow(tune::collect_extracts(result)) == nrow(multiclass_outer))
  }
}

if (all(vapply(c("vennabers", "tune", "yardstick", "rsample", "workflows", "recipes"),
               requireNamespace, logical(1), quietly = TRUE))) {
  regression_data <- data.frame(x = sin(1:100), truth = 10 + 3 * sin(1:100) + cos(1:100))
  regression_data$x[c(5, 20, 40)] <- NA_real_
  regression_base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x, regression_data) |>
      recipes::step_impute_median(x) |> recipes::step_normalize(x)) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  regression_outer <- withr::with_seed(301, rsample::vfold_cv(regression_data, v = 5))
  for (method in c("inductive", "cross")) {
    model <- vennabers::calibrated_workflow(regression_base, "truth", method = method,
                                           epsilon = .3, seed = 302)
    fitted <- generics::fit(model, regression_data)
    new_data <- data.frame(x = c(-.8, NA, .7))
    point <- predict(fitted, new_data, type = "numeric")
    details <- vennabers::venn_abers_details(fitted, new_data)
    stopifnot(identical(point, details$prediction))
    state <- workflows::extract_fit_engine(fitted)
    # state$folds retains positions, sizes and requested/effective moderation.
    # IVAR uses weighted points, including the zero-weight midpoint rule.
    # Its bounds concern regularized conditional expectations, not individual
    # outcome coverage. Cross averages fold points and each bound separately;
    # details$aggregate is descriptive, with no added coverage guarantee.
    # Calibration need not improve RMSE/MAE. Evaluate the complete procedure on
    # outer assessment data excluded from all its training and calibration.
    result <- withr::with_seed(303, tune::fit_resamples(model, regression_outer,
      metrics = yardstick::metric_set(yardstick::rmse, yardstick::mae),
      control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine,
        parallel_over = "resamples", pkgs = c("vennabers", "stats"))))
    checked <- check_evaluation(result, regression_outer, classification = "regression")
    print(tune::collect_metrics(result)) # standard yardstick RMSE/MAE in outcome units
    # Prepared local plans use the same resolver as classification workflows.
    prepared <- withr::with_seed(304, lapply(regression_outer$splits, function(s) {
      local <- rsample::analysis(s)
      if (method == "cross") rsample::vfold_cv(local, v = 5) else rsample::initial_split(local, .8)
    }))
    explicit_model <- vennabers::calibrated_workflow(regression_base, "truth", method = method,
      splits = vennabers::local_calibration_splits(prepared, names(regression_data)), epsilon = .3)
    explicit_result <- tune::fit_resamples(explicit_model, regression_outer,
      metrics = yardstick::metric_set(yardstick::rmse, yardstick::mae),
      control = tune::control_resamples(save_pred = TRUE))
    check_evaluation(explicit_result, regression_outer, classification = "regression")
  }
}

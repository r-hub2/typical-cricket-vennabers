testthat::test_that("matched-last-fit", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
packages <- c("tune", "yardstick", "rsample", "workflows", "recipes", "nnet")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  for (path in c("binary", "ovo", "ovr", "regression")) for (method in c("inductive", "cross")) {
    regression <- path == "regression"
    raw <- data.frame(x = sin(1:120), truth = if (regression) cos(1:120) else
      factor(rep(if (path == "binary") c("class", "id") else c("class", "id", ".pred"), length.out = 120)))
    schema <- do.call(venn_abers_schema, c(list(data = raw, outcome = "truth"),
      if (!regression) list(classification = path)))
    data <- venn_abers_encode(raw, schema)
    outcome <- schema$encoded_outcome
    spec <- if (regression) parsnip::linear_reg(engine = "lm") else
      if (path == "binary") parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(penalty = 0) |> parsnip::set_engine("nnet", trace = FALSE)
    base <- workflows::workflow() |> workflows::add_recipe(
      recipes::recipe(reformulate(names(data)[1L], outcome), data) |>
        recipes::step_normalize(recipes::all_numeric_predictors())) |> workflows::add_model(spec)
    make_model <- function(splits = NULL) calibrated_workflow(base, outcome, schema = schema,
      method = method, prediction_model = "matched", seed = 42, splits = splits)
    split <- rsample::make_splits(list(analysis = 1:90, assessment = 91:120), data)
    metrics <- if (regression) yardstick::metric_set(yardstick::rmse, yardstick::mae) else
      yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
    check <- function(result) {
      native <- tune::collect_predictions(result)
      decoded <- venn_abers_decode(native, schema)
      measured <- tune::collect_metrics(result)
      assert_conditions(nrow(tune::collect_notes(result)) == 0L, nrow(native) == 30L,
        identical(native$.row, 91:120), !anyDuplicated(native[c("id", ".row", ".config")]),
        all(native$id == result$id), all(native$.config == "pre0_mod0_post0"),
        identical(decoded$truth, raw$truth[91:120]), nrow(measured) == 2L,
        !anyDuplicated(measured[c(".metric", ".config")]), all(measured$.config == "pre0_mod0_post0"),
        all(is.finite(measured$.estimate)))
      if (regression) {
        expected <- c(rmse = sqrt(mean((decoded$truth - decoded$prediction)^2)),
          mae = mean(abs(decoded$truth - decoded$prediction)))
      } else {
        probability <- if (path == "binary") decoded$probability[, 1L] else decoded$probability
        expected <- c(brier_class = yardstick::brier_class_vec(decoded$truth, probability),
          mn_log_loss = yardstick::mn_log_loss_vec(decoded$truth, probability))
      }
      stopifnot(setequal(measured$.metric, names(expected)),
        all(abs(measured$.estimate - expected[measured$.metric]) <=
          1e-10 + if (regression) 1e-10 * abs(expected[measured$.metric]) else 0))
      list(decoded = decoded)
    }
    result <- tune::last_fit(make_model(), split, metrics = metrics)
    checked <- check(result)
    fit <- tune::extract_workflow(result)
    state <- workflows::extract_fit_engine(fit)
    direct <- generics::fit(make_model(), rsample::analysis(split))
    held <- rsample::assessment(split)
    stopifnot(identical(venn_abers_details(fit, held), venn_abers_details(direct, held)),
      identical(checked$decoded$identity$.row, 91:120),
      identical(state$schema, schema), identical(names(state$models), names(state$folds)))
    predicted <- venn_abers_decode(venn_abers_details(direct, held), schema)
    stopifnot(identical(checked$decoded$prediction, predicted$prediction),
      identical(unname(checked$decoded$probability), unname(predicted$probability)))
    local <- rsample::analysis(split)
    plans <- lapply(state$folds, function(f) rsample::make_splits(f[c("analysis", "assessment")], local))
    for (id in names(plans)) {
      recipe <- workflows::extract_recipe(state$models[[id]])
      stopifnot(abs(recipe$steps[[1L]]$means[[1L]] - mean(local[[1L]][state$folds[[id]]$analysis])) <= 1e-10)
    }
    plan <- if (method == "cross") rsample::manual_rset(plans, names(plans)) else plans[[1L]]
    prepared <- make_model(local_calibration_splits(list(plan), names(data)))
    replay <- tune::last_fit(prepared, split, metrics = metrics)
    check(replay)
    stopifnot(identical(tune::collect_predictions(result), tune::collect_predictions(replay)),
      identical(tune::collect_metrics(result), tune::collect_metrics(replay)))
    # Held-out targets cannot influence fitted models/calibration or native/decoded predictions.
    changed <- data
    changed[[outcome]][91:120] <- rev(changed[[outcome]][91:120])
    altered <- rsample::make_splits(list(analysis = 1:90, assessment = 91:120), changed)
    other <- tune::extract_workflow(tune::last_fit(make_model(), altered, metrics = metrics))
    stopifnot(identical(venn_abers_details(fit, held), venn_abers_details(other, held)),
      identical(state$calibration, workflows::extract_fit_engine(other)$calibration))
  }
}

  testthat::expect_true(TRUE)
})

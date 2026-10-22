testthat::test_that("schema-resampling", {
library(vennabers)
if (all(vapply(c("tune", "yardstick", "rsample", "workflows", "recipes"),
               requireNamespace, logical(1), quietly = TRUE))) {
  example <- new.env()
  sys.source(system.file("examples", "encoded-binary.R", package = "vennabers"), example)
  check <- example$check_encoded_evaluation
  original <- example$development
  error <- function(expr, pattern = "") {
    result <- tryCatch(suppressWarnings(force(expr)), error = identity)
    stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
  }
  for (method in c("inductive", "cross")) for (event in c("first", "second")) for (prediction_model in c("matched", "full")) {
    schema <- venn_abers_schema(original, ".pred_class", event_level = event)
    data <- venn_abers_encode(original, schema)
    outcome <- schema$encoded_outcome
    base <- workflows::workflow() |>
      workflows::add_recipe(recipes::recipe(reformulate(names(data)[1], outcome), data) |>
                              recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(parsnip::logistic_reg(engine = "glm"))
    old <- example$original_outer
    outer <- rsample::manual_rset(lapply(old$splits, function(s) {
      rsample::make_splits(list(analysis = as.integer(s, data = "analysis"),
                               assessment = as.integer(s, data = "assessment")), data)
    }), old[c("id", "id2")])
    stopifnot(identical(outer[c("id", "id2")], old[c("id", "id2")]))
    control <- tune::control_resamples(save_pred = TRUE, extract = identity, event_level = event,
                                       allow_par = FALSE, parallel_over = "resamples", pkgs = c("vennabers", "stats"))
    scoring <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
    run <- function(model, metrics = scoring, resamples = outer) {
      withr::with_seed(111, tune::fit_resamples(model, resamples, metrics = metrics, control = control))
    }
    model <- calibrated_workflow(base, outcome, schema = schema, method = method, prediction_model = prediction_model, seed = 42, loss = "brier", precision = 2)
    automatic <- run(model)
    checked <- check(automatic, outer, schema)
    full <- tune::collect_extracts(automatic)$.extracts
    stopifnot(nrow(checked$decoded$identity) == 2L * nrow(original))
    prepared <- lapply(seq_along(full), function(i) {
      state <- workflows::extract_fit_engine(full[[i]])
      local_data <- rsample::analysis(outer$splits[[i]])
      held <- rsample::assessment(outer$splits[[i]])
      direct <- venn_abers_decode(venn_abers_details(full[[i]], held), schema)
      stopifnot(all(direct$details$prediction_model == prediction_model),
                identical(direct$details$model_source, if (prediction_model == "full") rep("Full1", nrow(direct$details)) else direct$details$fold))
      rows <- checked$decoded$identity$id == outer$id[i] & checked$decoded$identity$id2 == outer$id2[i]
      stopifnot(identical(state$schema, schema),
                identical(unname(direct$probability), unname(checked$decoded$probability[rows, ])))
      split <- lapply(names(state$folds), function(id) {
        f <- state$folds[[id]]
        source <- if (prediction_model == "full") "Full1" else id
        learned <- workflows::extract_recipe(state$models[[source]])
        train <- if (prediction_model == "full") seq_len(nrow(local_data)) else f$analysis
        stopifnot(abs(learned$steps[[1]]$means[[1]] - mean(local_data[[1]][train])) <= 1e-10)
        rsample::make_splits(f[c("analysis", "assessment")], local_data)
      })
      if (method == "cross") rsample::manual_rset(split, names(state$folds)) else split[[1]]
    })
    resolver <- local_calibration_splits(prepared, names(data))
    seen <- list()
    callback <- function(local) {
      stopifnot(setequal(names(local), names(data)), identical(levels(local[[outcome]]), unname(schema$classes)))
      seen[[length(seen) + 1L]] <<- local[[1]]
      resolver(local)
    }
    planned <- calibrated_workflow(base, outcome, schema = schema, method = method, prediction_model = prediction_model, splits = callback,
                                    seed = 42, loss = "brier", precision = 2)
    explicit <- run(planned)
    check(explicit, outer, schema)
    stopifnot(identical(seen, lapply(outer$splits, function(s) rsample::analysis(s)[[1]])),
              identical(tune::collect_predictions(automatic, summarize = FALSE), tune::collect_predictions(explicit, summarize = FALSE)),
              identical(tune::collect_metrics(automatic, summarize = FALSE), tune::collect_metrics(explicit, summarize = FALSE)))
    # All-fit failure and metric failures remain visible to ordinary collectors.
    global <- if (method == "cross") rsample::vfold_cv(data, v = 2) else rsample::initial_split(data)
    bad <- suppressWarnings(run(calibrated_workflow(base, outcome, schema = schema, method = method, prediction_model = prediction_model, splits = global)))
    error(check(bad, outer, schema))
    stopifnot(nrow(tune::collect_notes(bad)) == nrow(outer))
    fail_metric <- yardstick::new_prob_metric(function(data, ...) stop("metric failed"), direction = "minimize")
    bad <- suppressWarnings(run(model, yardstick::metric_set(yardstick::brier_class, fail_metric)))
    error(check(bad, outer, schema))
    stopifnot(nrow(tune::collect_notes(bad)) > 0L)
    # Partial failure cannot become a successful mean over surviving assessments.
    first <- rsample::analysis(outer$splits[[1]])[[1]]
    partial <- function(local) {
      if (identical(local[[1]], first)) stop("local failure")
      resolver(local)
    }
    bad <- suppressWarnings(run(calibrated_workflow(base, outcome, schema = schema, method = method, prediction_model = prediction_model, splits = partial)))
    error(check(bad, outer, schema))
    stopifnot(nrow(tune::collect_notes(bad)) > 0L, nrow(tune::collect_predictions(bad)) < nrow(checked$decoded$identity))
    # Missing or duplicated rows, candidate identity and metric removal are caught.
    for (change in c("missing", "duplicate", "config", "metric", "truth")) {
      bad <- automatic
      if (change == "missing") bad$.predictions[[1]] <- bad$.predictions[[1]][-1, ]
      if (change == "duplicate") bad$.predictions[[1]] <- rbind(bad$.predictions[[1]], bad$.predictions[[1]][1, ])
      if (change == "config") bad$.predictions[[1]]$.config <- "unexpected"
      if (change == "metric") bad$.metrics[[1]] <- bad$.metrics[[1]][-1, ]
      if (change == "truth") bad$.predictions[[1]][[outcome]][1] <- NA
      error(check(bad, outer, schema))
    }
  }
}

  testthat::expect_true(TRUE)
})

testthat::test_that("schema-workflow", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")
if (all(vapply(c("parsnip", "workflows", "rsample", "recipes", "tune", "yardstick"),
               requireNamespace, logical(1), quietly = TRUE))) {
  error <- function(expr, pattern = "") {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
  }
  original <- data.frame(x = sin(1:80), y = factor(rep(c("class", ".pred_yes", ".pred_yes", "class"), 20),
                                                 levels = c("class", ".pred_yes")))
  names(original) <- c("vaV0001", ".pred_class")
  schema <- venn_abers_schema(original, ".pred_class", event_level = "second")
  data <- venn_abers_encode(original, schema)
  outcome <- schema$encoded_outcome
  base <- workflows::workflow() |>
    workflows::add_formula(reformulate(setdiff(names(data), outcome), outcome)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  model <- calibrated_workflow(base, outcome, schema = schema, seed = 42)
  fit <- generics::fit(model, data)
  state <- workflows::extract_fit_engine(fit)
  prob <- predict(fit, data, type = "prob")
  decoded <- venn_abers_decode(prob, schema)
  assert_conditions(identical(state$schema, schema), identical(state$event_class, unname(schema$classes[[schema$event_class]])),
            identical(state$class_levels, unname(schema$classes)),
            identical(venn_abers_details(fit, data)$prediction, prob),
            identical(decoded$probability[, 1], prob[[1]]))
  error(calibrated_workflow(base, outcome, schema = schema, event_level = "first"), "schema")
  error(calibrated_workflow(base, outcome, schema = schema, classification = "ovo"), "schema")
  error(calibrated_workflow(base, ".pred_class", schema = schema), "schema")
  bad <- data; bad[[outcome]] <- factor(bad[[outcome]], levels = rev(levels(bad[[outcome]])))
  error(generics::fit(model, bad), "schema")
  # The encoded target must retain both factor type and the schema's ordered status.
  for (ordered_target in c(FALSE, TRUE)) {
    contract_data <- original
    contract_data[[schema$outcome]] <- factor(contract_data[[schema$outcome]], ordered = ordered_target)
    contract_schema <- venn_abers_schema(contract_data, schema$outcome)
    contract_data <- venn_abers_encode(contract_data, contract_schema)
    contract_model <- calibrated_workflow(base, outcome, schema = contract_schema, seed = 42)
    contract_fit <- generics::fit(contract_model, contract_data)
    trained <- workflows::extract_fit_engine(contract_fit)$models[[1L]]
    stopifnot(identical(is.ordered(workflows::extract_mold(trained)$outcomes[[1L]]), ordered_target))
    bad <- contract_data
    bad[[outcome]] <- factor(bad[[outcome]], ordered = !ordered_target)
    error(generics::fit(contract_model, bad), "schema")
    error(predict(contract_fit, bad, type = "prob"), "schema")
    error(venn_abers_details(contract_fit, bad), "schema")
  }
  bad <- data; bad[[outcome]] <- unclass(bad[[outcome]])
  error(predict(fit, bad, type = "prob"), "schema")
  error(generics::fit(model, data[1]), "truth|outcome|column")
  foreign <- schema; foreign$event_level <- "first"; foreign$event_class <- "class"
  bad <- data; attr(bad, "vennabers_schema") <- foreign
  error(generics::fit(model, bad), "different schema")
  error(predict(fit, bad, type = "prob"), "different schema")
  error(venn_abers_details(fit, bad), "different schema")
  foreign_plan <- rsample::make_splits(list(analysis = 1:60, assessment = 61:80), bad)
  error(generics::fit(calibrated_workflow(base, outcome, schema = schema, splits = foreign_plan), data), "different schema")
  foreign_resolver <- local_calibration_splits(list(foreign_plan), names(data))
  error(generics::fit(calibrated_workflow(base, outcome, schema = schema, splits = foreign_resolver), data), "different schema")
  regression <- workflows::update_model(base, parsnip::linear_reg(engine = "lm"))
  error(calibrated_workflow(regression, outcome, schema = schema), "mode")
  error(venn_abers_decode(prob, foreign), "different schema")
  error(venn_abers_decode(venn_abers_details(fit, data), foreign), "different schema")
}
if (exists("model")) {
  # A public deterministic engine makes probability swaps and leakage observable.
  fits <- 0L
  schema_fit <- function(x, y) {
    fits <<- fits + 1L
    list(x = x, y = y, levels = levels(y))
  }
  schema_prob <- function(object, new_data) {
    p <- new_data[[1L]]
    stats::setNames(data.frame(p, 1 - p), paste0(".pred_", object$levels))[2:1]
  }
  register_reference_engine("va_schema_reference", "classification", "schema_fit", "schema_prob", "prob")
  reference <- parsnip::new_model_spec("va_schema_reference", args = list(), eng_args = NULL,
                                      mode = "classification", method = NULL, engine = "reference")
  class_names <- c("class", ".pred_yes", "a b", "a-b", "a.b", "1yes", "ümlaut", "vaL0001")
  outcomes <- c(".pred", ".pred_class", ".row", ".config", "id", "id2")
  for (i in seq_along(class_names)) for (method in c("inductive", "cross")) for (event in c("first", "second")) {
    original <- data.frame(x = seq(.013, .987, length.out = 40), row = 1:40,
                          y = factor(rep(c(1, 2, 2, 1), 10), levels = 1:2,
                                     labels = class_names[c(i, if (i == 8L) 1L else i + 1L)]))
    names(original) <- c("vaV0001", "row identity", outcomes[(i - 1L) %% length(outcomes) + 1L])
    schema <- venn_abers_schema(original, names(original)[3], event_level = event)
    data <- venn_abers_encode(original, schema)
    outcome <- schema$encoded_outcome
    base <- workflows::workflow() |>
      workflows::add_formula(reformulate(setdiff(names(data), outcome), outcome)) |>
      workflows::add_model(reference)
    # Every precision is exercised with both losses, event choices and methods.
    precisions <- if (i == 1L) list(NULL, 0, 2, -1) else list(2)
    for (precision in precisions) {
      model <- calibrated_workflow(base, outcome, schema = schema, method = method, seed = 42, precision = precision)
      fit <- generics::fit(model, data)
      state <- workflows::extract_fit_engine(fit)
      stopifnot(!state$parameters$explicit[["event_level"]], identical(state$schema, schema))
      new <- data[c(1, 9, 27), 1:2]
      # Compare to unchanged manual cores using identical model outputs and plans.
      ps <- ys <- fs <- list()
      new_p <- list()
      for (id in names(state$folds)) {
        f <- state$folds[[id]]
        engine <- workflows::extract_fit_engine(state$models[[id]])
        stopifnot(identical(as.numeric(engine$x[[2]]), as.numeric(data[[2]][f$analysis])),
                  !any(engine$x[[2]] %in% data[[2]][f$assessment]))
        ps[[id]] <- predict(state$models[[id]], data[f$assessment, ], type = "prob")[[paste0(".pred_", state$event_class)]]
        ys[[id]] <- data[[outcome]][f$assessment] == state$event_class
        fs[[id]] <- rep(id, length(f$assessment))
        new_p[[id]] <- predict(state$models[[id]], new, type = "prob")[[paste0(".pred_", state$event_class)]]
      }
      core <- if (method == "cross") cross_venn_abers(unlist(ys), unlist(ps), unlist(fs), precision = precision) else
        venn_abers(ys[[1]], ps[[1]], precision = precision)
      for (loss in c("log", "brier")) {
        expected <- predict(core, if (method == "cross") do.call(cbind, new_p) else new_p[[1]], loss = loss)$probability
        native <- predict(fit, new, type = "prob", loss = loss)
        decoded <- venn_abers_decode(native, schema)
        details <- venn_abers_details(fit, new, loss = loss)
        original_details <- venn_abers_decode(details, schema)
        stopifnot(max(abs(decoded$probability[, schema$event_class] - expected)) <= 1e-10,
                  identical(original_details$details[c("row", "fold", "lower", "upper", "probability")],
                            details$details[c("row", "fold", "lower", "upper", "probability")]),
                  all(original_details$details$event_class == schema$event_class),
                  identical(decoded$probability, original_details$probability))
        for (n in 0:1) {
          small <- new[seq_len(n), ]
          p <- venn_abers_decode(predict(fit, small, type = "prob", loss = loss), schema)
          cl <- venn_abers_decode(predict(fit, small, type = "class", loss = loss), schema)
          de <- venn_abers_decode(venn_abers_details(fit, small, loss = loss), schema)
          assert_conditions(identical(dim(p$probability), c(n, 2L)), is.double(p$probability),
                    identical(levels(cl$class), schema$class_levels), length(cl$class) == n,
                    identical(p$probability, de$probability))
        }
      }
      splits <- lapply(state$folds, function(f) rsample::make_splits(f[c("analysis", "assessment")], data))
      plan <- if (method == "cross") rsample::manual_rset(splits, names(state$folds)) else splits[[1]]
      explicit <- generics::fit(calibrated_workflow(base, outcome, schema = schema, method = method,
                                                   splits = plan, precision = precision), data)
      stopifnot(identical(predict(explicit, new, type = "prob"), predict(fit, new, type = "prob")))
    }
    # Balanced ties choose the event, including after decoding.
    tied <- data; tied[[1]] <- .5
    fit <- generics::fit(calibrated_workflow(base, outcome, schema = schema, method = method, seed = 42), tied)
    for (loss in c("log", "brier")) {
      stopifnot(all(predict(fit, tied, type = "prob", loss = loss) == .5),
                all(venn_abers_decode(predict(fit, tied, type = "class", loss = loss), schema)$class == schema$event_class))
    }
  }
}

  testthat::expect_true(TRUE)
})

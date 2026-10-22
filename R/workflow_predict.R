# `sources` maps every fold slot to the model it predicts with, so full-model
# reuse is just a repeated lookup instead of a `full` branch.
.source_predictions <- function(models, sources, new_data, predict_one, empty_one) {
  by_model <- lapply(names(models), function(id) {
    if (nrow(new_data)) predict_one(models[[id]], id) else empty_one()
  })
  names(by_model) <- names(models)
  stats::setNames(by_model[sources], names(sources))
}

.base_probability <- function(model, data, classes, event, id = "Inductive1") {
  p <- tryCatch(stats::predict(model, data, type = "prob"),
                error = function(e) stop(id, " base workflow probability prediction for original classes ",
                  paste(sQuote(classes), collapse = ", "), " and input columns ",
                  paste(sQuote(names(data)), collapse = ", "), ": ", conditionMessage(e), call. = FALSE))
  .validate_base_probabilities(p, nrow(data), classes, id)
  if (!is.null(event)) return(as.double(p[[paste0(".pred_", event)]]))
  p <- as.matrix(p[paste0(".pred_", classes)])
  colnames(p) <- classes
  # The matrix core normalizes once; repeating it can perturb tied scores.
  p
}

.validate_base_probabilities <- function(p, n, classes, id) {
  .validate_names(names(p), paste(id, "base probability column"))
  expected <- paste0(".pred_", classes)
  if (is.data.frame(p) && nrow(p) != n) {
    stop(id, " base workflow probability prediction returned ", nrow(p),
         " rows for ", n, " input rows.", call. = FALSE)
  }
  if (!is.data.frame(p) || !setequal(names(p), expected) ||
      !all(vapply(p, .valid_probability, logical(1))) ||
      any(abs(rowSums(p) - 1) > 1e-8)) {
    stop(id, " base workflow must return finite probabilities for original classes ",
         paste(sQuote(classes), collapse = ", "), " in columns ",
         paste(sQuote(expected), collapse = ", "), ".", call. = FALSE)
  }
  invisible(NULL)
}

.base_numeric <- function(model, data, id) {
  p <- tryCatch(stats::predict(model, data, type = "numeric"),
                error = function(e) stop(id, " base workflow numeric prediction: ", conditionMessage(e), call. = FALSE))
  if (!is.data.frame(p) || !identical(names(p), ".pred") || nrow(p) != nrow(data)) {
    stop(id, " base workflow must return one .pred column with one value per row.", call. = FALSE)
  }
  .validate_regression_vector(p$.pred, paste(id, "base regression prediction"))
  as.double(p$.pred)
}

.integration_prediction <- function(object, new_data, include_details = FALSE) {
  requested <- object$parameters$requested
  post <- requested$postprocessing
  if (!is.null(post) && post$method %in% c("tailor_logistic", "tailor_multinomial", "tailor_linear") &&
      !requireNamespace("tailor", quietly = TRUE)) {
    stop("Postprocessing prediction requires optional package 'tailor'.", call. = FALSE)
  }
  pairwise <- identical(requested$model_training, "pairwise")
  if (pairwise) .validate_pair_components(object)
  parts <- lapply(.fit_components(object), function(part) .predict_component(object, part, new_data, include_details))
  result <- if (pairwise) .couple_pairs(object, parts) else parts[[1L]]
  if (include_details) result$details <- tibble::as_tibble(result$details)
  result
}

# The read-side view of a fit as its components, matching .fit_domains(): a
# shared fit is its own single component, a pairwise fit stores one per class
# pair. Both shapes are saved since v0.9.0.
.fit_components <- function(state) {
  if (identical(state$parameters$requested$model_training, "pairwise")) state$pair_models else list(state)
}

# One component's calibrated prediction from its own source models; a pair
# component is binary over its two classes and checks each model's pair.
.predict_component <- function(object, part, new_data, include_details) {
  requested <- object$parameters$requested
  pair <- !is.null(part$pair_id)
  task <- if (pair) .binary_task() else .task_for(requested)
  classes <- part$class_levels
  event <- part$event_class
  p <- .source_predictions(part$models, part$sources, new_data, function(model, id) {
    if (pair) .validate_pair_model(model, classes, .domain_label(part, id))
    task$pred(model, new_data, .domain_label(part, id), classes, event)
  }, function() task$empty(classes))
  result <- task$predict(part$calibrator, p, include_details, classes, event, object$loss)
  if (!include_details) return(result)
  details <- result$details
  if (!is.null(requested$model_training)) {
    details$model_training <- rep(requested$model_training, nrow(details))
  }
  mode <- requested$prediction_model
  if (!is.null(mode)) {
    details$prediction_model <- rep(mode, nrow(details))
    details$model_source <- unname(part$sources[details$fold])
    if (!is.null(result$aggregate)) {
      result$aggregate$prediction_model <- rep(mode, nrow(new_data))
      result$aggregate$model_source <- rep(list(names(part$models)), nrow(new_data))
    }
  }
  result$details <- if (pair) .pair_details(part, details, mode) else details
  result
}

#' @export
print.vennabers_fit <- function(x, ...) {
  settings <- x$parameters$requested
  cat("Venn-Abers", settings$method, "calibrated predictor\n")
  if (!is.null(settings$model_training)) cat("Model training:", settings$model_training, "\n")
  if (!is.null(settings$prediction_model)) {
    cat("Prediction model:", settings$prediction_model,
        if (settings$prediction_model == "full") "(full-model transfer)" else "", "\n")
  }
  components <- .fit_components(x)
  source_fits <- sum(vapply(components, function(part) {
    if (is.null(part$postprocessing)) length(part$folds) + as.integer(identical(settings$prediction_model, "full")) else
      sum(vapply(part$postprocessing, `[[`, 1L, "observed_base_fits"))
  }, 1L))
  retained <- sum(vapply(components, function(part) sum(vapply(part$models, function(model) {
    if (inherits(model, "vennabers_ensemble")) length(model$models) else 1L
  }, 1L)), 1L))
  calibration_rows <- sum(vapply(components, function(part) {
    sum(vapply(part$calibration, function(values) length(values$row), 1L))
  }, 0))
  cat("Base fits:", source_fits, "\n")
  cat("Models per nonempty prediction:", retained, "\n")
  cat("Stored calibration rows:", calibration_rows, "(held-out truth and model predictions)\n")
  # Only size-requested plans store U; older fits and calibration_prop plans have none.
  for (part in components) {
    sized <- part$folds$Inductive1
    if (!is.null(sized$unassigned)) {
      cat(c(part$pair_id, "Inductive rows: training"), length(sized$analysis), "| calibration", length(sized$assessment),
          "| unassigned", length(sized$unassigned), "\n")
    }
  }
  for (part in components) for (fold in names(part$sources)) {
    cat(" ", part$pair_id, fold, "->", part$pair_id, part$sources[[fold]], "\n")
  }
  invisible(x)
}

# A fit's predictions in parsnip's shape: `.pred` for regression, otherwise
# one `.pred_<class>` column per stored class in their order.
.predict_fit <- function(object, new_data) UseMethod(".predict_fit")

.predict_fit.vennabers_fit <- function(object, new_data) .integration_prediction(object, new_data)$prediction

#' @export
predict_vennabers_numeric <- function(object, new_data) .predict_fit(object, new_data)$.pred

#' @export
predict_vennabers_prob <- function(object, new_data) .predict_fit(object, new_data)

#' @export
predict_vennabers_class <- function(object, new_data) {
  .class_decision(.predict_fit(object, new_data), object$class_levels, object$event_class)
}

#' @export
predict._vennabers_fit <- function(object, new_data, type = NULL, opts = list(), ..., loss = NULL) {
  .validate_dots(...)
  if (!is.null(object$fit$schema)) .schema_encoded_data(new_data, object$fit$schema)
  if (!missing(loss)) object$fit <- .with_loss(object$fit, loss)
  if (!is.null(type)) .validate_choice(type,
    if (object$fit$parameters$requested$mode == "regression") "numeric" else c("prob", "class"), "type")
  .bind_schema(parsnip::predict.model_fit(object, new_data, type = type, opts = opts), object$fit$schema)
}

# A prediction-time loss replaces the fitted one; regression has none.
.with_loss <- function(state, loss) {
  if (state$parameters$requested$mode == "regression") {
    stop("Regression does not accept `loss`.", call. = FALSE)
  }
  .validate_choice(loss, c("log", "brier"), "loss")
  state$loss <- loss
  state
}

#' @export
venn_abers_details <- function(object, new_data, ..., loss = NULL) {
  .validate_dots(...)
  .integration_packages()
  frozen <- inherits(object, "vennabers_frozen")
  state <- if (frozen) object else workflows::extract_fit_engine(object)
  if (!inherits(state, "vennabers_fit")) {
    stop("`object` must be a fitted calibrated_workflow() or calibrate_workflow().", call. = FALSE)
  }
  if (!frozen && !is.null(state$schema)) .schema_encoded_data(new_data, state$schema)
  if (!missing(loss)) state <- .with_loss(state, loss)
  predictors <- if (frozen) .frozen_data(state, new_data) else
    hardhat::forge(new_data, blueprint = workflows::extract_mold(object)$blueprint)$predictors
  result <- .integration_prediction(state, predictors, include_details = TRUE)
  # Postprocessing records are keyed by their source model, prefixed "<pair_id>/" in a pair.
  components <- .fit_components(state)
  key <- function(part, id) paste(c(part$pair_id, id), collapse = "/")
  postprocessing <- do.call(c, unname(lapply(components, function(part) {
    if (is.null(part$postprocessing)) return(NULL)
    stats::setNames(part$postprocessing, vapply(names(part$postprocessing), key, "", part = part, USE.NAMES = FALSE))
  })))
  if (!is.null(postprocessing)) {
    result$postprocessing <- postprocessing
    # Each detail row reads the record of its component's source for its fold.
    owners <- result$details[["pair_id"]]
    if (is.null(owners)) owners <- rep(1L, nrow(result$details))
    ids <- vapply(seq_along(owners), function(r) {
      part <- components[[owners[[r]]]]
      key(part, part$sources[[result$details$fold[[r]]]])
    }, "")
    result$details$postprocessing_method <- rep(state$parameters$requested$postprocessing$method, length(ids))
    result$details$postprocessing_input <- vapply(ids, function(id) postprocessing[[id]]$input, "")
    # Ensemble records keep their pairs under inner_folds, as saved since v0.9.0.
    if (state$parameters$requested$postprocessing$strategy == "ensemble") {
      result$details$inner_fold <- lapply(ids, function(id) names(postprocessing[[id]]$inner_folds))
      for (size in c("B_size", "P_size")) result$details[[size]] <- lapply(ids, function(id) {
        vapply(postprocessing[[id]]$inner_folds, `[[`, 1L, size)
      })
    } else {
      result$details$inner_fold <- rep("Holdout1", length(ids))
      result$details$B_size <- vapply(ids, function(id) postprocessing[[id]]$B_size, 1L)
      result$details$P_size <- vapply(ids, function(id) postprocessing[[id]]$P_size, 1L)
    }
  }
  result$prediction <- .bind_schema(result$prediction, state$schema)
  .bind_schema(result, state$schema)
}

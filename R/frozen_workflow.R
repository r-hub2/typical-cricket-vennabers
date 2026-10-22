#' @export
calibrate_workflow <- function(fitted_workflow, data_calibration, outcome, ...,
                               classification = "binary", event_level = "first",
                               precision = NULL, loss = "log", m = NULL, epsilon = NULL,
                               schema = NULL) {
  .exact_argument_names()
  .validate_dots(...)
  .integration_packages()
  explicit <- c(classification = !missing(classification), event_level = !missing(event_level),
                precision = !missing(precision), loss = !missing(loss), m = !missing(m), epsilon = !missing(epsilon))
  if (!is.null(schema)) {
    options <- .schema_options(schema, outcome, classification, event_level, explicit)
    classification <- options$classification
    event_level <- options$event_level
  }
  .validate_choice(classification, c("binary", "ovo", "ovr"), "classification")
  .validate_choice(event_level, c("first", "second"), "event_level")
  .validate_choice(loss, c("log", "brier"), "loss")
  .validate_names(outcome, "outcome")
  if (length(outcome) != 1L) stop("`outcome` must name one column.", call. = FALSE)
  if (!inherits(fitted_workflow, "workflow") || !workflows::is_trained_workflow(fitted_workflow)) {
    stop("`fitted_workflow` must be a fully fitted workflow.", call. = FALSE)
  }
  fit <- workflows::extract_fit_parsnip(fitted_workflow)
  mold <- workflows::extract_mold(fitted_workflow)
  if (!inherits(fit, "model_fit") || is.null(fit$fit) || inherits(fit$fit, "try-error") ||
      !fit$spec$mode %in% c("classification", "regression")) {
    stop("`fitted_workflow` must contain a fitted classification or regression model.", call. = FALSE)
  }
  mode <- fit$spec$mode
  regression <- mode == "regression"
  multiclass <- !regression && classification != "binary"
  .validate_task_options(regression, classification, explicit)
  targets <- mold$blueprint$ptypes$outcomes
  if (!identical(names(targets), outcome)) {
    stop("Frozen workflow outcome must be ", paste(sQuote(names(targets)), collapse = ", "),
         "; received ", sQuote(outcome), ".", call. = FALSE)
  }
  classes <- if (!regression) levels(targets[[outcome]]) else NULL
  if (!regression) {
    .validate_names(classes, "frozen class")
    if (length(classes) < 2L || (!multiclass && length(classes) != 2L) || !identical(fit$lvl, classes)) {
      stop("Frozen workflow must retain ", if (multiclass) "at least" else "exactly",
           " two classes in the trained outcome order.", call. = FALSE)
    }
  }
  pre <- workflows::extract_preprocessor(fitted_workflow)
  .schema_binding(fitted_workflow, schema)
  if (inherits(pre, "recipe")) .schema_binding(pre$template, schema)
  if (!is.null(schema)) {
    learned_columns <- if (inherits(pre, "recipe")) pre$var_info$variable else
      c(names(mold$blueprint$ptypes$predictors), names(targets))
    if (!identical(mode, schema$mode) || !setequal(learned_columns, unname(schema$columns)) ||
        !identical(classes, if (!regression) unname(schema$classes) else NULL) ||
        !identical(is.ordered(targets[[outcome]]), schema$ordered)) {
      stop("Frozen workflow mode, input columns and trained class order must agree with `schema`; ",
           "original outcome ", sQuote(schema$outcome), ". No learned state is renamed.", call. = FALSE)
    }
  }
  if (inherits(pre, "recipe") && !recipes::fully_trained(workflows::extract_recipe(fitted_workflow))) {
    stop("Frozen workflow recipe must be fully trained.", call. = FALSE)
  }
  # A trained workflow flag alone does not prove an attached postprocessor is fitted.
  post <- tryCatch(workflows::extract_postprocessor(fitted_workflow), error = function(e) NULL)
  if (!is.null(post)) {
    post <- workflows::extract_postprocessor(fitted_workflow, estimated = TRUE)
    if (nrow(generics::tidy(post)) && !tailor::tailor_fully_trained(post)) {
      stop("Frozen workflow postprocessor must be fully trained.", call. = FALSE)
    }
  }
  settings <- list(mode = mode, method = "inductive", classification = classification,
                   event_level = event_level, precision = precision, loss = loss, m = m, epsilon = epsilon)
  effective <- settings
  if (multiclass) effective$event_level <- NULL
  if (regression) effective[c("classification", "event_level", "precision", "loss")] <- NULL
  event <- if (regression || multiclass) NULL else .event_class(classes, event_level)
  state <- structure(list(parameters = list(requested = settings, effective = effective),
                          schema = schema, outcome = outcome,
                          predictors = names(mold$blueprint$ptypes$predictors),
                          class_levels = classes,
                          event_class = event,
                          loss = effective$loss, models = list(Inductive1 = fitted_workflow),
                          # A frozen fit has exactly one source model; the identity map keeps it
                          # on the same `sources` contract as matched/full calibrated_workflow() fits.
                          sources = c(Inductive1 = "Inductive1")),
                     class = c("vennabers_frozen", "vennabers_fit"))
  input <- .frozen_data(state, data_calibration, "calibration")
  if (!nrow(input)) stop("Frozen calibration data must be nonempty.", call. = FALSE)
  truth <- data_calibration[[outcome]]
  if (regression) {
    .validate_regression_vector(truth, "frozen regression truth")
  } else if (!is.factor(truth) || anyNA(truth) || !setequal(levels(truth), classes)) {
    stop("Frozen calibration outcome ", sQuote(outcome), " must be a nonmissing factor with trained classes ",
         paste(sQuote(classes), collapse = ", "), ".", call. = FALSE)
  }
  if (regression) {
    resolved <- .regression_moderation(nrow(input), m, epsilon)
    prediction <- .base_numeric(fitted_workflow, input, "Frozen calibration")
    state$calibrator <- .venn_abers_regression(truth, prediction, m, epsilon)
    state$parameters$effective$m <- resolved$m
    state$parameters$effective$epsilon <- resolved$epsilon
  } else {
    truth <- factor(truth, levels = classes)
    if (multiclass) .multiclass_eligibility(truth, classification)
    probability <- .base_probability(fitted_workflow, input, classes, event, "Frozen calibration")
    state$calibrator <- if (multiclass) venn_abers(truth, probability, precision = precision, method = classification) else
      venn_abers(truth == event, probability, precision = precision)
  }
  state$columns <- c(state$predictors, outcome)
  state$calibration_size <- nrow(input)
  packages <- unique(c("vennabers", names(.integration_versions), generics::required_pkgs(fitted_workflow)))
  state$versions <- c(R = as.character(getRversion()),
                      vapply(packages, function(p) as.character(utils::packageVersion(p)), ""))
  state
}

.frozen_data <- function(object, data, stage = "prediction") {
  .schema_binding(data, object$schema)
  if (!is.null(object$schema)) {
    # Match incoming levels by identity; the learned schema still fixes their order.
    if (is.data.frame(data) && object$outcome %in% names(data) &&
        is.factor(data[[object$outcome]]) && setequal(levels(data[[object$outcome]]), object$class_levels)) {
      data[[object$outcome]] <- factor(data[[object$outcome]], levels = object$class_levels,
                                       ordered = object$schema$ordered)
    }
    .schema_encoded_data(data, object$schema, truth = stage == "calibration")
  }
  if (!is.data.frame(data)) stop("Frozen ", stage, " data must be a data frame.", call. = FALSE)
  if (ncol(data)) .validate_names(names(data), paste("frozen", stage, "input column"))
  required <- c(object$predictors, if (stage == "calibration") object$outcome)
  absent <- setdiff(required, names(data))
  if (length(absent)) stop("Frozen ", stage, " is missing original column ",
                           sQuote(absent[1L]), ".", call. = FALSE)
  # Keep raw inputs, including additional recipe roles; only the truth is separated.
  data[setdiff(names(data), object$outcome)]
}

#' @export
print.vennabers_frozen <- function(x, ...) {
  cat("Venn-Abers inductive calibrated predictor\n")
  cat("Frozen workflow adopted unchanged\n")
  cat("Base fits:", 0L, "\n")
  cat("Models per nonempty prediction:", length(x$models), "\n")
  cat("Calibration rows:", x$calibration_size, "(new data; only the calibrator learned)", "\n")
  for (id in names(x$sources)) cat(" ", id, "->", x$sources[[id]], "\n")
  invisible(x)
}

#' @export
predict.vennabers_frozen <- function(object, new_data, type = NULL, ..., loss = NULL) {
  .validate_dots(...)
  .integration_packages()
  regression <- object$parameters$requested$mode == "regression"
  if (is.null(type)) type <- if (regression) "numeric" else "prob"
  .validate_choice(type, if (regression) "numeric" else c("prob", "class"), "type")
  if (!missing(loss)) object <- .with_loss(object, loss)
  new_data <- .frozen_data(object, new_data)
  result <- if (regression) tibble::tibble(.pred = predict_vennabers_numeric(object, new_data)) else
    if (type == "class") {
      # Keep the ordered status of the stored workflow's trained outcome.
      target <- workflows::extract_mold(object$models[[1L]])$blueprint$ptypes$outcomes[[object$outcome]]
      decision <- predict_vennabers_class(object, new_data)
      tibble::tibble(.pred_class = factor(decision, levels = levels(decision), ordered = is.ordered(target)))
    } else predict_vennabers_prob(object, new_data)
  .bind_schema(result, object$schema)
}

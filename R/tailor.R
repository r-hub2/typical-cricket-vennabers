# Tailor compatibility stays here. Never change its adjustments or learned engines.
.configure_tailor_workflow <- function(workflow, method) {
  canonical <- if (method == "linear") tailor::tailor() |>
    tailor::adjust_numeric_calibration(method = method, smooth = FALSE) else
    tailor::tailor() |> tailor::adjust_probability_calibration(method = method, smooth = FALSE)
  post <- tryCatch(workflows::extract_postprocessor(workflow), error = function(e) NULL)
  if (is.null(post)) return(workflows::add_tailor(workflow, canonical))
  if (!inherits(post, "tailor") || length(post$adjustments) != 1L ||
      !identical(post$type, canonical$type) || !is.null(post$columns) || ncol(post$ptype)) {
    stop("Postprocessing requires exactly one untrained canonical Tailor ", method, " step with smooth=FALSE.", call. = FALSE)
  }
  step <- post$adjustments[[1L]]
  expected <- canonical$adjustments[[1L]]
  if (!setequal(names(step$arguments), names(expected$arguments))) {
    stop("Postprocessing Tailor arguments must be exactly method='", method, "', smooth=FALSE.", call. = FALSE)
  }
  # Normalize only a comparison copy; preserve the supplied Tailor unchanged.
  step$arguments <- step$arguments[names(expected$arguments)]
  if (!identical(step, expected)) stop("Postprocessing Tailor step must match the untrained canonical ", method, " step.", call. = FALSE)
  workflow
}

.tailor_input <- function(workflow, data, classes) {
  numeric <- is.null(classes)
  prediction <- stats::predict(workflow, data, type = if (numeric) "numeric" else "prob")
  if (numeric) {
    if (!is.data.frame(prediction) || !identical(names(prediction), ".pred") || nrow(prediction) != nrow(data)) {
      stop("Tailor numeric input requires one .pred value per row.", call. = FALSE)
    }
    .validate_regression_vector(prediction$.pred, "Tailor numeric input")
  } else {
    .validate_base_probabilities(prediction, nrow(data), classes, "Tailor input")
    prediction <- prediction[paste0(".pred_", classes)]
    prediction$.pred_class <- .class_decision(prediction, classes)
    prediction <- prediction[c(".pred_class", paste0(".pred_", classes))]
  }
  prediction
}

# Tailor has no public learned-model extractor. Keep the read-only layout check
# at this boundary, and fail explicitly if an update changes the supported state.
.tailor_fit_info <- function(post, method, classes, size) {
  if (!inherits(post, "tailor") || !tailor::tailor_fully_trained(post) || length(post$adjustments) != 1L) {
    stop("P has no learned ", method, " state; failed fits and switches to none are not allowed.", call. = FALSE)
  }
  step <- post$adjustments[[1L]]
  fit <- step$results$fit
  if (!identical(step$arguments$method, method) || !inherits(fit, paste0("cal_estimate_", method)) ||
      length(fit$estimates) != 1L || !is.list(fit$estimates[[1L]])) {
    stop("P has no learned ", method, " state; unsupported Tailor state, failed fits and switches to none are not allowed.", call. = FALSE)
  }
  estimate <- fit$estimates[[1L]]$estimate
  rank <- if (method == "linear") 2L else length(classes)
  if (!inherits(estimate, if (method == "multinomial") "multinom" else "glm") ||
      !identical(estimate$rank, rank) || size <= rank) {
    stop("P requires full design rank ", rank, " and more rows than rank.", call. = FALSE)
  }
  converged <- if (method == "multinomial") isTRUE(estimate$convergence == 0L) else isTRUE(estimate$converged)
  coefficients <- stats::coef(estimate)
  if (!converged || !length(coefficients) || !all(is.finite(coefficients))) {
    stop("P ", method, " fit did not converge to finite coefficients.", call. = FALSE)
  }
  list(rank = rank, coefficients = coefficients)
}

.fit_tailor_workflow <- function(workflow, data, data_calibration, domain, method, context) {
  post <- workflows::extract_postprocessor(workflow)
  base <- generics::fit(.without_postprocessor(workflow), data)
  prediction <- .tailor_input(base, data_calibration, domain$classes)
  prediction[[domain$outcome]] <- data_calibration[[domain$outcome]]
  numeric <- is.null(domain$classes)
  fitted <- generics::fit(post, prediction, outcome = tidyselect::all_of(domain$outcome),
    estimate = tidyselect::all_of(if (numeric) ".pred" else ".pred_class"),
    probabilities = tidyselect::all_of(if (numeric) character() else paste0(".pred_", domain$classes)))
  info <- .tailor_fit_info(fitted, method, domain$classes, nrow(data_calibration))
  model <- structure(list(workflow = base, tailor = fitted, specification = post,
    method = method, classes = domain$classes, context = context), class = "vennabers_tailor_workflow")
  c(list(model = model), info)
}

#' @export
predict.vennabers_tailor_workflow <- function(object, new_data, type = if (is.null(object$classes)) "numeric" else "class", ...) {
  .validate_dots(...)
  numeric <- is.null(object$classes)
  .validate_choice(type, if (numeric) "numeric" else c("prob", "class"), "type")
  columns <- if (numeric) ".pred" else if (type == "class") ".pred_class" else paste0(".pred_", object$classes)
  if (!nrow(new_data)) {
    result <- if (numeric) data.frame(.pred = numeric()) else
      as.data.frame(stats::setNames(rep(list(numeric()), length(object$classes)), paste0(".pred_", object$classes)))
    if (!numeric) result$.pred_class <- factor(character(), levels = object$classes)
    return(tibble::as_tibble(result[columns]))
  }
  tryCatch({
    prediction <- .tailor_input(object$workflow, new_data, object$classes)
    # probably expects a matrix, but nnet drops dimensions for one row. Pad
    # only the prediction batch; invoke the base workflow and Tailor once each.
    if (object$method == "multinomial" && nrow(prediction) == 1L) prediction <- prediction[c(1L, 1L), , drop = FALSE]
    result <- stats::predict(object$tailor, prediction)
    if (!is.data.frame(result) || nrow(result) != nrow(prediction)) {
      stop("Tailor output requires one prediction per input row.", call. = FALSE)
    }
    result <- result[seq_len(nrow(new_data)), , drop = FALSE]
    if (numeric) .validate_regression_vector(result$.pred, "Tailor numeric output") else
      .validate_base_probabilities(result[paste0(".pred_", object$classes)], nrow(new_data), object$classes, "Tailor output")
    result[columns]
  }, error = function(e) stop(object$context, ": ", conditionMessage(e), call. = FALSE))
}

#' @exportS3Method hardhat::extract_fit_engine
extract_fit_engine.vennabers_tailor_workflow <- function(x, ...) workflows::extract_fit_engine(x$workflow, ...)

#' @exportS3Method hardhat::extract_fit_parsnip
extract_fit_parsnip.vennabers_tailor_workflow <- function(x, ...) workflows::extract_fit_parsnip(x$workflow, ...)

#' @exportS3Method hardhat::extract_recipe
extract_recipe.vennabers_tailor_workflow <- function(x, ...) workflows::extract_recipe(x$workflow, ...)

#' @exportS3Method hardhat::extract_preprocessor
extract_preprocessor.vennabers_tailor_workflow <- function(x, ...) workflows::extract_preprocessor(x$workflow, ...)

#' @exportS3Method hardhat::extract_tailor
extract_tailor.vennabers_tailor_workflow <- function(x, ..., estimated = TRUE) {
  .validate_dots(...)
  if (!is.logical(estimated) || length(estimated) != 1L || is.na(estimated)) stop("`estimated` must be TRUE or FALSE.", call. = FALSE)
  if (estimated) x$tailor else x$specification
}

#' @exportS3Method hardhat::extract_postprocessor
extract_postprocessor.vennabers_tailor_workflow <- function(x, ..., estimated = FALSE) {
  extract_tailor.vennabers_tailor_workflow(x, ..., estimated = estimated)
}

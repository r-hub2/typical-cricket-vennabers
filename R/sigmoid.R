.validate_sigmoid_workflow <- function(workflow, input) {
  post <- tryCatch(workflows::extract_postprocessor(workflow), error = function(e) NULL)
  if (!is.null(post)) stop("sigmoid B requires a base workflow without a postprocessor.", call. = FALSE)
  spec <- workflows::extract_spec_parsnip(workflow)
  if (spec$mode != "classification") stop("sigmoid B requires binary classification.", call. = FALSE)
  if (input == "glm_link") {
    if (!inherits(spec, "logistic_reg") || !identical(spec$engine, "glm")) {
      stop("sigmoid glm_link B requires a binomial GLM workflow with logit link.", call. = FALSE)
    }
    family <- parsnip::translate(spec)$method$fit$args$family
    family <- rlang::eval_tidy(family)
    if (is.character(family)) family <- get(family, envir = asNamespace("stats"), mode = "function")
    if (is.function(family)) family <- family()
    if (!identical(family$family, "binomial") || !identical(family$link, "logit")) {
      stop("sigmoid glm_link B requires a binomial GLM workflow with logit link.", call. = FALSE)
    }
  }
  invisible(workflow)
}

.sigmoid_input <- function(workflow, data, input, classes, role) {
  if (input == "probability") {
    p <- stats::predict(workflow, data, type = "prob")
    .validate_base_probabilities(p, nrow(data), classes, paste("sigmoid", role))
    return(as.double(p[[paste0(".pred_", classes[2L])]]))
  }
  # The public workflow path forges the recipe before requesting raw GLM links.
  x <- stats::predict(workflow, data, type = "raw", opts = list(type = "link"))
  .validate_regression_vector(x, paste("sigmoid glm_link", role))
  if (length(x) != nrow(data)) stop("sigmoid ", role, " input must have one value per row.", call. = FALSE)
  as.double(x)
}

.fit_sigmoid_workflow <- function(workflow, data, plan, domain, input, context) {
  model <- generics::fit(workflow, data[plan$B, , drop = FALSE])
  x <- .sigmoid_input(model, data[plan$P, , drop = FALSE], input, domain$classes, "P")
  positive <- data[[domain$outcome]][plan$P] == domain$classes[2L]
  n1 <- sum(positive); n0 <- sum(!positive)
  target <- ifelse(positive, (n1 + 1) / (n1 + 2), 1 / (n0 + 2))
  scale <- max(abs(x))
  if (scale < 30) scale <- 1
  x <- x / scale
  # Platt's unregularized logistic objective with smoothed binary targets.
  loss <- function(ab) {
    z <- -(ab[1L] * x + ab[2L])
    sum(pmax(z, 0) - target * z + log1p(exp(-abs(z))))
  }
  gradient <- function(ab) {
    residual <- stats::plogis(-(ab[1L] * x + ab[2L])) - target
    -c(sum(residual * x), sum(residual))
  }
  rank <- if (all(x == x[1L])) 1L else 2L
  if (rank == 1L) {
    ab <- c(0, -stats::qlogis(mean(target)))
  } else {
    pgtol <- 1e-6
    # sklearn's tolerances stop short of the exact optimum, so its result depends
    # on the optimizer path: lmm = 10 follows SciPy's L-BFGS-B memory (maxcor).
    fit <- stats::optim(c(0, log((n0 + 1) / (n1 + 1))), loss, gradient,
      method = "L-BFGS-B", control = list(factr = 64, pgtol = pgtol, maxit = 1000, lmm = 10))
    finite <- all(is.finite(fit$par)) && is.finite(fit$value)
    # With sklearn's ftol (factr = 64) the line search can stop at the optimum
    # because no representable step decreases the loss (code 52); sklearn keeps
    # that point. Accept it only when stationary: mean gradient per P row <= pgtol.
    stationary <- finite && fit$convergence == 52L &&
      isTRUE(max(abs(gradient(fit$par))) <= pgtol * length(x))
    if (!finite || !(fit$convergence == 0L || stationary)) {
      stop("sigmoid P fit did not converge to finite coefficients.", call. = FALSE)
    }
    ab <- fit$par
  }
  coefficients <- stats::setNames(ab, c("a", "b"))
  composite <- structure(list(workflow = model, input = input, classes = domain$classes,
    coefficients = coefficients, scale = scale, context = context), class = "vennabers_sigmoid")
  list(model = composite, rank = rank, coefficients = coefficients, scale = scale)
}

#' @export
predict.vennabers_sigmoid <- function(object, new_data, type = "prob", ...) {
  .validate_dots(...)
  .validate_choice(type, "prob", "type")
  tryCatch({
    x <- if (nrow(new_data)) .sigmoid_input(object$workflow, new_data, object$input, object$classes, "prediction") else numeric()
    p <- stats::plogis(-(object$coefficients[1L] * (x / object$scale) + object$coefficients[2L]))
    result <- stats::setNames(data.frame(1 - p, p), paste0(".pred_", object$classes))
    .validate_base_probabilities(result, nrow(new_data), object$classes, "sigmoid output")
    result
  }, error = function(e) stop(object$context, ": ", conditionMessage(e), call. = FALSE))
}

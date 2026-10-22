#' @export
cross_venn_abers_regression <- function(truth, prediction, fold, ..., m = 1, epsilon = NULL) {
  if (length(list(...))) stop("`...` must be empty; unsupported fit argument.", call. = FALSE)
  .validate_regression_vector(truth, "truth")
  .validate_regression_vector(prediction, "prediction")
  if (!length(truth) || length(prediction) != length(truth)) {
    stop("`truth` and `prediction` must have the same nonzero length.", call. = FALSE)
  }
  fold <- .calibration_folds(fold, length(truth))
  # Explicit options follow the venn_abers_regression() contract.
  given <- c(m = !missing(m), epsilon = !missing(epsilon))
  .validate_moderation(m, epsilon, given)
  .cross_venn_abers_regression(truth, prediction, fold, if (given[["m"]]) m, if (given[["epsilon"]]) epsilon)
}

.cross_venn_abers_regression <- function(truth, prediction, fold, m = NULL, epsilon = NULL) {
  identities <- unique(fold)
  calibrators <- lapply(identities, function(id) {
    rows <- fold == id
    tryCatch(.venn_abers_regression(truth[rows], prediction[rows], m, epsilon),
             error = function(error) stop("Fold ", sQuote(id), ": ", conditionMessage(error), call. = FALSE))
  })
  names(calibrators) <- identities
  structure(list(folds = identities, calibrators = calibrators,
                 observations = length(truth)), class = "cross_venn_abers_regression")
}

# Mean of a nonempty finite numeric vector, using only double arithmetic.
# Cancel opposite signs before scaling and divide after compensated summation
# to retain small residuals and equal subnormals. Scale only if the remaining
# same-sign sum overflows; the mean stays finite. A mean below the subnormal
# range can round to zero; this is not a correctly-rounded summation guarantee.
.regression_mean <- function(values) {
  values <- sort(values)
  left <- 1L
  right <- length(values)
  total <- correction <- 0
  while (left <= right) {
    # Choose the opposite sign when available, so overflow cannot precede
    # cancellation. Sorting also makes the result independent of fold order.
    if (total > 0) {
      value <- values[left]
      left <- left + 1L
    } else {
      value <- values[right]
      right <- right - 1L
    }
    next_total <- total + value
    if (!is.finite(next_total)) {
      total <- next_total
      break
    }
    correction <- correction + if (abs(total) >= abs(value)) {
      (total - next_total) + value
    } else {
      (value - next_total) + total
    }
    total <- next_total
  }
  result <- (total + correction) / length(values)
  if (is.finite(result)) return(result)
  scale <- max(abs(values))
  .regression_mean(values / scale) * scale
}

#' @export
predict.cross_venn_abers_regression <- function(object, newdata, ..., type = "response") {
  if (length(list(...))) stop("`...` must be empty; unsupported prediction argument.", call. = FALSE)
  .validate_choice(type, c("response", "details"), "type")
  if (!is.matrix(newdata) || !is.numeric(newdata) || any(!is.finite(newdata))) {
    stop("`newdata` must be a numeric matrix of finite values.", call. = FALSE)
  }
  columns <- colnames(newdata)
  if (is.null(columns) || anyNA(columns) || any(!nzchar(columns)) ||
      anyDuplicated(columns) || !setequal(columns, object$folds)) {
    stop("`newdata` column names must be non-empty, unique, and match the fitted fold IDs.", call. = FALSE)
  }
  folds <- lapply(object$folds, function(id) {
    stats::predict(object$calibrators[[id]], newdata[, id])
  })
  names(folds) <- object$folds
  prediction <- lapply(c("prediction", "lower", "upper"), function(column) {
    values <- do.call(cbind, lapply(folds, `[[`, column))
    vapply(seq_len(nrow(newdata)), function(i) .regression_mean(values[i, ]), 0)
  })
  names(prediction) <- c("prediction", "lower", "upper")
  prediction <- as.data.frame(prediction)
  if (type == "response") return(prediction)
  list(prediction = prediction, folds = folds)
}

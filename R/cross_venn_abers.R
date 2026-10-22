#' @export
cross_venn_abers <- function(truth, probability, fold, precision = NULL, method = "ovo") {
  .validate_choice(method, c("ovo", "ovr"), "method")
  if (is.factor(truth)) {
    return(.fit_multiclass(truth, probability, precision, method,
                           fold = .calibration_folds(fold, length(truth))))
  }
  if (!missing(method)) stop("`method` requires factor truth and a class matrix.", call. = FALSE)
  .validate_binary_input(truth, probability)
  fold <- .calibration_folds(fold, length(truth))
  identities <- unique(fold)
  calibrators <- lapply(identities, function(id) {
    rows <- fold == id
    venn_abers(truth[rows], probability[rows], precision = precision)
  })
  names(calibrators) <- identities
  structure(
    list(folds = identities, calibrators = calibrators, observations = length(truth)),
    class = "cross_venn_abers"
  )
}

.calibration_folds <- function(fold, observations) {
  if (!is.vector(fold) ||
      !(is.character(fold) || is.numeric(fold)) || anyNA(fold)) {
    stop("`fold` must be a character or numeric vector without missing IDs.", call. = FALSE)
  }
  if (is.character(fold)) {
    if (any(!nzchar(fold))) {
      stop("`fold` labels must be non-empty.", call. = FALSE)
    }
  } else if (any(!is.finite(fold) | fold <= 0 | fold != floor(fold))) {
    stop("`fold` numeric IDs must be finite positive whole numbers.", call. = FALSE)
  }
  if (length(fold) != observations) {
    stop("`fold` must have the same length as `truth`.", call. = FALSE)
  }
  fold <- as.character(fold)
  if (length(unique(fold)) < 2L) {
    stop("`fold` must contain at least two non-empty folds.", call. = FALSE)
  }
  fold
}

.cross_probability <- function(lower, upper) {
  stats::plogis(rowMeans(log(upper)) - rowMeans(log1p(-lower)))
}

#' @export
predict.cross_venn_abers <- function(object, newdata, ..., loss = "log", type = "response") {
  if (length(list(...))) {
    stop("`...` must be empty; unsupported prediction argument.", call. = FALSE)
  }
  .validate_choice(loss, c("log", "brier"), "loss")
  .validate_choice(type, c("response", "class", "details"), "type")
  if (!is.null(object$classes)) return(.predict_multiclass(object, newdata, loss, type))
  if (!is.matrix(newdata) || !is.numeric(newdata) ||
      any(!is.finite(newdata) | newdata < 0 | newdata > 1)) {
    stop("`newdata` must be a numeric matrix of finite values between 0 and 1.", call. = FALSE)
  }
  columns <- colnames(newdata)
  if (is.null(columns) || anyNA(columns) || any(!nzchar(columns)) ||
      anyDuplicated(columns) || !setequal(columns, object$folds)) {
    stop("`newdata` column names must be non-empty, unique, and match the fitted fold IDs.", call. = FALSE)
  }
  endpoints <- lapply(seq_along(object$folds), function(k) {
    stats::predict(object$calibrators[[k]], newdata[, object$folds[k]], loss = loss)
  })
  lower <- do.call(cbind, lapply(endpoints, `[[`, "lower"))
  upper <- do.call(cbind, lapply(endpoints, `[[`, "upper"))
  probability <- if (loss == "brier") {
    rowMeans(do.call(cbind, lapply(endpoints, `[[`, "probability")))
  } else {
    .cross_probability(lower, upper)
  }
  if (type == "class") return(.event_chosen(probability))
  prediction <- data.frame(probability = probability)
  if (type == "response") return(prediction)
  list(
    prediction = prediction,
    folds = data.frame(
      row = rep(seq_len(nrow(newdata)), times = length(object$folds)),
      fold = rep(object$folds, each = nrow(newdata)),
      lower = as.vector(lower),
      upper = as.vector(upper),
      probability = unlist(lapply(endpoints, `[[`, "probability"), use.names = FALSE)
    )
  )
}

#' @export
print.cross_venn_abers <- function(x, ...) {
  if (!is.null(x$classes)) {
    cat("<cross_venn_abers: ", x$method, ">\nClasses: ",
        paste(x$classes, collapse = ", "), "\nFolds: ", length(x$folds),
        "\nCalibration observations: ", x$observations, "\n", sep = "")
    return(invisible(x))
  }
  cat(
    "<cross_venn_abers>\n",
    "Folds: ", length(x$folds), "\n",
    "Calibration observations: ", x$observations, "\n",
    sep = ""
  )
  invisible(x)
}

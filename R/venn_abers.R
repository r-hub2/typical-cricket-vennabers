.cross <- function(x, y, a, b, c) {
  (x[b] - x[a]) * (y[c] - y[b]) - (y[b] - y[a]) * (x[c] - x[b])
}

.upper_scan <- function(weight, events, upper = 1) {
  k <- length(weight)
  x <- c(-1, 0, cumsum(weight))
  y <- c(-upper, 0, cumsum(events))
  stack <- integer(k + 2L)
  stack[1:2] <- 1:2
  size <- 2L

  for (i in seq_len(k)) {
    point <- i + 2L
    while (size > 1L && .cross(x, y, stack[size - 1L], stack[size], point) <= 0) {
      size <- size - 1L
    }
    size <- size + 1L
    stack[size] <- point
  }

  stack[seq_len(size)] <- rev(stack[seq_len(size)])
  result <- numeric(k)
  for (i in seq_len(k)) {
    top <- stack[size]
    next_top <- stack[size - 1L]
    result[i] <- (y[next_top] - y[top]) / (x[next_top] - x[top])

    point <- i + 1L
    x[point] <- x[i] + x[i + 2L] - x[point]
    y[point] <- y[i] + y[i + 2L] - y[point]
    if (.cross(x, y, top, next_top, point) >= 0) next

    size <- size - 1L
    while (size > 1L && .cross(x, y, point, stack[size], stack[size - 1L]) <= 0) {
      size <- size - 1L
    }
    size <- size + 1L
    stack[size] <- point
  }
  result
}

#' @export
venn_abers <- function(truth, probability, precision = NULL, method = "ovo") {
  .validate_choice(method, c("ovo", "ovr"), "method")
  if (is.factor(truth)) {
    return(.fit_multiclass(truth, probability, precision, method))
  }
  if (!missing(method)) stop("`method` requires factor truth and a class matrix.", call. = FALSE)
  .validate_binary_input(truth, probability)
  .validate_precision(precision)
  if (!is.null(precision)) probability <- round(probability, digits = precision)

  ordered <- order(probability, method = "radix")
  probability <- probability[ordered]
  truth <- as.numeric(truth[ordered])
  runs <- rle(probability)
  ends <- cumsum(runs$lengths)
  starts <- c(1L, ends[-length(ends)] + 1L)
  events <- vapply(seq_along(starts), function(i) sum(truth[starts[i]:ends[i]]), 0)
  weight <- as.numeric(runs$lengths)

  upper <- .upper_scan(weight, events)
  lower <- 1 - rev(.upper_scan(rev(weight), rev(weight - events)))
  structure(
    list(
      knots = runs$values,
      lower = c(0, lower),
      upper = c(upper, 1),
      observations = length(truth),
      events = sum(truth)
    ),
    class = "venn_abers"
  )
}

.regularized_probability <- function(lower, upper, loss) {
  if (loss == "brier") upper + (lower^2 - upper^2) / 2 else upper / (1 - lower + upper)
}

# The Klassenentscheidung on `prob` columns in `classes` order: with an `event`,
# the binary event at a probability of at least 0.5; without one, the first
# stored class among the highest probabilities.
.class_decision <- function(prob, classes, event = NULL) {
  prob <- as.matrix(prob)
  chosen <- if (is.null(event)) max.col(prob, ties.method = "first") else {
    at <- match(event, classes)
    ifelse(.event_chosen(prob[, at]), at, 3L - at)
  }
  factor(classes[chosen], levels = classes)
}

# The binary event class: the first or second of the two stored classes.
.event_class <- function(classes, event_level) classes[if (event_level == "first") 1L else 2L]

# The binary half of the Klassenentscheidung on bare event probabilities.
.event_chosen <- function(probability) probability >= .5

#' @export
predict.venn_abers <- function(object, newdata, ..., loss = "log", type = "response") {
  if (length(list(...))) {
    stop("`...` must be empty; unsupported prediction argument.", call. = FALSE)
  }
  .validate_choice(loss, c("log", "brier"), "loss")
  if (!is.null(object$classes)) return(.predict_multiclass(object, newdata, loss, type))
  .validate_choice(type, c("response", "class"), "type")
  if (!.valid_probability(newdata)) {
    stop("`newdata` must be a numeric vector of finite values between 0 and 1.", call. = FALSE)
  }

  lower <- object$lower[findInterval(newdata, object$knots) + 1L]
  upper <- object$upper[findInterval(newdata, object$knots, left.open = TRUE) + 1L]
  probability <- .regularized_probability(lower, upper, loss)
  if (type == "class") return(.event_chosen(probability))
  data.frame(
    lower = lower,
    upper = upper,
    probability = probability
  )
}

#' @export
print.venn_abers <- function(x, ...) {
  if (!is.null(x$classes)) {
    method <- if (identical(x$method, "ovr")) "ovr" else "ovo"
    cat("<venn_abers: ", method, ">\nCalibration observations: ", x$observations,
        "\nClasses: ", paste(x$classes, collapse = ", "),
        if (method == "ovr") "\nClasses against rest: " else "\nClass pairs: ",
        length(x$calibrators), "\n", sep = "")
    return(invisible(x))
  }
  cat(
    "<venn_abers>\n",
    "Calibration observations: ", x$observations, "\n",
    "Event observations: ", x$events, "\n",
    "Unique probability knots: ", length(x$knots), "\n",
    sep = ""
  )
  invisible(x)
}

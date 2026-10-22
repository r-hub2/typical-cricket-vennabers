#' @export
brier_score <- function(truth, probability) {
  if (is.factor(truth)) {
    probability <- .validate_multiclass_input(truth, probability)
    probability[cbind(seq_along(truth), as.integer(truth))] <-
      probability[cbind(seq_along(truth), as.integer(truth))] - 1
    return(mean(rowSums(probability^2)))
  }
  .validate_binary_input(truth, probability)
  mean((truth - probability)^2)
}

#' @export
log_loss <- function(truth, probability) {
  if (is.factor(truth)) {
    probability <- .validate_multiclass_input(truth, probability)
    return(-mean(log(probability[cbind(seq_along(truth), as.integer(truth))])))
  }
  .validate_binary_input(truth, probability)
  -mean(ifelse(truth == 1, log(probability), log1p(-probability)))
}

#' @export
calibration_plot <- function(truth, probability, bins = 10, ...) {
  if (is.factor(truth)) {
    probability <- .validate_multiclass_input(truth, probability)
    classes <- levels(truth)
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old))
    rows <- ceiling(sqrt(length(classes)))
    graphics::par(mfrow = c(rows, ceiling(length(classes) / rows)))
    summaries <- lapply(classes, function(class) {
      options <- utils::modifyList(list(main = paste("Calibration:", class)), list(...))
      summary <- do.call(calibration_plot, c(list(
        truth = truth == class, probability = probability[, class], bins = bins
      ), options))
      data.frame(class = factor(rep(class, nrow(summary)), levels = classes), summary)
    })
    return(invisible(do.call(rbind, summaries)))
  }
  .validate_binary_input(truth, probability)
  # The 0.9.0 rule, looser than .is_count(): attributes and counts beyond the
  # integer range stay valid.
  if (!is.numeric(bins) || !is.null(dim(bins)) || length(bins) != 1L ||
    !is.finite(bins) || bins < 1 || bins != floor(bins)) {
    stop("`bins` must be one finite whole number of at least 1.", call. = FALSE)
  }

  upper_break <- ceiling(probability * bins)
  bin <- pmin(pmax(upper_break + (probability >= upper_break / bins), 1), bins)
  observations <- split(seq_along(probability), bin)
  summary <- data.frame(
    mean_probability = vapply(observations, function(i) mean(probability[i]), 0),
    event_rate = vapply(observations, function(i) mean(truth[i]), 0),
    count = vapply(observations, length, 0L),
    row.names = NULL
  )

  options <- utils::modifyList(list(
    type = "o",
    xlim = c(0, 1),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    xlab = "Mean predicted probability",
    ylab = "Observed event rate",
    main = "Calibration plot"
  ), list(...))
  options$x <- summary$mean_probability
  options$y <- summary$event_rate
  do.call(graphics::plot, options)
  graphics::abline(a = 0, b = 1, lty = 2)
  invisible(summary)
}

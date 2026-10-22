#' @export
regression_plot <- function(truth, prediction, which = "observed", ...) {
  .validate_regression_vector(truth, "truth")
  .validate_regression_vector(prediction, "prediction")
  if (!length(truth) || length(prediction) != length(truth)) {
    stop("`truth` and `prediction` must have the same nonzero length.", call. = FALSE)
  }
  which <- .plot_which(which, c("observed", "residual"))
  table <- data.frame(prediction = as.numeric(prediction), truth = as.numeric(truth))
  table$residual <- table$truth - table$prediction
  # Both diagrams stack in one column over the same predictions.
  if (length(which) > 1L) {
    old <- graphics::par(mfrow = c(2L, 1L))
    on.exit(graphics::par(old))
  }
  diagrams <- list(
    observed = list(y = "truth", ylab = "Observed value", main = "Observed vs predicted",
                    reference = list(a = 0, b = 1)),
    residual = list(y = "residual", ylab = "Residual (observed minus predicted)",
                    main = "Residuals vs predicted", reference = list(h = 0))
  )
  for (diagram in diagrams[which]) {
    options <- utils::modifyList(list(
      xlab = "Predicted value", ylab = diagram$ylab, main = diagram$main
    ), list(...))
    options$x <- table$prediction
    options$y <- table[[diagram$y]]
    do.call(graphics::plot, options)
    do.call(graphics::abline, c(diagram$reference, lty = 2))
  }
  invisible(table)
}

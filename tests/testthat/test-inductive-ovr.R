testthat::test_that("inductive-ovr", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

close <- function(x, y) stopifnot(identical(dim(x), dim(y)), all(abs(x - y) <= 1e-10))

# Each class is the event against all remaining observations. Constant scores
# give bounds n_class / (n + 1), (n_class + 1) / (n + 1).
classes <- c("z", "a", "rest")
truth <- factor(c("z", "a", "a", "rest", "rest", "rest"), levels = classes)
p <- matrix(1 / 3, 6, 3, dimnames = list(NULL, classes))
x <- p[1, , drop = FALSE]
fit <- venn_abers(truth, p, method = "ovr")
details <- predict(fit, x, type = "details")
assert_conditions(identical(names(details), c("prediction", "classes")),
          identical(names(details$classes), c("row", "event", "lower", "upper", "probability")),
          identical(details$classes$event, classes),
          identical(details$classes$row, rep(1L, 3)),
          identical(predict(fit, x, type = "class"), factor("rest", levels = classes)))
close(details$classes$lower, c(1, 2, 3) / 7)
close(details$classes$upper, c(2, 3, 4) / 7)
close(details$classes$probability, c(2, 3, 4) / 8)
close(details$prediction, matrix(c(2, 3, 4) / 9, 1))
brier <- predict(fit, x, loss = "brier", type = "details")
close(brier$classes$probability, c(25, 37, 49) / 98)
close(brier$prediction, matrix(c(25, 37, 49) / 111, 1))
stopifnot(identical(details$classes[3:4], brier$classes[3:4]))
stopifnot(identical(venn_abers(truth, p), venn_abers(truth, p, method = "ovo")),
          any(grepl("<venn_abers: ovr>", capture.output(print(fit)), fixed = TRUE)))

# Evaluate the actual OvR output against hand-worked scores and plot groups.
# The assessment subset can omit global classes while retaining their levels.
assessment <- factor("z", levels = classes)
close(brier_score(assessment, details$prediction), 74 / 81)
close(log_loss(assessment, details$prediction), log(9 / 2))
close(brier_score(assessment, brier$prediction), 11166 / 12321)
close(log_loss(assessment, brier$prediction), log(111 / 25))
local({
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) })
  for (probability in list(details$prediction, brier$prediction)) {
    summary <- calibration_plot(assessment, probability, bins = 2)
    stopifnot(identical(levels(summary$class), classes),
              identical(as.character(summary$class), classes))
    close(summary$mean_probability, as.numeric(probability))
    close(summary$event_rate, c(1, 0, 0))
    close(summary$count, rep(1, 3))
  }
})

# Stable schemas and the first-level tie rule for two through four classes.
# One calibration observation per class is sufficient; row names are retained.
for (k in 2:4) {
  classes <- c("z", "a", "rest", "a/z")[seq_len(k)]
  truth <- factor(classes, levels = classes)
  p <- diag(k)
  colnames(p) <- classes
  x <- matrix(1 / k, 1, k, dimnames = list("new", classes))
  fit <- venn_abers(truth, p, method = "ovr")
  for (loss in c("log", "brier")) {
    close(predict(fit, x, loss = loss), matrix(1 / k, 1, k))
    stopifnot(identical(predict(fit, x[, rev(classes), drop = FALSE], loss = loss, type = "class"),
                        factor("z", levels = classes)))
    for (n in 0:1) {
      small <- x[seq_len(n), , drop = FALSE]
      actual <- predict(fit, small, loss = loss, type = "details")
      assert_conditions(identical(dim(actual$prediction), c(n, k)),
                identical(dimnames(actual$prediction), dimnames(small)),
                identical(actual$prediction, predict(fit, small, loss = loss)),
                identical(names(actual$classes), c("row", "event", "lower", "upper", "probability")),
                identical(actual$classes$row, rep(seq_len(n), k)),
                identical(actual$classes$event, rep(classes, each = n)),
                all(vapply(actual$classes[3:5], is.double, FALSE)),
                identical(predict(fit, small, loss = loss, type = "class"),
                          factor(rep("z", n), levels = classes)))
    }
  }
}

# Native R can round 0.545 differently from NumPy at two digits. All rows enter
# each class fit; rounding follows input normalization and never touches newdata.
classes <- c("z", "a")
truth <- factor(c("z", "a", "z", "a"), levels = classes)
p <- cbind(z = c(.545, .544, .551, .559), a = 1 - c(.545, .544, .551, .559))
x <- cbind(z = c(0, .541, .545, .552, 1), a = 1 - c(0, .541, .545, .552, 1))
rounding_changes_prediction <- FALSE
for (precision in list(NULL, 0, 2, -1)) {
  fit <- venn_abers(truth, p * (1 - 5e-9), precision = precision, method = "ovr")
  normalized <- p * (1 - 5e-9)
  normalized <- normalized / rowSums(normalized)
  for (loss in c("log", "brier")) {
    actual <- predict(fit, x, loss = loss, type = "details")
    for (class in classes) {
      scores <- normalized[, class]
      if (!is.null(precision)) scores <- round(scores, precision)
      binary <- venn_abers(truth == class, scores)
      expected <- predict(binary, x[, class], loss = loss)
      parts <- actual$classes[actual$classes$event == class, ]
      close(as.matrix(parts[3:5]), as.matrix(expected))
      if (!is.null(precision)) {
        rounding_changes_prediction <- rounding_changes_prediction ||
          !identical(expected, predict(binary, round(x[, class], precision), loss = loss))
      }
    }
  }
}
stopifnot(rounding_changes_prediction)

  testthat::expect_true(TRUE)
})

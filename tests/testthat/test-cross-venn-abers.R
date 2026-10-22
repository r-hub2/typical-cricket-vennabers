testthat::test_that("cross-venn-abers", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

# Unequal folds with distinct model columns, from the pinned Python reference.
truth <- c(0, 1, 0, 0, 1, 0, 1, 1)
probability <- c(.1, .2, .3, .2, .2, .4, .4, .8)
fold <- c(rep("A", 3), rep("B", 5))
newdata <- cbind(A = c(0, .2, .25, 1), B = c(1, .4, .2, 0))
calibrator <- cross_venn_abers(truth, probability, fold)
actual <- stats::predict(calibrator, newdata)
assert_conditions(inherits(calibrator, "cross_venn_abers"), identical(names(actual), "probability"),
    is.double(actual$probability), nrow(actual) == nrow(newdata), max(abs(actual$probability - c(0.5,
        0.513167019494862, 0.486832980505138, 0.486832980505138))) <= 1e-10)

expect_error <- function(expr, argument) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"),
            grepl(argument, conditionMessage(error), fixed = TRUE))
}

# Alignment uses vector positions and literal column identities, never names
# attached to calibration vector elements or a prediction column's position.
rows <- c(8, 1, 5, 2, 7, 3, 6, 4)
permuted <- cross_venn_abers(truth[rows], probability[rows], fold[rows])
named <- cross_venn_abers(
  setNames(truth, rev(letters[1:8])),
  setNames(probability, letters[1:8]),
  setNames(fold, rep("same", 8))
)
assert_conditions(identical(stats::predict(calibrator, newdata), actual), identical(stats::predict(calibrator,
    newdata[, 2:1]), actual), max(abs(stats::predict(permuted, newdata)$probability - actual$probability)) <=
    1e-10, identical(stats::predict(named, newdata), actual), identical(stats::predict(calibrator, newdata[c(4,
    1, 3, 2), ])$probability, actual$probability[c(4, 1, 3, 2)]), max(abs(stats::predict(calibrator,
    newdata[, c("B", "A")])$probability - stats::predict(calibrator, `colnames<-`(newdata, c("B", "A")))$probability)) >
    1e-05, identical(stats::predict(calibrator, newdata[FALSE, ]), data.frame(probability = numeric())))

# Supported degenerate inputs: integer/logical truth, integer probabilities,
# numeric fold IDs, singleton and single-class folds.
for (y in list(c(FALSE, TRUE), c(0L, 1L), c(0, 1))) {
  for (ids in list(c(1L, 2L), c(1, 1e20), c(" A", "A"), c("a", "A"), c("1", "01"))) {
    fit <- cross_venn_abers(y, c(0L, 1L), ids)
    input <- matrix(c(0L, 1L, 0L, 1L), 2, dimnames = list(NULL, as.character(ids)))
    out <- stats::predict(fit, input)
    assert_conditions(identical(names(out), "probability"), is.double(out$probability), nrow(out) == 2L,
        all(is.finite(out$probability)), all(out$probability >= 0 & out$probability <= 1))
  }
}
literal <- cross_venn_abers(c(0, 1), c(0, 1), c(" A", "A"))
expect_error(stats::predict(literal, cbind(A = .5, a = .5)), "newdata")

invalid_truth <- list(NULL, numeric(), "0", factor(0), list(0), 0i,
                      NA, NaN, Inf, -Inf, -1, .5, 2, matrix(0), array(0, c(1, 1, 1)))
invalid_probability <- list(NULL, numeric(), TRUE, "0.5", factor(.5), list(.5), .5 + 0i,
                            NA_real_, NaN, Inf, -Inf, -.1, 1.1, matrix(.5), array(.5, c(1, 1, 1)))
invalid_fold <- list(NULL, character(), factor(c("A", "B")), c(TRUE, FALSE),
                     list("A", "B"), c(1i, 2i), c(NA, "B"), c("", "B"),
                     c(1, NA), c(1, NaN), c(1, Inf), c(1, -Inf), c(1, 0),
                     c(1, -1), c(1, 1.5), matrix(c(1, 2)), array(c(1, 2), c(1, 1, 2)),
                     "A", c("A", "A"), c("A", "B", "C"))
for (bad in invalid_truth) {
  expect_error(cross_venn_abers(bad, c(.1, .9), c("A", "B")), "truth")
}
for (bad in invalid_probability) {
  expect_error(cross_venn_abers(c(0, 1), bad, c("A", "B")), "probability")
}
for (bad in invalid_fold) {
  expect_error(cross_venn_abers(c(0, 1), c(.1, .9), bad), "fold")
}
expect_error(cross_venn_abers(truth, probability[-1], fold), "probability")
expect_error(cross_venn_abers(truth, probability, fold, loss = "log"), "loss")
stopifnot(identical(cross_venn_abers(truth, probability, fold, NULL), calibrator))

invalid_newdata <- list(NULL, numeric(), c(.1, .9), list(A = .1, B = .9),
                        as.data.frame(newdata), array(.5, c(1, 2, 1)),
                        matrix(TRUE, 1, 2), matrix("0.5", 1, 2), matrix(.5 + 0i, 1, 2),
                        unname(newdata), newdata[, "A", drop = FALSE],
                        cbind(newdata, extra = .5), matrix(numeric(), 4, 0))
for (bad in invalid_newdata) {
  expect_error(stats::predict(calibrator, bad), "newdata")
}
for (names in list(c("A", "A"), c("A", ""), c("A", NA), c("A", "C"), c("a", "B"))) {
  expect_error(stats::predict(calibrator, `colnames<-`(newdata, names)), "newdata")
}
for (value in c(NA_real_, NaN, Inf, -Inf, -.1, 1.1)) {
  bad <- newdata
  bad[3, 2] <- value
  expect_error(stats::predict(calibrator, bad), "newdata")
}
stopifnot(identical(stats::predict(calibrator, newdata, loss = "log"), actual))
expect_error(stats::predict(calibrator, newdata, loss = "unknown"), "loss")
expect_error(stats::predict(calibrator, newdata, NULL), "...")

three_folds <- cross_venn_abers(c(0, 1, 0, 1, 0, 1), c(.1, .8, .2, .9, .3, .7),
                               rep(c("A", "B", "C"), each = 2))
for (fit in list(calibrator, three_folds)) {
  printed <- capture.output(returned <- withVisible(print(fit)))
  counts <- if (identical(fit, calibrator)) c(2, 8) else c(3, 6)
  stopifnot(identical(printed, c("<cross_venn_abers>", paste("Folds:", counts[1]),
                               paste("Calibration observations:", counts[2]))),
            identical(returned$value, fit), !returned$visible)
}

# Exercise the production merger, not a detached copy of the formula.
# Exact identical-fold reductions catch plausible but wrong finite products.
merge <- getFromNamespace(".cross_probability", "vennabers")
for (case in list(c(1024, .9, .91, 91 / 101), c(4096, .4, .7, 7 / 13),
                  c(2048, .1, .2, 2 / 11), c(2, 0, 1, .5))) {
  for (k in c(1, case[1])) {
    result <- merge(matrix(case[2], 2, k), matrix(case[3], 2, k))
    stopifnot(length(result) == 2L, is.double(result), all(is.finite(result)),
              max(abs(result - case[4])) <= 1e-10)
  }
}

# Existing public consumers accept the regularized output unchanged.
test_truth <- c(0, 1, 1, 0)
stopifnot(is.finite(brier_score(test_truth, actual$probability)),
          is.finite(log_loss(test_truth, actual$probability)))
plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file)
summary <- withVisible(calibration_plot(test_truth, actual$probability))
grDevices::dev.off()
unlink(plot_file)
stopifnot(!summary$visible, sum(summary$value$count) == nrow(newdata),
          identical(names(summary$value), c("mean_probability", "event_rate", "count")))

  testthat::expect_true(TRUE)
})

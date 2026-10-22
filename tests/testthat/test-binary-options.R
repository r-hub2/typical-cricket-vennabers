testthat::test_that("binary-options", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

expect_error <- function(expr, argument) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"),
            grepl(argument, conditionMessage(error), fixed = TRUE))
}

# Fit-time rounding uses R's ties-to-even rule; prediction values stay exact.
truth <- c(0, 1, 0, 1)
probability <- c(.501, .502, .509, .510)
newdata <- c(0, .499, .5, .501, .505, .51, .511, 1)
fold <- c("A", "A", "B", "B")
input <- cbind(A = newdata, B = rev(newdata))
for (precision in list(NULL, 0, 1L, 2, -1, 1e10, -1e10)) {
  rounded <- if (is.null(precision)) probability else round(probability, precision)
  stopifnot(
    identical(predict(venn_abers(truth, probability, precision), newdata),
              predict(venn_abers(truth, rounded), newdata)),
    identical(predict(cross_venn_abers(truth, probability, fold, precision), input),
              predict(cross_venn_abers(truth, rounded, fold), input))
  )
}
rounded_fit <- venn_abers(truth, probability, precision = 2)
stopifnot(
  !identical(predict(rounded_fit, newdata), predict(rounded_fit, round(newdata, 2))),
  !identical(predict(rounded_fit, newdata), predict(venn_abers(truth, probability), newdata))
)
for (bad in list(TRUE, FALSE, "2", NA_real_, NaN, Inf, -Inf, .5, -1.5,
                 numeric(), c(1, 2), list(2), 2i, matrix(2), factor(2))) {
  expect_error(venn_abers(truth, probability, precision = bad), "precision")
  expect_error(cross_venn_abers(truth, probability, fold, precision = bad), "precision")
}
for (bad in c(-.01, 1.01, NA_real_, NaN, Inf, -Inf)) {
  expect_error(venn_abers(0, bad, precision = 0), "probability")
  expect_error(cross_venn_abers(c(0, 1), c(bad, .5), c("A", "B"), 0), "probability")
}

# Binary-exact half-way inputs: .125 -> .12, .375 -> .38 in native R.
stopifnot(identical(
  predict(venn_abers(c(0, 1), c(.125, .375), 2), c(.12, .125, .375, .38)),
  predict(venn_abers(c(0, 1), c(.12, .38)), c(.12, .125, .375, .38))
))
# Decimal boundaries can depend on the R/platform version. Compare native R
# directly; the pinned NumPy reference rounds .545 to .55 and .575 to .57.
# These cases are kept separate from the Python parity tolerance.
boundary <- c(.545, .575)
stopifnot(identical(
  predict(venn_abers(c(0, 1), boundary, 2), boundary),
  predict(venn_abers(c(0, 1), round(boundary, 2)), boundary)
))

# Independent hand-worked Brier values from intervals (0, 1/2), (0, 1),
# and (1/2, 1), including the exact event tie.
fit <- venn_abers(c(0, 1), c(0, 1))
before <- serialize(fit, NULL)
log_prediction <- predict(fit, c(0, .5, 1))
brier_prediction <- predict(fit, c(0, .5, 1), loss = "brier")
stopifnot(
  identical(log_prediction, predict(fit, c(0, .5, 1), loss = "log")),
  identical(brier_prediction$probability, c(3 / 8, 1 / 2, 5 / 8)),
  identical(log_prediction[c("lower", "upper")], brier_prediction[c("lower", "upper")]),
  identical(serialize(fit, NULL), before)
)

# Unequal folds: (0, 1/2) and (2/3, 1) give Brier 3/8 and 13/18,
# hence 79/144 with equal weights, rather than the size-weighted 61/96.
cross_fit <- cross_venn_abers(c(0, 1, 1, 1), c(.5, .2, .5, .8), c("A", "B", "B", "B"))
cross_input <- cbind(A = .5, B = .5)
cross_before <- serialize(cross_fit, NULL)
stopifnot(
  abs(predict(cross_fit, cross_input, loss = "brier")$probability - 79 / 144) < 1e-15,
  abs(predict(cross_fit, cross_input)$probability - sqrt(3) / (sqrt(2) + sqrt(3))) < 1e-15,
  identical(serialize(cross_fit, NULL), cross_before)
)
for (bad in list(NULL, NA_character_, "Brier", "l", "", c("log", "brier"),
                 TRUE, 1, list("log"), factor("log"), matrix("log"))) {
  expect_error(predict(fit, .5, loss = bad), "loss")
  expect_error(predict(cross_fit, cross_input, loss = bad), "loss")
}

for (loss in c("log", "brier")) {
  stopifnot(identical(predict(fit, c(0, .5, 1), loss = loss, type = "class"),
                      c(FALSE, TRUE, TRUE)))
  # Balanced tied calibration gives .5; one changed truth lies just to either side.
  for (events in 999:1001) {
    y <- c(rep(1, events), rep(0, 2000 - events))
    tied <- venn_abers(y, rep(.5, 2000))
    tied_cross <- cross_venn_abers(rep(y, 2), rep(.5, 4000), rep(c("A", "B"), each = 2000))
    stopifnot(
      abs(predict(tied, .5, loss = loss)$probability - .5) < .001,
      identical(predict(tied, .5, loss = loss, type = "class"), events >= 1000),
      identical(predict(tied_cross, cross_input, loss = loss, type = "class"), events >= 1000)
    )
  }
  for (n in c(0L, 1L, 3L)) {
    rows <- seq_len(n)
    x <- c(0, .5, 1)[rows]
    matrix_input <- cbind(A = x, B = rev(x))
    rownames(matrix_input) <- rep("same", n)
    details <- predict(cross_fit, matrix_input, loss = loss, type = "details")
    assert_conditions(identical(names(details), c("prediction", "folds")), identical(details$prediction,
        predict(cross_fit, matrix_input, loss = loss)), identical(names(details$folds), c("row", "fold",
        "lower", "upper", "probability")), identical(details$folds$row, rep(rows, 2)), identical(details$folds$fold,
        rep(c("A", "B"), each = n)), identical(predict(cross_fit, matrix_input[, 2:1, drop = FALSE], loss = loss,
        type = "details"), details), identical(predict(cross_fit, matrix_input, loss = loss, type = "class"),
        details$prediction$probability >= 0.5), identical(predict(fit, x, loss = loss, type = "class"), predict(fit,
        x, loss = loss)$probability >= 0.5), identical(predict(fit, x, loss = loss, type = "response"), predict(fit,
        x, loss = loss)), identical(predict(cross_fit, matrix_input, loss = loss, type = "response"), details$prediction),
        all(vapply(details$folds[c("lower", "upper", "probability")], is.double, FALSE)))
    for (row in rows) {
      single <- predict(cross_fit, matrix_input[row, , drop = FALSE], loss = loss, type = "details")
      parts <- details$folds[details$folds$row == row, ]
      reconstructed <- if (loss == "log") {
        plogis(mean(log(parts$upper)) - mean(log1p(-parts$lower)))
      } else mean(parts$probability)
      stopifnot(abs(reconstructed - details$prediction$probability[row]) < 1e-15,
                identical(single$prediction$probability, details$prediction$probability[row]),
                identical(single$folds$probability, parts$probability))
    }
  }
  # Two non-event fold classes lose to the combined probability from all three.
  vote_fit <- cross_venn_abers(c(0, 0, rep(1, 9)), rep(.5, 11), c("A", "B", rep("C", 9)))
  vote_input <- cbind(A = .5, B = .5, C = .5)
  vote <- predict(vote_fit, vote_input, loss = loss, type = "details")
  stopifnot(sum(vote$folds$probability >= .5) == 1L,
            identical(predict(vote_fit, vote_input, loss = loss, type = "class"), TRUE))
}
log_details <- predict(cross_fit, cross_input, type = "details")
brier_details <- predict(cross_fit, cross_input, loss = "brier", type = "details")
stopifnot(
  identical(log_details$folds[1:4], brier_details$folds[1:4]),
  max(abs(brier_details$folds$probability - c(3 / 8, 13 / 18))) < 1e-15,
  identical(serialize(cross_fit, NULL), cross_before),
  identical(serialize(fit, NULL), before)
)
for (bad in list(NULL, NA_character_, "Class", "c", "", c("class", "response"),
                 TRUE, 1, list("class"), factor("class"), matrix("class"))) {
  expect_error(predict(fit, .5, type = bad), "type")
  expect_error(predict(cross_fit, cross_input, type = bad), "type")
}
expect_error(predict(fit, .5, type = "details"), "type")
expect_error(predict(fit, .5, precision = 2), "...")
expect_error(predict(cross_fit, cross_input, precision = 2), "...")

  testthat::expect_true(TRUE)
})

testthat::test_that("either binary class can be the event before calibration rounding", {
  truth <- c(0, 1, 0, 1)
  probability <- c(.2, .2, .4, .4)
  newdata <- c(0, .2, .3, .4, 1)
  # Unrounded endpoints from the pinned Python example in test-venn-abers.R.
  lower <- c(0, 1 / 3, 1 / 3, 2 / 5, 2 / 5)
  upper <- c(3 / 5, 3 / 5, 2 / 3, 2 / 3, 1)
  points <- list(log = c(3 / 8, 9 / 19, 1 / 2, 10 / 19, 5 / 8),
                 brier = c(21 / 50, 107 / 225, 1 / 2, 118 / 225, 29 / 50))
  for (precision in list(NULL, 0, 2, -1)) {
    for (event in c(1, 0)) {
      y <- truth == event
      p <- if (event == 1) probability else 1 - probability
      x <- if (event == 1) newdata else 1 - newdata
      expected_lower <- if (event == 1) lower else 1 - upper
      expected_upper <- if (event == 1) upper else 1 - lower
      expected_points <- if (event == 1) points else lapply(points, function(p) 1 - p)
      if (!is.null(precision) && precision <= 0) {
        # Rounding makes one balanced knot: zero, or one for event 0 at precision 0.
        at_knot <- if (precision == 0 && event == 0) x == 1 else x == 0
        below <- precision == 0 & event == 0 & !at_knot
        expected_lower <- ifelse(below, 0, 2 / 5)
        expected_upper <- ifelse(at_knot | below, 3 / 5, 1)
        expected_points <- list(log = ifelse(at_knot, 1 / 2, ifelse(below, 3 / 8, 5 / 8)),
                                brier = ifelse(at_knot, 1 / 2, ifelse(below, 21 / 50, 29 / 50)))
      }
      fit <- vennabers::venn_abers(y, p, precision)
      cross <- vennabers::cross_venn_abers(rep(y, 2), rep(p, 2),
                                         rep(c("A", "B"), each = 4), precision)
      for (loss in c("log", "brier")) {
        actual <- predict(fit, x, loss = loss)
        testthat::expect_lte(max(abs(actual$lower - expected_lower)), 1e-10)
        testthat::expect_lte(max(abs(actual$upper - expected_upper)), 1e-10)
        testthat::expect_lte(max(abs(actual$probability - expected_points[[loss]])), 1e-10)
        testthat::expect_identical(predict(fit, x, loss = loss, type = "class"),
                                   expected_points[[loss]] >= .5)
        details <- predict(cross, cbind(B = x, A = x), loss = loss, type = "details")
        testthat::expect_lte(max(abs(details$prediction$probability - expected_points[[loss]])), 1e-10)
        testthat::expect_lte(max(abs(details$folds$lower - rep(expected_lower, 2))), 1e-10)
        testthat::expect_lte(max(abs(details$folds$upper - rep(expected_upper, 2))), 1e-10)
        testthat::expect_identical(predict(cross, cbind(A = x, B = x), loss = loss, type = "class"),
                                   expected_points[[loss]] >= .5)
      }
    }
  }
})

testthat::test_that("many-fold public Cross prediction stays stable for either event", {
  ids <- paste0("fold_", seq_len(1024))
  for (event in c(TRUE, FALSE)) {
    fit <- vennabers::cross_venn_abers(rep(event, 9 * length(ids)), rep(.5, 9 * length(ids)),
                                      rep(ids, each = 9))
    input <- matrix(.5, 1, length(ids), dimnames = list(NULL, rev(ids)))
    # Identical folds have (9/10, 1) or (0, 1/10); products underflow.
    expected <- if (event) c(log = 10 / 11, brier = 181 / 200) else c(log = 1 / 11, brier = 19 / 200)
    for (loss in names(expected)) {
      testthat::expect_lte(abs(predict(fit, input, loss = loss)$probability - expected[[loss]]), 1e-10)
      testthat::expect_identical(predict(fit, input, loss = loss, type = "class"), event)
    }
  }
})

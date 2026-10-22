testthat::test_that("venn-abers", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

stopifnot(identical(
  sort(getNamespaceExports("vennabers")),
  sort(c("venn_abers", "venn_abers_regression", "cross_venn_abers_regression",
         "cross_venn_abers", "brier_score", "log_loss", "calibration_plot",
         "corp_plot", "corp_decomposition_plot", "regression_plot",
         "venn_abers_schema", "venn_abers_encode", "venn_abers_decode", "venn_abers_one_hot",
         "calibrated_workflow", "calibrate_workflow", "model_adapter", "local_calibration_splits", "venn_abers_details", "venn_abers_tuning", "fit_vennabers",
         "fit_vennabers_proxy",
         "predict_vennabers_prob", "predict_vennabers_class", "predict_vennabers_numeric"))
))

expect_error <- function(expr, argument) {
  error <- tryCatch({
    force(expr)
    NULL
  }, error = identity)
  stopifnot(
    inherits(error, "error"),
    grepl(paste0("`", argument, "`"), conditionMessage(error), fixed = TRUE)
  )
}

expect_no_warning <- function(expr) {
  warned <- FALSE
  value <- withCallingHandlers(expr, warning = function(warning) {
    warned <<- TRUE
    invokeRestart("muffleWarning")
  })
  stopifnot(!warned)
  value
}

fixture <- function(truth, probability, newdata, lower, upper, regularized) {
  list(
    truth = truth,
    probability = probability,
    newdata = newdata,
    expected = data.frame(lower, upper, probability = regularized)
  )
}

# Values captured from Python venn-abers 1.5.4 with precision=None.
fixtures <- list(
  fixture(
    c(0, 1, 0), c(.1, .2, .3), c(.05, .1, .15, .2, .25, .3, .4),
    c(0, 0, 0, 1 / 3, 1 / 3, 1 / 3, 1 / 3),
    c(1 / 2, 1 / 2, 2 / 3, 2 / 3, 2 / 3, 2 / 3, 1),
    c(1 / 3, 1 / 3, 2 / 5, 1 / 2, 1 / 2, 1 / 2, 3 / 5)
  ),
  fixture(
    c(0, 1, 0, 1), c(.2, .2, .4, .4), c(.1, .2, .3, .4, .5),
    c(0, 1 / 3, 1 / 3, 2 / 5, 2 / 5),
    c(3 / 5, 3 / 5, 2 / 3, 2 / 3, 1),
    c(3 / 8, 9 / 19, 1 / 2, 10 / 19, 5 / 8)
  ),
  fixture(
    c(0, 1, 0, 1), c(.501, .502, .509, .510), c(.501, .505, .510),
    c(0, 1 / 3, 1 / 2), c(1 / 2, 2 / 3, 1), c(1 / 3, 1 / 2, 2 / 3)
  ),
  fixture(
    c(0, 0, 0), c(.2, .5, .8), c(.1, .5, .9),
    c(0, 0, 0), c(1 / 4, 1 / 3, 1), c(1 / 5, 1 / 4, 1 / 2)
  ),
  fixture(
    c(1, 1, 1), c(.2, .5, .8), c(.1, .5, .9),
    c(0, 2 / 3, 3 / 4), c(1, 1, 1), c(1 / 2, 3 / 4, 4 / 5)
  ),
  fixture(
    0, .2, c(0, .2, 1),
    c(0, 0, 0), c(1 / 2, 1 / 2, 1), c(1 / 3, 1 / 3, 1 / 2)
  ),
  fixture(
    1, .2, c(0, .2, 1),
    c(0, 1 / 2, 1 / 2), c(1, 1, 1), c(1 / 2, 2 / 3, 2 / 3)
  ),
  fixture(
    c(0, 1, 0), c(.4, .4, .4), c(0, .4, 1),
    c(0, 1 / 4, 1 / 4), c(1 / 2, 1 / 2, 1), c(1 / 3, 2 / 5, 4 / 7)
  ),
  fixture(
    c(0, 1), c(0, 1), c(0, .5, 1),
    c(0, 0, 1 / 2), c(1 / 2, 1, 1), c(1 / 3, 1 / 2, 2 / 3)
  )
)

for (case in fixtures) {
  calibrator <- expect_no_warning(venn_abers(case$truth, case$probability))
  actual <- predict(calibrator, case$newdata)
  assert_conditions(
    inherits(calibrator, "venn_abers"),
    identical(names(actual), c("lower", "upper", "probability")),
    all(vapply(actual, is.double, FALSE)),
    max(abs(as.matrix(actual) - as.matrix(case$expected))) <= 1e-10,
    all(actual$lower <= actual$probability),
    all(actual$probability <= actual$upper),
    all(actual >= 0 & actual <= 1)
  )
}

calibrator <- venn_abers(c(0, 1, 0), c(.1, .2, .3))
ordered <- predict(calibrator, c(.4, .05, .2))
stopifnot(max(abs(as.matrix(ordered) - as.matrix(fixtures[[1L]]$expected[c(7, 1, 4), ]))) <= 1e-10)

logical_calibrator <- venn_abers(c(FALSE, TRUE), c(0L, 1L))
stopifnot(nrow(predict(logical_calibrator, c(1L, 0L))) == 2L)

empty <- predict(calibrator, numeric(0))
stopifnot(
  identical(names(empty), c("lower", "upper", "probability")),
  nrow(empty) == 0L,
  all(vapply(empty, is.double, FALSE))
)

invalid_truths <- list("0", factor(0), list(0), 0 + 0i, c(0, NA), c(0, Inf), c(0, .5), matrix(c(0, 1)))
invalid_probabilities <- list(
  TRUE, "0.5", factor(.5), list(.5), .5 + 0i,
  NA_real_, NaN, Inf, -.1, 1.1, matrix(.5)
)

numeric_truth <- c(0, 1, 1)
logical_truth <- as.logical(numeric_truth)
regularized_probability <- c(.1, .8, .6)
names(regularized_probability) <- letters[1:3]
for (truth in list(numeric_truth, logical_truth)) {
  brier <- brier_score(truth, regularized_probability)
  loss <- log_loss(truth, regularized_probability)
  assert_conditions(
    isTRUE(all.equal(brier, .07)),
    is.double(brier),
    length(brier) == 1L,
    is.null(names(brier)),
    isTRUE(all.equal(loss, 0.2797765635793423)),
    is.double(loss),
    length(loss) == 1L,
    is.null(names(loss))
  )
}

stopifnot(
  identical(brier_score(c(0, 1), c(0, 1)), 0),
  identical(brier_score(c(0, 1), c(1, 0)), 1),
  identical(log_loss(c(0, 1), c(0, 1)), 0),
  identical(log_loss(c(0, 1), c(1, 0)), Inf)
)

for (consumer in list(venn_abers, brier_score, log_loss)) {
  expect_error(consumer(numeric(0), numeric(0)), "truth")
  expect_error(consumer(c(0, 1), .5), "probability")
  for (bad_truth in invalid_truths) {
    expect_error(consumer(bad_truth, rep(.5, length(bad_truth))), "truth")
  }
  for (bad_probability in invalid_probabilities) {
    expect_error(consumer(0, bad_probability), "probability")
  }
}

for (bad_newdata in invalid_probabilities) {
  expect_error(predict(calibrator, bad_newdata), "newdata")
}
expect_error(predict(calibrator, .5, unused = TRUE), "...")
expect_error(predict(calibrator, .5, NULL), "...")

printed <- capture.output(returned <- print(calibrator))
stopifnot(
  identical(printed, c(
    "<venn_abers>",
    "Calibration observations: 3",
    "Event observations: 1",
    "Unique probability knots: 3"
  )),
  identical(returned, calibrator)
)

plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file)
calibration <- withVisible(calibration_plot(
  c(0, 1, 1, 0, 1),
  c(0, .25, .5, .75, 1),
  bins = 4
))
assert_conditions(
  !calibration$visible,
  identical(names(calibration$value), c("mean_probability", "event_rate", "count")),
  identical(calibration$value$count, c(1L, 1L, 1L, 2L)),
  isTRUE(all.equal(calibration$value$mean_probability, c(0, .25, .5, .875))),
  isTRUE(all.equal(calibration$value$event_rate, c(0, 1, 1, .5)))
)

decimal_boundary <- calibration_plot(c(0, 1), c(.579, .58), bins = 50)
stopifnot(
  identical(decimal_boundary$mean_probability, c(.579, .58)),
  identical(decimal_boundary$event_rate, c(0, 1)),
  identical(decimal_boundary$count, c(1L, 1L))
)

for (bad_bins in list(numeric(0), c(1, 2), NA_real_, NaN, Inf, -Inf, 0, -1, 1.5, TRUE, "2", list(2), matrix(2))) {
  expect_error(calibration_plot(0, .5, bins = bad_bins), "bins")
}

# 0.9.0 also accepted a count with extra attributes or beyond the integer range.
for (good_bins in list(structure(2, label = "two"), 3e9)) {
  expect_identical(
    calibration_plot(c(0, 1), c(.25, .75), bins = good_bins),
    data.frame(mean_probability = c(.25, .75), event_rate = c(0, 1), count = c(1L, 1L))
  )
}

degenerate_plots <- list(
  expect_no_warning(calibration_plot(1, .5)),
  expect_no_warning(calibration_plot(c(0, 0), c(.1, .9), bins = 1)),
  expect_no_warning(calibration_plot(c(1, 1, 1), rep(.5, 3), bins = 1e9)),
  expect_no_warning(calibration_plot(
    c(FALSE, TRUE), c(0, 1), col = "red", pch = 19, main = "Override", x = 1:2, y = 0
  ))
)
grDevices::dev.off()
unlink(plot_file)
stopifnot(
  identical(degenerate_plots[[1L]], data.frame(mean_probability = .5, event_rate = 1, count = 1L)),
  identical(degenerate_plots[[2L]], data.frame(mean_probability = .5, event_rate = 0, count = 2L)),
  identical(degenerate_plots[[3L]], data.frame(mean_probability = .5, event_rate = 1, count = 3L)),
  nrow(degenerate_plots[[4L]]) == 2L
)

for (bad_truth in invalid_truths) {
  expect_error(calibration_plot(bad_truth, rep(.5, length(bad_truth))), "truth")
}
for (bad_probability in invalid_probabilities) {
  expect_error(calibration_plot(0, bad_probability), "probability")
}
expect_error(calibration_plot(numeric(0), numeric(0)), "truth")
expect_error(calibration_plot(c(0, 1), .5), "probability")

  testthat::expect_true(TRUE)
})

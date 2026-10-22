testthat::test_that("randomized-parity", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

calibration_path <- file.path("extdata", "randomized-calibration.csv")
prediction_path <- file.path("extdata", "randomized-predictions.csv")
stopifnot(nzchar(calibration_path), nzchar(prediction_path))

calibration <- read.csv(calibration_path, comment.char = "#")
prediction <- read.csv(prediction_path, comment.char = "#")
case_ids <- unique(calibration$case_id)
sizes <- table(calibration$case_id)
truth_class_counts <- vapply(
  split(calibration$truth, calibration$case_id),
  function(x) length(unique(x)),
  0L
)

assert_conditions(
  length(case_ids) == 100L,
  min(sizes) == 1L,
  max(sizes) == 200L,
  identical(as.integer(table(prediction$case_id)), rep(25L, 100L)),
  any(duplicated(calibration[c("case_id", "probability")])),
  any(vapply(split(calibration, calibration$case_id), function(x) !anyDuplicated(x$probability), FALSE)),
  1L %in% truth_class_counts,
  2L %in% truth_class_counts
)

maximum_error <- 0
for (case_id in case_ids) {
  calibration_set <- calibration[calibration$case_id == case_id, ]
  expected <- prediction[prediction$case_id == case_id, ]
  actual <- predict(
    venn_abers(calibration_set$truth, calibration_set$probability),
    expected$probability
  )
  expected_output <- expected[c("lower", "upper", "regularized_probability")]
  names(expected_output)[3L] <- "probability"

  assert_conditions(
    identical(names(actual), c("lower", "upper", "probability")),
    nrow(actual) == 25L,
    all(vapply(actual, is.double, FALSE)),
    all(is.finite(as.matrix(actual))),
    all(actual >= 0 & actual <= 1),
    all(actual$lower <= actual$probability),
    all(actual$probability <= actual$upper)
  )
  maximum_error <- max(maximum_error, abs(as.matrix(actual) - as.matrix(expected_output)))

  stopifnot(
    0 %in% expected$probability,
    1 %in% expected$probability,
    any(expected$probability < min(calibration_set$probability)),
    any(expected$probability > max(calibration_set$probability))
  )
}

stopifnot(maximum_error <= 1e-10)

  testthat::expect_true(TRUE)
})

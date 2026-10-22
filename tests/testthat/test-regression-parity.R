testthat::test_that("regression-parity", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

read_fixture <- function(name, text_columns) {
  path <- file.path("extdata", name)
  result <- read.csv(path, colClasses = "character")
  for (column in setdiff(names(result), text_columns)) result[[column]] <- as.numeric(result[[column]])
  result
}
calibration <- read_fixture("regression-calibration.csv", "case_id")
expected <- read_fixture("regression-predictions.csv", "case_id")
options <- read_fixture("regression-options.csv", c("case_id", "parameter"))
close <- function(actual, wanted, scale) {
  stopifnot(all(is.finite(actual)), all(is.finite(wanted)),
            all(abs(actual - wanted) <= 1e-10 * scale + 1e-10 * abs(wanted)))
}

# Independent max-min characterization of isotonic regression: append the
# hypothetical observation, pool equal scores, then consider every contiguous
# block containing its score. This is deliberately not the production hull scan.
isotonic_at <- function(truth, prediction, new, outer) {
  truth <- c(truth, outer)
  prediction <- c(prediction, new)
  knots <- sort(unique(prediction))
  index <- match(new, knots)
  max(vapply(seq_len(index), function(left) {
    min(vapply(seq.int(index, length(knots)), function(right) {
      mean(truth[prediction >= knots[left] & prediction <= knots[right]])
    }, 0))
  }, 0))
}

stopifnot(nrow(options) == 116L, nrow(expected) == 1436L,
          sum(expected$direct_valid) == 1099L,
          sum(!expected$direct_valid) == 337L)
for (i in seq_len(nrow(options))) {
  option <- options[i, ]
  cal <- calibration[calibration$case_id == option$case_id, ]
  wanted <- expected[expected$case_id == option$case_id, ]
  fit <- do.call(venn_abers_regression,
                 c(list(truth = cal$truth, prediction = cal$prediction),
                   setNames(list(option$requested), option$parameter)))
  assert_conditions(fit$m == option$m, fit$epsilon == option$epsilon,
            fit$a == option$a, fit$b == option$b,
            fit$observations == nrow(cal),
            identical(cal$row, as.numeric(seq_len(nrow(cal)))),
            identical(wanted$row, as.numeric(seq_len(nrow(wanted)))))
  actual <- predict(fit, wanted$input)
  stopifnot(identical(names(actual), c("prediction", "lower", "upper")),
            nrow(actual) == nrow(wanted), all(vapply(actual, typeof, "") == "double"))
  for (column in names(actual)) {
    close(actual[[column]], wanted[[column]], option$scale)
    valid <- wanted$direct_valid == 1
    close(actual[[column]][valid], wanted[[paste0("direct_", column)]][valid], option$scale)
  }
  rounding <- 1e-10 * option$scale + 1e-10 * abs(actual$prediction)
  stopifnot(all(actual$lower <= actual$prediction + rounding),
            all(actual$prediction <= actual$upper + rounding))
  # Joint permutation preserves fit; arbitrary new-row permutations preserve output.
  permutation <- rev(seq_len(nrow(cal)))
  permuted <- do.call(venn_abers_regression,
                      c(list(truth = cal$truth[permutation], prediction = cal$prediction[permutation]),
                        setNames(list(option$requested), option$parameter)))
  close(as.matrix(predict(permuted, rev(wanted$input))),
        as.matrix(actual[nrow(actual):1L, ]), option$scale)
  if (endsWith(option$case_id, "_ordinary")) {
    lower <- vapply(wanted$input, function(new) isotonic_at(cal$lower_truth, cal$prediction, new, option$a), 0)
    upper <- vapply(wanted$input, function(new) isotonic_at(cal$upper_truth, cal$prediction, new, option$b), 0)
    close(actual$lower, lower, option$scale)
    close(actual$upper, upper, option$scale)
  }
}

  testthat::expect_true(TRUE)
})

testthat::test_that("cross-regression-parity", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

read_fixture <- function(name, text_columns) {
  result <- read.csv(file.path("extdata", name),
                     colClasses = "character")
  for (column in setdiff(names(result), text_columns)) result[[column]] <- as.numeric(result[[column]])
  result
}
calibration <- read_fixture("cross-regression-calibration.csv", c("case_id", "fold"))
components <- read_fixture("cross-regression-predictions.csv", c("case_id", "fold"))
options <- read_fixture("cross-regression-options.csv", c("case_id", "fold", "parameter"))
aggregate <- read_fixture("cross-regression-aggregate.csv", "case_id")
columns <- c("prediction", "lower", "upper")
close <- function(actual, wanted, scale) {
  stopifnot(identical(dim(actual), dim(wanted)),
            all(is.finite(actual)), all(is.finite(wanted)),
            all(abs(actual - wanted) <= 1e-10 * scale + 1e-10 * abs(wanted)))
}
schema <- function(actual, rows) {
  # The IVAR contract permits endpoint interpolation roundoff. Check ordering
  # separately from reference tolerances, at a few machine epsilons of the bounds.
  rounding <- 8 * .Machine$double.eps * pmax(abs(actual$lower), abs(actual$upper))
  assert_conditions(is.data.frame(actual), nrow(actual) == rows, identical(names(actual), columns), all(vapply(actual,
      typeof, "") == "double"), all(is.finite(as.matrix(actual))), all(actual$lower <= actual$prediction +
      rounding), all(actual$prediction <= actual$upper + rounding))
}
assert_conditions(length(unique(options$case_id)) == 80L, nrow(options) == 220L, nrow(components) ==
    1540L, nrow(aggregate) == 560L, sum(components$direct_valid) == 1243L, sum(!components$direct_valid) ==
    297L)
for (case in unique(options$case_id)) {
  option <- options[options$case_id == case, ]
  cal <- calibration[calibration$case_id == case, ]
  wanted <- aggregate[aggregate$case_id == case, ]
  expected_folds <- components[components$case_id == case, ]
  arguments <- list(truth = cal$truth, prediction = cal$prediction, fold = cal$fold)
  if (option$parameter[1] != "default") arguments[[option$parameter[1]]] <- option$requested[1]
  fit <- do.call(cross_venn_abers_regression, arguments)
  newdata <- vapply(option$fold, function(id) expected_folds$input[expected_folds$fold == id], numeric(7))
  # Deliberately reverse columns independently of the first-observed fold order.
  newdata <- newdata[, rev(seq_len(ncol(newdata))), drop = FALSE]
  actual <- predict(fit, newdata, type = "details")
  stopifnot(identical(fit$folds, unique(cal$fold)),
            identical(names(actual$folds), fit$folds),
            identical(cal$row, as.numeric(seq_len(nrow(cal)))),
            identical(wanted$row, as.numeric(seq_len(nrow(wanted)))))
  schema(actual$prediction, nrow(wanted))
  close(as.matrix(actual$prediction), as.matrix(wanted[columns]), option$scale[1])
  for (id in fit$folds) {
    opt <- option[option$fold == id, ]
    expected <- expected_folds[expected_folds$fold == id, ]
    component <- fit$calibrators[[id]]
    assert_conditions(identical(expected$row, as.numeric(seq_len(nrow(expected)))), component$observations ==
        opt$observations, component$m == opt$m, component$epsilon == opt$epsilon, component$a == opt$a, component$b ==
        opt$b)
    requested <- if (opt$parameter == "default") list() else setNames(list(opt$requested), opt$parameter)
    stopifnot(identical(component$requested, requested))
    schema(actual$folds[[id]], nrow(expected))
    close(as.matrix(actual$folds[[id]]), as.matrix(expected[columns]), opt$scale)
    valid <- expected$direct_valid == 1
    close(as.matrix(actual$folds[[id]][valid, ]),
          as.matrix(expected[valid, paste0("direct_", columns)]), opt$scale)
    # Fit the stored IVAR inputs independently through its public interface.
    rows <- cal$fold == id
    ivar <- do.call(venn_abers_regression, c(list(truth = cal$truth[rows], prediction = cal$prediction[rows]), requested))
    stopifnot(identical(actual$folds[[id]], predict(ivar, newdata[, id])))
  }
  # Independently reconstruct all three columns with ordinary arithmetic on
  # these predeclared <=1e120 fixtures (extreme overflow tests are separate).
  for (column in columns) {
    reconstructed <- Reduce(`+`, lapply(actual$folds, `[[`, column)) / length(fit$folds)
    close(actual$prediction[[column]], reconstructed, option$scale[1])
  }
  reversed <- lapply(arguments, rev)
  if (option$parameter[1] != "default") reversed[[option$parameter[1]]] <- arguments[[option$parameter[1]]]
  reordered <- do.call(cross_venn_abers_regression, reversed)
  close(as.matrix(predict(reordered, newdata[7:1, , drop = FALSE])),
        as.matrix(actual$prediction[7:1, ]), option$scale[1])
}

  testthat::expect_true(TRUE)
})

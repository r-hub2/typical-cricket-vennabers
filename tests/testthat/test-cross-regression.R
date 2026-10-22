testthat::test_that("cross-regression", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

close_regression <- function(actual, expected, scale = 1) {
  stopifnot(identical(dim(actual), dim(expected)),
            all(is.finite(actual)), all(is.finite(expected)),
            all(abs(actual - expected) <= 1e-10 * scale + 1e-10 * abs(expected)))
}

# A tied five-row fold has L=-1, U=8, point=29/7 at score zero.
# A constant three-row fold contributes 10 to each column. Equal fold means
# differ from size weighting, pooling, midpoint and merging averaged bounds.
truth <- c(-10, 0, 2, 6, 20, 10, 10, 10)
prediction <- rep(0, 8)
fold <- rep(c("tied", "constant"), c(5, 3))
fit <- cross_venn_abers_regression(truth, prediction, fold)
newdata <- cbind(constant = c(1, 0, -1), tied = c(-1, 0, 1))
details <- predict(fit, newdata, type = "details")
close_regression(as.matrix(details$prediction),
                 cbind(prediction = c(9, 99 / 14, 4.5),
                       lower = c(0, 4.5, 4.5), upper = c(9, 9, 15)))
stopifnot(identical(names(details$folds), c("tied", "constant")),
          identical(predict(fit, newdata), details$prediction))

expect_error <- function(expr, argument) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(argument, conditionMessage(error), fixed = TRUE))
}
check_schema <- function(result, rows) {
  assert_conditions(is.data.frame(result), nrow(result) == rows, identical(names(result), c("prediction",
      "lower", "upper")), all(vapply(result, typeof, "") == "double"), all(is.finite(as.matrix(result))),
      all(result$lower <= result$prediction), all(result$prediction <= result$upper))
}
for (column in names(details$prediction)) {
  close_regression(details$prediction[[column]],
                   (details$folds$tied[[column]] + details$folds$constant[[column]]) / 2)
}
stopifnot(abs(details$prediction$prediction[2] - (5 * 29 / 7 + 3 * 10) / 8) > .1,
          abs(details$prediction$prediction[2] - 6.75) > .1)
pooled <- predict(venn_abers_regression(truth, prediction), 0)
stopifnot(abs(details$prediction$prediction[2] - pooled$prediction) > .1)
# Averaged outer bounds are a=0, b=15; merging L=4.5, U=9 gives 99/14
# here by coincidence. Use two nonconstant folds below to distinguish ordering.
other <- cross_venn_abers_regression(c(truth[1:5], 0, 0, 0, 4, 10), rep(0, 10),
                                    rep(c("A", "B"), each = 5))
other_result <- predict(other, cbind(A = 0, B = 0))
close_regression(other_result$prediction, (29 / 7 + 116 / 33) / 2)
merged_after <- ((1 / 6 + 5) / 6 + (15 - 6) * 6) / ((1 / 6 + 5) + (15 - 6))
stopifnot(abs(other_result$prediction - merged_after) > .01)

# Metadata retains the requested access and effective per-fold moderation.
stopifnot(identical(fit$folds, c("tied", "constant")), fit$observations == 8,
          identical(names(fit$calibrators), fit$folds))
for (id in fit$folds) {
  component <- fit$calibrators[[id]]
  assert_conditions(inherits(component, "venn_abers_regression"), identical(component$requested, list()),
      component$m == 1, component$observations == sum(fold == id), component$epsilon == 2/(sum(fold ==
          id) + 1))
  stopifnot(identical(details$folds[[id]], predict(component, newdata[, id])))
}
moderation_folds <- rep(c("five", "eight"), c(5, 8))
for (epsilon in c(1 / 3, 2 / 3 - 1e-12, 2 / 3, 2 / 3 + 1e-12)) {
  moderated <- cross_venn_abers_regression(1:13, 13:1, moderation_folds, epsilon = epsilon)
  for (id in moderated$folds) {
    component <- moderated$calibrators[[id]]
    k <- sum(moderation_folds == id)
    stopifnot(identical(component$requested, list(epsilon = epsilon)),
              component$m == floor(epsilon * (k + 1) / 2),
              component$epsilon == 2 * component$m / (k + 1))
  }
  if (epsilon < 2 / 3) stopifnot(moderated$calibrators$five$epsilon != moderated$calibrators$eight$epsilon)
}
explicit <- cross_venn_abers_regression(1:13, 13:1, moderation_folds, m = 2)
stopifnot(all(vapply(explicit$calibrators, function(x) identical(x$requested, list(m = 2)), TRUE)))
# At tied score zero, epsilon=.5 gives m=1 for five rows, m=2 for eight.
# The latter's moderated sums are 24 (lower) and 32 (upper); adding a=-2
# or b=10 gives L=22/9, U=14/3 and weighted point=362/99.
unequal <- cross_venn_abers_regression(c(-10, 0, 2, 6, 20, -8, -2, 0, 1, 3, 8, 10, 30),
                                     rep(0, 13), moderation_folds, epsilon = .5)
unequal_details <- predict(unequal, cbind(eight = 0, five = 0), type = "details")
close_regression(as.matrix(unequal_details$folds$five), matrix(c(29 / 7, -1, 8), 1))
close_regression(as.matrix(unequal_details$folds$eight), matrix(c(362 / 99, 22 / 9, 14 / 3), 1))
close_regression(as.matrix(unequal_details$prediction), matrix(c(5405 / 1386, 13 / 18, 19 / 3), 1))
expect_error(cross_venn_abers_regression(truth, prediction, fold, m = 2), "constant")
expect_error(cross_venn_abers_regression(1:7, 1:7, rep(c("valid", "too small"), c(5, 2))), "too small")
expect_error(cross_venn_abers_regression(truth, prediction, fold, epsilon = .4), "constant")
expect_error(cross_venn_abers_regression(1:9, 1:9, rep(c("five", "four"), c(5, 4)), epsilon = .9), "four")
expect_error(cross_venn_abers_regression(truth, prediction, fold, m = 1, epsilon = .5), "only one")
expect_error(cross_venn_abers_regression(truth, prediction, fold, m = 1, epsilon = NULL), "only one")
for (bad in list(NULL, TRUE, "1", NA_real_, NaN, Inf, -Inf, 1i, factor(1),
                 list(1), matrix(1), structure(1, class = "Date"), numeric(), c(1, 2), 0, -1, 1.5, 3)) {
  expect_error(cross_venn_abers_regression(truth, prediction, fold, m = bad), "m")
}
for (bad in list(NULL, TRUE, "0.5", NA_real_, NaN, Inf, -Inf, .5i,
                 factor(.5), list(.5), matrix(.5), numeric(), c(.4, .5), 0, -1, 1, 1.1)) {
  expect_error(cross_venn_abers_regression(truth, prediction, fold, epsilon = bad), "epsilon")
}
for (bad in list(NULL, TRUE, "1", factor(1:8), list(1:8), matrix(1:8),
                 data.frame(x = 1:8), 1i, structure(1:8, class = "Date"),
                 c(NA, 1:7), c(NaN, 1:7), c(Inf, 1:7), c(-Inf, 1:7))) {
  expect_error(cross_venn_abers_regression(bad, prediction, fold), "truth")
  expect_error(cross_venn_abers_regression(truth, bad, fold), "prediction")
}
expect_error(cross_venn_abers_regression(numeric(), numeric(), character()), "nonzero")
expect_error(cross_venn_abers_regression(truth, prediction[-1], fold), "same")
for (bad in list(NULL, character(), rep("A", 8), fold[-1], c(fold, "A"),
                 c("", fold[-1]), c(NA, fold[-1]), factor(fold), as.list(fold),
                 matrix(fold), rep(TRUE, 8), c(0, 1:7), c(-1, 1:7), c(.5, 1:7),
                 c(Inf, 1:7), c(-Inf, 1:7), c(NaN, 1:7), c(NA, 1:7), rep(1i, 8))) {
  expect_error(cross_venn_abers_regression(truth, prediction, bad), "fold")
}
for (bad in list(NULL, 1:2, data.frame(newdata), list(newdata), matrix(TRUE, 3, 2),
                 matrix("1", 3, 2), matrix(1i, 3, 2), newdata + NA_real_,
                 newdata + NaN, newdata * Inf, newdata * -Inf)) {
  expect_error(predict(fit, bad), "newdata")
}
for (columns in list(NULL, c("tied", "tied"), c("tied", "other"),
                     c("tied", ""), c("tied", NA_character_))) {
  bad <- newdata
  colnames(bad) <- columns
  expect_error(predict(fit, bad), "column names")
}
expect_error(predict(fit, newdata[, 1, drop = FALSE]), "column names")
expect_error(predict(fit, cbind(newdata, extra = 0)), "column names")
for (bad in list(NULL, NA_character_, "class", "raw", "det", c("response", "details"), TRUE)) {
  expect_error(predict(fit, newdata, type = bad), "type")
}
expect_error(cross_venn_abers_regression(truth, prediction, fold, precision = 2), "unsupported")
expect_error(predict(fit, newdata, loss = "brier"), "unsupported")
expect_error(predict(fit, newdata, m = 2), "unsupported")
expect_error(predict(fit, newdata, "details"), "unsupported")

# IDs are literal; all row association is positional, even with conflicting names.
for (ids in list(c(20, 3), c("fold / B", "fold:A"), c(" A", "a"))) {
  named <- cross_venn_abers_regression(setNames(truth, letters[1:8]),
                                      setNames(prediction, rev(letters[1:8])),
                                      setNames(rep(ids, c(5, 3)), rev(LETTERS[1:8])))
  input <- newdata[, c("tied", "constant")]
  colnames(input) <- as.character(ids)
  rownames(input) <- c("third", "first", "second")
  stopifnot(identical(predict(named, input), details$prediction))
}
# Numeric 0/1 outcomes retain the regression bounds and weighted merger.
# At tied score zero, three rows give (.5, 0, .5); five give (17/30, .5, 2/3).
binary_fit <- cross_venn_abers_regression(c(0, 1, 0, 0, 0, 1, 1, 1), rep(0, 8),
                                        rep(c("three", "five"), c(3, 5)))
binary_details <- predict(binary_fit, cbind(five = 0, three = 0), type = "details")
close_regression(as.matrix(binary_details$folds$three), matrix(c(.5, 0, .5), 1))
close_regression(as.matrix(binary_details$folds$five), matrix(c(17 / 30, .5, 2 / 3), 1))
close_regression(as.matrix(binary_details$prediction), matrix(c(8 / 15, .25, 7 / 12), 1))
check_schema(binary_details$prediction, 1)
permutation <- c(8, 1, 6, 3, 5, 2, 7, 4)
permuted <- cross_venn_abers_regression(truth[permutation], prediction[permutation], fold[permutation])
close_regression(as.matrix(predict(permuted, newdata[3:1, 2:1])),
                 as.matrix(details$prediction[3:1, ]))
before <- serialize(fit, NULL)
for (rows in list(integer(), 2L, c(3L, 1L, 2L))) {
  input <- newdata[rows, , drop = FALSE]
  output <- predict(fit, input, type = "details")
  check_schema(output$prediction, length(rows))
  stopifnot(identical(names(output$folds), fit$folds))
  for (component in output$folds) check_schema(component, length(rows))
  close_regression(as.matrix(output$prediction), as.matrix(details$prediction[rows, ]))
}
stopifnot(identical(serialize(fit, NULL), before))

# Exact fitted endpoints and their zero-weight midpoint propagate through Cross.
endpoint_fit <- cross_venn_abers_regression(rep(c(-1e16, -3, -3), 2),
                                           rep(1:3, 2), rep(c("A", "B"), each = 3))
endpoint_result <- predict(endpoint_fit, cbind(B = 0, A = 0), type = "details")
endpoint_expected <- data.frame(prediction = -1e16 / 2 - 3 / 2, lower = -1e16, upper = -3)
stopifnot(identical(endpoint_result$folds$A, endpoint_expected),
          identical(endpoint_result$folds$B, endpoint_expected),
          identical(endpoint_result$prediction, endpoint_expected),
          identical(predict(endpoint_fit, cbind(B = 0, A = 0)), endpoint_expected))

# Nonconstant zero weights, small positive weights and changed target units.
for (scale in c(1e-200, 1e-12, 1, 1e120)) {
  zero <- cross_venn_abers_regression(c(0, 0, 2, 2, 0, 0, 0, 4, 10) * scale,
                                     c(1:4, rep(0, 5)), rep(c("zero", "tiny"), c(4, 5)))
  output <- predict(zero, cbind(zero = 2.5, tiny = 0), type = "details")
  close_regression(as.matrix(output$prediction),
                   matrix(c((1 + 116 / 33) / 2, 2 / 3, 3), 1) * scale, scale)
  for (component in c(list(output$prediction), output$folds)) check_schema(component, 1)
}
for (constants in list(c(0, 0), c(-7, -7), c(5e-324, 5e-324),
                       c(1e308, 1e308), c(-1e308, -1e308),
                       rep(.Machine$double.xmax, 3),
                       c(-1.7e308, 1.7e308, 1.7e308))) {
  ids <- as.character(seq_along(constants))
  constant_fit <- cross_venn_abers_regression(rep(constants, each = 3),
                                             rep(1:3, length(constants)), rep(ids, each = 3))
  input <- matrix(rep(c(-100, 2, 100), length(ids)), ncol = length(ids), dimnames = list(NULL, ids))
  output <- predict(constant_fit, input, type = "details")
  check_schema(output$prediction, 3)
  wanted <- if (constants[1] != constants[2]) 1.7e308 / 3 else constants[1]
  close_regression(as.matrix(output$prediction), matrix(wanted, 3, 3), max(abs(constants), 5e-324))
  if (length(constants) == 2L) stopifnot(all(as.matrix(output$prediction) == constants[1]))
}

# Cancellation must retain the small residual, regardless of fold order.
for (magnitude in c(1e16, 1e308)) {
  for (residual in c(1, 1e-16, 1e-300, 3 * 5e-324)) {
    constants <- c(magnitude, residual, -magnitude)
    for (permutation in list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
                             c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))) {
      cancellation <- cross_venn_abers_regression(rep(constants[permutation], each = 3),
                                                 rep(1:3, 3), rep(letters[1:3], each = 3))
      output <- predict(cancellation, cbind(c = 0, a = 0, b = 0))
      stopifnot(all(is.finite(as.matrix(output))),
                all(abs(as.matrix(output) / (residual / 3) - 1) < 1e-14))
    }
  }
}

# Overflow in input order must not trigger scaling before cancellation.
fold_mean <- vennabers:::.regression_mean
for (residual in c(1e-16, 1e-300, 5 * 5e-324)) {
  values <- c(1e308, 1e308, residual, -1e308, -1e308)
  for (direction in c(-1, 1)) {
    for (x in list(values, rev(values))) {
      stopifnot(abs(fold_mean(direction * x) / (direction * residual / 5) - 1) < 1e-14)
    }
  }
}
assert_conditions(identical(fold_mean(c(4.94065645841247e-324, 4.94065645841247e-324)), 4.94065645841247e-324),
    identical(fold_mean(c(1e+16, 1, 1)), (1e+16 + 2)/3), identical(fold_mean(rep(.Machine$double.xmax,
        3)), .Machine$double.xmax), identical(fold_mean(rep(-.Machine$double.xmax, 3)), -.Machine$double.xmax),
    identical(fold_mean(c(4.94065645841247e-324, 0, 0)), 0))

local({
  input <- tempfile(fileext = ".rds")
  output <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  on.exit(unlink(c(input, output, script)))
  saveRDS(list(fit = fit, newdata = newdata), input)
  writeLines(c(paste0(".libPaths(", deparse(.libPaths()[1L]), ")"),
               "library(vennabers)", "args <- commandArgs(TRUE)",
               "x <- readRDS(args[1L])", "saveRDS(predict(x$fit, x$newdata, type = 'details'), args[2L])"), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script), shQuote(input), shQuote(output)))
  stopifnot(status == 0, identical(readRDS(output), details))
})

  testthat::expect_true(TRUE)
})

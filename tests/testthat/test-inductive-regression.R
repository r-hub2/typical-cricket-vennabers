testthat::test_that("inductive-regression", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

close_regression <- function(actual, expected, scale = 1) {
  stopifnot(identical(dim(actual), dim(expected)),
            all(is.finite(actual)), all(is.finite(expected)),
            all(abs(actual - expected) <= 1e-10 * scale + 1e-10 * abs(expected)))
}

# All predictions tied: each isotonic endpoint is the mean of the five
# moderated truths and the added outer truth. Lower: (-10,0,2,6,6;-10),
# upper: (0,0,2,6,20;20), giving -1 and 8. Merger: (-9+96)/21.
fit <- venn_abers_regression(c(-10, 0, 2, 6, 20), rep(0, 5))
close_regression(as.matrix(predict(fit, c(-1, 0, 1))),
                 cbind(prediction = c(8, 29 / 7, -1),
                       lower = c(-10, -1, -1), upper = c(8, 8, 20)))

expect_error <- function(expr, argument) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(argument, conditionMessage(error), fixed = TRUE))
}

assert_conditions(identical(fit$requested, list()), fit$m == 1, fit$epsilon == 1 / 3,
          fit$a == -10, fit$b == 20, fit$observations == 5)
explicit <- venn_abers_regression(1:5, 1:5, m = 2)
stopifnot(identical(explicit$requested, list(m = 2)), explicit$m == 2,
          explicit$epsilon == 2 / 3)
# Maximum m=2 replaces extremes rather than discarding rows: lower truths
# (0,0,2,2,2;0) average to 1, upper (2,2,2,6,6;6) to 4.
# At the tie, weights 1 and 2 give point 3; the midpoint would be 2.5.
maximum <- venn_abers_regression(c(-10, 0, 2, 6, 20), rep(0, 5), m = 2)
close_regression(as.matrix(predict(maximum, c(-1, 0, 1))),
                 cbind(prediction = c(4, 3, 1),
                       lower = c(0, 1, 1), upper = c(4, 4, 6)))
for (epsilon in c(1 / 3, 2 / 3 - 1e-12, 2 / 3, 2 / 3 + 1e-12)) {
  fit_e <- venn_abers_regression(1:5, 1:5, epsilon = epsilon)
  stopifnot(identical(fit_e$requested, list(epsilon = epsilon)),
            fit_e$m == if (epsilon < 2 / 3) 1 else 2,
            identical(predict(fit_e, 0:6), predict(venn_abers_regression(1:5, 1:5, m = fit_e$m), 0:6)))
}
for (bad in list(NULL, TRUE, FALSE, "1", NA_real_, NaN, Inf, -Inf, 1i,
                 factor(1), list(1), matrix(1), structure(1, class = "Date"),
                 numeric(), c(1, 2), 0, -1, 1.5, 3)) {
  expect_error(venn_abers_regression(1:5, 1:5, m = bad), "m")
}
for (bad in list(NULL, TRUE, "0.5", NA_real_, NaN, Inf, -Inf, .5i,
                 factor(.5), list(.5), matrix(.5), numeric(), c(.4, .5),
                 0, -1, 1, 1.1, 1 / 3 - 1e-12)) {
  expect_error(venn_abers_regression(1:5, 1:5, epsilon = bad),
               if (identical(bad, 1 / 3 - 1e-12)) "m" else "epsilon")
}
expect_error(venn_abers_regression(1:4, 1:4, epsilon = .9), "m")
expect_error(venn_abers_regression(1:5, 1:5, m = 1, epsilon = .5), "only one")
expect_error(venn_abers_regression(1:5, 1:5, m = 1, epsilon = NULL), "only one")
expect_error(venn_abers_regression(1:2, 1:2), "three")
expect_error(venn_abers_regression(numeric(), numeric()), "nonzero")
expect_error(venn_abers_regression(1:3, 1:4), "same")
expect_error(venn_abers_regression(1:3, 1:3, precision = 2), "unsupported")
expect_error(predict(fit, 0, type = "class"), "unsupported")
for (bad in list(NULL, TRUE, "1", factor(1:3), list(1:3), matrix(1:3),
                 data.frame(x = 1:3), 1i, structure(1:3, class = "Date"),
                 c(NA, 1, 2), c(NaN, 1, 2), c(Inf, 1, 2), c(-Inf, 1, 2))) {
  expect_error(venn_abers_regression(bad, 1:3), "truth")
  expect_error(venn_abers_regression(1:3, bad), "prediction")
  expect_error(predict(fit, bad), "newdata")
}

# Zero weights with nonconstant truth: the gap separates two flat isotonic blocks.
zero <- venn_abers_regression(c(0, 0, 2, 2), 1:4)
stopifnot(identical(predict(zero, 2.5), data.frame(prediction = 1, lower = 0, upper = 2)))
# Interpolation must preserve endpoints for scalar and vector bounds, including
# opposite signs whose difference overflows, and stay accurate near either end.
interpolate <- vennabers:::.regression_interpolate
endpoints <- list(c(-1e16, -3), c(3, 1e16), c(-1e308, 1e308),
                  c(-1e308, 0), c(0, 1e308), c(5e-324, 5e-324))
for (bounds in endpoints) {
  stopifnot(identical(interpolate(bounds[1], bounds[2], c(0, .5, 1))[c(1, 3)], bounds))
}
lower <- vapply(endpoints, `[`, 0, 1)
upper <- vapply(endpoints, `[`, 0, 2)
stopifnot(identical(interpolate(lower, upper, rep(0, length(lower))), lower),
          identical(interpolate(lower, upper, rep(1, length(upper))), upper),
          identical(interpolate(c(1, -2^54), c(2^54, -1), c(2^-54, 1 - 2^-53)), c(2, -3)))
# The constant upper-moderated truth must remain exactly b, so both exterior
# merger weights vanish and the point is the midpoint, rather than upper.
endpoint_fit <- venn_abers_regression(c(-1e16, -3, -3), 1:3)
stopifnot(all(endpoint_fit$upper == -3),
          identical(predict(endpoint_fit, 0),
                    data.frame(prediction = -1e16 / 2 - 3 / 2, lower = -1e16, upper = -3)))
# Tiny positive weights must retain the weighted answer, not the midpoint.
for (scale in c(1e-200, 1e-12, 1, 1e150)) {
  tiny <- venn_abers_regression(c(0, 0, 0, 4, 10) * scale, rep(0, 5))
  # L=4/3, U=4 in base units: wL=4/3, wU=6 -> prediction=116/33.
  close_regression(as.matrix(predict(tiny, 0)), matrix(c(116 / 33, 4 / 3, 4), 1) * scale, scale)
}
for (constant in c(0, -7, 1e-300, 1e308, -1e308)) {
  constant_fit <- venn_abers_regression(rep(constant, 3), c(-9, 0, 6))
  result <- predict(constant_fit, c(-100, -9, 0, 6, 100))
  stopifnot(all(as.matrix(result) == constant), all(is.finite(as.matrix(result))))
}
# Finite extremes would overflow unscaled cumulative sums and merger products.
huge <- venn_abers_regression(c(-1e308, -1e308, 1e308, 1e308), 1:4)
stopifnot(identical(predict(huge, 2.5), data.frame(prediction = 0, lower = -1e308, upper = 1e308)))
for (truth in list(c(-1e308, 0, 1e308), c(1e308, 1.2e308, 1.7e308))) {
  result <- predict(venn_abers_regression(truth, 1:3), 0:4)
  stopifnot(all(is.finite(as.matrix(result))),
            all(result$lower <= result$prediction), all(result$prediction <= result$upper))
}

# Double schema for zero/one/many rows; names never realign observations.
named <- venn_abers_regression(setNames(c(-10, 0, 2, 6, 20), letters[1:5]),
                               setNames(rep(0, 5), rev(letters[1:5])))
stopifnot(identical(predict(fit, c(-1, 0, 1)), predict(named, c(z = -1, x = 0, y = 1))),
          inherits(venn_abers_regression(c(0, 1, 0), 1:3), "venn_abers_regression"))
# Numeric 0/1 outcomes still use the regression merger: tied lower truths
# (0,0,0;0) give 0, upper (0,1,0;1) give .5, hence point .5.
binary_truth <- venn_abers_regression(c(0, 1, 0), rep(0, 3))
stopifnot(identical(predict(binary_truth, c(-1, 0, 1)),
                    data.frame(prediction = rep(.5, 3), lower = rep(0, 3),
                               upper = c(.5, .5, 1))))
before <- serialize(fit, NULL)
for (input in list(numeric(), 0, c(1, -1, 0))) {
  result <- predict(fit, input)
  stopifnot(is.data.frame(result), nrow(result) == length(input),
            identical(names(result), c("prediction", "lower", "upper")),
            all(vapply(result, typeof, "") == "double"))
}
stopifnot(identical(serialize(fit, NULL), before))

# A fresh R process must need only the installed package and serialized fit.
local({
  input <- tempfile(fileext = ".rds")
  output <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  on.exit(unlink(c(input, output, script)))
  saveRDS(list(fit = fit, newdata = c(1, 0, -1)), input)
  writeLines(c(paste0(".libPaths(", deparse(.libPaths()[1L]), ")"),
               "library(vennabers)", "args <- commandArgs(TRUE)",
               "x <- readRDS(args[1L])", "saveRDS(predict(x$fit, x$newdata), args[2L])"), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script), shQuote(input), shQuote(output)))
  stopifnot(status == 0, identical(readRDS(output), predict(fit, c(1, 0, -1))))
})

  testthat::expect_true(TRUE)
})

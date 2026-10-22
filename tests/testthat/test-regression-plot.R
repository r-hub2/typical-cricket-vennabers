testthat::test_that("regression plots of held-out pairs", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
error <- function(expr, pattern) {
  result <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result), fixed = TRUE))
}
path <- tempfile(fileext = ".pdf")
grDevices::pdf(path)
on.exit({ grDevices::dev.off(); unlink(path) })
layout <- graphics::par("mfrow")
columns <- c("prediction", "truth", "residual")

# One row per pair in input order, residuals observed minus predicted,
# returned invisibly in the same form whatever `which` holds.
truth <- c(3, -1, 2.5, 10, 0)
prediction <- c(2, 0, 2.5, 7, -4)
drawn <- withVisible(regression_plot(truth, prediction))
table <- drawn$value
assert_conditions(!drawn$visible, is.data.frame(table), identical(names(table), columns),
                  identical(table$prediction, prediction), identical(table$truth, truth),
                  identical(table$residual, c(1, -1, 0, 3, 4)),
                  identical(rownames(table), as.character(1:5)),
                  identical(graphics::par("mfrow"), layout))
for (which in list("residual", c("observed", "residual"), c("residual", "observed"))) {
  assert_conditions(identical(regression_plot(truth, prediction, which = which), table),
                    identical(graphics::par("mfrow"), layout))
}

# Integer and named inputs give plain double columns; a single pair is valid.
named <- regression_plot(c(a = 1L, b = 4L), c(a = 2L, b = 2L))
assert_conditions(is.double(named$prediction), is.double(named$truth), is.double(named$residual),
                  identical(named$residual, c(-1, 2)), identical(rownames(named), c("1", "2")))
assert_conditions(identical(regression_plot(5, 5, which = "residual")$residual, 0))

# Graphical options override the defaults: a single panel leaves its
# coordinates, so `xlim` shows through `usr` on an unpadded axis.
regression_plot(truth, prediction, xlim = c(-5, 8), xaxs = "i", main = "Own title")
assert_conditions(identical(graphics::par("usr")[1:2], c(-5, 8)))

# Both diagrams share one grid, which a failing draw restores as well.
error(regression_plot(truth, prediction, which = c("observed", "residual"), ylim = "bad"), "ylim")
assert_conditions(identical(graphics::par("mfrow"), layout))

# `which` names distinct diagrams.
for (bad in list(character(), NA_character_, c("observed", NA), "map", c("observed", "observed"),
                 1, list("observed"), factor("observed"))) {
  error(regression_plot(truth, prediction, which = bad),
        "`which` must be one or more distinct values of: observed, residual.")
}

# Inputs are finite numeric vectors of equal nonzero length, without `na.rm`.
error(regression_plot(c(1, NA), c(1, 2)), "`truth` must be a numeric vector of finite values.")
error(regression_plot(c(1, Inf), c(1, 2)), "`truth` must be a numeric vector of finite values.")
error(regression_plot(c("1", "2"), c(1, 2)), "`truth` must be a numeric vector of finite values.")
error(regression_plot(factor(1:2), c(1, 2)), "`truth` must be a numeric vector of finite values.")
error(regression_plot(c(1, 2), c(1, NaN)), "`prediction` must be a numeric vector of finite values.")
error(regression_plot(c(1, 2), matrix(1:2)), "`prediction` must be a numeric vector of finite values.")
error(regression_plot(numeric(), numeric()), "`truth` and `prediction` must have the same nonzero length.")
error(regression_plot(c(1, 2), 1), "`truth` and `prediction` must have the same nonzero length.")
assert_conditions(identical(graphics::par("mfrow"), layout))

  testthat::expect_true(TRUE)
})

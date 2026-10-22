testthat::test_that("fit diagnoses draw binary and Cross calibration maps", {
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
values <- c("lower", "upper", "probability")

# Binary: one row per knot plus both range ends, valued exactly by predict().
truth <- c(0, 1, 0, 1, 1, 0, 1, 1)
probability <- c(.1, .2, .3, .4, .4, .6, .7, .9)
fit <- venn_abers(truth, probability)
drawn <- withVisible(plot(fit))
table <- drawn$value
assert_conditions(!drawn$visible, is.data.frame(table),
                  identical(names(table), c("score", values, "width")),
                  identical(table$score, c(0, fit$knots, 1)),
                  identical(table[values], predict(fit, table$score)),
                  identical(table$width, table$upper - table$lower),
                  identical(graphics::par("mfrow"), layout))

# The table does not depend on `which`; both diagrams draw map before width
# in their own grid and leave the caller's layout in place.
# The width panel, drawn last, keeps its own axis: the widest gap between
# knots (upper of the right knot minus lower of the left one) plus 5 %.
assert_conditions(identical(graphics::par("usr")[3:4], c(0, 1)))
gap <- max(table$upper[-1L] - table$lower[-nrow(table)])
for (which in list("width", c("map", "width"), c("width", "map"))) {
  assert_conditions(identical(plot(fit, which = which), table),
                    identical(graphics::par("mfrow"), layout),
                    isTRUE(all.equal(graphics::par("usr")[3:4], c(0, 1.05 * gap))))
}
brier <- plot(fit, "map", loss = "brier")
assert_conditions(identical(brier[values], predict(fit, brier$score, loss = "brier")))

# Knots at the range ends are not repeated.
edges <- plot(venn_abers(c(0, 1, 1), c(0, .5, 1)))
assert_conditions(identical(edges$score, c(0, .5, 1)))

# Graphical options override the defaults: a single panel leaves its
# coordinates, so `xlim` shows through `usr`.
plot(fit, xlim = c(.2, .8), main = "Own title")
assert_conditions(identical(graphics::par("usr")[1:2], c(.2, .8)))

# Cross: stored fold IDs lead, folds keep their stored order, and each fold
# map equals its own calibrator.
fold <- c(2, 1, 2, 1, 2, 1, 2, 1)
cross <- cross_venn_abers(truth, probability, fold = fold)
for (loss in c("log", "brier")) {
  maps <- plot(cross, c("width", "map"), loss = loss)
  assert_conditions(identical(names(maps), c("fold", "score", values, "width")),
                    identical(unique(maps$fold), cross$folds), identical(cross$folds, c("2", "1")),
                    identical(rownames(maps), as.character(seq_len(nrow(maps)))),
                    identical(maps$width, maps$upper - maps$lower),
                    identical(graphics::par("mfrow"), layout))
  for (id in cross$folds) {
    calibrator <- cross$calibrators[[id]]
    own <- maps[maps$fold == id, ]
    assert_conditions(identical(own$score, c(0, calibrator$knots, 1)),
                      identical(as.list(own[values]), as.list(predict(calibrator, own$score, loss = loss))))
  }
}

# A failing draw restores the grid as well.
error(plot(fit, c("map", "width"), ylim = "bad"), "ylim")
assert_conditions(identical(graphics::par("mfrow"), layout))

# Choices are validated like other options.
for (bad in list(character(), c("map", "map"), "fit", NA_character_, c("map", NA), 1, list("map"))) {
  error(plot(fit, which = bad), "`which` must be one or more distinct values of: map, width.")
  error(plot(cross, which = bad), "`which` must be one or more distinct values of: map, width.")
}
for (bad in list("hinge", c("log", "brier"), NA_character_, 1)) {
  error(plot(fit, loss = bad), "`loss` must be one of: log, brier.")
  error(plot(cross, loss = bad), "`loss` must be one of: log, brier.")
}

# Fits without a stored calibration map stop clearly.
stored <- "`x` has no stored calibration map (`knots`, `lower`, `upper`); refit it with vennabers 1.0.0."
for (field in c("knots", "lower", "upper")) {
  old <- fit
  old[[field]] <- NULL
  error(plot(old), stored)
  old_cross <- cross
  old_cross$calibrators[[2L]][[field]] <- NULL
  error(plot(old_cross), stored)
}
old_cross <- cross
old_cross$calibrators <- NULL
error(plot(old_cross), stored)

testthat::expect_true(TRUE)
})

testthat::test_that("fit diagnoses draw one panel per multiclass subproblem", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
error <- function(expr, pattern) {
  result <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result), fixed = TRUE))
}
path <- tempfile(fileext = ".pdf")
grDevices::pdf(path)
layout <- graphics::par("mfrow")
# Every new panel records the grid it is drawn in.
grids <- list()
setHook("plot.new", function() grids[[length(grids) + 1L]] <<- graphics::par("mfrow"))
on.exit({ setHook("plot.new", NULL, "replace"); grDevices::dev.off(); unlink(path) })
panels <- function(expr) {
  grids <<- list()
  force(expr)
  grids
}
values <- c("score", "lower", "upper", "probability", "width")

classes <- c("a", "b", "c")
truth <- factor(rep(classes, 4), levels = classes)
probability <- matrix(c(.6, .3, .1, .2, .5, .3, .1, .2, .7, .5, .1, .4,
                        .3, .4, .3, .2, .2, .6, .4, .4, .2, .1, .6, .3,
                        .3, .3, .4, .5, .3, .2, .2, .7, .1, .2, .2, .6), 12, 3, byrow = TRUE,
                      dimnames = list(NULL, classes))
fold <- rep(c(2, 1), 6)

for (method in c("ovr", "ovo")) {
  labels <- if (method == "ovr") "class" else c("class_1", "class_2")
  fits <- list(inductive = venn_abers(truth, probability, method = method),
               cross = cross_venn_abers(truth, probability, fold = fold, method = method))
  for (kind in names(fits)) {
    fit <- fits[[kind]]
    problems <- if (method == "ovr") as.list(classes) else fit$pairs
    cross <- kind == "cross"
    for (loss in c("log", "brier")) {
      drawn <- panels(maps <- plot(fit, loss = loss))
      assert_conditions(identical(names(maps), c(if (cross) "fold", labels, values)),
                        identical(rownames(maps), as.character(seq_len(nrow(maps)))),
                        identical(maps$width, maps$upper - maps$lower),
                        length(drawn) == 3L, identical(drawn[[1L]], c(2L, 2L)),
                        identical(graphics::par("mfrow"), layout))
      # Each subproblem holds the map of its own calibrator, Cross folds in
      # the stored order of the fit.
      for (i in seq_along(problems)) {
        problem <- problems[[i]]
        own <- if (method == "ovr") maps[maps$class == problem, ] else
          maps[maps$class_1 == problem[1L] & maps$class_2 == problem[2L], ]
        calibrator <- fit$calibrators[[i]]
        if (cross) {
          assert_conditions(identical(unique(own$fold), fit$folds))
          for (id in fit$folds) {
            map <- own[own$fold == id, values]
            assert_conditions(identical(map$score, c(0, calibrator$calibrators[[id]]$knots, 1)),
                              identical(as.list(map[c("lower", "upper", "probability")]),
                                        as.list(predict(calibrator$calibrators[[id]], map$score, loss = loss))))
          }
        } else {
          assert_conditions(identical(own$score, c(0, calibrator$knots, 1)),
                            identical(as.list(own[c("lower", "upper", "probability")]),
                                      as.list(predict(calibrator, own$score, loss = loss))))
        }
      }
    }
    # Both diagrams double the grid rows: each row of maps is followed by the
    # widths of the same subproblems, the short last row padded by one empty
    # panel. The table does not depend on `which`.
    drawn <- panels(both <- plot(fit, c("width", "map")))
    assert_conditions(identical(both, plot(fit)), length(drawn) == 7L,
                      all(vapply(drawn, identical, TRUE, c(4L, 2L))),
                      identical(graphics::par("mfrow"), layout))
  }
}

# Folds of a subproblem follow the fit even where they occur in another
# order among its rows.
cross <- cross_venn_abers(truth, probability, fold = c(1, 1, 3, 2, 2, 1, 3, 3, 2, 1, 2, 3), method = "ovo")
assert_conditions(!identical(cross$calibrators[[1L]]$folds, cross$folds))
maps <- plot(cross)
for (i in seq_along(cross$pairs)) {
  own <- maps[maps$class_1 == cross$pairs[[i]][1L] & maps$class_2 == cross$pairs[[i]][2L], ]
  assert_conditions(identical(unique(own$fold), cross$folds))
}

# Two classes in OvO are a single pair, drawn in the current layout.
pair <- factor(rep(c("a", "b"), 3), levels = c("a", "b"))
fit <- venn_abers(pair, cbind(a = c(.7, .4, .6, .3, .5, .2), b = c(.3, .6, .4, .7, .5, .8)), method = "ovo")
drawn <- panels(maps <- plot(fit))
assert_conditions(identical(unique(maps[c("class_1", "class_2")]), data.frame(class_1 = "a", class_2 = "b")),
                  length(drawn) == 1L, identical(drawn[[1L]], layout))

# Graphical options reach every panel, and failures restore the grid.
fit <- venn_abers(truth, probability, method = "ovr")
plot(fit, xlim = c(.2, .8), main = "Own title")
assert_conditions(identical(graphics::par("usr")[1:2], c(.2, .8)))
error(plot(fit, ylim = "bad"), "ylim")
assert_conditions(identical(graphics::par("mfrow"), layout))

# Fits without stored subproblem maps stop clearly.
stored <- "`x` has no stored calibration map (`knots`, `lower`, `upper`); refit it with vennabers 1.0.0."
old <- fit
old$calibrators[[2L]]$knots <- NULL
error(plot(old), stored)
old$calibrators <- NULL
error(plot(old), stored)
old <- cross_venn_abers(truth, probability, fold = fold, method = "ovo")
old$calibrators[[3L]]$calibrators[[1L]]$upper <- NULL
error(plot(old, "width"), stored)
assert_conditions(identical(graphics::par("mfrow"), layout))

testthat::expect_true(TRUE)
})

testthat::test_that("fit diagnoses draw regression maps over the knot range", {
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
values <- c("lower", "upper", "prediction")

# Inductive: one row per knot plus the ends of the knot range widened by 5 %
# on each side, valued exactly by predict().
truth <- c(-4, 1, 0, 3, 2, 6, 5, 9, 7, 12)
prediction <- c(-2, 0, 0, 1, 2, 3, 3, 5, 6, 8)
fit <- venn_abers_regression(truth, prediction)
drawn <- withVisible(plot(fit))
table <- drawn$value
assert_conditions(!drawn$visible, is.data.frame(table),
                  identical(names(table), c("score", values, "width")),
                  identical(table$score, c(-2.5, fit$knots, 8.5)),
                  identical(table[values], predict(fit, table$score)[values]),
                  identical(table$width, table$upper - table$lower),
                  identical(graphics::par("mfrow"), layout),
                  identical(graphics::par("usr")[1:2], c(-2.5, 8.5)))

# The table does not depend on `which`, and both diagrams draw in their own
# grid. The width panel, drawn last, reaches 5 % above the widest gap.
gap <- max(table$upper[-1L] - table$lower[-nrow(table)])
for (which in list("width", c("map", "width"), c("width", "map"))) {
  assert_conditions(identical(plot(fit, which = which), table),
                    identical(graphics::par("mfrow"), layout),
                    isTRUE(all.equal(graphics::par("usr")[3:4], c(0, 1.05 * gap))))
}

# A given `xlim` replaces the range: its ends bound the table, which keeps
# only the knots within, and the drawn coordinates.
narrow <- plot(fit, xlim = c(.5, 4), main = "Own title")
assert_conditions(identical(narrow$score, c(.5, 1, 2, 3, 4)),
                  identical(narrow[values], predict(fit, narrow$score)[values]),
                  identical(graphics::par("usr")[1:2], c(.5, 4)))
wide <- plot(fit, xlim = c(20, -10))
assert_conditions(identical(wide$score, c(-10, fit$knots, 20)),
                  identical(graphics::par("usr")[1:2], c(20, -10)))

# A single knot widens by 5 % of its absolute value, or 0.05 at zero.
tied <- venn_abers_regression(c(-10, 0, 2, 6, 20), rep(4, 5))
assert_conditions(identical(plot(tied)$score, c(3.8, 4, 4.2)))
tied <- venn_abers_regression(c(-10, 0, 2, 6, 20), rep(0, 5))
assert_conditions(identical(plot(tied)$score, c(-.05, 0, .05)))
# Knots near the largest double keep finite range ends.
huge <- venn_abers_regression(1:3, c(-.Machine$double.xmax, 0, .Machine$double.xmax))
assert_conditions(identical(range(plot(huge)$score), c(-.Machine$double.xmax, .Machine$double.xmax)))

# Cross: stored fold IDs lead, all folds share one range, and each fold map
# equals its own calibrator.
fold <- c(2, 1, 2, 1, 2, 1, 2, 1, 2, 1)
cross <- cross_venn_abers_regression(truth, prediction, fold = fold)
maps <- plot(cross, c("width", "map"))
ends <- c(-2.5, 8.5)
assert_conditions(identical(names(maps), c("fold", "score", values, "width")),
                  identical(unique(maps$fold), cross$folds), identical(cross$folds, c("2", "1")),
                  identical(rownames(maps), as.character(seq_len(nrow(maps)))),
                  identical(maps$width, maps$upper - maps$lower),
                  identical(graphics::par("mfrow"), layout))
for (id in cross$folds) {
  calibrator <- cross$calibrators[[id]]
  own <- maps[maps$fold == id, ]
  assert_conditions(identical(own$score, c(ends[1L], calibrator$knots, ends[2L])),
                    identical(as.list(own[values]), as.list(predict(calibrator, own$score)[values])))
}
narrow <- plot(cross, xlim = c(0, 5))
for (id in cross$folds) assert_conditions(identical(range(narrow$score[narrow$fold == id]), c(0, 5)))

# `loss` belongs to classification.
for (bad in list("log", "brier", NULL)) {
  error(plot(fit, loss = bad), "Regression does not accept `loss`.")
  error(plot(cross, loss = bad), "Regression does not accept `loss`.")
}

# A failing draw restores the grid as well, and `which` is validated.
error(plot(fit, c("map", "width"), ylim = "bad"), "ylim")
error(plot(cross, c("map", "width"), ylim = "bad"), "ylim")
assert_conditions(identical(graphics::par("mfrow"), layout))
for (bad in list(character(), c("map", "map"), "fit", NA_character_)) {
  error(plot(fit, which = bad), "`which` must be one or more distinct values of: map, width.")
  error(plot(cross, which = bad), "`which` must be one or more distinct values of: map, width.")
}

# Fits without a stored calibration map stop clearly.
stored <- "`x` has no stored calibration map (`knots`, `lower`, `upper`); refit it with vennabers 1.0.0."
for (field in c("knots", "lower", "upper")) {
  old <- fit
  old[[field]] <- NULL
  error(plot(old), stored)
  old_cross <- cross
  old_cross$calibrators[[2L]][[field]] <- NULL
  error(plot(old_cross), stored)
}
old_cross <- cross
old_cross$calibrators <- NULL
error(plot(old_cross), stored)

testthat::expect_true(TRUE)
})

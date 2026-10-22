testthat::test_that("inductive-multiclass", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

expect_error <- function(expr, message) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(message, conditionMessage(error), fixed = TRUE))
}
close <- function(x, y) stopifnot(identical(dim(x), dim(y)), all(abs(x - y) <= 1e-10))
fixture <- function(name) read.csv(file.path("extdata", name))
calibration <- fixture("ovo-calibration.csv")
inputs <- fixture("ovo-input.csv")
case_ids <- unique(calibration$case_id)
stopifnot(length(case_ids) == 33L, sum(startsWith(case_ids, "random")) == 30L)
# fixed_3_00 at precision 2 also protects binary64 row normalization: extended
# row sums on x86 otherwise move unrounded predictions across calibration knots.

# Independent coupling from the public table, using odds rather than the
# implementation's reciprocal accumulator. Class identity never relies on names
# joined with a separator (class names may themselves contain separators).
reconstruct <- function(details) {
  result <- details$prediction
  for (row in seq_len(nrow(result))) {
    parts <- details$pairs[details$pairs$row == row, ]
    for (class in colnames(result)) {
      q <- c(1 - parts$probability[parts$non_event == class],
             parts$probability[parts$event == class])
      result[row, class] <- 1 / (1 + sum((1 - q) / q))
    }
  }
  result / rowSums(result)
}

for (method in c("ovo", "ovr")) {
  reference <- fixture(paste0(method, "-output.csv"))
  binary_reference <- fixture(if (method == "ovo") "ovo-pairs.csv" else "ovr-classes.csv")
  maximum_error <- 0
  for (case_id in case_ids) {
    cal <- calibration[calibration$case_id == case_id, ]
    new <- inputs[inputs$case_id == case_id, ]
    k <- cal$classes[1L]
    classes <- c("zebra", "a b", "class/3", "4")[seq_len(k)]
    decode <- function(data) {
      result <- matrix(as.numeric(as.matrix(data[paste0("p", seq_len(k))])), nrow(data), k)
      colnames(result) <- classes
      result
    }
    p <- decode(cal)
    x <- decode(new)
    truth <- factor(classes[cal$truth], levels = classes)
    for (setting in c("none", "0", "2", "-1")) {
      precision <- if (setting == "none") NULL else as.numeric(setting)
      fit <- venn_abers(truth, p, precision = precision, method = method)
      before <- serialize(fit, NULL)
      stopifnot(identical(fit, venn_abers(truth, p[, rev(classes), drop = FALSE], precision, method = method)))
      endpoints <- NULL
      for (loss in c("log", "brier")) {
        expected <- reference[reference$case_id == case_id & reference$precision == setting & reference$loss == loss, ]
        binary_expected <- binary_reference[binary_reference$case_id == case_id & binary_reference$precision == setting &
                                             binary_reference$loss == loss, ]
        actual <- predict(fit, x, loss = loss, type = "details")
        parts <- actual[[2L]]
        assert_conditions(identical(names(actual), c("prediction", if (method == "ovo") "pairs" else "classes")),
                  identical(names(parts), c("row", if (method == "ovo") "non_event", "event", "lower", "upper", "probability")),
                  identical(colnames(actual$prediction), classes),
                  identical(actual$prediction, predict(fit, x, loss = loss)),
                  identical(actual, predict(fit, x[, rev(classes), drop = FALSE], loss = loss, type = "details")),
                  identical(parts$row, binary_expected$row),
                  identical(parts$event, classes[binary_expected$event]),
                  identical(predict(fit, x, loss = loss, type = "class"),
                            factor(classes[max.col(actual$prediction, ties.method = "first")], levels = classes)),
                  all(is.finite(actual$prediction)), all(actual$prediction >= 0 & actual$prediction <= 1),
                  all(abs(rowSums(actual$prediction) - 1) < 1e-14),
                  all(parts$lower <= parts$upper))
        if (method == "ovo") stopifnot(identical(parts$non_event, classes[binary_expected$non_event]))
        if (!is.null(endpoints)) stopifnot(identical(parts[c("lower", "upper")], endpoints))
        endpoints <- parts[c("lower", "upper")]
        maximum_error <- max(maximum_error, abs(actual$prediction - decode(expected)),
                             abs(as.matrix(parts[c("lower", "upper", "probability")]) -
                                   as.matrix(binary_expected[c("lower", "upper", "probability")])))
        if (maximum_error > 1e-10) stop(method, " parity failed: ", case_id, "/", setting, "/", loss)
        # Pinned Python probabilities choose different classes for both losses.
        # Protect against dispatching OvR class prediction back through OvO.
        if (case_id == "random_3_02" && setting == "none") {
          stopifnot(as.character(predict(fit, x, loss = loss, type = "class"))[10L] ==
                      if (method == "ovo") "zebra" else "a b")
        }
        if (method == "ovo") {
          close(reconstruct(actual), actual$prediction)
        } else {
          # Rebuild by observation and literal event identity, independent of
          # diagnostic row order and matrix column position.
          recovered <- vapply(seq_len(nrow(x)), function(row) {
            values <- setNames(parts$probability[parts$row == row], parts$event[parts$row == row])
            values[classes] / sum(values)
          }, numeric(k))
          close(t(recovered), actual$prediction)
        }
      }
      stopifnot(identical(serialize(fit, NULL), before))
    }
  }
  cat("Inductive", method, "maximum absolute reference error:", maximum_error, "\n")
}

classes <- c("z", "a", "a/z")
truth <- factor(rep(classes, 2), levels = classes)
p <- rbind(c(1, 0, 0), c(1, 0, 0), c(0, .4, .6), c(.1, .5, .4), c(.2, .3, .5), c(.6, .1, .3))
colnames(p) <- classes
x <- rbind(c(1, 0, 0), c(.28, .32, .4), c(.2, .2, .6), c(.1, .445, .455))
colnames(x) <- classes

# R's defined zero-pair and rounding behavior is checked against separately
# prepared binary fits. It is intentionally outside invalid Python zero-pair cases.
for (precision in list(NULL, 0, 2, -1)) {
  fit <- venn_abers(truth, p, precision)
  for (loss in c("log", "brier")) {
    actual <- predict(fit, x, loss = loss, type = "details")
    for (pair in combn(classes, 2, simplify = FALSE)) {
      rows <- truth == pair[1] | truth == pair[2]
      a <- p[rows, pair[1]]
      b <- p[rows, pair[2]]
      score <- ifelse(a + b == 0, .5, b / (a + b))
      if (!is.null(precision)) score <- round(score, precision)
      a <- x[, pair[1]]
      b <- x[, pair[2]]
      new_score <- ifelse(a + b == 0, .5, b / (a + b))
      expected <- predict(venn_abers(truth[rows] == pair[2], score), new_score, loss = loss)
      parts <- actual$pairs[actual$pairs$non_event == pair[1] & actual$pairs$event == pair[2], ]
      close(as.matrix(parts[4:6]), as.matrix(expected))
    }
    close(actual$prediction, reconstruct(actual))
    for (n in 0:1) {
      small <- x[seq_len(n), , drop = FALSE]
      details <- predict(fit, small, loss = loss, type = "details")
      assert_conditions(identical(dim(details$prediction), c(n, 3L)), nrow(details$pairs) == n * 3L,
                identical(names(details$pairs), names(actual$pairs)),
                identical(levels(predict(fit, small, type = "class", loss = loss)), classes),
                length(predict(fit, small, type = "class", loss = loss)) == n,
                all(vapply(details$pairs[4:6], is.double, FALSE)))
      close(details$prediction, actual$prediction[seq_len(n), , drop = FALSE])
    }
  }
}
# This case distinguishes rounding after pair normalization from rounding the
# original class matrix; new scores distinguish rounded/unrounded prediction.
rounded <- venn_abers(truth, p, precision = 0)
zero_pair <- predict(rounded, x, type = "details")$pairs
az <- zero_pair[zero_pair$non_event == "a" & zero_pair$event == "a/z", ]
expected <- predict(venn_abers(c(0, 1, 0, 1), c(0, 1, 1, 1)), c(.5, .4 / .72, .75, .455 / .9))
close(as.matrix(az[4:6]), as.matrix(expected))
stopifnot(!identical(expected, predict(venn_abers(c(0, 1, 0, 1), c(0, 1, 1, 1)),
                                      round(c(.5, .4 / .72, .75, .455 / .9), 0))))

# Two-class factor ties use the first level; binary ties still choose the event.
two_truth <- factor(c("z", "a"), levels = c("z", "a"))
two_p <- cbind(z = c(1, 0), a = c(0, 1))
two_x <- cbind(a = .5, z = .5)
for (loss in c("log", "brier")) {
  two_fit <- venn_abers(two_truth, two_p)
  close(predict(two_fit, two_x, loss = loss), matrix(.5, 1, 2))
  stopifnot(identical(predict(two_fit, two_x, loss = loss, type = "class"), factor("z", levels = c("z", "a"))),
            identical(predict(venn_abers(c(0, 1), c(0, 1)), .5, loss = loss, type = "class"), TRUE))
  all_ties <- venn_abers(factor(classes, levels = classes), matrix(1 / 3, 3, 3, dimnames = list(NULL, classes)))
  stopifnot(identical(predict(all_ties, x[1, , drop = FALSE] * 0 + 1 / 3, type = "class", loss = loss),
                      factor("z", levels = classes)))
}

for (method in c("ovo", "ovr")) {
  fit <- venn_abers(truth, p, method = method)
  for (bad in list(TRUE, "2", NA_real_, Inf, .5, c(1, 2), matrix(2))) {
    expect_error(venn_abers(truth, p, precision = bad, method = method), "precision")
  }
  for (bad in list("unknown", "OVO", "o", NA_character_, c("ovo", "ovo"), TRUE, NULL)) {
    expect_error(venn_abers(truth, p, method = bad), "method")
  }
  expect_error(venn_abers(c(0, 1), c(.2, .8), method = "ovo"), "method")
  for (bad in list("bad", NA_character_, c("log", "brier"), NULL, 1)) {
    expect_error(predict(fit, x, loss = bad), "loss")
    expect_error(predict(fit, x, type = bad), "type")
  }
  expect_error(predict(fit, x, precision = 2), "...")
  for (bad in list(factor(character(), levels = classes), factor(c("z", NA), levels = classes),
                   factor(c("", "a")), factor("a"), structure(truth, dim = c(2L, 3L)))) {
    expect_error(venn_abers(bad, p, method = method), "truth")
  }
  expect_error(venn_abers(truth[1:2], p[1:2, ], method = method), "Every class")
  expect_error(venn_abers(truth, p[-1, ], method = method), "one row")
  bad_matrices <- list(unname(p), p[, -1], cbind(p, extra = 0), as.data.frame(p),
                       p[, 1], matrix(TRUE, 6, 3), p * 0, p * .9)
  for (names in list(c("z", "z", "a/z"), c("z", "", "a/z"), c("z", NA, "a/z"))) {
    bad <- p
    colnames(bad) <- names
    bad_matrices[[length(bad_matrices) + 1L]] <- bad
  }
  for (value in c(-.01, 1.01, Inf, NaN, NA_real_)) {
    bad <- p
    bad[1, 1] <- value
    bad_matrices[[length(bad_matrices) + 1L]] <- bad
  }
  for (bad in bad_matrices) {
    expect_error(venn_abers(truth, bad, precision = 0, method = method), "probability")
    expect_error(predict(fit, bad), "newdata")
    expect_error(brier_score(truth, bad), "probability")
    expect_error(log_loss(truth, bad), "probability")
    expect_error(calibration_plot(truth, bad), "probability")
  }
  near <- x * (1 - 5e-9)
  # Explicit binary64 normalization also keeps this oracle stable on x86.
  normalized <- near / (near[, 1L] + near[, 2L] + near[, 3L])
  close(predict(fit, near), predict(fit, normalized))
  expect_error(predict(fit, x * (1 - 2e-8)), "row sums")
  tiny <- x[1, , drop = FALSE]
  tiny[1, ] <- c(1, 1e-310, 2e-310)
  stopifnot(all(is.finite(predict(fit, tiny))))
}

# Assessment subsets retain unobserved levels. Hand-worked scores include
# zero observed probability and the two-column Brier scale.
assessment <- factor(c("z", "a"), levels = classes)
prob <- rbind(c(.5, .25, .25), c(0, 1, 0))
colnames(prob) <- classes
close(brier_score(assessment, prob), .1875)
close(log_loss(assessment, prob), -log(.5) / 2)
stopifnot(identical(log_loss(assessment, prob[2:1, ]), Inf))
close(brier_score(two_truth, cbind(z = c(.8, .3), a = c(.2, .7))),
      2 * brier_score(c(0, 1), c(.2, .7)))
close(brier_score(assessment, prob[, rev(classes)]), brier_score(assessment, prob))
close(log_loss(assessment, prob[, rev(classes)]), log_loss(assessment, prob))
local({
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) })
  layout <- graphics::par("mfrow")
  summary <- withVisible(calibration_plot(assessment, prob, bins = 2))
  stopifnot(!summary$visible, identical(graphics::par("mfrow"), layout),
            identical(names(summary$value), c("class", "mean_probability", "event_rate", "count")),
            identical(levels(summary$value$class), classes))
  # Literal groups: z has 0/.5, a has .25/1, and the unobserved third
  # class has one group averaging 0 and .25. No binary plotting oracle.
  stopifnot(identical(as.character(summary$value$class), c("z", "z", "a", "a", "a/z")))
  close(summary$value$mean_probability, c(0, .5, .25, 1, .125))
  close(summary$value$event_rate, c(0, 1, 0, 1, 0))
  close(summary$value$count, c(1, 1, 1, 1, 2))
  expect_error(calibration_plot(assessment, prob, bins = 0), "bins")
  stopifnot(identical(graphics::par("mfrow"), layout))
})

  testthat::expect_true(TRUE)
})

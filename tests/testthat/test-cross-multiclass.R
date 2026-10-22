testthat::test_that("cross-multiclass", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

close <- function(x, y) stopifnot(identical(dim(x), dim(y)), all(abs(x - y) <= 1e-10))
expect_error <- function(expr, text) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}

# Unequal single-class folds, including a singleton. Constant scores give
# exact endpoints; equal fold weights differ from pooling or size weighting.
classes <- c("z", "a")
truth <- factor(c("z", "a", "a", "a"), levels = classes)
p <- matrix(.5, 4, 2, dimnames = list(NULL, classes))
fold <- c("A", "B", "B", "B")
x <- list(B = p[1, , drop = FALSE], A = p[1, , drop = FALSE])
for (method in c("ovo", "ovr")) {
  fit <- cross_venn_abers(truth, p, fold, method = method)
  close(predict(fit, x), matrix(c(sqrt(2) - 1, 2 - sqrt(2)), 1))
  close(predict(fit, x, loss = "brier"), matrix(c(27, 37) / 64, 1))
}

# Scores and plots consume each finished matrix, including two-class
# multiclass's unhalved Brier scale. These constants come from the counts
# above: pooled bounds (1/5, 2/5) versus singleton/three-row fold bounds.
local({
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) })
  expected <- list(
    inductive = list(log = c(1/3, 2/3, 7/18, log(3) - 3 * log(2)/4),
                     brier = c(17/50, 33/50, 489/1250, -(log(17/50) + 3 * log(33/50))/4)),
    cross = list(log = c(sqrt(2) - 1, 2 - sqrt(2), (15 - 10 * sqrt(2))/2,
                        -(log(sqrt(2) - 1) + 3 * log(2 - sqrt(2)))/4),
                 brier = c(27/64, 37/64, 889/2048, -(log(27/64) + 3 * log(37/64))/4))
  )
  for (method in c("ovo", "ovr")) {
    fits <- list(inductive = venn_abers(truth, p, method = method),
                 cross = cross_venn_abers(truth, p, fold, method = method))
    newdata <- list(inductive = p[, rev(classes)], cross = list(B = p[, rev(classes)], A = p))
    for (family in names(fits)) for (loss in c("log", "brier")) {
      oracle <- expected[[family]][[loss]]
      probability <- predict(fits[[family]], newdata[[family]], loss = loss)
      close(probability, matrix(rep(oracle[1:2], each = 4), 4))
      close(brier_score(truth, probability), oracle[3L])
      close(log_loss(truth, probability), oracle[4L])
      groups <- calibration_plot(truth, probability[, rev(classes)], bins = 2)
      stopifnot(identical(as.character(groups$class), classes))
      close(groups$mean_probability, oracle[1:2])
      close(groups$event_rate, c(1/4, 3/4))
      close(groups$count, c(4, 4))
    }
  }
})

# Independent reconstruction: public binary IVAP fits (already checked against
# pinned Python fixtures), geometric means, and odds-form OvO coupling. This is
# component parity, not a claim of a public manual Python Cross-multiclass API.
combine <- function(points, problems, classes, method) {
  result <- matrix(0, nrow(points), length(classes), dimnames = list(NULL, classes))
  for (j in seq_along(classes)) {
    if (method == "ovr") {
      result[, j] <- points[, j]
    } else {
      relevant <- which(vapply(problems, function(pair) classes[j] %in% pair, FALSE))
      odds <- lapply(relevant, function(i) {
        q <- if (classes[j] == problems[[i]][2L]) points[, i] else 1 - points[, i]
        (1 - q) / q
      })
      result[, j] <- 1 / (1 + Reduce(`+`, odds))
    }
  }
  result / rowSums(result)
}
score <- function(p, problem) {
  if (length(problem) == 1L) return(p[, problem])
  denominator <- p[, problem[1L]] + p[, problem[2L]]
  ifelse(denominator == 0, .5, p[, problem[2L]] / denominator)
}
set.seed(305)
differs_from_fold_average <- differs_from_pooling <- differs_from_weighting <- FALSE
rounding_order_matters <- prediction_rounding_matters <- FALSE
for (k in 2:4) {
  classes <- c("z", "a", "rest", "a/z")[seq_len(k)]
  # Every fold contains k-1 classes; some pair calibrations are singletons.
  truth <- factor(c(classes[-k], rep(classes[-1L], 3), classes), levels = classes)
  fold <- rep(c(" B", "A", "01"), c(k - 1L, 3L * (k - 1L), k))
  p <- matrix(runif(length(truth) * k), ncol = k)
  p <- p / rowSums(p)
  colnames(p) <- classes
  # For k >= 3, the observed z row has zero z/a pair mass at fit time.
  p[1, ] <- c(rep(0, k - 1L), 1)
  newdata <- setNames(lapply(unique(fold), function(id) {
    x <- matrix(runif(7 * k), 7, k)
    x <- x / rowSums(x)
    x[1, ] <- c(1, rep(0, k - 1L))
    dimnames(x) <- list(paste0("new", 1:7), classes)
    x
  }), unique(fold))
  for (method in c("ovo", "ovr")) {
    problems <- if (method == "ovo") combn(classes, 2, simplify = FALSE) else as.list(classes)
    for (precision in list(NULL, 2, 0, -1)) {
      fit <- cross_venn_abers(truth, p, fold, precision, method = method)
      before <- serialize(fit, NULL)
      endpoints <- NULL
      for (loss in c("log", "brier")) {
        actual <- predict(fit, newdata, loss = loss, type = "details")
        parts <- actual[[2L]]
        assert_conditions(identical(names(actual), c("prediction", if (method == "ovo") "pairs" else "classes")),
            identical(names(parts), c("row", if (method == "ovo") "non_event", "event", "fold", "lower", "upper",
                "probability")), identical(dimnames(actual$prediction), dimnames(newdata[[1L]])), identical(actual$prediction,
                predict(fit, newdata, loss = loss)), identical(predict(fit, newdata, loss = loss, type = "class"),
                factor(classes[max.col(actual$prediction, ties.method = "first")], levels = classes)), !anyDuplicated(parts[c("row",
                "fold", if (method == "ovo") "non_event", "event")]), nrow(parts) == 7L * length(unique(fold)) *
                length(problems), is.integer(parts$row), is.character(parts$fold), is.character(parts$event),
            all(vapply(parts[c("lower", "upper", "probability")], is.double, FALSE)), all(is.finite(actual$prediction)),
            all(actual$prediction >= 0 & actual$prediction <= 1))
        close(rowSums(actual$prediction), rep(1, 7))
        if (!is.null(endpoints)) stopifnot(identical(parts[c("lower", "upper")], endpoints))
        endpoints <- parts[c("lower", "upper")]
        points <- weighted <- pooled <- matrix(0, 7, length(problems))
        fold_points <- array(0, c(7, length(problems), 3))
        for (i in seq_along(problems)) {
          problem <- problems[[i]]
          rows <- if (method == "ovo") truth %in% problem else rep(TRUE, length(truth))
          normalized <- p / rowSums(p)
          scores <- score(normalized, problem)
          if (!is.null(precision)) scores <- round(scores, precision)
          lower <- upper <- matrix(0, 7, 3)
          for (f in seq_along(newdata)) {
            id <- names(newdata)[f]
            cal_rows <- rows & fold == id
            event <- truth[cal_rows] == tail(problem, 1)
            binary <- venn_abers(event, scores[cal_rows])
            input <- newdata[[id]] / rowSums(newdata[[id]])
            query <- score(input, problem)
            expected <- predict(binary, query, loss = loss)
            selected <- parts$event == tail(problem, 1) & parts$fold == id
            if (method == "ovo") selected <- selected & parts$non_event == problem[1L]
            d <- parts[selected, ]
            close(as.matrix(d[c("lower", "upper", "probability")]), as.matrix(expected))
            lower[, f] <- d$lower
            upper[, f] <- d$upper
            fold_points[, i, f] <- d$probability
            if (!is.null(precision)) {
              prediction_rounding_matters <- prediction_rounding_matters ||
                any(abs(predict(binary, round(query, precision), loss = loss)$probability - expected$probability) > 1e-10)
              if (method == "ovo") {
                wrong <- score(round(normalized, precision), problem)
                rounding_order_matters <- rounding_order_matters || any(abs(wrong - scores) > 1e-10)
              }
            }
          }
          # Geometric means are safe for these small reference cases. Production
          # uses log space; existing binary underflow checks protect that merger.
          a <- apply(1 - lower, 1, prod)^(1 / 3)
          b <- apply(upper, 1, prod)^(1 / 3)
          points[, i] <- if (loss == "log") b / (a + b) else rowMeans(upper + (lower^2 - upper^2) / 2)
          weights <- as.numeric(table(factor(fold[rows], levels = names(newdata))))
          if (loss == "log") {
            wa <- exp(as.vector(log(1 - lower) %*% weights) / sum(weights))
            wb <- exp(as.vector(log(upper) %*% weights) / sum(weights))
            weighted[, i] <- wb / (wa + wb)
          } else weighted[, i] <- as.vector(fold_points[, i, ] %*% weights) / sum(weights)
          pooled[, i] <- predict(venn_abers(truth[rows] == tail(problem, 1), scores[rows]),
                                  score(newdata[[1L]], problem), loss = loss)$probability
        }
        close(combine(points, problems, classes, method), actual$prediction)
        wrong_average <- Reduce(`+`, lapply(1:3, function(f) {
          combine(matrix(fold_points[, , f], 7), problems, classes, method)
        })) / 3
        differs_from_fold_average <- differs_from_fold_average || any(abs(wrong_average - actual$prediction) > 1e-5)
        differs_from_pooling <- differs_from_pooling || any(abs(combine(pooled, problems, classes, method) - actual$prediction) > 1e-5)
        differs_from_weighting <- differs_from_weighting || any(abs(combine(weighted, problems, classes, method) - actual$prediction) > 1e-5)
        reordered <- lapply(rev(newdata), function(x) x[, rev(classes), drop = FALSE])
        stopifnot(identical(actual, predict(fit, reordered, loss = loss, type = "details")))
        permuted_fit <- cross_venn_abers(truth, p[, rev(classes), drop = FALSE], fold, precision, method)
        stopifnot(identical(actual, predict(permuted_fit, newdata, loss = loss, type = "details")))
        for (n in 0:1) {
          small <- lapply(newdata, function(x) x[seq_len(n), , drop = FALSE])
          out <- predict(fit, small, loss = loss, type = "details")
          assert_conditions(identical(dim(out$prediction), c(n, k)), identical(out$prediction, actual$prediction[seq_len(n),
              , drop = FALSE]), identical(names(out[[2L]]), names(parts)), identical(vapply(out[[2L]], typeof,
              ""), vapply(parts, typeof, "")), nrow(out[[2L]]) == n * length(problems) * 3L, identical(predict(fit,
              small, loss = loss, type = "class"), predict(fit, newdata, loss = loss, type = "class")[seq_len(n)]))
        }
      }
      stopifnot(identical(serialize(fit, NULL), before))
    }
  }
}
assert_conditions(differs_from_fold_average, differs_from_pooling, differs_from_weighting, rounding_order_matters,
    prediction_rounding_matters)

# Validation is shared with binary Cross and Inductive multiclass; exercise
# the new dispatch and nested list boundary as well as contextual empty pairs.
classes <- c("z", "a", "rest")
truth <- factor(c("z", "a", "a", "rest"), levels = classes)
p <- matrix(1 / 3, 4, 3, dimnames = list(NULL, classes))
fold <- c("B", "B", "A", "A")
x <- list(A = p[1:2, ], B = p[1:2, ])
for (method in c("ovo", "ovr")) {
  fit <- cross_venn_abers(truth, p, fold, method = method)
  for (bad in list(NULL, character(), rep("B", 4), c("B", "B", "A"),
                   c("B", "B", "A", NA), c("B", "B", "A", ""),
                   c(1, 1, 2, Inf), c(1, 1, 2, 0), factor(fold))) {
    expect_error(cross_venn_abers(truth, p, bad, method = method), "fold")
  }
  for (bad in list(NULL, p, unname(x), x[1], c(x, C = x[1]),
                   setNames(x, c("A", "A")), setNames(x, c("A", NA)),
                   setNames(x, c("A", "")), list(A = p, B = p[1, , drop = FALSE]),
                   list(A = as.data.frame(p), B = p))) {
    expect_error(predict(fit, bad), "newdata")
  }
  for (bad in list(unname(p), p[, -1], cbind(p, extra = 0),
                   `colnames<-`(p, c("a", "a", "rest")),
                   `colnames<-`(p, c("z", "", "rest")),
                   `colnames<-`(p, c("z", NA, "rest")), p * 1.001)) {
    expect_error(cross_venn_abers(truth, bad, fold, method = method), "probability")
    expect_error(predict(fit, list(A = bad, B = p)), "newdata")
  }
  for (value in c(NA, NaN, Inf, -Inf, -.1, 1.1)) {
    bad <- p
    bad[1, 1] <- value
    expect_error(cross_venn_abers(truth, bad, fold, method = method), "probability")
    expect_error(predict(fit, list(A = bad, B = p)), "newdata")
  }
  for (bad in list(TRUE, NA_real_, Inf, .5, "2", c(1, 2))) {
    expect_error(cross_venn_abers(truth, p, fold, precision = bad, method = method), "precision")
  }
  expect_error(cross_venn_abers(truth[0], p[0, ], fold[0], method = method), "truth")
  expect_error(cross_venn_abers(truth[-4], p[-4, ], fold[-4], method = method), "Every class")
  expect_error(predict(fit, x, loss = "LOG"), "loss")
  expect_error(predict(fit, x, type = "interval"), "type")
  expect_error(predict(fit, x, precision = 2), "...")
  # Normalization within input tolerance is applied separately to every fold.
  close(predict(cross_venn_abers(truth, p * (1 - 5e-9), fold, method = method),
                lapply(x, function(p) p * (1 + 5e-9))), predict(fit, x))
  numeric_fit <- cross_venn_abers(truth, p, c(20, 20, 1, 1), method = method)
  close(predict(numeric_fit, setNames(x, c("1", "20"))), predict(fit, x))
  printed <- withVisible(print(fit))
  stopifnot(identical(printed$value, fit), !printed$visible)
}
expect_error(cross_venn_abers(truth, p, fold, method = "OVO"), "method")
expect_error(cross_venn_abers(c(0, 1), c(.1, .9), c(1, 2), method = "ovo"), "method")
stopifnot(identical(cross_venn_abers(truth, p, fold), cross_venn_abers(truth, p, fold, method = "ovo")))
empty_fold <- c("only z", "others", "others", "others")
for (context in c("a, rest", "only z")) {
  expect_error(cross_venn_abers(truth, p, empty_fold), context)
}
# OvR permits entire singleton/single-class folds even when OvO cannot.
fit <- cross_venn_abers(truth, p, empty_fold, method = "ovr")
stopifnot(all(is.finite(predict(fit, setNames(x, c("others", "only z"))))))

# Ties use the first global class, even for two classes and reordered names.
for (k in 2:4) {
  classes <- c("z", "a", "rest", "a/z")[seq_len(k)]
  truth <- factor(rep(classes, 2), levels = classes)
  p <- matrix(1 / k, 2 * k, k, dimnames = list(NULL, classes))
  x <- list(B = p[1, rev(classes), drop = FALSE], A = p[1, , drop = FALSE])
  for (method in c("ovo", "ovr")) {
    fit <- cross_venn_abers(truth, p, rep(c("A", "B"), each = k), method = method)
    for (loss in c("log", "brier")) {
      close(predict(fit, x, loss = loss), matrix(1 / k, 1, k))
      stopifnot(identical(predict(fit, x, loss = loss, type = "class"), factor("z", levels = classes)))
    }
  }
}

  testthat::expect_true(TRUE)
})

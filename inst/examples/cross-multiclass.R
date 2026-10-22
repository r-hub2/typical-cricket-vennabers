# After package installation: Rscript inst/examples/cross-multiclass.R
library(vennabers)

# Precomputed out-of-fold probabilities from two retained multiclass models.
# Model A excluded rows 1:2; model B excluded rows 3:6, including preprocessing
# and tuning. Retain both models yourself. No model is trained by this package.
classes <- c("cat", "dog", "bird")
truth <- factor(c("cat", "dog", "cat", "bird", "dog", "bird"), levels = classes)
fold <- c("A", "A", "B", "B", "B", "B")
p <- rbind(c(.7, .2, .1), c(.2, .6, .2), c(.5, .3, .2),
           c(.1, .2, .7), c(.3, .5, .2), c(.2, .3, .5))
colnames(p) <- classes
# The same three untouched assessment observations, in the same row order,
# predicted by each corresponding retained model. List order is immaterial.
new_p <- list(B = rbind(c(.6, .3, .1), c(.1, .4, .5), c(.1, .8, .1)),
              A = rbind(c(.5, .4, .1), c(.2, .2, .6), c(.2, .7, .1)))
new_p <- lapply(new_p, function(x) { colnames(x) <- classes; x })
assessment <- factor(c("cat", "bird", "dog"), levels = classes)

local({
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({ grDevices::dev.off(); unlink(plot_file) })
  for (method in c("ovo", "ovr")) {
    fit <- cross_venn_abers(truth, p, fold, method = method, precision = 2)
    print(fit)
    for (loss in c("log", "brier")) {
      details <- predict(fit, new_p, loss = loss, type = "details")
      probability <- predict(fit, new_p, loss = loss)
      predicted_class <- predict(fit, new_p, loss = loss, type = "class")
      # Reconstruct from the fold endpoints/points, then combine subproblems.
      parts <- details[[2L]]
      problems <- unique(parts[c(if (method == "ovo") "non_event", "event")])
      merged <- vapply(seq_len(nrow(problems)), function(i) {
        selected <- parts$event == problems$event[i]
        if (method == "ovo") selected <- selected & parts$non_event == problems$non_event[i]
        vapply(seq_len(nrow(probability)), function(row) {
          f <- parts[selected & parts$row == row, ]
          if (loss == "log") {
            plogis(mean(log(f$upper)) - mean(log1p(-f$lower)))
          } else mean(f$probability)
        }, 0)
      }, numeric(nrow(probability)))
      recovered <- if (method == "ovr") merged else vapply(classes, function(class) {
        q <- cbind(1 - merged[, problems$non_event == class, drop = FALSE],
                   merged[, problems$event == class, drop = FALSE])
        1 / (1 + rowSums((1 - q) / q))
      }, numeric(nrow(probability)))
      recovered <- recovered / rowSums(recovered)
      stopifnot(max(abs(recovered - probability)) <= 1e-10,
                identical(probability, details$prediction),
                identical(predicted_class,
                          factor(classes[max.col(probability, ties.method = "first")],
                                 levels = classes)))
      print(details)
      print(predicted_class)
      print(c(brier = brier_score(assessment, probability),
              log = log_loss(assessment, probability)))
      groups <- calibration_plot(assessment, probability, bins = 1)
      stopifnot(max(abs(groups$mean_probability - colMeans(probability))) <= 1e-10,
                all(groups$count == length(assessment)), all(groups$event_rate == 1 / 3))
      print(groups)
    }
  }
})
# One shared fold plan serves all pairs/classes. Fold A lacks bird but every
# OvO pair remains nonempty. Empty pairs fail with the pair and fold identified.
# Endpoints describe binary subproblems, with no aggregate confidence interval.
# Manual Cross multiclass is an R extension; Python has no matching public
# manual whole-procedure API. Numerical tests reconstruct binary components.

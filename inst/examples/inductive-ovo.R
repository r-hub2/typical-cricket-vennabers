# After package installation: Rscript inst/examples/inductive-ovo.R
library(vennabers)

# Illustrative caller-supplied outputs from one already fixed multiclass model.
# Calibration and assessment observations are separate from that model's
# training, preprocessing, and tuning; the package does not train the model.
classes <- c("cat", "dog", "bird")
truth <- factor(rep(classes, 2), levels = classes)
p <- rbind(c(.7, .2, .1), c(.2, .6, .2), c(.1, .2, .7),
           c(.5, .3, .2), c(.3, .5, .2), c(.2, .3, .5))
colnames(p) <- classes
fit <- venn_abers(truth, p, precision = 2, method = "ovo")
print(fit)
new_p <- rbind(c(.6, .3, .1), c(.1, .4, .5), c(.1, .8, .1))
colnames(new_p) <- classes
assessment <- factor(c("cat", "bird", "dog"), levels = classes)

local({
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({ grDevices::dev.off(); unlink(plot_file) })
  for (loss in c("log", "brier")) {
    details <- predict(fit, new_p[, rev(classes)], loss = loss, type = "details")
    probability <- predict(fit, new_p, loss = loss)
    predicted_class <- predict(fit, new_p, loss = loss, type = "class")
    reconstructed <- probability
    for (row in seq_len(nrow(new_p))) {
      pairs <- details$pairs[details$pairs$row == row, ]
      for (class in classes) {
        q <- c(pairs$probability[pairs$event == class],
               1 - pairs$probability[pairs$non_event == class])
        reconstructed[row, class] <- 1 / (sum(1 / q) - (length(classes) - 2))
      }
    }
    reconstructed <- reconstructed / rowSums(reconstructed)
    stopifnot(max(abs(reconstructed - probability)) <= 1e-10,
              identical(probability, details$prediction),
              identical(predicted_class,
                        factor(classes[max.col(probability, ties.method = "first")], levels = classes)))
    print(details)
    print(predicted_class)
    print(c(brier = brier_score(assessment, probability), log = log_loss(assessment, probability)))
    groups <- calibration_plot(assessment, probability)
    stopifnot(all(tapply(groups$count, groups$class, sum) == length(assessment)))
    print(groups)
  }
})
# Pair endpoints have only their binary interpretation. Assessment scores use
# the final matrix; Brier sums class errors without dividing by class count.

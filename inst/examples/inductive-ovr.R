# After package installation: Rscript inst/examples/inductive-ovr.R
library(vennabers)

# Illustrative probabilities from an already fixed multiclass model.
# Calibration and assessment are separate from its training and tuning.
classes <- c("cat", "dog", "bird")
truth <- factor(rep(classes, 2), levels = classes)
p <- rbind(c(.7, .2, .1), c(.2, .6, .2), c(.1, .2, .7),
           c(.5, .3, .2), c(.3, .5, .2), c(.2, .3, .5))
colnames(p) <- classes
fit <- venn_abers(truth, p, method = "ovr", precision = 2)
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
    points <- vapply(classes, function(class) {
      parts <- details$classes[details$classes$event == class, ]
      parts$probability[match(seq_len(nrow(new_p)), parts$row)]
    }, numeric(nrow(new_p)))
    stopifnot(max(abs(points / rowSums(points) - probability)) <= 1e-10,
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
# Endpoints concern each binary class-against-rest problem. Scores and plots
# use the final normalized matrix. Manual Brier is an R extension over Python.

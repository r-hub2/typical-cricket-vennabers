# After installation: source(system.file("examples", "binary-options.R", package = "vennabers"))
library(vennabers)

# Precomputed probabilities from a fixed model on its separate calibration set.
# No calibration or test truths were used to train that model or preprocessing.
truth <- c(0, 1, 0, 1)
calibration_probability <- c(.501, .502, .509, .510)
new_probability <- c(.499, .505, .511)
inductive <- venn_abers(truth, calibration_probability, precision = 2)
for (loss in c("log", "brier")) {
  prediction <- predict(inductive, new_probability, loss = loss)
  classes <- predict(inductive, new_probability, loss = loss, type = "class")
  stopifnot(identical(classes, prediction$probability >= .5))
  print(prediction)
}
# Precision affects calibration only; choosing a new precision requires a fit.
# The new values above are deliberately not rounded to two decimal places.

# To select class 0 as the event, reverse truth AND the supplied probabilities.
# Round only after selecting the event; rounded results need not be complements.
other_event <- venn_abers(truth == 0, 1 - calibration_probability, precision = 2)
print(predict(other_event, 1 - new_probability))
# Inductive bounds are VA endpoints, not confidence intervals.

# Cross uses held-out probabilities from each corresponding retained fold model.
# These illustrative folds have unequal sizes; each still receives equal weight.
# See cross-venn-abers.R for complete model fitting and held-out data generation.
cross <- cross_venn_abers(
  truth = c(0, 1, 1, 1), probability = c(.501, .201, .501, .801),
  fold = c("A", "B", "B", "B"), precision = 2
)
# Each column comes from the model that excluded its named calibration fold.
# All columns refer to the same new observations, in the same order.
new_probability <- cbind(B = c(.201, .501, .901), A = c(.199, .501, .801))
test_truth <- c(0, 1, 1)  # Used only for assessment, after fitting.

local({
  plot_file <- tempfile(fileext = ".pdf")
  on.exit(unlink(plot_file), add = TRUE)
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  for (loss in c("log", "brier")) {
    details <- predict(cross, new_probability, loss = loss, type = "details")
    reconstructed <- vapply(seq_len(nrow(new_probability)), function(row) {
      f <- details$folds[details$folds$row == row, ]
      if (loss == "log") {
        plogis(mean(log(f$upper)) - mean(log1p(-f$lower)))
      } else mean(f$probability)
    }, 0)
    stopifnot(max(abs(reconstructed - details$prediction$probability)) < 1e-10,
              identical(predict(cross, new_probability, loss = loss, type = "class"),
                        details$prediction$probability >= .5))
    print(details)
    print(c(brier = brier_score(test_truth, reconstructed),
            log = log_loss(test_truth, reconstructed)))
    calibration_plot(test_truth, reconstructed)
  }
})
# Classes use the combined probability with TRUE at .5, never a fold vote.
# Fold details concern new observations; there is no aggregate Cross interval.
# No improvement in test scores or perfect Cross calibration is promised.

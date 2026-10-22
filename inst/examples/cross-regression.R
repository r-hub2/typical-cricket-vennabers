# After installation: source(system.file("examples", "cross-regression.R",
#                                       package = "vennabers"))
library(vennabers)

# IID observations; fold membership and the separate assessment role are fixed
# before examining outcomes. Unequal folds still receive equal aggregate weight.
set.seed(20260908)
x <- runif(240, -3, 3)
data <- data.frame(id = seq_along(x), x = x, truth = 2 + x + .4 * x^2 + rnorm(240))
development <- data[1:180, ]
assessment <- data[181:240, ]
fold <- sample(rep(c("A", "B", "C"), c(50, 60, 70)))
ids <- unique(fold)
models <- setNames(vector("list", length(ids)), ids)
out_of_fold <- rep(NA_real_, nrow(development))
newdata <- matrix(NA_real_, nrow(assessment), length(ids), dimnames = list(NULL, ids))
for (id in ids) {
  held_out <- fold == id
  training <- development[!held_out, ]
  calibration <- development[held_out, ]
  stopifnot(!any(training$id %in% c(calibration$id, assessment$id)),
            !any(calibration$id %in% assessment$id))
  # No learned preprocessing or tuning outside this training set. The caller
  # retains this exact model for both its calibration and assessment predictions.
  models[[id]] <- lm(truth ~ x, data = training)
  out_of_fold[held_out] <- predict(models[[id]], calibration)
  newdata[, id] <- predict(models[[id]], assessment)
}
fit <- cross_venn_abers_regression(development$truth, out_of_fold, fold, epsilon = .1)
details <- predict(fit, newdata, type = "details")
print(lapply(fit$calibrators, function(x) x[c("observations", "requested", "m", "epsilon")]))
result <- details$prediction
for (column in names(result)) {
  reconstructed <- Reduce(`+`, lapply(details$folds, `[[`, column)) / length(ids)
  stopifnot(all(abs(result[[column]] - reconstructed) <= 1e-10 + 1e-10 * abs(reconstructed)))
}

# Matching names cannot prove model provenance or semantic row identity. Each
# column above uses the same assessment rows in the same order. A model refitted
# on all development data cannot replace these fold models. Predict trains nothing.
# Assessment truths enter only in scoring. Aggregated bounds are descriptive
# fold summaries, with no additional coverage or selector guarantee; improvement
# in point-prediction scores is not promised.
scores <- vapply(list(uncalibrated = rowMeans(newdata), calibrated = result$prediction),
                 function(prediction) {
                   c(RMSE = sqrt(mean((assessment$truth - prediction)^2)),
                     MAE = mean(abs(assessment$truth - prediction)))
                 }, c(RMSE = 0, MAE = 0))
stopifnot(all(is.finite(scores)), nrow(result) == nrow(assessment),
          all(result$lower <= result$prediction), all(result$prediction <= result$upper))
for (column in colnames(scores)) {
  residual <- assessment$truth - if (column == "calibrated") result$prediction else rowMeans(newdata)
  stopifnot(abs(scores["RMSE", column]^2 - sum(residual^2) / length(residual)) < 1e-12,
            abs(scores["MAE", column] - sum(abs(residual)) / length(residual)) < 1e-12)
}
print(scores)

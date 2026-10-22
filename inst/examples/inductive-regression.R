# After installation: source(system.file("examples", "inductive-regression.R",
#                                       package = "vennabers"))
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

# Independent IID observations; roles are fixed before examining their outcomes.
set.seed(20260907)
x <- runif(240, -3, 3)
data <- data.frame(id = seq_along(x), x = x, truth = 2 + x + .4 * x^2 + rnorm(240))
training <- data[1:120, ]
calibration <- data[121:180, ]
assessment <- data[181:240, ]
stopifnot(!anyDuplicated(c(training$id, calibration$id, assessment$id)))

# Only training data determine this model; no preprocessing or tuning uses the
# calibration or assessment observations. The caller retains the fixed model.
model <- lm(truth ~ x, data = training)
calibration_prediction <- predict(model, calibration)
fit <- venn_abers_regression(calibration$truth, calibration_prediction, epsilon = .1)
uncalibrated <- predict(model, assessment)
result <- predict(fit, uncalibrated)
print(list(requested = fit$requested, m = fit$m, epsilon = fit$epsilon))

# Assessment truths enter only here. These are point-prediction scores in target
# units; the VA bounds concern regularized conditional expectations, not coverage
# of individual outcomes. Neither improvement nor a general selector guarantee
# for this particular point merger is promised.
scores <- vapply(list(uncalibrated = uncalibrated, calibrated = result$prediction),
                 function(prediction) {
                   c(RMSE = sqrt(mean((assessment$truth - prediction)^2)),
                     MAE = mean(abs(assessment$truth - prediction)))
                 }, c(RMSE = 0, MAE = 0))
assert_conditions(all(is.finite(scores)), nrow(result) == nrow(assessment),
          all(result$lower <= result$prediction), all(result$prediction <= result$upper),
          abs(scores["RMSE", "calibrated"]^2 -
              sum((assessment$truth - result$prediction)^2) / nrow(assessment)) < 1e-12,
          abs(scores["MAE", "calibrated"] -
              sum(abs(assessment$truth - result$prediction)) / nrow(assessment)) < 1e-12)
print(scores)

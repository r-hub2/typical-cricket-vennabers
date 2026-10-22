# Run from the checkout root after installation:
# Rscript inst/examples/cross-venn-abers.R
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

set.seed(20260906)
observations <- data.frame(x1 = rnorm(900L), x2 = rnorm(900L))
observations$truth <- rbinom(900L, 1L, stats::plogis(
  -0.4 + observations$x1 - 0.6 * observations$x2
))  # The event is 1 (equivalently TRUE).

# Set aside test observations before any model fitting or calibration.
test_rows <- sample.int(nrow(observations), 300L)
test <- observations[test_rows, ]
development <- observations[-test_rows, ]
fold_ids <- c("A", "B", "C")
fold <- sample(rep(fold_ids, length.out = nrow(development)))
models <- setNames(vector("list", length(fold_ids)), fold_ids)
out_of_fold_probability <- numeric(nrow(development))

for (id in fold_ids) {
  held_out <- fold == id
  # Exclusion must cover the complete learned pipeline, including any learned
  # preprocessing or tuning. This example needs neither.
  models[[id]] <- stats::glm(truth ~ x1 + x2,
    data = development[!held_out, ], family = stats::binomial()
  )
  # Restore held-out probabilities to their original development positions.
  out_of_fold_probability[held_out] <- predict(models[[id]],
    newdata = development[held_out, ], type = "response"
  )
}

# Retain these models: each predicts the same test rows in the same order.
# A single all-data refit cannot replace them. Names cannot prove model
# provenance, training exclusion, unique observations, or semantic alignment.
test_probability <- vapply(models, function(model) {
  predict(model, newdata = test[c("x1", "x2")], type = "response")
}, numeric(nrow(test)))

# The package owns per-fold calibration and equal-weight log-loss aggregation.
# Model fitting and fold creation remain caller responsibilities.
calibrator <- cross_venn_abers(
  development$truth, out_of_fold_probability, fold
)
prediction <- predict(calibrator, test_probability)
assert_conditions(
  is.data.frame(prediction), nrow(prediction) == nrow(test),
  identical(names(prediction), "probability"),
  is.double(prediction$probability), all(is.finite(prediction$probability)),
  all(prediction$probability >= 0 & prediction$probability <= 1)
)

# Test truth is used only for evaluation. No score improvement is promised.
scores <- c(
  brier_score = brier_score(test$truth, prediction$probability),
  log_loss = log_loss(test$truth, prediction$probability)
)
stopifnot(all(is.finite(scores)))
print(scores)

# A temporary PDF also works without an interactive display. Close and remove
# it even on error. For an interactive plot, call calibration_plot() directly.
local({
  plot_file <- tempfile(fileext = ".pdf")
  on.exit(unlink(plot_file), add = TRUE)
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  summary <- calibration_plot(test$truth, prediction$probability)
  stopifnot(sum(summary$count) == nrow(test))
})

# Cross output is a regularized event probability, with no public aggregate
# interval. Individual fold diagnostics are available with type = "details".
# The combined probability inherits no IVAP perfect-calibration
# guarantee and makes no confidence-interval coverage claim.

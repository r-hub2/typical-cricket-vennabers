# Run after installation: Rscript inst/examples/shared-multiclass-grid.R
library(vennabers)

data <- iris
names(data)[5L] <- "truth"
base <- workflows::workflow() |>
  workflows::add_recipe(recipes::recipe(truth ~ ., data) |>
    recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
  workflows::add_model(parsnip::multinom_reg(engine = "nnet"))

for (classification in c("ovo", "ovr")) for (strategy in c("joint", "two_stage")) {
  tuning <- list(strategy = strategy, selection_prop = .4,
    model = list(method = "grid", grid = data.frame(degree = c(1L, 2L)),
      resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)),
    calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")),
      resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
  fitted <- calibrated_workflow(base, "truth", classification = classification,
    method = "cross", v = 2L, prediction_model = "full",
    loss = tune::tune("va_loss"), tuning = tuning, seed = 19) |>
    generics::fit(data)
  selection <- venn_abers_tuning(fitted)
  probability <- predict(fitted, data[1:3, ], type = "prob")
  stopifnot(identical(names(probability), paste0(".pred_", levels(data$truth))),
            max(abs(rowSums(probability) - 1)) < 1e-8,
            "Full1" %in% names(selection$models))
  print(c(classification = classification, strategy = strategy,
          calibration = selection$calibration$chosen,
          full_model = selection$models$Full1$chosen))
  print(probability)
}

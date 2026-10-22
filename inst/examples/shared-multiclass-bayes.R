# Run after installation: Rscript inst/examples/shared-multiclass-bayes.R
library(vennabers)

data <- iris
names(data)[5L] <- "truth"
base <- workflows::workflow() |>
  workflows::add_formula(truth ~ Sepal.Length + Petal.Length) |>
  workflows::add_model(parsnip::decision_tree(tree_depth = tune::tune("depth"),
    mode = "classification", engine = "rpart"))
folds <- function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)
model_params <- dials::parameters(dials::new_quant_param(type = "integer",
  range = c(1L, 8L), inclusive = c(TRUE, TRUE), label = c(depth = "Depth")))
calibration_params <- dials::parameters(dials::new_qual_param(type = "character",
  values = c("log", "brier"), label = c(va_loss = "Loss")))

for (classification in c("ovo", "ovr")) for (strategy in c("joint", "two_stage")) {
  tuning <- list(strategy = strategy, selection_prop = .4,
    model = list(method = "bayes", initial = 4L, iter = 1L,
      param_info = model_params, resamples = folds,
      control = tune::control_bayes(seed = 4, save_pred = TRUE)),
    calibration = list(method = "bayes", initial = 2L, iter = 1L,
      param_info = calibration_params, resamples = folds,
      control = tune::control_bayes(seed = 5, save_pred = TRUE)))
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

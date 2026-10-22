# Run after installation: Rscript inst/examples/pairwise-grid.R
library(vennabers)

# Three overlapping classes: unlike iris, no class pair is separable in its
# small local splits, so the base GLMs converge.
data <- withr::with_seed(5, {
  truth <- factor(rep(c("setosa", "versicolor", "virginica"), 50))
  data.frame(Sepal.Length = stats::rnorm(150, c(5.4, 5.9, 6.4)[truth], .7), truth = truth)
})
base <- workflows::workflow() |>
  workflows::add_recipe(recipes::recipe(truth ~ Sepal.Length, data) |>
    recipes::step_poly(Sepal.Length, degree = tune::tune("degree"))) |>
  workflows::add_model(parsnip::logistic_reg(engine = "glm"))

for (strategy in c("joint", "two_stage")) {
  tuning <- list(strategy = strategy, selection_prop = .6,
    model = list(method = "grid", grid = data.frame(degree = c(1L, 2L))),
    calibration = list(method = "grid",
      grid = data.frame(va_loss = c("log", "brier"),
                        rounding = c("none", "-2")),
      resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)))
  if (strategy == "two_stage") tuning$model$resamples <- function(d)
    rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)
  outer <- if (strategy == "joint")
    calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
      loss = tune::tune("va_loss"), precision = tune::tune("rounding"),
      tuning = tuning, seed = 19) else
    calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
      method = "cross", v = 2L, prediction_model = "full",
      loss = tune::tune("va_loss"), precision = tune::tune("rounding"),
      tuning = tuning, seed = 19)
  fitted <- generics::fit(outer, data)
  selection <- venn_abers_tuning(fitted)
  probability <- predict(fitted, data[1:3, ], type = "prob")
  stopifnot(identical(names(probability), paste0(".pred_", levels(data$truth))),
            max(abs(rowSums(probability) - 1)) < 1e-8,
            identical(names(selection$models), paste0("pair_", 1:3)),
            if (strategy == "two_stage")
              all(vapply(selection$models, function(pair) "Full1" %in% names(pair), TRUE)) else
              nrow(selection$joint$candidates) == 16L)
  print(c(strategy = strategy, calibration = selection$calibration$chosen))
  print(if (strategy == "joint") selection$joint$selected$pairs else
    lapply(selection$models, function(pair) vapply(pair, `[[`, "", "chosen")))
  print(probability)
}

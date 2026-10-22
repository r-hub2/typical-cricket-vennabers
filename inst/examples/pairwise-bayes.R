# Run after installation: Rscript inst/examples/pairwise-bayes.R
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
folds <- function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)
tuning <- list(strategy = "joint", selection_prop = .6,
  model = list(method = "bayes", initial = 2L, iter = 1L,
    param_info = dials::parameters(dials::new_quant_param(type = "integer",
      range = c(1L, 2L), inclusive = c(TRUE, TRUE), label = c(degree = "Degree"))),
    resamples = folds, control = tune::control_bayes(seed = 9, save_pred = TRUE)),
  calibration = list(method = "bayes", initial = 2L, iter = 1L,
    param_info = dials::parameters(dials::new_qual_param(type = "character",
      values = c("log", "brier"), label = c(va_loss = "Loss"))),
    resamples = folds, control = tune::control_bayes(seed = 10, save_pred = TRUE)))
fitted <- calibrated_workflow(base, "truth", classification = "ovo",
  model_training = "pairwise", prediction_model = "full",
  loss = tune::tune("va_loss"), tuning = tuning, seed = 11) |>
  generics::fit(data)
selection <- venn_abers_tuning(fitted)
stopifnot(identical(names(selection$models), paste0("pair_", 1:3)),
  selection$calibration$chosen %in% selection$calibration$candidates$id)
print(vapply(selection$models, function(pair) pair$Full1$chosen, ""))
print(predict(fitted, data[1:3, ], type = "prob"))

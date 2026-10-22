# IVAR/CVAR grid selection: local learning, independent assessment.
library(vennabers)
set.seed(19)
data <- data.frame(x = stats::rnorm(360))
data$truth <- data$x^2 + stats::rnorm(360, sd = .3)
development <- data[1:300, ]
assessment <- data[301:360, ]
base <- workflows::workflow() |>
  workflows::add_recipe(recipes::recipe(truth ~ x, development) |>
    recipes::step_poly(x, degree = tune::tune("degree"))) |>
  workflows::add_model(parsnip::linear_reg(engine = "lm"))

for (strategy in c("joint", "two_stage")) for (method in c("inductive", "cross")) {
  config <- list(strategy = strategy, selection_prop = .5,
    model = list(method = "grid", grid = data.frame(degree = c(1L, 2L)),
      resamples = function(d) rsample::vfold_cv(d, v = 2)),
    calibration = list(method = "grid", grid = data.frame(moderation = c(.4, .6)),
      resamples = function(d) rsample::vfold_cv(d, v = 2),
      report_metrics = list(yardstick::mae)))
  args <- list(workflow = base, outcome = "truth", method = method,
    prediction_model = "full", epsilon = tune::tune("moderation"),
    tuning = config, seed = 31)
  if (method == "cross") args$v <- 2L
  fitted <- generics::fit(do.call(calibrated_workflow, args, quote = TRUE), development)
  choices <- venn_abers_tuning(fitted)
  stopifnot(!length(intersect(choices$calibration$selection_rows,
    choices$calibration$final_rows)), "Full1" %in% names(choices$models))
  print(c(strategy = strategy, method = method,
    epsilon = choices$calibration$selected$moderation,
    full_degree = choices$models$Full1$selected$degree))
  prediction <- predict(fitted, assessment)$.pred
  stopifnot(all(is.finite(prediction)))
  print(c(rmse = sqrt(mean((assessment$truth - prediction)^2)),
    mae = mean(abs(assessment$truth - prediction))))
  print(choices$fit_roles)
}

# To search m instead, use m = tune::tune("moderation") and a positive integer
# grid, for example c(1L, 2L), omitting epsilon. Every inner calibration size k
# must support 1 <= m <= floor((k-1)/2); invalid candidates are recorded/excluded.
# S selects moderation, F fits the final predictor; assessment stays untouched.
# Joint scores complete calibrated predictions, two_stage scores raw fits first.
# RMSE is the default selection metric; MAE above is an independent report only.
# Full1 chooses its own model parameters; it does not inherit a fold's choice.
# CVAR averages fold points and each bound separately, with equal weights.
# These bounds describe regularized conditional means, not prediction intervals.

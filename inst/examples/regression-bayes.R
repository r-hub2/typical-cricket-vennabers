# IVAR/CVAR Bayes selection on development data; assessment remains untouched.
library(vennabers)
set.seed(19)
data <- data.frame(x = stats::rnorm(240))
data$truth <- data$x^2 + stats::rnorm(240, sd = .3)
development <- data[1:200, ]
assessment <- data[201:240, ]
base <- workflows::workflow() |>
  workflows::add_recipe(recipes::recipe(truth ~ x, development) |>
    recipes::step_poly(x, degree = tune::tune("degree"))) |>
  workflows::add_model(parsnip::linear_reg(engine = "lm"))
model_params <- dials::parameters(dials::new_quant_param(
  type = "integer", range = c(1L, 3L), inclusive = c(TRUE, TRUE),
  label = c(degree = "Degree")))
moderation_params <- dials::parameters(dials::new_quant_param(
  type = "double", range = c(.3, .4), inclusive = c(TRUE, TRUE),
  label = c(moderation = "Moderation")))

for (strategy in c("joint", "two_stage")) for (method in c("inductive", "cross")) {
  search <- function(params, seed) list(method = "bayes", initial = 3L, iter = 1L,
    param_info = params, resamples = function(d) rsample::vfold_cv(d, v = 2),
    control = tune::control_bayes(seed = seed, save_pred = TRUE, no_improve = 2L))
  config <- list(strategy = strategy, selection_prop = .5,
    model = search(model_params, 3), calibration = search(moderation_params, 4))
  args <- list(workflow = base, outcome = "truth", method = method,
    prediction_model = "full", epsilon = tune::tune("moderation"),
    tuning = config, seed = 31)
  if (method == "cross") args$v <- 2L else args$calibration_prop <- .5
  fitted <- generics::fit(do.call(calibrated_workflow, args, quote = TRUE), development)
  choices <- venn_abers_tuning(fitted)
  stopifnot(!length(intersect(choices$calibration$selection_rows,
    choices$calibration$final_rows)), "Full1" %in% names(choices$models))
  prediction <- predict(fitted, assessment)$.pred
  print(c(strategy = strategy, method = method,
    epsilon = choices$calibration$selected$moderation,
    full_degree = choices$models$Full1$selected$degree,
    rmse = sqrt(mean((assessment$truth - prediction)^2))))
}

# To search m instead, use m = tune::tune("moderation") and an integer
# dials range such as 1:2, omitting epsilon. Every inner calibration size
# must support the proposed moderation value. Joint scores complete calibrated
# predictions; two_stage selects the raw fit first. Full1 chooses locally.

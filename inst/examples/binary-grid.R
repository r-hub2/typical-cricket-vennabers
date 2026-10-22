# Binary grid selection. Every callback sees only its current local data.
library(vennabers)
set.seed(19)
data <- data.frame(x = stats::rnorm(600))
data$truth <- factor(ifelse(stats::runif(600) < stats::plogis(data$x), "yes", "no"),
                     levels = c("yes", "no"))
development <- data[1:400, ]
validation <- data[401:500, ]
assessment <- data[501:600, ]
recipe <- recipes::recipe(truth ~ x, development) |>
  recipes::step_poly(x, degree = tune::tune("degree"))
base <- workflows::workflow() |>
  workflows::add_recipe(recipe) |>
  workflows::add_model(parsnip::logistic_reg(engine = "glm"))

for (strategy in c("joint", "two_stage")) {
  config <- list(strategy = strategy, selection_prop = .3,
    model = list(method = "grid", grid = data.frame(degree = c(1L, 2L)),
      resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)),
    calibration = list(method = "grid",
      grid = data.frame(va_loss = c("log", "brier"), va_rounding = c("none", "-2")),
      resamples = function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0),
      primary_metric = yardstick::brier_class))
  outer <- calibrated_workflow(base, "truth", seed = 31,
    loss = tune::tune("va_loss"), precision = tune::tune("va_rounding"), tuning = config)
  fitted <- generics::fit(outer, development)
  choices <- venn_abers_tuning(fitted)
  stopifnot(!length(intersect(choices$calibration$selection_rows,
                              choices$calibration$final_rows)))
  print(c(strategy = strategy, calibration = choices$calibration$chosen,
          model = choices$models$Inductive1$chosen))
  print(head(predict(fitted, assessment, type = "prob")))
}

# Narrow loss-only choice after an already fixed fit: no refitting.
fixed <- workflows::workflow() |>
  workflows::add_formula(truth ~ x) |>
  workflows::add_model(parsnip::logistic_reg(engine = "glm"))
fitted_fixed <- generics::fit(calibrated_workflow(fixed, "truth", seed = 31), development)
losses <- vapply(c("log", "brier"), function(loss) {
  probability <- predict(fitted_fixed, validation, type = "prob", loss = loss)$.pred_yes
  yardstick::mn_log_loss_vec(validation$truth, probability, event_level = "first")
}, numeric(1))
chosen_loss <- c("log", "brier")[which.min(losses)]
print(head(predict(fitted_fixed, assessment, type = "prob", loss = chosen_loss)))

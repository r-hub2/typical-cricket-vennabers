## ----setup, include=FALSE-----------------------------------------------------
# The code needs optional packages from Suggests; without them it is shown
# but not run.
required <- c("parsnip", "workflows", "hardhat", "rsample", "recipes", "rlang",
              "tibble", "withr", "generics", "tidyselect", "tune", "dials",
              "yardstick")
unavailable <- required[!vapply(required, requireNamespace, logical(1),
                                quietly = TRUE)]
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", error = FALSE,
                      eval = length(unavailable) == 0)

## ----missing-packages, eval=length(unavailable) > 0, echo=FALSE, results="asis"----
# cat("> **Note:** The code in this vignette was not run because these",
#     "suggested packages are not installed:",
#     paste0("`", unavailable, "`", collapse = ", "), "\n")

## ----data-and-recipe----------------------------------------------------------
library(vennabers)
set.seed(19)
data <- data.frame(x = rnorm(800))
data$truth <- factor(ifelse(runif(800) < plogis(.5 * data$x),
                            "yes", "no"), levels = c("yes", "no"))
set.seed(20)
split <- rsample::initial_split(data, prop = .8, strata = truth)
development <- rsample::analysis(split)

recipe <- recipes::recipe(truth ~ x, data = development) |>
  recipes::step_poly(x, degree = tune::tune("degree"))
base <- workflows::workflow() |>
  workflows::add_recipe(recipe) |>
  workflows::add_model(parsnip::logistic_reg(engine = "glm"))

# Always construct resamples from the supplied local data.
folds <- function(data) {
  rsample::vfold_cv(data, v = 2, strata = truth, pool = 0)
}

## ----joint-grid---------------------------------------------------------------
config <- list(
  strategy = "joint",
  selection_prop = .3,
  model = list(method = "grid", grid = data.frame(degree = c(1L, 2L)),
               resamples = folds, primary_metric = yardstick::brier_class),
  calibration = list(method = "grid",
    grid = data.frame(va_loss = c("log", "brier")),
    resamples = folds, primary_metric = yardstick::brier_class)
)
tuned_spec <- calibrated_workflow(base, "truth", seed = 31,
  loss = tune::tune("va_loss"), tuning = config)
joint_fit <- generics::fit(tuned_spec, development)
choices <- venn_abers_tuning(joint_fit)
choices$calibration$selected
choices$models$Inductive1$selected
c(selection_rows = length(choices$calibration$selection_rows),
  final_fit_rows = length(choices$calibration$final_rows))
stopifnot(length(intersect(choices$calibration$selection_rows,
                           choices$calibration$final_rows)) == 0L)

## ----two-stage-grid-----------------------------------------------------------
two_stage <- config
two_stage$strategy <- "two_stage"
two_stage_spec <- calibrated_workflow(base, "truth", seed = 31,
  loss = tune::tune("va_loss"), tuning = two_stage)
two_stage_fit <- generics::fit(two_stage_spec, development)
venn_abers_tuning(two_stage_fit)$models$Inductive1$selected
venn_abers_tuning(two_stage_fit)$calibration$selected

## ----bayes--------------------------------------------------------------------
bayes_config <- list(
  strategy = "joint",
  model = list(
    method = "bayes", initial = 4L, iter = 1L,
    param_info = dials::parameters(list(degree = dials::degree_int(c(1L, 6L)))),
    resamples = folds,
    control = tune::control_bayes(seed = 4, save_pred = TRUE, uncertain = Inf)
  )
)
bayes_fit <- calibrated_workflow(base, "truth", tuning = bayes_config,
                                 seed = 31) |>
  generics::fit(development)
venn_abers_tuning(bayes_fit)$models$Inductive1$selected

## ----outer-resampling---------------------------------------------------------
set.seed(61)
outer_folds <- rsample::vfold_cv(development, v = 2, strata = truth, pool = 0)
evaluation <- tune::fit_resamples(
  tuned_spec, resamples = outer_folds,
  metrics = yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss),
  control = tune::control_resamples(save_pred = TRUE,
                                    pkgs = c("vennabers", "stats"))
)
tune::collect_metrics(evaluation)
tune::collect_notes(evaluation)

## ----final-assessment---------------------------------------------------------
final <- tune::last_fit(
  tuned_spec, split = split,
  metrics = yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
)
tune::collect_metrics(final)
final_workflow <- tune::extract_workflow(final)
head(predict(final_workflow, rsample::assessment(split)["x"], type = "prob"))

## ----further-examples---------------------------------------------------------
system.file("examples", c("regression-grid.R", "shared-multiclass-grid.R",
  "pairwise-grid.R", "postprocessing-grid.R", "postprocessing-bayes.R"),
  package = "vennabers")


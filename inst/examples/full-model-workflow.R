# Fixed binary settings; the same permitted splits compare two procedures.
source(system.file("examples", "assertions.R", package = "vennabers"))
packages <- c("parsnip", "workflows", "hardhat", "rsample", "recipes", "rlang",
              "tibble", "withr", "generics", "tidyselect")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  library(vennabers)
  set.seed(12)
  data <- data.frame(x = rnorm(160))
  data$truth <- factor(ifelse(runif(160) < plogis(data$x), "yes", "no"), levels = c("yes", "no"))
  final_split <- rsample::make_splits(list(analysis = 1:120, assessment = 121:160), data)
  development <- rsample::analysis(final_split)
  assessment <- rsample::assessment(final_split)
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x, development) |> recipes::step_normalize(x)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (method in c("inductive", "cross")) {
    # This function receives only the current local fit domain, also in last_fit.
    plan <- function(data) {
      if (method == "cross") rsample::vfold_cv(data, v = 3, strata = truth) else
        rsample::make_splits(list(analysis = 1:90, assessment = 91:120), data)
    }
    shared_plan <- withr::with_seed(42, plan(development))
    fits <- lapply(c("matched", "full"), function(mode) generics::fit(
      calibrated_workflow(base, "truth", method = method, prediction_model = mode,
                          splits = shared_plan, seed = 43), development))
    matched <- workflows::extract_fit_engine(fits[[1L]])
    full <- workflows::extract_fit_engine(fits[[2L]])
    stopifnot(identical(matched$calibrator, full$calibrator),
              identical(matched$calibration, full$calibration),
              identical(names(full$models), "Full1"))
    p <- lapply(fits, predict, new_data = assessment, type = "prob")
    print(data.frame(method, prediction_model = c("matched", "full"),
      brier = vapply(p, function(x) brier_score(assessment$truth == "yes", x$.pred_yes), 0)))
    print(full)
    print(head(venn_abers_details(fits[[2L]], assessment)$details))
    # All source calibration rows/values remain inspectable after source models
    # are discarded. The full recipe is learned anew on all 120 local rows.
    print(head(full$calibration[[1L]]))
    print(workflows::extract_recipe(full$models$Full1))
    if (all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
      model <- calibrated_workflow(base, "truth", method = method, prediction_model = "full",
                                    splits = shared_plan, seed = 43)
      final <- tune::last_fit(model, final_split,
        metrics = yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss))
      prediction <- tune::collect_predictions(final)
      assert_conditions(nrow(tune::collect_notes(final)) == 0L,
                identical(prediction$.row, 121:160), !anyDuplicated(prediction$.row),
                all(prediction$.config == "pre0_mod0_post0"),
                identical(as.matrix(prediction[names(p[[2L]])]), as.matrix(p[[2L]])))
    }
  }
  # Full transfer costs 2 (Inductive) or K+1 (Cross) base fits. Prediction uses
  # one full model per batch; every held-out calibrator still participates.
  # This model change does not transfer calibration/coverage guarantees and need
  # not improve quality. Comparing procedures is selection: use fresh untouched
  # data if the comparison above determines the procedure eventually deployed.
})

# A previously learned recipe, model and Tailor remain frozen. The user must
# keep new VA calibration data outside ALL earlier training/model selection.
# This adopts the existing workflow unchanged; calibrated_workflow() instead
# trains new base workflows on its own split.
packages <- c("workflows", "parsnip", "recipes", "tailor", "probably", "rsample",
              "hardhat", "rlang", "tibble", "withr", "generics", "tidyselect", "yardstick")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  set.seed(613)
  data <- data.frame(x = rnorm(450), z = rnorm(450))
  data$truth <- factor(ifelse(runif(450) < plogis(data$x - .5 * data$z), "yes", "no"),
                       levels = c("yes", "no"))
  data$x[c(2, 155, 255, 355)] <- NA_real_
  training <- data[1:150, ]
  earlier_postprocessing <- data[151:250, ]
  calibration <- data[251:350, ]
  validation <- data[351:400, ]
  assessment <- data[401:450, ]

  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ ., training) |>
      recipes::step_impute_median(x) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm")) |>
    workflows::add_tailor(tailor::tailor() |>
      tailor::adjust_probability_calibration(method = "logistic", smooth = FALSE))
  # These fits happened earlier, on their own data.
  frozen <- generics::fit(base, training, data_calibration = earlier_postprocessing)

  # Only the VA calibrator learns now. No bake() and no subsequent fit().
  calibrated <- vennabers::calibrate_workflow(frozen, calibration, "truth")
  print(calibrated)  # Base fits: 0; the adopted workflow is the only model

  # Narrow Log/Brier choice on separate validation rows with existing predict()
  # and yardstick calls. Nothing is refitted; pass the chosen rule explicitly.
  # Validation rows are selection data, so assess on still other rows.
  losses <- vapply(c(log = "log", brier = "brier"), function(loss) {
    p <- predict(calibrated, validation[c("x", "z")], type = "prob", loss = loss)
    yardstick::mn_log_loss_vec(validation$truth, p$.pred_yes)
  }, numeric(1))
  chosen <- names(losses)[which.min(losses)]
  print(losses)

  new_data <- assessment[c("x", "z")]
  probability <- predict(calibrated, new_data, type = "prob", loss = chosen)
  classes <- predict(calibrated, new_data, type = "class", loss = chosen)
  details <- vennabers::venn_abers_details(calibrated, new_data, loss = chosen)
  print(c(log_loss = vennabers::log_loss(assessment$truth == "yes", probability$.pred_yes),
          brier = vennabers::brier_score(assessment$truth == "yes", probability$.pred_yes)))
  stopifnot(identical(levels(classes$.pred_class), levels(training$truth)),
            identical(details$prediction, probability),
            identical(calibrated$models[[1L]], frozen), calibrated$loss == "log")

  path <- tempfile(fileext = ".rds")
  saveRDS(calibrated, path)
  stopifnot(identical(probability, predict(readRDS(path), new_data, loss = chosen)))
  unlink(path)
  # Retain the required packages. Engines using external resources may need
  # their own persistence procedure. RDS is not a universal portability promise.
})

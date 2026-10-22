# Separate earlier model learning, earlier postprocessing, new VA calibration,
# and independent assessment. All calls use the complete learned workflow.
packages <- c("workflows", "parsnip", "recipes", "tailor", "probably", "nnet", "mgcv",
              "rsample", "hardhat", "rlang", "tibble", "withr", "generics", "tidyselect")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  set.seed(614)
  data <- data.frame(x = rnorm(400), z = rnorm(400))
  data$truth <- factor(rep(c("third", "first", "second"), length.out = 400),
                       levels = c("third", "first", "second"))
  training <- data[1:150, ]
  earlier_postprocessing <- data[151:250, ]
  calibration <- data[251:350, ]
  assessment <- data[351:400, ]
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ ., training) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::multinom_reg(penalty = .1) |>
      parsnip::set_engine("nnet", trace = FALSE)) |>
    workflows::add_tailor(tailor::tailor() |>
      tailor::adjust_probability_calibration(method = "multinomial", smooth = TRUE))
  frozen <- generics::fit(base, training, data_calibration = earlier_postprocessing)
  for (method in c("ovo", "ovr")) {
    calibrated <- vennabers::calibrate_workflow(frozen, calibration, "truth", classification = method, precision = 2)
    new <- assessment[c("x", "z")]
    probability <- predict(calibrated, new)
    classes <- predict(calibrated, new, type = "class")
    details <- vennabers::venn_abers_details(calibrated, new, loss = "brier")
    stopifnot(identical(levels(classes$.pred_class), levels(training$truth)),
              identical(details$prediction, predict(calibrated, new, loss = "brier")),
              identical(calibrated$models[[1L]], frozen))
    print(head(probability))
  }

  # Numeric truth (including 0/1) stays regression because the workflow says so.
  data$truth <- 2 + data$x + cos(1:400)
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ ., data[1:150, ]) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm")) |>
    workflows::add_tailor(tailor::tailor() |>
      tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE))
  frozen <- generics::fit(base, data[1:150, ], data_calibration = data[151:250, ])
  calibrated <- vennabers::calibrate_workflow(frozen, data[251:350, ], "truth", epsilon = .1)
  points <- predict(calibrated, data[351:400, ], type = "numeric")
  details <- vennabers::venn_abers_details(calibrated, data[351:400, ])
  stopifnot(identical(points, details$prediction), identical(calibrated$models[[1L]], frozen),
            calibrated$parameters$effective$m == 5)
  print(c(rmse = sqrt(mean((points$.pred - data$truth[351:400])^2))))
  # details$details contains row, fold, prediction, lower and upper. These IVAR
  # bounds describe regularized conditional expectations, not prediction intervals.
})

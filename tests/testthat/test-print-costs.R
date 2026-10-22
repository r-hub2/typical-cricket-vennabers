testthat::test_that("print-costs", {
library(vennabers)
packages <- c("parsnip", "workflows", "rsample", "withr", "nnet")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  # Public print output describes the complete fixed procedure, including discarded fits.
  expected <- data.frame(method = c("inductive", "inductive", "cross", "cross"),
    mode = c("matched", "full", "matched", "full"), fits = c(1, 2, 5, 6),
    retained = c(1, 1, 5, 1), calibration = c(24, 24, 120, 120))
  for (path in c("binary", "ovo", "ovr", "pairwise", "regression")) {
    data <- data.frame(x = sin(1:120), truth = if (path == "regression") cos(1:120) else
      factor(rep(if (path == "binary") c("a", "b") else c("a", "b", "c"), length.out = 120)))
    spec <- if (path == "regression") parsnip::linear_reg(engine = "lm") else
      if (path %in% c("binary", "pairwise")) parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(penalty = 0) |> parsnip::set_engine("nnet", trace = FALSE)
    base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |> workflows::add_model(spec)
    for (i in 1:4) {
      model <- do.call(calibrated_workflow, c(list(workflow = base, outcome = "truth",
        method = expected$method[i], prediction_model = expected$mode[i], strata = FALSE, seed = 42),
        if (path != "regression") list(classification = if (path == "pairwise") "ovo" else path),
        if (path == "pairwise") list(model_training = "pairwise"),
        if (path == "regression") list(m = 1)))
      state <- workflows::extract_fit_engine(generics::fit(model, data))
      output <- capture.output(print(state))
      multiple <- if (path == "pairwise") 3 else 1
      stopifnot(any(grepl(paste("Base fits:", multiple * expected$fits[i]), output, fixed = TRUE)),
        any(grepl(paste("Models per nonempty prediction:", multiple * expected$retained[i]), output, fixed = TRUE)),
        any(grepl(paste("Stored calibration rows:", (if (path == "pairwise") 2 else 1) * expected$calibration[i]), output, fixed = TRUE)))
    }
    # A frozen fit adopts one complete workflow; only its calibrator learns.
    if (path == "pairwise") next
    frozen <- generics::fit(base, data[1:60, ])
    adopted <- do.call(calibrate_workflow, c(list(frozen, data[61:120, ], "truth"),
      if (path %in% c("ovo", "ovr")) list(classification = path)))
    output <- capture.output(shown <- withVisible(print(adopted)))
    stopifnot(identical(output, c("Venn-Abers inductive calibrated predictor",
      "Frozen workflow adopted unchanged", "Base fits: 0 ", "Models per nonempty prediction: 1 ",
      "Calibration rows: 60 (new data; only the calibrator learned) ", "  Inductive1 -> Inductive1 ")),
      identical(shown, list(value = adopted, visible = FALSE)))
  }
}

  testthat::expect_true(TRUE)
})

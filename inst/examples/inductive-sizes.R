# Independent Inductive training and calibration sizes, as in sklearn's
# train_size/test_size. Rows in neither set form the unassigned set U: the
# matched source model and the Venn-Abers calibrator never see them, while a
# full-model refit (Full1) learns on all local rows, including C and U.
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

if (all(vapply(c("parsnip", "workflows", "rsample", "withr"),
               requireNamespace, logical(1), quietly = TRUE))) {
  data <- data.frame(x = sin(1:120) + rep(c(1, 0, 0, 1), 30),
                     truth = factor(rep(c("yes", "no", "no", "yes"), 30)))
  development <- data[1:100, ]
  assessment <- data[101:120, ]
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))

  # 100 local rows: 50 train the source model, 20 calibrate, 30 stay unassigned.
  # Whole numbers are exact row counts; a value in (0, 1) is a proportion
  # (training rounds down, calibration rounds up). Classification draws are
  # stratified by outcome with exact totals unless strata = FALSE.
  fits <- list()
  for (mode in c("matched", "full")) {
    fits[[mode]] <- generics::fit(calibrated_workflow(base, "truth", seed = 42,
      training_size = 50, calibration_size = 20, prediction_model = mode), development)
    state <- workflows::extract_fit_engine(fits[[mode]])
    print(state)
    roles <- state$folds$Inductive1
    assert_conditions(length(roles$analysis) == 50L, length(roles$assessment) == 20L,
      length(roles$unassigned) == 30L, identical(state$calibration$Inductive1$row, roles$assessment))
    stopifnot(identical(sort(c(roles$analysis, roles$assessment, roles$unassigned)), 1:100),
      identical(unname(state$sources), if (mode == "full") "Full1" else "Inductive1"))
  }
  # Both modes calibrate on the same 20 rows; only the prediction model differs.
  stopifnot(identical(workflows::extract_fit_engine(fits$matched)$folds,
                      workflows::extract_fit_engine(fits$full)$folds))
  for (mode in names(fits)) {
    p <- predict(fits[[mode]], assessment, type = "prob")
    probability <- as.matrix(p)
    colnames(probability) <- sub("^\\.pred_", "", colnames(probability))
    print(c(mode = mode, brier = round(brier_score(assessment$truth, probability), 4)))
  }

  # A single size takes the complement, so U is empty: 70 training rows, 30 C.
  complement <- workflows::extract_fit_engine(generics::fit(calibrated_workflow(base, "truth",
    seed = 42, training_size = .7), development))
  stopifnot(lengths(complement$folds$Inductive1[c("analysis", "assessment", "unassigned")]) == c(70L, 30L, 0L))

  # Sizes are checked strictly against the local rows and never shrunk.
  message <- tryCatch(generics::fit(calibrated_workflow(base, "truth", seed = 42,
    training_size = 80, calibration_size = 30), development), error = conditionMessage)
  stopifnot(grepl("request 110 rows, but only the 100 local rows are available", message, fixed = TRUE))

  # Pairwise multiclass resolves the sizes on each class pair's rows, with U
  # per pair. Classes of 40, 30 and 20 rows give pairs of 70, 60 and 50 rows:
  # half of each pair trains, 10 rows calibrate and the rest stay unassigned.
  three <- data.frame(x = sin(1:90) + rep(c(0, .5, 1), c(40, 30, 20)),
                      truth = factor(rep(c("a", "b", "c"), c(40, 30, 20))))
  pairwise <- workflows::extract_fit_engine(generics::fit(calibrated_workflow(base, "truth",
    classification = "ovo", model_training = "pairwise", seed = 42,
    training_size = .5, calibration_size = 10), three))
  print(pairwise)
  stopifnot(identical(pairwise$parameters$effective$unassigned_size,
                      c(pair_1 = 25L, pair_2 = 20L, pair_3 = 15L)))
  # Searches whose candidates are calibrated workflows (calibration or joint
  # model tuning) resolve the sizes again on resample rows and therefore
  # accept proportions only; two-stage model searches also take row counts.
}

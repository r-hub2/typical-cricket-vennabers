# Fixed engine options and model-internal validation. Options set with
# parsnip::set_engine() travel with the base workflow, so every local fit
# (source models, class pairs, Full1 and postprocessing B) uses them on its own
# rows only. calibrated_workflow() and fit() take no extra fit arguments.
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

if (all(vapply(c("parsnip", "workflows", "rsample", "withr", "rpart"),
               requireNamespace, logical(1), quietly = TRUE))) {
  data <- data.frame(x = sin(1:200) + rep(c(1, 0, 0, 1), 50), z = cos(1:200),
                     truth = factor(rep(c("yes", "no", "no", "yes"), 50)))
  development <- data[1:160, ]
  assessment <- data[161:200, ]
  # xval = 0 skips rpart's internal cross-validation; maxdepth is a fixed
  # main argument of the model specification.
  tree <- parsnip::decision_tree(mode = "classification", tree_depth = 3L) |>
    parsnip::set_engine("rpart", xval = 0L, maxsurrogate = 0L)
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ .) |>
    workflows::add_model(tree)

  fit <- generics::fit(calibrated_workflow(base, "truth", method = "cross", v = 3L,
    prediction_model = "full", seed = 42), development)
  state <- workflows::extract_fit_engine(fit)
  print(state)
  # Full mode keeps one fresh model on all 160 rows; it carries the options.
  engine <- workflows::extract_fit_engine(state$models$Full1)
  assert_conditions(engine$control$xval == 0L, engine$control$maxsurrogate == 0L,
    engine$control$maxdepth == 3L, length(engine$where) == 160L)
  print(head(predict(fit, assessment, type = "prob")))

  # Fit options belong in set_engine(); extra arguments are rejected.
  message <- tryCatch(calibrated_workflow(base, "truth", xval = 0L), error = conditionMessage)
  stopifnot(grepl("unsupported argument: xval", message, fixed = TRUE))

  # Native early stopping works the same way. Each local fit splits its
  # validation share from its own prepared rows; calibration rows never reach
  # it. Not run here because xgboost is not a vennabers dependency:
  if (FALSE) {
    spec <- parsnip::boost_tree(mode = "classification", trees = 100L, stop_iter = 10L) |>
      parsnip::set_engine("xgboost", validation = 0.2)
    boosted <- workflows::workflow() |>
      workflows::add_formula(truth ~ .) |>
      workflows::add_model(spec)
    fitted <- generics::fit(calibrated_workflow(boosted, "truth", seed = 42), development)
  }
}

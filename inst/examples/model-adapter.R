# model_adapter(): your own fit(x, y) and predict(object, new_data) functions in
# a calibrated workflow, without registering a parsnip engine. The workflow keeps
# preprocessing, splits and calibration: every local fit trains a fresh model on
# its own prepared rows, and the calibrator never sees those training rows.
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

examples <- list()
if (all(vapply(c("parsnip", "workflows", "recipes", "rsample", "withr", "nnet"),
               requireNamespace, logical(1), quietly = TRUE))) {
  # Each factory binds its helpers and fixed settings locally, so a saved fit
  # carries the actual functions into a fresh R process. Package functions use
  # pkg::fun; nothing is looked up in the global environment later.
  glm_adapter <- function() model_adapter("classification",
    fit = function(x, y) list(model = stats::glm(y ~ ., data = cbind(x, y = y), family = stats::binomial()),
                              classes = levels(y)),
    predict = function(object, new_data) {
      # A binomial GLM models the second level. Return both named columns.
      p <- stats::predict(object$model, new_data, type = "response")
      matrix(c(1 - p, p), ncol = 2L, dimnames = list(NULL, object$classes))
    })
  multinom_adapter <- function(decay = 0) {
    force(decay)  # bind the value now, not a promise to the caller's variable
    model_adapter("classification",
      fit = function(x, y) nnet::multinom(y ~ ., data = cbind(x, y = y), decay = decay, trace = FALSE),
      predict = function(object, new_data) {
        # multinom drops a single row to a vector; keep one column per class.
        p <- stats::predict(object, new_data, type = "probs")
        matrix(p, nrow(new_data), length(object$lev), dimnames = list(NULL, object$lev))
      })
  }
  lm_adapter <- function() model_adapter("regression",
    fit = function(x, y) stats::lm(y ~ ., data = cbind(x, y = y)),
    predict = function(object, new_data) unname(stats::predict(object, new_data)))

  base <- function(spec, data) workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ ., data) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(spec)

  # Binary: the same GLM as the native engine, so both calibrate identically.
  binary <- data.frame(x = sin(1:200) + rep(c(1, 0, 0, 1), 50), z = cos(1:200),
                       truth = factor(rep(c("yes", "no", "no", "yes"), 50), levels = c("yes", "no")))
  development <- binary[1:160, ]
  assessment <- binary[161:200, c("x", "z")]
  adapter <- calibrated_workflow(base(glm_adapter(), development), "truth", seed = 42)
  native <- calibrated_workflow(base(parsnip::logistic_reg(engine = "glm"), development), "truth", seed = 42)
  fit <- generics::fit(adapter, development)
  p <- predict(fit, assessment, type = "prob")
  print(head(p))
  stopifnot(all.equal(p, predict(generics::fit(native, development), assessment, type = "prob")))

  # Shared multiclass: the callback learns and predicts all local classes.
  classes <- factor(rep(c("a", "b", "c"), 70))
  multiclass <- data.frame(x = sin(1:210) + as.integer(classes), z = cos(1:210), truth = classes)
  shared <- calibrated_workflow(base(multinom_adapter(decay = .1), multiclass[1:180, ]), "truth",
                                classification = "ovr", method = "cross", v = 3, seed = 42)
  shared_fit <- generics::fit(shared, multiclass[1:180, ])
  print(head(predict(shared_fit, multiclass[181:210, c("x", "z")], type = "prob")))

  # Regression: one finite number per row; Full1 refits on all local rows.
  regression <- data.frame(x = sin(1:150), z = cos(1:150))
  regression$truth <- 10 + 3 * regression$x + regression$z^2
  numeric <- calibrated_workflow(base(lm_adapter(), regression[1:120, ]), "truth",
                                 prediction_model = "full", seed = 42)
  numeric_fit <- generics::fit(numeric, regression[1:120, ])
  print(head(venn_abers_details(numeric_fit, regression[121:150, c("x", "z")])$details))

  # The stored state is the object your fit() returned, inside the workflow.
  source_model <- workflows::extract_fit_engine(workflows::extract_fit_engine(fit)$models$Inductive1)
  assert_conditions(inherits(source_model, "vennabers_adapter_fit"),
    inherits(source_model$state$model, "glm"),
    length(source_model$state$model$y) == length(workflows::extract_fit_engine(fit)$folds$Inductive1$analysis))

  # Class pairs: the binary GLM callback carries pairwise OvO. Each call sees
  # one pair's rows, with its two original levels. Overlapping classes keep
  # every pair's GLM away from complete separation.
  overlap <- transform(multiclass, x = x + .5 * sin(1:210))
  pairwise <- calibrated_workflow(base(glm_adapter(), overlap[1:180, ]), "truth",
                                  classification = "ovo", model_training = "pairwise", seed = 42)
  pairwise_fit <- generics::fit(pairwise, overlap[1:180, ])
  print(head(predict(pairwise_fit, overlap[181:210, c("x", "z")], type = "prob")))

  # Postprocessing: B fits the callback, disjoint P learns Sigmoid on its
  # probabilities, C learns Venn-Abers. The glm_link input needs the native GLM.
  sigmoid <- list(method = "sigmoid", sigmoid_input = "probability", strategy = "holdout", prop = .3)
  post_fit <- generics::fit(calibrated_workflow(base(glm_adapter(), development), "truth", seed = 42,
                                                postprocessing = sigmoid), development)
  stopifnot(all.equal(predict(post_fit, assessment, type = "prob"),
    predict(generics::fit(calibrated_workflow(base(parsnip::logistic_reg(engine = "glm"), development), "truth",
                                              seed = 42, postprocessing = sigmoid), development), assessment, type = "prob")))

  examples <- list(
    binary = list(data = binary, model = adapter, fit = fit, new = assessment),
    multiclass = list(data = multiclass, model = shared, fit = shared_fit, new = multiclass[181:210, c("x", "z")]),
    regression = list(data = regression, model = numeric, fit = numeric_fit, new = regression[121:150, c("x", "z")]))

  # Searches tune calibration and recipe parameters. Callback settings are fixed
  # in the functions, so the adapter itself declares no tune() parameters.
  if (all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
    poly <- workflows::workflow() |>
      workflows::add_recipe(recipes::recipe(truth ~ x + z, development) |>
        recipes::step_poly(x, degree = tune::tune("degree"))) |>
      workflows::add_model(glm_adapter())
    folds <- function(data) rsample::vfold_cv(data, v = 2)
    searched <- calibrated_workflow(poly, "truth", loss = tune::tune("va_loss"), seed = 42,
      tuning = list(strategy = "two_stage", selection_prop = .5,
        model = list(method = "grid", grid = data.frame(degree = 1:2), resamples = folds),
        calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")), resamples = folds)))
    searched_fit <- generics::fit(searched, development)
    print(venn_abers_tuning(searched_fit)$calibration$selected)
    examples$search <- list(data = binary, model = searched, fit = searched_fit, new = assessment)
  }
}

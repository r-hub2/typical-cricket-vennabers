# Reference-style row partitions; no downloads or new splitter API.
# Inspired by the binary, multiclass and regression applications of the
# Venn-Abers papers. These small synthetic examples illustrate roles, not
# reproductions of paper benchmarks.
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

# Calibration gets ceiling(n * calibration_prop) rows; training gets the
# preceding n - ceiling(n * calibration_prop). Never sort, shuffle or stratify.
ordered_holdout <- function(data, calibration_prop = .25) {
  assert_conditions(is.numeric(calibration_prop), length(calibration_prop) == 1L,
            is.finite(calibration_prop), calibration_prop > 0, calibration_prop < 1)
  n <- nrow(data)
  cut <- n - ceiling(n * calibration_prop)
  stopifnot(cut >= 1L, cut < n)
  rsample::make_splits(list(analysis = seq_len(cut),
                           assessment = seq.int(cut + 1L, n)), data)
}

if (all(vapply(c("parsnip", "workflows", "rsample", "recipes", "nnet", "withr"),
               requireNamespace, logical(1), quietly = TRUE))) {
  # A small known membership, including rounding of 11 * .25 to 3 held rows.
  small <- data.frame(x = 11:1)
  stopifnot(identical(as.integer(ordered_holdout(small), data = "analysis"), 1:8),
            identical(as.integer(ordered_holdout(small), data = "assessment"), 9:11))

  # Non-shuffled, stratified three-fold reference assignment. Class order is
  # first appearance (b, a, c); class counts are 3, 4, 4. Remainders are spread
  # across folds (sizes 4, 4, 3). Positions within each class keep their order.
  small$truth <- factor(c("b", "b", "a", "a", "a", "c", "c", "c", "c", "b", "a"))
  stratified_held <- list(c(1L, 3L, 4L, 6L), c(2L, 5L, 7L, 8L), 9:11)
  stratified <- rsample::manual_rset(lapply(stratified_held, function(h) {
    rsample::make_splits(list(analysis = setdiff(1:11, h), assessment = h), small)
  }), c("First", "Second", "Third"))
  stopifnot(identical(sort(unlist(stratified_held)), 1:11))

  # All eight trainable paths: binary, shared OvO/OvR, and moderated regression,
  # each with an ordered holdout and ordinary non-shuffled Cross folds.
  reference_results <- list()
  for (kind in c("binary", "ovo", "ovr", "regression")) {
    data <- data.frame(x = sin(1:72), z = cos(1:72))
    data$truth <- if (kind == "regression") 10 + 2 * data$x + data$z else
      factor(rep(if (kind == "binary") c("working", "fault") else c("soil", "road", "water"),
                 length.out = nrow(data)))
    # Final assessment is excluded from model learning AND calibration.
    development <- data[1:60, ]
    assessment <- data[61:72, ]
    spec <- switch(kind, regression = parsnip::linear_reg(engine = "lm"),
                   binary = parsnip::logistic_reg(engine = "glm"),
                   parsnip::multinom_reg(engine = "nnet", penalty = .1))
    base <- workflows::workflow() |>
      workflows::add_recipe(recipes::recipe(truth ~ ., development) |>
                              recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(spec)
    # These are literal local positions, not forward validation. In particular,
    # the first Cross model trains on later rows 21:60.
    cross <- rsample::manual_rset(lapply(list(1:20, 21:40, 41:60), function(h) {
      rsample::make_splits(list(analysis = setdiff(1:60, h), assessment = h), development)
    }), c("Early", "Middle", "Late"))
    for (method in c("inductive", "cross")) {
      plan <- if (method == "inductive") ordered_holdout(development) else cross
      resolve <- local_calibration_splits(list(plan), names(development))
      options <- if (kind == "regression") list(m = 1L) else list(classification = kind)
      model <- do.call(calibrated_workflow, c(list(workflow = base, outcome = "truth",
        method = method, splits = resolve, seed = 42), options))
      fit <- generics::fit(model, development)
      p <- predict(fit, assessment)
      score <- if (kind == "regression") {
        c(rmse = sqrt(mean((p$.pred - assessment$truth)^2)))
      } else {
        p <- predict(fit, assessment, type = "prob")
        probability <- as.matrix(p)
        colnames(probability) <- sub("^\\.pred_", "", colnames(probability))
        c(brier = brier_score(assessment$truth, probability),
          log_loss = log_loss(assessment$truth, probability))
      }
      state <- workflows::extract_fit_engine(fit)
      reference_results[[paste(kind, method, sep = "_")]] <- list(
        model = model, fit = fit, development = development, assessment = assessment,
        prediction = p, score = score, folds = state$folds)
      stopifnot(all(is.finite(score)), identical(names(state$folds), names(state$models)))
    }
  }
  # For mixed plans, use splits = NULL with seed, or a seeded initial_split()/
  # vfold_cv() callback. Automatic classification stratifies by outcome unless
  # strata = FALSE; regression does not stratify. Same seed requires the same
  # R RNGkind, package/engine versions and model settings for reproducibility.
  # R/Python seeds alone do not prove equal assignments or model outputs.
  #
  # Machine-failure detection is classification, not a failure-time forecast.
  # Temperature forecast correction is regression; its subject does not make
  # these partitions time-series validation. No rolling window, horizon, embargo,
  # grouping or online update is introduced. Cross uses exact fold complements.
  # Scores may get worse. VA bounds are not individual confidence/prediction
  # intervals; neither dependent data nor Cross aggregates gain universal coverage.
}

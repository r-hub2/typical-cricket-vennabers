# Separate pair/fold workflows with matching Cross calibration at fixed settings.
source(system.file("examples", "assertions.R", package = "vennabers"))
packages <- c("vennabers", "parsnip", "workflows", "recipes", "rsample")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  library(vennabers)
  data <- data.frame(x = sin(1:120), z = cos(1:120),
    truth = factor(rep(c("class", "a b", "a.b"), 40), levels = c("class", "a b", "a.b")))
  names(data) <- c("a b", "a.b", ".pred_class")
  schema <- venn_abers_schema(data, ".pred_class", classification = "ovo")
  encoded <- venn_abers_encode(data, schema)
  outer_split <- rsample::make_splits(list(analysis = 1:90, assessment = 91:120), encoded)
  development <- rsample::analysis(outer_split)
  assessment <- rsample::assessment(outer_split)
  formula <- reformulate(names(encoded)[1:2], schema$encoded_outcome)
  base <- workflows::workflow() |> workflows::add_recipe(
    recipes::recipe(formula, development) |> recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  pairs <- combn(levels(development[[schema$encoded_outcome]]), 2, simplify = FALSE)
  # A local function also works inside fit_resamples()/last_fit(). Here the
  # class pair selects 3, 4 or 5 folds. Fold1 in different pairs is unrelated.
  pair_plan <- function(data) {
    j <- which(vapply(pairs, identical, logical(1), y = levels(data[[schema$encoded_outcome]])))
    v <- 2L + j
    held <- split(seq_len(nrow(data)), cut(seq_len(nrow(data)), v, labels = FALSE))
    rsample::manual_rset(lapply(held, function(h) rsample::make_splits(
      list(analysis = setdiff(seq_len(nrow(data)), h), assessment = h), data)), paste0("Fold", seq_len(v)))
  }
  plans <- setNames(lapply(pairs, function(pair) {
    local <- development[development[[schema$encoded_outcome]] %in% pair, , drop = FALSE]
    local[[schema$encoded_outcome]] <- factor(local[[schema$encoded_outcome]], levels = pair)
    pair_plan(local)
  }), paste0("pair_", seq_along(pairs)))
  explicit <- calibrated_workflow(base, schema$encoded_outcome, schema = schema,
    model_training = "pairwise", method = "cross", splits = plans, seed = 42)
  fit <- generics::fit(explicit, development)
  state <- workflows::extract_fit_engine(fit)
  print(state)
  pair_table <- data.frame(pair_id = names(state$pair_models),
    class_1 = vapply(state$pair_models, function(s) s$class_levels[1L], ""),
    class_2 = vapply(state$pair_models, function(s) s$class_levels[2L], ""))
  print(pair_table)
  native <- venn_abers_details(fit, assessment)
  decoded <- venn_abers_decode(native, schema)
  print(head(native$details))
  print(head(decoded$probability))
  print(head(decoded$details))
  model <- calibrated_workflow(base, schema$encoded_outcome, schema = schema,
    model_training = "pairwise", method = "cross", splits = pair_plan, seed = 42)
  stopifnot(identical(native, venn_abers_details(generics::fit(model, development), assessment)))
  # Reconstruct each pair from its binary fold endpoints, then couple globally.
  # There is no global fold matrix and no averaging of coupled fold predictions.
  for (loss in c("log", "brier")) {
    details <- venn_abers_details(fit, assessment, loss = loss)
    q <- vapply(names(plans), function(id) {
      d <- details$details[details$details$pair_id == id, ]
      lower <- matrix(d$lower, nrow(assessment))
      upper <- matrix(d$upper, nrow(assessment))
      if (loss == "brier") rowMeans(matrix(d$probability, nrow(assessment))) else {
        a <- apply(upper, 1, prod)^(1 / ncol(upper))
        b <- apply(1 - lower, 1, prod)^(1 / ncol(lower))
        a / (a + b)
      }
    }, numeric(nrow(assessment)))
    p <- vapply(state$class_levels, function(cl) {
      relevant <- which(pair_table$class_1 == cl | pair_table$class_2 == cl)
      odds <- lapply(relevant, function(j) {
        probability <- if (pair_table$class_2[j] == cl) q[, j] else 1 - q[, j]
        (1 - probability) / probability
      })
      1 / (1 + Reduce(`+`, odds))
    }, numeric(nrow(assessment)))
    stopifnot(all(abs(p / rowSums(p) - as.matrix(details$prediction)) <= 1e-10))
  }
  cost <- sum(vapply(state$pair_models, function(p) length(p$models), 1L))
  stopifnot(cost == 12L)
  cat("Base fits per complete local fit:", cost, "(3 + 4 + 5)\n")
  if (all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
    scoring <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
    outer <- withr::with_seed(41, rsample::vfold_cv(development, v = 2, repeats = 2))
    result <- tune::fit_resamples(model, outer, metrics = scoring,
      control = tune::control_resamples(save_pred = TRUE, extract = identity, allow_par = FALSE))
    checker <- new.env()
    sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), checker)
    checked <- checker$check_encoded_evaluation(result, outer, schema)
    print(checked$metrics)
    final <- tune::last_fit(model, outer_split, metrics = scoring)
    predictions <- tune::collect_predictions(final)
    metrics <- tune::collect_metrics(final)
    assert_conditions(nrow(tune::collect_notes(final)) == 0L, identical(predictions$.row, 91:120),
      !anyDuplicated(predictions$.row), all(predictions$.config == "pre0_mod0_post0"),
      identical(predictions[[schema$encoded_outcome]], assessment[[schema$encoded_outcome]]),
      nrow(metrics) == 2L, all(is.finite(metrics$.estimate)),
      all(abs(as.matrix(predictions[names(native$prediction)]) - as.matrix(native$prediction)) <= 1e-10))
    cat("Outer resampling base fits:", nrow(outer) * cost, "; last_fit base fits:", cost, "\n")
  }
  # Matched: all 12 complete workflows/calibrators remain stored. Nonempty
  # predict calls each model exactly once and never learns; empty inputs call
  # no models.
  # Outer assessment outcomes enter only evaluation, never the local pair plans.
  # One-sided calibration folds are valid when their complements train and return
  # valid two-class probabilities. Every row is assessed once within each pair.
  # Full transfer keeps the held-out calibrators but fits one extra workflow
  # on each complete local pair, including its calibration rows.
  full_model <- calibrated_workflow(base, schema$encoded_outcome, schema = schema,
    model_training = "pairwise", method = "cross", prediction_model = "full",
    splits = pair_plan, seed = 42)
  full_fit <- generics::fit(full_model, development)
  full_state <- workflows::extract_fit_engine(full_fit)
  for (id in names(state$pair_models)) stopifnot(
    identical(state$pair_models[[id]]$calibration, full_state$pair_models[[id]]$calibration),
    identical(state$pair_models[[id]]$calibrator, full_state$pair_models[[id]]$calibrator),
    identical(names(full_state$pair_models[[id]]$models), "Full1"))
  full_native <- venn_abers_details(full_fit, assessment)
  full_decoded <- venn_abers_decode(full_native, schema)
  print(head(full_decoded$details))
  for (mode in c("matched", "full")) {
    probability <- if (mode == "full") full_decoded$probability else decoded$probability
    truth <- venn_abers_decode(assessment[schema$encoded_outcome], schema)$truth
    cat(mode, "Brier:", brier_score(truth, probability),
        "log loss:", log_loss(truth, probability), "\n")
  }
  # Fifteen base fits (3 + 4 + 5 sources, three Full1); only the Full1 workflows
  # remain, each predicting once per nonempty batch. The printed slot map sends
  # every held-out calibrator to its pair's Full1.
  print(full_state)
  # Each pair source is named pair_N/Full1; pair_id + fold identifies its
  # calibrator. Source workflows are discarded after all full fits succeed.
  # No in-sample recalibration or pooling occurs. Outer assessment stays excluded.
  # Matching-model guarantees do not automatically transfer to the new model;
  # neither coverage nor improved quality is promised. Settings are fixed here;
  # pairwise-grid.R and pairwise-bayes.R show automatic local selection.
}

# Separate binary training with matching Inductive calibration, at fixed settings.
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
  # Receives only this pair's current local raw data, with exactly two levels.
  # This is local row order, not a shared global time cutoff. Both classes must
  # occur in training and calibration. Mixed holdouts are also supported.
  pair_plan <- function(data) {
    n <- nrow(data)
    held <- seq.int(n - ceiling(n / 4) + 1L, n)
    rsample::make_splits(list(analysis = setdiff(seq_len(n), held), assessment = held), data)
  }
  model <- calibrated_workflow(base, schema$encoded_outcome, schema = schema,
    model_training = "pairwise", splits = pair_plan, seed = 42)
  fit <- generics::fit(model, development)
  state <- workflows::extract_fit_engine(fit)
  print(state)
  pair_table <- data.frame(pair_id = names(state$pair_models),
    class_1 = vapply(state$pair_models, function(s) s$class_levels[1L], ""),
    class_2 = vapply(state$pair_models, function(s) s$class_levels[2L], ""))
  print(pair_table)
  native <- venn_abers_details(fit, assessment)
  decoded <- venn_abers_decode(native, schema)
  print(head(decoded$probability))
  print(head(decoded$details))
  # Direct fits may instead supply exactly pair_1, pair_2, ... as a plan list.
  # The table above maps these stable IDs to classes; class names are never joined.
  plans <- lapply(state$pair_models, function(component) {
    local <- development[component$rows, , drop = FALSE]
    local[[schema$encoded_outcome]] <- factor(local[[schema$encoded_outcome]], levels = component$class_levels)
    pair_plan(local)
  })
  explicit <- generics::fit(calibrated_workflow(base, schema$encoded_outcome, schema = schema,
    model_training = "pairwise", splits = plans, seed = 42), development)
  stopifnot(identical(native, venn_abers_details(explicit, assessment)))
  # Independent complete binary workflows: preprocessing sees only each pair's
  # training rows, and held-out model probabilities fit the existing binary core.
  q <- matrix(0, nrow(assessment), nrow(pair_table))
  for (j in seq_along(plans)) {
    training <- rsample::analysis(plans[[j]])
    held <- rsample::assessment(plans[[j]])
    reference <- workflows::workflow() |> workflows::add_recipe(
      recipes::recipe(formula, training) |> recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(parsnip::logistic_reg(engine = "glm"))
    binary <- generics::fit(reference, training)
    event <- pair_table$class_2[j]
    column <- paste0(".pred_", event)
    cal <- venn_abers(held[[schema$encoded_outcome]] == event, predict(binary, held, type = "prob")[[column]])
    expected <- predict(cal, predict(binary, assessment, type = "prob")[[column]])
    parts <- native$details[native$details$pair_id == pair_table$pair_id[j], c("lower", "upper", "probability")]
    stopifnot(all(abs(as.matrix(parts) - as.matrix(expected)) <= 1e-10))
    q[, j] <- expected$probability
  }
  # Reciprocal OvO coupling, written independently as sums of binary odds.
  p <- vapply(state$class_levels, function(cl) {
    relevant <- which(pair_table$class_1 == cl | pair_table$class_2 == cl)
    odds <- lapply(relevant, function(j) {
      probability <- if (pair_table$class_2[j] == cl) q[, j] else 1 - q[, j]
      (1 - probability) / probability
    })
    1 / (1 + Reduce(`+`, odds))
  }, numeric(nrow(assessment)))
  stopifnot(all(abs(p / rowSums(p) - as.matrix(native$prediction)) <= 1e-10))
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
  }
  # Matched: P = K*(K-1)/2 pairs cost P base fits and P predictions per nonempty
  # batch. No extra shared model or OvR ensemble is trained. All P recipes and
  # calibrators remain available; predict never learns. A row may belong to
  # several pairs.
  # Outer assessment labels never reach these local fits or calibration plans.
  # Prepared outer plans can use local_calibration_splits() over all local pairs;
  # a development-wide plan cannot be reused on a smaller outer analysis domain.
  # Matching-model assumptions still apply; no universal accuracy gain is promised.
  # Full transfer keeps the held-out calibrators but fits one extra workflow
  # on each complete local pair, including its calibration rows.
  full_model <- calibrated_workflow(base, schema$encoded_outcome, schema = schema,
    model_training = "pairwise", method = "inductive", prediction_model = "full",
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
  # Six base fits (three sources, three Full1); only the Full1 workflows remain,
  # each predicting once per nonempty batch. The printed slot map sends every
  # held-out calibrator to its pair's Full1.
  print(full_state)
  # Each pair source is named pair_N/Full1; pair_id + fold identifies its
  # calibrator. Source workflows are discarded after all full fits succeed.
  # No in-sample recalibration or pooling occurs. Outer assessment stays excluded.
  # Matching-model guarantees do not automatically transfer to the new model;
  # neither coverage nor improved quality is promised. Settings are fixed here;
  # pairwise-grid.R and pairwise-bayes.R show automatic local selection.
}

# Class pairs (pair_N) and inner B/P model/postprocessor pairs are different:
# each class pair owns its VA folds, and each VA fold retains its own B/P paths.
source(system.file("examples", "assertions.R", package = "vennabers"))
packages <- c("vennabers", "workflows", "recipes", "rsample", "tune", "yardstick", "tailor", "probably")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  library(vennabers)
  raw <- data.frame(x = sin(1:360), z = cos(1:360),
    truth = factor(rep(c("class", "a b", "class", "a.b", "a b", "class"), 60),
      levels = c("class", "a b", "a.b")))
  names(raw) <- c("input value", "a.b", ".pred_class")
  schema <- venn_abers_schema(raw, ".pred_class", classification = "ovo")
  d <- venn_abers_encode(raw, schema); outcome <- schema$encoded_outcome
  predictors <- setdiff(names(d), outcome)
  base <- workflows::workflow() |> workflows::add_recipe(
    recipes::recipe(reformulate(predictors, outcome), d) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  outer <- withr::with_seed(41, rlang::inject(rsample::vfold_cv(d, v = 2, repeats = 2, strata = !!outcome, pool = 0)))
  split <- rsample::make_splits(list(analysis = 1:270, assessment = 271:360), d)
  metrics <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
  checker <- new.env()
  sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), checker)
  cases <- list()
  for (kind in c("tailor_logistic", "probability", "glm_link")) for (strategy in c("holdout", "ensemble"))
    for (method in c("inductive", "cross")) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = strategy)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    if (strategy == "holdout") post$prop <- .3 else post$v <- 2L
    args <- list(base, outcome, schema = schema, model_training = "pairwise", method = method,
      seed = 42, postprocessing = post)
    if (method == "cross") args$v <- 2L
    model <- do.call(calibrated_workflow, args)
    result <- tune::fit_resamples(model, outer, metrics = metrics,
      control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, allow_par = FALSE))
    checker$check_encoded_evaluation(result, outer, schema)
    states <- tune::collect_extracts(result)$.extracts
    for (i in seq_along(states)) {
      local <- rsample::analysis(outer$splits[[i]])
      analysis_rows <- as.integer(outer$splits[[i]], "analysis")
      A <- as.integer(outer$splits[[i]], "assessment")
      for (component in states[[i]]$pair_models) for (diagnostic in component$postprocessing) {
        members <- if (strategy == "ensemble") diagnostic$inner_folds else list(diagnostic)
        for (p in members) assert_conditions(
          all(as.character(local[[outcome]][component$rows]) %in% component$class_levels),
          identical(p$pair_id, component$pair_id), identical(p$classes, component$class_levels),
          !any(analysis_rows[component$rows[c(p$B, p$P, p$C)]] %in% A),
          !length(intersect(p$B, p$P)), !length(intersect(p$Q, p$C)))
      }
    }
    last <- tune::last_fit(model, split, metrics = metrics)
    native <- tune::collect_predictions(last); scores <- tune::collect_metrics(last)
    decoded <- venn_abers_decode(native, schema)
    values <- c(brier_class = yardstick::brier_class_vec(decoded$truth, decoded$probability),
      mn_log_loss = yardstick::mn_log_loss_vec(decoded$truth, decoded$probability))
    fit <- tune::extract_workflow(last)
    new <- d[271:360, predictors]
    assert_conditions(nrow(tune::collect_notes(last)) == 0L, identical(native$.row, 271:360),
      !anyDuplicated(native$.row), all(native$.config == "pre0_mod0_post0"),
      identical(native[[outcome]], d[[outcome]][271:360]), nrow(scores) == 2L,
      setequal(scores$.metric, names(values)), all(is.finite(scores$.estimate)),
      all(abs(scores$.estimate - values[scores$.metric]) <= 1e-10),
      all(abs(as.matrix(native[paste0(".pred_", schema$classes)]) - as.matrix(predict(fit, new, type = "prob"))) <= 1e-10))
    state <- workflows::extract_fit_engine(fit)
    details <- venn_abers_details(fit, new)
    cost <- sum(vapply(state$pair_models, function(p) sum(vapply(p$postprocessing, `[[`, 1L, "observed_base_fits")), 1L))
    stopifnot(cost == 3L * (if (method == "cross") 2L else 1L) * (if (strategy == "ensemble") 2L else 1L),
      identical(details$prediction, predict(fit, new, type = "prob")))
    # A truth and predictor changes affect assessment only, never any class pair.
    for (column in c(predictors[1], outcome)) {
      changed <- d
      if (column == outcome) changed[[column]][271:360] <- schema$classes[1] else changed[[column]][271:360] <- 1000
      other <- tune::last_fit(model, rsample::make_splits(list(analysis = 1:270, assessment = 271:360), changed), metrics = metrics)
      other_fit <- tune::extract_workflow(other)
      stopifnot(identical(venn_abers_details(fit, new), venn_abers_details(other_fit, new)))
    }
    cases[[paste(kind, strategy, method)]] <- list(fit = fit, new = new, schema = schema)
    cat(kind, strategy, method, ": 3 class pairs;", cost, "retained model/postprocessor pairs.\n")
  }
}

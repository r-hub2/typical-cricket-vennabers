# Installed example: shared OvO/OvR retains full workflows trained on B, with
# one Tailor step learned on disjoint P, before VA calibration on held-out C.
# Choose logistic explicitly for two classes, multinomial for three or more.
# Multinomial P uses probably/nnet's default 100 iterations, independently of
# the B engine. Weakly informative probabilities may fail with fold/pair/B/P
# context, without retry or method switch; postprocessing has no maxit option.
packages <- c("vennabers", "workflows", "tailor", "probably", "recipes", "nnet", "tune", "yardstick")
cases <- list()
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  source(system.file("examples", "assertions.R", package = "vennabers"), local = TRUE)
  audit <- new.env()
  sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), audit)
  for (k in 2:4) for (method in c("inductive", "cross")) for (classification in c("ovo", "ovr")) {
    original <- withr::with_seed(123, {
      truth <- factor(rep(c("z z", "a.a", "a a", ".pred_class")[seq_len(k)], length.out = 480),
        levels = c("z z", "a.a", "a a", ".pred_class")[seq_len(k)])
      data.frame(x = rnorm(480, sd = 2) + c(-1, 1, -1, 1)[truth],
        z = rnorm(480, sd = 2) + c(-1, -1, 1, 1)[truth],
        w = rnorm(480, sd = 2) + c(1, -1, -1, 1)[truth], truth = truth)
    })
    names(original) <- c("input value", "input.value", "other value", ".pred_class")
    schema <- vennabers::venn_abers_schema(original, ".pred_class", classification = classification)
    d <- vennabers::venn_abers_encode(original, schema)
    outcome <- schema$encoded_outcome
    base <- workflows::workflow() |>
      workflows::add_recipe(recipes::recipe(reformulate(setdiff(names(d), outcome), outcome), d) |>
        recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(if (k == 2L) parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
    # Both B and P need all k classes, P needs at least k+1 rows and the actual
    # Tailor design must have rank k. Failed convergence is an error, not fallback.
    post <- list(method = if (k == 2L) "tailor_logistic" else "tailor_multinomial", strategy = "holdout", prop = .4)
    spec <- vennabers::calibrated_workflow(base, outcome, schema = schema, method = method,
      seed = 24, postprocessing = post)
    split <- rsample::make_splits(list(analysis = 1:384, assessment = 385:480), d)
    metrics <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
    result <- tune::last_fit(spec, split, metrics = metrics)
    fit <- result$.workflow[[1]]
    native <- vennabers::venn_abers_details(fit, d[385:480, setdiff(names(d), outcome)])
    decoded <- vennabers::venn_abers_decode(native, schema)
    predictions <- tune::collect_predictions(result)
    scores <- tune::collect_metrics(result)
    assert_conditions(nrow(tune::collect_notes(result)) == 0L, identical(predictions$.row, 385:480),
      !anyDuplicated(predictions$.row), all(predictions$.config == "pre0_mod0_post0"),
      nrow(scores) == 2L, all(is.finite(scores$.estimate)),
      identical(as.matrix(predictions[names(native$prediction)]), as.matrix(native$prediction)),
      identical(colnames(decoded$probability), schema$class_levels),
      identical(decoded$postprocessing[[1]]$classes, schema$class_levels),
      identical(decoded$postprocessing[[1]]$columns, names(original)),
      all(native$details$postprocessing_method == post$method),
      abs(scores$.estimate[scores$.metric == "brier_class"] -
        vennabers::brier_score(original[[schema$outcome]][385:480], decoded$probability) / 2) <= 1e-10)
    if (k > 2L) {
      coefficients <- decoded$postprocessing[[1]]$coefficients
      assert_conditions(identical(rownames(coefficients), schema$class_levels[-1L]),
        all(colnames(coefficients)[-1L] %in% paste0(".pred_", schema$class_levels)))
    }
    outer <- withr::with_seed(74, rsample::vfold_cv(d, v = 2, repeats = 2))
    resampled <- tune::fit_resamples(spec, outer, metrics = metrics,
      control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, allow_par = FALSE))
    audit$check_encoded_evaluation(resampled, outer, schema)
    states <- tune::collect_extracts(resampled)$.extracts
    for (i in seq_along(states)) for (plan in states[[i]]$postprocessing) {
      Q <- as.integer(outer$splits[[i]], "analysis")
      A <- as.integer(outer$splits[[i]], "assessment")
      assert_conditions(!any(Q[c(plan$B, plan$P, plan$C)] %in% A),
        !length(intersect(plan$B, plan$P)), !length(intersect(plan$Q, plan$C)),
        plan$observed_base_fits == 1L, plan$observed_postprocessor_fits == 1L)
    }
    # Changing outer A never changes the learned state.
    changed <- d; changed[385:480, 1L] <- 1000
    changed[[outcome]][385:480] <- levels(changed[[outcome]])[1L]
    other <- tune::last_fit(spec, rsample::make_splits(list(analysis = 1:384, assessment = 385:480), changed), metrics = metrics)
    assert_conditions(identical(workflows::extract_fit_engine(fit)$postprocessing,
      workflows::extract_fit_engine(other$.workflow[[1]])$postprocessing),
      identical(predict(fit, d, type = "prob"), predict(other$.workflow[[1]], d, type = "prob")))
    cases[[paste(k, method, classification)]] <- list(fit = fit, schema = schema,
      new = d[setdiff(names(d), outcome)])
  }
}

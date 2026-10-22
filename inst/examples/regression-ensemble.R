# Each inner pair learns its complete recipe/engine on B_j and its numeric
# postprocessor on disjoint P_j. P folds partition Q; B_j is Q minus P_j.
# Average finished numeric pair predictions equally, then calibrate VA on C.
# CVAR averages outer fold points and bounds separately: these bounds are
# descriptive VA bounds, not confidence or prediction intervals.
packages <- c("vennabers", "workflows", "tailor", "probably", "recipes", "rsample", "tune", "yardstick")
cases <- list()
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  source(system.file("examples", "assertions.R", package = "vennabers"), local = TRUE)
  audit <- new.env()
  sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), audit)
  original <- data.frame("input value" = sin(1:120),
    "target y" = 10 + 3 * sin(1:120) + cos(1:120) + sin((1:120) / 4), check.names = FALSE)
  schema <- vennabers::venn_abers_schema(original, "target y")
  d <- vennabers::venn_abers_encode(original, schema)
  outcome <- schema$encoded_outcome
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(reformulate(setdiff(names(d), outcome), outcome), d) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  metrics <- yardstick::metric_set(yardstick::rmse, yardstick::mae)
  outer <- withr::with_seed(74, rsample::vfold_cv(d, v = 2, repeats = 2))
  split <- rsample::make_splits(list(analysis = 1:90, assessment = 91:120), d)
  for (post_method in c("tailor_linear", "isotonic_numeric")) for (method in c("inductive", "cross")) {
    # Linear requires >=3 P rows and rank two; isotonic requires >=2 finite
    # P rows and permits constant inputs/truths. Failed pairs abort the fit.
    post <- list(method = post_method, strategy = "ensemble", v = 3L)
    options <- if (method == "inductive") list(m = 1) else list(epsilon = .5, v = 3L)
    spec <- do.call(vennabers::calibrated_workflow, c(list(workflow = base,
      outcome = outcome, schema = schema, method = method, seed = 42,
      postprocessing = post), options))
    resampled <- tune::fit_resamples(spec, outer, metrics = metrics,
      control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, allow_par = FALSE))
    audit$check_encoded_evaluation(resampled, outer, schema)
    states <- tune::collect_extracts(resampled)$.extracts
    for (i in seq_along(states)) for (ensemble in states[[i]]$postprocessing) {
      assert_conditions(ensemble$observed_base_fits == 3L, ensemble$observed_postprocessor_fits == 3L,
        identical(sort(unlist(lapply(ensemble$inner_folds, `[[`, "P"), use.names = FALSE)), sort(ensemble$Q)))
      for (plan in ensemble$inner_folds) {
        Q <- as.integer(outer$splits[[i]], "analysis")
        A <- as.integer(outer$splits[[i]], "assessment")
        assert_conditions(!any(Q[c(plan$B, plan$P, plan$C)] %in% A),
          !length(intersect(plan$B, plan$P)), !length(intersect(plan$Q, plan$C)),
          plan$observed_base_fits == 1L, plan$observed_postprocessor_fits == 1L)
      }
    }
    last <- tune::last_fit(spec, split, metrics = metrics)
    fit <- last$.workflow[[1L]]
    new <- d[91:120, setdiff(names(d), outcome), drop = FALSE]
    native <- vennabers::venn_abers_details(fit, new)
    decoded <- vennabers::venn_abers_decode(native, schema)
    collected <- tune::collect_predictions(last, summarize = FALSE)
    scores <- tune::collect_metrics(last)
    expected <- c(rmse = yardstick::rmse_vec(original[[schema$outcome]][91:120], decoded$prediction),
      mae = yardstick::mae_vec(original[[schema$outcome]][91:120], decoded$prediction))
    assert_conditions(nrow(tune::collect_notes(last)) == 0L, identical(collected$.row, 91:120),
      !anyDuplicated(collected$.row), all(collected$.config == "pre0_mod0_post0"),
      nrow(scores) == 2L, all(abs(scores$.estimate - expected[scores$.metric]) <= 1e-10),
      identical(collected$.pred, native$prediction$.pred), all(is.finite(decoded$prediction)),
      all(native$details$postprocessing_method == post_method),
      all(native$details$postprocessing_input == "numeric"), all(lengths(native$details$inner_fold) == 3L),
      identical(decoded$postprocessing[[1]]$inner_folds[[1]]$columns, names(original)))
    # Mutate A predictors and truth separately: neither may affect any learning.
    before <- workflows::extract_fit_engine(fit)
    for (column in c(names(d)[1L], outcome)) {
      changed <- d; changed[91:120, column] <- changed[91:120, column] + 100
      other <- tune::last_fit(spec, rsample::make_splits(list(analysis = 1:90, assessment = 91:120), changed), metrics = metrics)
      assert_conditions(nrow(tune::collect_notes(other)) == 0L)
      other_fit <- other$.workflow[[1L]]; after <- workflows::extract_fit_engine(other_fit)
      assert_conditions(identical(before$postprocessing, after$postprocessing),
        identical(before$calibrator, after$calibrator), identical(predict(fit, new), predict(other_fit, new)))
      for (id in names(before$models)) for (j in names(before$models[[id]]$models)) {
        left <- before$models[[id]]$models[[j]]; right <- after$models[[id]]$models[[j]]
        if (post_method == "isotonic_numeric") {
          assert_conditions(identical(left$isotonic, right$isotonic))
          left <- left$workflow; right <- right$workflow
        }
        assert_conditions(identical(coef(workflows::extract_fit_engine(left)), coef(workflows::extract_fit_engine(right))),
          identical(workflows::extract_recipe(left)$steps, workflows::extract_recipe(right)$steps))
      }
    }
    cases[[paste(post_method, method)]] <- list(fit = fit, new = new, schema = schema)
  }
}

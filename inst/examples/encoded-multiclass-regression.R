# Installed examples: shared OvO/OvR and regression, Inductive and Cross.
source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), local = TRUE)
if (all(vapply(c("parsnip", "workflows", "hardhat", "rsample", "recipes", "rlang", "tibble",
                 "withr", "generics", "tidyselect", "tune", "yardstick", "nnet"),
               requireNamespace, logical(1), quietly = TRUE))) {
  examples <- list()
  for (classification in c("ovo", "ovr", "regression")) {
    regression <- classification == "regression"
    development <- data.frame(x = sin(1:120), y = if (regression) 3 + sin(1:120) + cos(1:120) else
      factor(rep(c("class", ".pred_yes", "ümlaut"), 40), levels = c("ümlaut", "class", ".pred_yes")))
    names(development) <- c("vaV0001", if (regression) ".pred" else ".pred_class")
    schema <- if (regression) vennabers::venn_abers_schema(development, ".pred") else
      vennabers::venn_abers_schema(development, ".pred_class", classification = classification)
    encoded <- vennabers::venn_abers_encode(development, schema)
    formula <- reformulate(names(encoded)[1], schema$encoded_outcome)
    base <- workflows::workflow() |>
      workflows::add_recipe(recipes::recipe(formula, encoded) |>
                              recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(if (regression) parsnip::linear_reg(engine = "lm") else
        parsnip::multinom_reg(engine = "nnet", penalty = .01))
    # Rebuild prescribed outer memberships and IDs on the encoded data.
    original_outer <- withr::with_seed(41, rsample::vfold_cv(development, v = 2, repeats = 2))
    outer <- rsample::manual_rset(lapply(original_outer$splits, function(s) {
      rsample::make_splits(list(analysis = as.integer(s, data = "analysis"),
                               assessment = as.integer(s, data = "assessment")), encoded)
    }), original_outer[c("id", "id2")])
    scoring <- if (regression) yardstick::metric_set(yardstick::rmse, yardstick::mae) else
      yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
    control <- tune::control_resamples(save_pred = TRUE, extract = identity, allow_par = FALSE,
      parallel_over = "resamples", pkgs = c("vennabers", "stats", "nnet"))
    for (method in c("inductive", "cross")) {
      model <- vennabers::calibrated_workflow(base, schema$encoded_outcome, schema = schema,
                                             method = method, seed = 42)
      fitted <- generics::fit(model, encoded)
      new <- vennabers::venn_abers_encode(development[1:3, "vaV0001", drop = FALSE], schema)
      native <- predict(fitted, new, type = if (regression) "numeric" else "prob")
      decoded <- vennabers::venn_abers_decode(native, schema)
      details <- vennabers::venn_abers_decode(vennabers::venn_abers_details(fitted, new), schema)
      if (!regression) {
        class <- vennabers::venn_abers_decode(predict(fitted, new, type = "class"), schema)$class
        stopifnot(identical(colnames(decoded$probability), schema$class_levels), identical(levels(class), schema$class_levels))
      } else stopifnot(is.double(decoded$prediction), is.null(decoded$probability), is.null(decoded$class))
      evaluation <- withr::with_seed(43, tune::fit_resamples(model, outer, metrics = scoring, control = control))
      checked <- check_encoded_evaluation(evaluation, outer, schema)
      stopifnot(nrow(checked$decoded$identity) == 2L * nrow(development), is.null(checked$decoded$details))
      examples[[paste(classification, method)]] <- list(schema = schema, data = encoded, base = base,
        model = model, fit = fitted, outer = outer, evaluation = evaluation, decoded = decoded, details = details)
    }
  }
}

# Full means an additional complete B/P composition on F (each local class
# pair for pairwise OvO), not an engine fitted on every F row. F gets its own
# holdout/inner folds. Held C predictions and VA calibrators stay unchanged;
# after success only F workflows remain, alongside all source diagnostics.
# Transfer to F changes the predictor distribution: split-calibration guarantees
# do not automatically transfer. CVAR bounds remain descriptive VA bounds.
source(system.file("examples", "assertions.R", package = "vennabers"))
packages <- c("vennabers", "workflows", "recipes", "rsample", "tune", "yardstick", "tailor", "probably", "nnet")
cases <- list()
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  library(vennabers)
  checker <- new.env()
  sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), checker)
  # Ten domain/method choices x two strategies x two VA methods x two sources.
  domains <- data.frame(path = c("binary", "binary", "pairwise", "pairwise",
    "ovo3", "ovr3", "ovo2", "ovr2", "regression", "regression"),
    post = c("sigmoid", "tailor_logistic", "sigmoid", "tailor_logistic",
      "tailor_multinomial", "tailor_multinomial", "tailor_logistic", "tailor_logistic",
      "isotonic_numeric", "tailor_linear"))
  for (i in seq_len(nrow(domains))) {
    path <- domains$path[i]; kind <- domains$post[i]
    regression <- path == "regression"; pairwise <- path == "pairwise"
    classification <- if (pairwise) "ovo" else if (path %in% c("binary", "regression")) "binary" else substr(path, 1, 3)
    k <- if (pairwise || path %in% c("ovo3", "ovr3")) 3L else 2L
    raw <- withr::with_seed(123, {
      cl <- factor(rep(c("a a", "a.a", "class")[seq_len(k)], length.out = 480),
        levels = c("a a", "a.a", "class")[seq_len(k)])
      x <- rnorm(480, sd = 2) + c(-1, 1, -1)[cl]
      z <- rnorm(480, sd = 2) + c(-1, -1, 1)[cl]
      w <- rnorm(480, sd = 2) + c(1, -1, -1)[cl]
      data.frame(x = x, z = z, w = w, truth = if (regression) 10 + 2*x + z + rnorm(480) else cl)
    })
    names(raw) <- c("input value", "input.value", "other value", ".pred_class")
    schema <- if (regression) venn_abers_schema(raw, ".pred_class") else
      venn_abers_schema(raw, ".pred_class", classification = classification)
    d <- venn_abers_encode(raw, schema); outcome <- schema$encoded_outcome
    predictors <- setdiff(names(d), outcome)
    engine <- if (regression) parsnip::linear_reg(engine = "lm") else
      if (pairwise || k == 2L) parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE, maxit = 1000)
    # maxit = 1000 controls only B. Tailor multinomial P still uses probably's
    # default 100 iterations. Any unconverged pair fails the whole composition
    # with fold/pair/B/P context, without retry, method switch or partial ensemble.
    base <- workflows::workflow() |> workflows::add_recipe(
      recipes::recipe(reformulate(predictors, outcome), d) |>
        recipes::step_normalize(recipes::all_numeric_predictors())) |> workflows::add_model(engine)
    # Repeated outer folds: identities are checked on (id, id2, .row, .config).
    # Small nested data: with seed 44 one pairwise Tailor P of 25 rows stops
    # unconverged and that resample fails (no retry, no switch to none); seeds
    # 43 and 45 keep a converged separable pairwise B engine with its retained
    # glm warning, a tune note. This example uses 42, whose evaluation has none.
    outer <- withr::with_seed(42, rsample::vfold_cv(d, v = 2, repeats = 2))
    split <- rsample::make_splits(list(analysis = 1:384, assessment = 385:480), d)
    metrics <- if (regression) yardstick::metric_set(yardstick::rmse, yardstick::mae) else
      yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
    for (strategy in c("holdout", "ensemble")) for (method in c("inductive", "cross"))
      for (mode in c("matched", "full")) {
      post <- list(method = kind, strategy = strategy)
      if (kind == "sigmoid") post$sigmoid_input <- "probability"
      if (strategy == "holdout") post$prop <- .3 else post$v <- 3L
      args <- list(base, outcome, schema = schema, method = method,
        prediction_model = mode, seed = 24, postprocessing = post)
      if (method == "cross") args$v <- 2L
      if (pairwise) args$model_training <- "pairwise"
      model <- do.call(calibrated_workflow, args)
      result <- tune::fit_resamples(model, outer, metrics = metrics,
        control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, allow_par = FALSE))
      checker$check_encoded_evaluation(result, outer, schema)
      states <- checker$check_outer_identities(result, outer, metrics = 2L)$extracts
      for (h in seq_along(states)) {
        components <- if (pairwise) states[[h]]$pair_models else list(states[[h]])
        F_rows <- as.integer(outer$splits[[h]], "analysis"); A_rows <- as.integer(outer$splits[[h]], "assessment")
        for (component in components) {
          rows <- if (pairwise) component$rows else seq_along(F_rows)
          for (diagnostic in component$postprocessing) {
            members <- if (strategy == "ensemble") diagnostic$inner_folds else list(diagnostic)
            for (p in members) assert_conditions(!any(F_rows[rows[c(p$B, p$P, p$C)]] %in% A_rows),
              !length(intersect(p$B, p$P)), !length(intersect(p$Q, p$C)),
              p$observed_base_fits == 1L, p$observed_postprocessor_fits == 1L)
          }
          if (mode == "full") stopifnot(identical(names(component$models), "Full1"),
            identical(component$postprocessing$Full1$Q, seq_along(rows)))
        }
      }
      last <- tune::last_fit(model, split, metrics = metrics)
      fit <- tune::extract_workflow(last); new <- d[385:480, predictors]
      native <- tune::collect_predictions(last)
      decoded <- venn_abers_decode(native, schema)
      values <- if (regression) c(rmse = yardstick::rmse_vec(decoded$truth, decoded$prediction),
        mae = yardstick::mae_vec(decoded$truth, decoded$prediction)) else {
        p <- decoded$probability
        if (ncol(p) == 2L) p <- p[, 1L]
        c(brier_class = yardstick::brier_class_vec(decoded$truth, p),
          mn_log_loss = yardstick::mn_log_loss_vec(decoded$truth, p))
      }
      scores <- tune::collect_metrics(last)
      assert_conditions(nrow(tune::collect_notes(last)) == 0L, identical(native$.row, 385:480),
        !anyDuplicated(native$.row), all(native$.config == "pre0_mod0_post0"),
        identical(native[[outcome]], d[[outcome]][385:480]), nrow(scores) == 2L,
        setequal(scores$.metric, names(values)), all(is.finite(scores$.estimate)),
        all(abs(scores$.estimate - values[scores$.metric]) <= 1e-10))
      prediction <- predict(fit, new, type = if (regression) "numeric" else "prob")
      stopifnot(identical(as.matrix(native[names(prediction)]), as.matrix(prediction)))
      if (mode == "full") for (column in c(predictors[1L], outcome)) {
        changed <- d
        if (column != outcome || regression) changed[[column]][385:480] <- changed[[column]][385:480] + 100 else
          changed[[column]][385:480] <- rev(changed[[column]][385:480])
        other <- tune::last_fit(model, rsample::make_splits(list(analysis = 1:384, assessment = 385:480), changed), metrics = metrics)
        stopifnot(nrow(tune::collect_notes(other)) == 0L,
          identical(venn_abers_details(fit, new), venn_abers_details(tune::extract_workflow(other), new)))
      }
      cases[[paste(path, kind, strategy, method, mode)]] <- list(fit = fit, new = new, schema = schema)
      cat(path, kind, strategy, method, mode, ": native evaluation complete.\n")
    }
  }
  stopifnot(length(cases) == 80L)
}

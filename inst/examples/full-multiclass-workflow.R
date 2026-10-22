# Shared multiclass full-model transfer is an R extension. OvO still means
source(system.file("examples", "assertions.R", package = "vennabers"))
# pairwise calibration of one shared class matrix, not pairwise model training.
packages <- c("vennabers", "parsnip", "workflows", "recipes", "rsample", "nnet")
examples <- list()
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  library(vennabers)
  data <- data.frame(x = sin(1:120), z = cos(1:120),
    truth = factor(rep(c("class", ".pred_class", "a b"), 40),
      levels = c("class", ".pred_class", "a b")))
  names(data) <- c("a b", "a.b", ".pred_class")
  for (method in c("inductive", "cross")) for (classification in c("ovo", "ovr")) {
    schema <- venn_abers_schema(data, ".pred_class", classification = classification)
    encoded <- venn_abers_encode(data, schema)
    outer_split <- rsample::make_splits(list(analysis = 1:90, assessment = 91:120), encoded)
    development <- rsample::analysis(outer_split)
    assessment <- rsample::assessment(outer_split)
    base <- workflows::workflow() |>
      workflows::add_recipe(recipes::recipe(reformulate(names(encoded)[1:2], schema$encoded_outcome), development) |>
        recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
    plan <- local({
      cross <- method == "cross"
      function(data) {
        n <- nrow(data)
        held <- if (cross) split(seq_len(n), cut(seq_len(n), 3, labels = FALSE)) else
          list(seq.int(n - ceiling(n * .25) + 1L, n))
        parts <- lapply(held, function(h) rsample::make_splits(
          list(analysis = setdiff(seq_len(n), h), assessment = h), data))
        if (cross) rsample::manual_rset(parts, paste0("Part", seq_along(parts))) else parts[[1L]]
      }
    })
    models <- lapply(c("matched", "full"), function(mode) calibrated_workflow(base,
      schema$encoded_outcome, schema = schema, method = method, model_training = "shared",
      prediction_model = mode, splits = plan, seed = 42, precision = 2))
    fits <- lapply(models, generics::fit, data = development)
    source <- workflows::extract_fit_engine(fits[[1L]])
    state <- workflows::extract_fit_engine(fits[[2L]])
    stopifnot(identical(source$calibrator, state$calibrator), identical(source$calibration, state$calibration))
    native <- venn_abers_details(fits[[2L]], assessment)
    decoded <- venn_abers_decode(native, schema)
    print(state)
    print(head(native$details))
    print(head(decoded$details))
    print(head(decoded$probability))
    print(head(venn_abers_decode(predict(fits[[2L]], assessment, type = "class"), schema)$class))
    case <- list(schema = schema, data = encoded, base = base, model = models[[2L]],
      fit = fits[[2L]], matched = fits[[1L]], assessment = assessment)
    if (all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
      audit <- new.env()
      sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), audit)
      outer <- withr::with_seed(41, rsample::vfold_cv(development, v = 2, repeats = 2))
      scoring <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
      result <- tune::fit_resamples(models[[2L]], outer, metrics = scoring,
        control = tune::control_resamples(save_pred = TRUE, extract = identity, allow_par = FALSE))
      checked <- audit$check_encoded_evaluation(result, outer, schema)
      print(checked$metrics)
      final <- tune::last_fit(models[[2L]], outer_split, metrics = scoring)
      final_p <- tune::collect_predictions(final)
      final_m <- tune::collect_metrics(final)
      assert_conditions(nrow(tune::collect_notes(final)) == 0L, identical(final_p$.row, 91:120),
        !anyDuplicated(final_p$.row), all(final_p$.config == "pre0_mod0_post0"),
        identical(as.matrix(final_p[names(native$prediction)]), as.matrix(native$prediction)),
        nrow(final_m) == 2L, all(is.finite(final_m$.estimate)))
      # yardstick uses half the manual multiclass matrix Brier score, including
      # the two-class multiclass contract. Log clipping affects metrics only.
      truth <- venn_abers_decode(final_p, schema)$truth
      stopifnot(abs(final_m$.estimate[final_m$.metric == "brier_class"] -
        brier_score(truth, decoded$probability) / 2) <= 1e-10)
      case$evaluation <- result
      case$outer <- outer
    }
    examples[[paste(method, classification)]] <- case
  }
  # Full transfer costs two Inductive or K+1 Cross fits and one prediction per
  # nonempty batch. The same full matrix reaches all held-out calibrators.
  # This R extension makes no new calibration, coverage or quality guarantee.
}

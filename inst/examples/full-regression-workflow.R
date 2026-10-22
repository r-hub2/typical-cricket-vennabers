# Fixed-setting IVAR/CVAR full-model transfer. Python 1.5.4 exposes this option
source(system.file("examples", "assertions.R", package = "vennabers"))
# in the exported lower VennAbersCV class, not its regression facade.
packages <- c("vennabers", "parsnip", "workflows", "recipes", "rsample")
examples <- list()
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  library(vennabers)
  data <- data.frame(x = sin(1:120), z = cos(1:120), truth = 10 + 3 * sin(1:120) + cos(1:120)^2)
  names(data) <- c("a b", "a.b", ".pred")
  schema <- venn_abers_schema(data, ".pred")
  encoded <- venn_abers_encode(data, schema)
  outer_split <- rsample::make_splits(list(analysis = 1:90, assessment = 91:120), encoded)
  development <- rsample::analysis(outer_split)
  assessment <- rsample::assessment(outer_split)
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(reformulate(names(encoded)[1:2], schema$encoded_outcome), development) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  for (method in c("inductive", "cross")) {
    # Ordered local partitions are replayed for both modes and rebuilt inside
    # each outer analysis domain. This synthetic example assumes IID data.
    plan <- local({
      cross <- method == "cross"
      function(data) {
        n <- nrow(data)
        held <- if (cross) split(seq_len(n), cut(seq_len(n), 3, labels = FALSE)) else list(seq.int(n - 9L, n))
        parts <- lapply(held, function(h) rsample::make_splits(
          list(analysis = setdiff(seq_len(n), h), assessment = h), data))
        if (cross) rsample::manual_rset(parts, paste0("Part", seq_along(parts))) else parts[[1L]]
      }
    })
    models <- lapply(c("matched", "full"), function(mode) calibrated_workflow(base,
      schema$encoded_outcome, schema = schema, method = method, prediction_model = mode,
      splits = plan, m = 2, seed = 42))
    fits <- lapply(models, generics::fit, data = development)
    native <- venn_abers_details(fits[[2L]], assessment)
    decoded <- venn_abers_decode(native, schema)
    print(data.frame(method = method, matched = predict(fits[[1L]], assessment)$.pred, full = decoded$prediction))
    print(head(decoded$details))
    # Reconstruct the final points and descriptive bounds independently from
    # complete per-fold results. Never merge the averaged bounds into a point.
    reconstructed <- stats::aggregate(cbind(prediction, lower, upper) ~ row, decoded$details, mean)
    stopifnot(max(abs(reconstructed$prediction - decoded$prediction)) <= 1e-10)
    if (method == "cross") {
      stopifnot(max(abs(as.matrix(reconstructed[c("lower", "upper")]) -
        as.matrix(decoded$aggregate[c("lower", "upper")]))) <= 1e-10)
      print(head(decoded$aggregate))
    }
    case <- list(schema = schema, data = encoded, base = base, model = models[[2L]],
      fit = fits[[2L]], matched = fits[[1L]], assessment = assessment)
    if (all(vapply(c("tune", "yardstick"), requireNamespace, logical(1), quietly = TRUE))) {
      audit <- new.env()
      sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), audit)
      scoring <- yardstick::metric_set(yardstick::rmse, yardstick::mae)
      truth <- assessment[[schema$encoded_outcome]]
      print(data.frame(method = method, mode = c("matched", "full"),
        rmse = vapply(fits, function(f) yardstick::rmse_vec(truth, predict(f, assessment)$.pred), 0),
        mae = vapply(fits, function(f) yardstick::mae_vec(truth, predict(f, assessment)$.pred), 0)))
      outer <- withr::with_seed(41, rsample::vfold_cv(development, v = 2, repeats = 2))
      result <- tune::fit_resamples(models[[2L]], outer, metrics = scoring,
        control = tune::control_resamples(save_pred = TRUE, extract = identity, allow_par = FALSE))
      checked <- audit$check_encoded_evaluation(result, outer, schema)
      print(checked$metrics)
      final <- tune::last_fit(models[[2L]], outer_split, metrics = scoring)
      p <- tune::collect_predictions(final); metrics <- tune::collect_metrics(final)
      assert_conditions(nrow(tune::collect_notes(final)) == 0L, identical(p$.row, 91:120), !anyDuplicated(p$.row),
        all(p$.config == "pre0_mod0_post0"), identical(p[[schema$encoded_outcome]], truth),
        identical(p$.pred, native$prediction$.pred), nrow(metrics) == 2L,
        setequal(metrics$.metric, c("rmse", "mae")), all(is.finite(metrics$.estimate)),
        all(metrics$.config == "pre0_mod0_post0"))
      stopifnot(abs(metrics$.estimate[metrics$.metric == "rmse"] - yardstick::rmse_vec(truth, decoded$prediction)) <= 1e-10,
        abs(metrics$.estimate[metrics$.metric == "mae"] - yardstick::mae_vec(truth, decoded$prediction)) <= 1e-10)
      case$evaluation <- result; case$outer <- outer
    }
    examples[[method]] <- case
  }
  # Full transfer fits 2 / K+1 models, keeps only the full workflow plus held-out
  # calibration values/calibrators, and predicts once per nonempty batch.
  # No automatic coverage, calibration or quality guarantee follows. Aggregate
  # regression bounds are descriptive, not confidence or prediction intervals.
}

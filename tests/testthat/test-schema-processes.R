testthat::test_that("schema-processes", {
library(vennabers)
if (all(vapply(c("parsnip", "workflows", "recipes", "tune", "yardstick", "future", "rsample"),
               requireNamespace, logical(1), quietly = TRUE))) local({
  saved <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  output <- tempfile(fileext = ".log")
  on.exit(unlink(c(saved, script, output)), add = TRUE)
  data <- data.frame(x = sin(1:80), y = factor(rep(c("class", ".pred_yes", ".pred_yes", "class"), 20),
                                             levels = c("class", ".pred_yes")))
  names(data) <- c("vaV0001", ".pred_class")
  schema <- venn_abers_schema(data, ".pred_class", event_level = "second")
  encoded <- venn_abers_encode(data, schema)
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(reformulate(names(encoded)[1], schema$encoded_outcome), encoded) |>
                            recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  outer <- withr::with_seed(123, rsample::vfold_cv(encoded, v = 2, repeats = 2))
  for (method in c("inductive", "cross")) for (prediction_model in c("matched", "full")) {
    model <- calibrated_workflow(base, schema$encoded_outcome, schema = schema, method = method, seed = 42, prediction_model = prediction_model)
    fit <- generics::fit(model, encoded)
    expected <- lapply(c(0L, 1L, 3L), function(n) {
      new <- encoded[seq_len(n), 1, drop = FALSE]
      list(probability = venn_abers_decode(predict(fit, new, type = "prob"), schema),
           class = venn_abers_decode(predict(fit, new, type = "class"), schema),
           details = venn_abers_decode(venn_abers_details(fit, new), schema))
    })
    saveRDS(list(schema = schema, model = model, fit = fit, data = encoded, outer = outer, expected = expected), saved)
    writeLines(c(
      paste0('.libPaths(', paste(deparse(.libPaths()), collapse = ''), ')'),
      'library(vennabers); loadNamespace("workflows"); loadNamespace("stats")',
      paste0('s <- readRDS(', deparse(saved), ')'),
      'stopifnot(identical(workflows::extract_fit_engine(s$fit)$schema, s$schema))',
      'trace("fit.workflow", where = asNamespace("workflows"), tracer = quote(stop("Unexpected learning")), print = FALSE)',
      'trace("prep.recipe", where = asNamespace("recipes"), tracer = quote(stop("Unexpected recipe learning")), print = FALSE)',
      'trace("venn_abers", where = asNamespace("vennabers"), tracer = quote(stop("Unexpected calibration learning")), print = FALSE)',
      'actual <- lapply(c(0L, 1L, 3L), function(n) {',
      '  new <- s$data[seq_len(n), 1, drop = FALSE]',
      '  list(probability = venn_abers_decode(predict(s$fit, new, type = "prob"), s$schema),',
      '       class = venn_abers_decode(predict(s$fit, new, type = "class"), s$schema),',
      '       details = venn_abers_decode(venn_abers_details(s$fit, new), s$schema))',
      '})',
      'stopifnot(identical(actual, s$expected))',
      'untrace("fit.workflow", where = asNamespace("workflows"))',
      'untrace("prep.recipe", where = asNamespace("recipes"))',
      'untrace("venn_abers", where = asNamespace("vennabers"))',
      'Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")',
      'control <- tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine,',
      '  event_level = s$schema$event_level, parallel_over = "resamples", pkgs = c("vennabers", "stats"))',
      'run <- function() withr::with_seed(124, tune::fit_resamples(s$model, s$outer,',
      '  metrics = yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss), control = control))',
      'future::plan(future::multisession, workers = 2)',
      'one <- run(); two <- run()',
      'future::plan(future::sequential)',
      'stopifnot(nrow(tune::collect_notes(one)) == 0L, nrow(tune::collect_notes(two)) == 0L,',
      '  identical(tune::collect_predictions(one), tune::collect_predictions(two)),',
      '  identical(tune::collect_metrics(one), tune::collect_metrics(two)))',
      'decoded <- venn_abers_decode(tune::collect_predictions(one, summarize = FALSE), s$schema)',
      'stopifnot(nrow(decoded$identity) == 2L * nrow(s$data),',
      '  all(vapply(tune::collect_extracts(one)$.extracts, function(x) identical(x$schema, s$schema), logical(1))))'
    ), script)
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = output, stderr = output)
    if (status != 0L) stop(paste(readLines(output), collapse = "\n"))
  }
})

  testthat::expect_true(TRUE)
})

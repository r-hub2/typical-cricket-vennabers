testthat::test_that("full-multiclass-processes", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")
packages <- c("workflows", "parsnip", "recipes", "rsample", "nnet")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  example <- new.env()
  sys.source(system.file("examples", "full-multiclass-workflow.R", package = "vennabers"), example)
  error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
  snapshot <- function(fit, new, schema) lapply(c("log", "brier"), function(loss) {
    lapply(c(0L, 1L, 7L), function(n) {
      x <- new[seq_len(n), , drop = FALSE]
      native <- list(prob = predict(fit, x, type = "prob", loss = loss),
        class = predict(fit, x, type = "class", loss = loss), details = venn_abers_details(fit, x, loss = loss))
      list(native = native, decoded = lapply(native, venn_abers_decode, schema = schema))
    })
  })
  saved <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  log <- tempfile(fileext = ".log")
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  workers <- worker_replay_available()
  for (case in example$examples) {
    schema <- case$schema
    fit <- case$fit
    new <- case$assessment[setdiff(names(case$assessment), schema$encoded_outcome)]
    native <- venn_abers_details(fit, new)
    decoded <- venn_abers_decode(native, schema)
    assert_conditions(identical(colnames(decoded$probability), schema$class_levels),
      all(decoded$details$model_training == "shared"), all(decoded$details$prediction_model == "full"),
      all(decoded$details$model_source == "Full1"),
      identical(decoded$details$fold, native$details$fold),
      all(decoded$details$event_class %in% schema$class_levels),
      identical(attr(native, "vennabers_schema"), schema))
    if (schema$classification == "ovo") stopifnot(all(decoded$details$event_class == decoded$details$class_2)) else
      stopifnot(all(decoded$details$event_class == decoded$details$class))
    foreign <- venn_abers_schema(data.frame(x = 1:6, truth = factor(rep(c("a", "b", "c"), 2))), "truth",
      classification = schema$classification)
    error(venn_abers_decode(native, foreign))
    stopifnot(identical(native, venn_abers_details(fit, new[rev(names(new))])))
    wrong_levels <- case$assessment
    wrong_levels[[schema$encoded_outcome]] <- factor(wrong_levels[[schema$encoded_outcome]], levels = rev(unname(schema$classes)))
    wrong_schema <- new
    attr(wrong_schema, "vennabers_schema") <- foreign
    for (bad in list(wrong_levels, wrong_schema, new[-1L], setNames(new, c("unknown", names(new)[2L])))) {
      error(predict(fit, bad, type = "prob"))
      error(venn_abers_details(fit, bad))
    }
    error(calibrated_workflow(case$base, schema$encoded_outcome, schema = schema,
      classification = if (schema$classification == "ovo") "ovr" else "ovo", prediction_model = "full"))
    saveRDS(list(fit = fit, schema = schema, new = new, snapshot = snapshot,
      expected = snapshot(fit, new, schema), replay = if (workers) list(model = case$model,
        outer = case$outer, evaluation = case$evaluation, schema = schema, new = new)), saved)
    writeLines(c(
      paste0('.libPaths(', paste(deparse(.libPaths()), collapse = ''), ')'),
      'library(vennabers)',
      paste0('s <- readRDS(', deparse(saved), ')'),
      'trace("fit.workflow", where = asNamespace("workflows"), tracer = quote(stop("Unexpected model fit")), print = FALSE)',
      'trace("prep.recipe", where = asNamespace("recipes"), tracer = quote(stop("Unexpected preprocessing")), print = FALSE)',
      'trace("venn_abers", where = asNamespace("vennabers"), tracer = quote(stop("Unexpected calibrator fit")), print = FALSE)',
      'stopifnot(identical(s$expected, s$snapshot(s$fit, s$new, s$schema)),',
      '  identical(names(workflows::extract_fit_engine(s$fit)$models), "Full1"))',
      if (workers) worker_replay("yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)",
        c("vennabers", "stats", "nnet"), c(workflows = "fit.workflow", recipes = "prep.recipe", vennabers = "venn_abers"))
    ), script)
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
    if (status != 0L) stop(paste(readLines(log), collapse = "\n"))
  }
})

  testthat::expect_true(TRUE)
})

testthat::test_that("adapter fits restore in fresh processes and replay on outer workers", {
  skip_if_no_integration("nnet", "tune", "yardstick")
  # The example's factories bind helpers and settings locally, as a user script would.
  example <- new.env(parent = globalenv())
  sys.source(system.file("examples", "model-adapter.R", package = "vennabers"), example)
  testthat::expect_named(example$examples, c("binary", "multiclass", "regression", "search"))
  snapshot <- function(fit, new, schema = NULL) lapply(c(0L, 1L, 7L), function(n) {
    batch <- new[seq_len(n), , drop = FALSE]
    regression <- workflows::extract_fit_engine(fit)$parameters$requested$mode == "regression"
    list(prediction = predict(fit, batch), probability = if (!regression) predict(fit, batch, type = "prob"),
         details = venn_abers_details(fit, batch))
  })
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  for (case in names(example$examples)) {
    e <- example$examples[[case]]
    # Fresh process with only vennabers attached: no learning, every output identical.
    saveRDS(list(fit = e$fit, new = e$new, expected = snapshot(e$fit, e$new), snapshot = snapshot), saved)
    writeLines(c(installed_vennabers_lines(),
      paste0("s <- readRDS(", deparse(saved), ")"),
      "guard <- function(...) stop('Unexpected learning')",
      "registerS3method('fit', 'workflow', guard, envir = asNamespace('generics'))",
      "for (f in c('venn_abers', 'venn_abers_regression', 'cross_venn_abers', 'cross_venn_abers_regression')) trace(f, where = asNamespace('vennabers'), tracer = quote(guard()), print = FALSE)",
      "trace('prep.recipe', where = asNamespace('recipes'), tracer = quote(guard()), print = FALSE)",
      "stopifnot(identical(s$snapshot(s$fit, s$new), s$expected))"), script)
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
    testthat::expect_equal(status, 0L, info = paste(case, paste(readLines(log), collapse = "\n")))
    if (!worker_replay_available()) next
    # Two multisession workers refit the calibrated adapter workflow per outer fold.
    outer <- withr::with_seed(93, rsample::vfold_cv(e$data, v = 2))
    metrics <- if (case == "regression") "yardstick::metric_set(yardstick::rmse)" else "yardstick::metric_set(yardstick::mn_log_loss)"
    saveRDS(list(snapshot = snapshot, replay = list(model = e$model, outer = outer, new = e$new, schema = NULL)), saved)
    run_guarded_replay(saved, metrics, c("vennabers", "recipes", "nnet"), info = case)
  }
})

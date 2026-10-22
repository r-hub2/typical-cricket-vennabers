testthat::test_that("pairwise postprocessing supports native evaluation and fresh RDS prediction", {
  for (pkg in c("tailor", "probably", "recipes", "tune", "yardstick")) testthat::skip_if_not_installed(pkg)
  example <- new.env()
  sys.source(system.file("examples", "pairwise-postprocessing.R", package = "vennabers"), example)
  testthat::expect_length(example$cases, 12L)
  snapshot <- function(fit, new, schema) lapply(c("log", "brier"), function(loss) {
    lapply(c(0L, 1L, 7L), function(n) {
      batch <- new[seq_len(n), , drop = FALSE]
      details <- venn_abers_details(fit, batch, loss = loss)
      list(prob = predict(fit, batch, type = "prob", loss = loss),
        class = predict(fit, batch, type = "class", loss = loss),
        details = details, decoded = venn_abers_decode(details, schema))
    })
  })
  records <- lapply(example$cases, function(r) {
    state <- workflows::extract_fit_engine(r$fit)
    testthat::expect_true(all(c("R", "vennabers", "workflows") %in% names(state$versions)))
    if (state$parameters$requested$postprocessing$method == "tailor_logistic")
      testthat::expect_true(all(c("tailor", "probably") %in% names(state$versions)))
    r$expected <- snapshot(r$fit, r$new, r$schema)
    r
  })
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  saveRDS(records, saved)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "library(vennabers); loadNamespace('workflows')",
    "guard <- function(...) stop('Unexpected learning')",
    "registerS3method('fit', 'workflow', guard, envir=asNamespace('generics'))",
    "registerS3method('fit', 'tailor', guard, envir=asNamespace('generics'))",
    "trace('prep.recipe', where=asNamespace('recipes'), tracer=quote(stop('Unexpected learning')), print=FALSE)",
    "trace('optim', where=asNamespace('stats'), tracer=quote(stop('Unexpected learning')), print=FALSE)",
    "trace('glm.fit', where=asNamespace('stats'), tracer=quote(stop('Unexpected learning')), print=FALSE)",
    paste0("snapshot <- ", paste(deparse(snapshot), collapse = "\n")),
    paste0("records <- readRDS(", deparse(saved), ")"),
    "for (r in records) stopifnot(identical(r$expected, snapshot(r$fit, r$new, r$schema)))"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

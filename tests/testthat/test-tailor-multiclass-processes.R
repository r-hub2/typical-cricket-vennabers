testthat::test_that("multiclass Tailor survives native evaluation and fresh R without learning", {
  for (pkg in c("tailor", "probably", "recipes", "nnet", "tune", "yardstick")) testthat::skip_if_not_installed(pkg)
  example <- new.env()
  sys.source(system.file("examples", "tailor-multiclass-holdout.R", package = "vennabers"), example)
  snapshot <- function(fit, new, schema) lapply(c("log", "brier"), function(loss) {
    lapply(c(0L, 1L, 7L), function(n) {
      x <- new[seq_len(n), , drop = FALSE]
      native <- list(prob = predict(fit, x, type = "prob", loss = loss),
        class = predict(fit, x, type = "class", loss = loss), details = venn_abers_details(fit, x, loss = loss))
      list(native = native, decoded = lapply(native, venn_abers_decode, schema = schema))
    })
  })
  records <- lapply(example$cases, function(r) {
    r$expected <- snapshot(r$fit, r$new, r$schema)
    r
  })
  testthat::expect_length(records, 12L)
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  saveRDS(records, saved)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "library(vennabers); loadNamespace('workflows'); loadNamespace('tailor')",
    "guard <- function(...) stop('Unexpected learning')",
    "registerS3method('fit', 'workflow', guard, envir=asNamespace('generics'))",
    "registerS3method('fit', 'tailor', guard, envir=asNamespace('generics'))",
    "trace('prep.recipe', where=asNamespace('recipes'), tracer=quote(stop('Unexpected learning')), print=FALSE)",
    paste0("snapshot <- ", paste(deparse(snapshot), collapse = "\n")),
    paste0("records <- readRDS(", deparse(saved), ")"),
    "for (r in records) stopifnot(identical(r$expected, snapshot(r$fit, r$new, r$schema)))"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

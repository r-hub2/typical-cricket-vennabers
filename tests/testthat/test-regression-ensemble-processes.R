testthat::test_that("numeric ensembles survive native evaluation and fresh R without learning", {
  for (pkg in c("tailor", "probably", "recipes", "tune", "yardstick")) testthat::skip_if_not_installed(pkg)
  example <- new.env()
  sys.source(system.file("examples", "regression-ensemble.R", package = "vennabers"), example)
  testthat::expect_length(example$cases, 4L)
  snapshot <- function(fit, new, schema) lapply(c(0L, 1L, 7L), function(n) {
    batch <- new[seq_len(n), , drop = FALSE]
    native <- list(prediction = predict(fit, batch), details = venn_abers_details(fit, batch))
    list(native = native, decoded = lapply(native, venn_abers_decode, schema = schema))
  })
  records <- lapply(example$cases, function(r) {
    native <- venn_abers_details(r$fit, r$new)
    decoded <- venn_abers_decode(native, r$schema)
    for (id in names(native$postprocessing)) for (j in names(native$postprocessing[[id]]$inner_folds)) {
      left <- native$postprocessing[[id]]$inner_folds[[j]]
      right <- decoded$postprocessing[[id]]$inner_folds[[j]]
      testthat::expect_identical(left$inner_fold, right$inner_fold)
      testthat::expect_identical(left$schema, r$schema)
      testthat::expect_identical(left$B, right$B)
      testthat::expect_identical(left$P, right$P)
      testthat::expect_identical(left$input, "numeric")
      testthat::expect_identical(left$format_version, 1L)
    }
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
    "trace('isoreg', where=asNamespace('stats'), tracer=quote(stop('Unexpected learning')), print=FALSE)",
    paste0("snapshot <- ", paste(deparse(snapshot), collapse = "\n")),
    paste0("records <- readRDS(", deparse(saved), ")"),
    "for (r in records) stopifnot(identical(r$expected, snapshot(r$fit, r$new, r$schema)))"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

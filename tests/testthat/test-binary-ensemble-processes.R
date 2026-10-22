testthat::test_that("Binary ensembles produce complete native outer results", {
  for (pkg in c("tune", "yardstick", "recipes", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  source(system.file("examples", "binary-ensemble.R", package = "vennabers"), local = new.env())
  testthat::expect_true(TRUE)
})

testthat::test_that("Ensemble schemas and full predictive state survive fresh R without learning", {
  for (pkg in c("workflows", "recipes", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  original <- data.frame(x = sin(1:160), truth = factor(rep(c("a b", "a.b", "a.b", "a b"), 40)))
  names(original) <- c("input value", ".pred_class")
  schema <- venn_abers_schema(original, ".pred_class", event_level = "second")
  d <- venn_abers_encode(original, schema)
  outcome <- schema$encoded_outcome
  base <- workflows::workflow() |> workflows::add_formula(reformulate(setdiff(names(d), outcome), outcome)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  records <- list()
  outputs <- function(fit, new) lapply(c("log", "brier"), function(loss) {
    lapply(c(0L, 1L, nrow(new)), function(n) {
      x <- new[seq_len(n), , drop = FALSE]
      list(prob = predict(fit, x, type = "prob", loss = loss),
           class = predict(fit, x, type = "class", loss = loss),
           details = venn_abers_details(fit, x, loss = loss))
    })
  })
  for (kind in c("tailor_logistic", "probability", "glm_link")) for (method in c("inductive", "cross")) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = "ensemble", v = 3L)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    fit <- generics::fit(calibrated_workflow(base, outcome, schema = schema, method = method, seed = 24,
      postprocessing = post), d)
    new <- d[setdiff(names(d), outcome)]
    details <- venn_abers_details(fit, new)
    decoded <- venn_abers_decode(details, schema)
    testthat::expect_identical(colnames(decoded$probability), levels(original[[schema$outcome]]))
    testthat::expect_identical(decoded$postprocessing[[1]]$inner_folds[[1]]$classes, schema$class_levels)
    testthat::expect_identical(decoded$postprocessing[[1]]$inner_folds[[1]]$columns, names(original))
    testthat::expect_identical(details$postprocessing[[1]]$inner_folds[[1]]$schema, schema)
    if (kind != "tailor_logistic") testthat::expect_identical(decoded$postprocessing[[1]]$inner_folds[[1]]$sigmoid_class, schema$class_levels[2])
    records[[paste(kind, method)]] <- list(fit = fit, new = new, expected = outputs(fit, new))
  }
  rds <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(rds, script, log)), add = TRUE)
  saveRDS(records, rds)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "library(vennabers)",
    "guard <- function(...) stop('Unexpected learning or plan resolution')",
    "loadNamespace('workflows')",
    "registerS3method('fit', 'workflow', guard, envir=asNamespace('generics'))",
    "registerS3method('fit', 'tailor', guard, envir=asNamespace('generics'))",
    "registerS3method('prep', 'recipe', guard, envir=asNamespace('recipes'))",
    "trace('optim', where=asNamespace('stats'), tracer=quote(stop('Unexpected learning')), print=FALSE)",
    "trace('glm.fit', where=asNamespace('stats'), tracer=quote(stop('Unexpected learning')), print=FALSE)",
    paste0("outputs <- ", paste(deparse(outputs), collapse = "\n")),
    paste0("records <- readRDS(", deparse(rds), ")"),
    "for (r in records) stopifnot(identical(r$expected, outputs(r$fit, r$new)))"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

testthat::test_that("pairwise-processes", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")
packages <- c("workflows", "parsnip", "recipes", "rsample")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) for (example_name in c("pairwise-inductive.R", "pairwise-cross.R")) local({
  example <- new.env()
  sys.source(system.file("examples", example_name, package = "vennabers"), example)
  objects <- as.list(example)
  engine <- function(x) tryCatch(workflows::extract_fit_engine(x), error = function(...) NULL)
  fits <- function(mode) Filter(function(x) {
    state <- engine(x)
    inherits(x, "workflow") && !is.null(state) &&
      identical(state$parameters$requested$model_training, "pairwise") &&
      identical(state$parameters$requested$prediction_model, mode)
  }, objects)
  for (mode in c("matched", "full")) {
    fit <- fits(mode)[[1L]]
    state <- engine(fit)
    schema <- state$schema
    assessment <- Filter(function(x) is.data.frame(x) && schema$encoded_outcome %in% names(x), objects)[[1L]]
    new <- assessment[setdiff(names(assessment), schema$encoded_outcome)]
    snapshot <- function(fit, new, schema) lapply(c("log", "brier"), function(loss) {
      lapply(c(0L, 1L, 7L), function(n) {
        x <- new[seq_len(n), , drop = FALSE]
        native <- list(prob = predict(fit, x, type = "prob", loss = loss),
          class = predict(fit, x, type = "class", loss = loss), details = venn_abers_details(fit, x, loss = loss))
        list(native = native, decoded = lapply(native, venn_abers_decode, schema = schema))
      })
    })
    native <- venn_abers_details(fit, new)
    decoded <- venn_abers_decode(native, schema)
    assert_conditions(identical(colnames(decoded$probability), schema$class_levels),
      all(decoded$details$model_training == "pairwise"), all(decoded$details$prediction_model == mode),
      identical(decoded$details$model_source, native$details$model_source),
      identical(decoded$details$pair_id, native$details$pair_id),
      identical(decoded$details$fold, native$details$fold),
      all(decoded$details$event_class == decoded$details$class_2),
      all(decoded$details$event_class %in% schema$class_levels),
      identical(attr(native, "vennabers_schema"), schema))
    state <- workflows::extract_fit_engine(fit)
    stopifnot(identical(state$schema, schema))
    for (component in state$pair_models) for (model in component$models) {
      mold <- workflows::extract_mold(model)
      stopifnot(nlevels(mold$outcomes[[1L]]) == 2L,
        identical(levels(mold$outcomes[[1L]]), component$class_levels))
    }
    error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
    wrong <- new
    attr(wrong, "vennabers_schema") <- venn_abers_schema(iris, "Species", classification = "ovo")
    error(predict(fit, wrong, type = "prob"))
    error(venn_abers_details(fit, wrong))
    error(venn_abers_decode(native, attr(wrong, "vennabers_schema")))
    saved <- tempfile(fileext = ".rds")
    script <- tempfile(fileext = ".R")
    log <- tempfile(fileext = ".log")
    on.exit(unlink(c(saved, script, log)), add = TRUE)
    workers <- !is.null(objects$outer) && worker_replay_available()
    saveRDS(list(fit = fit, schema = schema, new = new, snapshot = snapshot,
      expected = snapshot(fit, new, schema), replay = if (workers) list(
        model = objects[[if (mode == "full") "full_model" else "model"]], outer = objects$outer,
        evaluation = if (mode == "matched") objects$result, schema = schema, new = new)), saved)
    writeLines(c(
      paste0('.libPaths(', paste(deparse(.libPaths()), collapse = ''), ')'),
      'library(vennabers)',
      paste0('s <- readRDS(', deparse(saved), ')'),
      'trace("fit.workflow", where = asNamespace("workflows"), tracer = quote(stop("Unexpected model fit")), print = FALSE)',
      'trace("prep.recipe", where = asNamespace("recipes"), tracer = quote(stop("Unexpected preprocessing")), print = FALSE)',
      'trace("venn_abers", where = asNamespace("vennabers"), tracer = quote(stop("Unexpected calibrator fit")), print = FALSE)',
      'stopifnot(identical(s$expected, s$snapshot(s$fit, s$new, s$schema)))',
      if (workers) worker_replay("yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)",
        c("vennabers", "stats"), c(workflows = "fit.workflow", recipes = "prep.recipe", vennabers = "venn_abers"))
    ), script)
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
    if (status != 0L) stop(paste(readLines(log), collapse = "\n"))
  }
})

  testthat::expect_true(TRUE)
})

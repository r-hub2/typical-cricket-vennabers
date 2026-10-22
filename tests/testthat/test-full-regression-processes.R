testthat::test_that("full-regression-processes", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")
packages <- c("workflows", "parsnip", "recipes", "rsample")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  example <- new.env()
  sys.source(system.file("examples", "full-regression-workflow.R", package = "vennabers"), example)
  error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
  snapshot <- function(fit, new, schema) lapply(c(0L, 1L, 7L), function(n) {
    x <- new[seq_len(n), , drop = FALSE]
    native <- list(prediction = predict(fit, x), details = venn_abers_details(fit, x))
    list(native = native, decoded = lapply(native, venn_abers_decode, schema = schema))
  })
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile(fileext = ".log")
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  workers <- worker_replay_available()
  for (method in c("inductive", "cross")) for (outcome in c(".pred", ".row", "id")) {
    case <- example$examples[[method]]
    raw <- example$data; names(raw)[3L] <- outcome
    schema <- venn_abers_schema(raw, outcome)
    encoded <- venn_abers_encode(raw, schema)
    options <- workflows::extract_fit_engine(case$fit)$parameters$requested
    base <- workflows::update_recipe(case$base,
      recipes::recipe(reformulate(names(encoded)[1:2], schema$encoded_outcome), encoded[1:90, ]) |>
        recipes::step_normalize(recipes::all_numeric_predictors()))
    model <- calibrated_workflow(base, schema$encoded_outcome, schema = schema,
      method = method, prediction_model = "full", m = 2, seed = 42, splits = options$splits)
    fit <- generics::fit(model, encoded[1:90, ])
    new <- encoded[91:120, setdiff(names(encoded), schema$encoded_outcome)]
    native <- venn_abers_details(fit, new)
    decoded <- venn_abers_decode(native, schema)
    assert_conditions(identical(decoded$outcome, outcome), identical(decoded$prediction, native$prediction$.pred),
      identical(decoded$details, native$details), all(decoded$details$model_source == "Full1"),
      identical(attr(native, "vennabers_schema"), schema))
    if (method == "cross") stopifnot(identical(decoded$aggregate, native$aggregate))
    stopifnot(identical(native, venn_abers_details(fit, new[rev(names(new))])))
    rows <- native$prediction
    rows[[schema$encoded_outcome]] <- encoded[[schema$encoded_outcome]][91:120]
    rows$.row <- 91:120; rows$id <- "Assessment"; rows$.config <- "pre0_mod0_post0"
    decoded_rows <- venn_abers_decode(rows, schema)
    stopifnot(identical(decoded_rows$truth, raw[[outcome]][91:120]),
      identical(decoded_rows$identity$.row, 91:120), identical(decoded_rows$outcome, outcome))
    foreign <- venn_abers_schema(data.frame(x = 1:6, other = as.double(1:6)), "other")
    error(venn_abers_decode(native, foreign))
    bad <- new; attr(bad, "vennabers_schema") <- foreign
    for (input in list(bad, new[-1L], setNames(new, c("unknown", names(new)[2L])))) {
      error(predict(fit, input)); error(venn_abers_details(fit, input))
    }
    # The shipped evaluation belongs to the example's own outcome name.
    replay <- workers && outcome == ".pred"
    saveRDS(list(fit = fit, new = new, schema = schema, snapshot = snapshot,
      expected = snapshot(fit, new, schema), replay = if (replay) list(model = case$model, outer = case$outer,
        evaluation = case$evaluation, schema = case$schema,
        new = case$assessment[setdiff(names(case$assessment), case$schema$encoded_outcome)])), saved)
    writeLines(c(
      paste0('.libPaths(', paste(deparse(.libPaths()), collapse = ''), ')'), 'library(vennabers)',
      paste0('s <- readRDS(', deparse(saved), ')'),
      'trace("fit.workflow", where = asNamespace("workflows"), tracer = quote(stop("Unexpected model fit")), print = FALSE)',
      'trace("prep.recipe", where = asNamespace("recipes"), tracer = quote(stop("Unexpected preprocessing")), print = FALSE)',
      'trace("venn_abers_regression", where = asNamespace("vennabers"), tracer = quote(stop("Unexpected calibration")), print = FALSE)',
      'stopifnot(identical(s$expected, s$snapshot(s$fit, s$new, s$schema)),',
      '  identical(names(workflows::extract_fit_engine(s$fit)$models), "Full1"))',
      if (replay) worker_replay("yardstick::metric_set(yardstick::rmse, yardstick::mae)", c("vennabers", "stats"),
        c(workflows = "fit.workflow", recipes = "prep.recipe", vennabers = "venn_abers_regression"))
    ), script)
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
    if (status != 0L) stop(paste(readLines(log), collapse = "\n"))
  }
})

  testthat::expect_true(TRUE)
})

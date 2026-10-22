testthat::test_that("frozen-processes", {
library(vennabers)

packages <- c("workflows", "parsnip", "recipes", "rsample", "hardhat", "rlang",
              "tibble", "withr", "generics", "tidyselect")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  run <- function(code) {
    script <- tempfile(fileext = ".R")
    log <- tempfile(fileext = ".log")
    on.exit(unlink(c(script, log)), add = TRUE)
    writeLines(c(paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"), code), script)
    status <- system2(file.path(R.home("bin"), "Rscript"),
                       c("--vanilla", shQuote(script)), stdout = log, stderr = log)
    if (status != 0L) stop(paste(readLines(log), collapse = "\n"))
  }
  post_available <- all(vapply(c("tailor", "probably"), requireNamespace, logical(1), quietly = TRUE))
  for (post in if (post_available) c(FALSE, TRUE) else FALSE) {
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    prefix <- c('library(vennabers)', paste0('post <- ', post))
    guards <- c(
      'trace("glm.fit", where = asNamespace("stats"), tracer = quote(stop("Unexpected model learning")), print = FALSE)',
      'trace("fit.workflow", where = asNamespace("workflows"), tracer = quote(stop("Unexpected workflow learning")), print = FALSE)',
      'trace("prep.recipe", where = asNamespace("recipes"), tracer = quote(stop("Unexpected recipe learning")), print = FALSE)',
      'if (post) trace("fit.tailor", where = asNamespace("tailor"), tracer = quote(stop("Unexpected postprocessor learning")), print = FALSE)'
    )
    outputs <- c(
      'outputs <- function(fit, new) lapply(c("log", "brier"), function(loss) list(prob = predict(fit, new, loss = loss), class = predict(fit, new, type = "class", loss = loss), details = venn_abers_details(fit, new, loss = loss)))'
    )
    run(c(prefix,
      'set.seed(830)',
      'd <- data.frame(x = rnorm(240), z = rnorm(240))',
      'd$truth <- factor(ifelse(runif(240) < plogis(d$x), "yes", "no"), levels = c("yes", "no"))',
      'd$x[c(1, 121, 201)] <- NA_real_',
      'base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ ., d[1:80, ]) |> recipes::step_impute_median(x) |> recipes::step_normalize(recipes::all_numeric_predictors())) |> workflows::add_model(parsnip::logistic_reg(engine = "glm"))',
      'if (post) base <- workflows::add_tailor(base, tailor::tailor() |> tailor::adjust_probability_calibration(method = "logistic", smooth = FALSE))',
      'frozen <- if (post) generics::fit(base, d[1:80, ], data_calibration = d[81:120, ]) else generics::fit(base, d[1:80, ])',
      guards,
      'fit <- calibrate_workflow(frozen, d[121:200, ], "truth", loss = "brier", event_level = "second", precision = 2)',
      'stopifnot(identical(fit$models[[1L]], frozen), all(c("R", "vennabers", "stats", "recipes", "workflows") %in% names(fit$versions)))',
      'if (post) stopifnot(all(c("tailor", "probably") %in% names(fit$versions)))',
      'stopifnot(!any(vapply(fit[setdiff(names(fit), "models")], is.data.frame, logical(1))))',
      'new <- d[201:240, c("x", "z")]', outputs,
      paste0('saveRDS(list(fit = fit, new = new, expected = outputs(fit, new), empty = outputs(fit, new[0, ]), single = outputs(fit, new[1, , drop = FALSE])), ', deparse(path), ')')
    ))
    run(c(prefix, guards, outputs,
      'trace("venn_abers", where = asNamespace("vennabers"), tracer = quote(stop("Unexpected calibration learning")), print = FALSE)',
      paste0('saved <- readRDS(', deparse(path), ')'),
      'stopifnot(identical(saved$expected, outputs(saved$fit, saved$new)), identical(saved$empty, outputs(saved$fit, saved$new[0, ])), identical(saved$single, outputs(saved$fit, saved$new[1, , drop = FALSE])))',
      'stopifnot(saved$fit$loss == "brier", saved$fit$event_class == "no", saved$fit$parameters$effective$precision == 2)'
    ))
    unlink(path)
  }
})

  testthat::expect_true(TRUE)
})

testthat::test_that("frozen-multiclass-processes", {
library(vennabers)

packages <- c("workflows", "parsnip", "recipes", "rsample", "hardhat", "rlang", "nnet",
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
  post_available <- all(vapply(c("tailor", "probably", "mgcv"), requireNamespace, logical(1), quietly = TRUE))
  for (method in c("ovo", "ovr", "regression")) for (k in if (method == "regression") 1L else 2:3) {
    for (post in if (post_available) c(FALSE, TRUE) else FALSE) {
      path <- tempfile(fileext = ".rds")
      on.exit(unlink(path), add = TRUE)
      prefix <- c('library(vennabers)', paste0('method <- ', deparse(method)),
                  paste0('k <- ', k), paste0('post <- ', post))
      guards <- c(
        'for (name in c("glm.fit", "lm.fit")) trace(name, where = asNamespace("stats"), tracer = quote(stop("Unexpected model learning")), print = FALSE)',
        'trace("multinom", where = asNamespace("nnet"), tracer = quote(stop("Unexpected multinomial learning")), print = FALSE)',
        'trace("fit.workflow", where = asNamespace("workflows"), tracer = quote(stop("Unexpected workflow learning")), print = FALSE)',
        'trace("prep.recipe", where = asNamespace("recipes"), tracer = quote(stop("Unexpected recipe learning")), print = FALSE)',
        'if (post) trace("fit.tailor", where = asNamespace("tailor"), tracer = quote(stop("Unexpected postprocessor learning")), print = FALSE)',
        'if (post) trace("gam", where = asNamespace("mgcv"), tracer = quote(stop("Unexpected smooth learning")), print = FALSE)'
      )
      outputs <- c(
        'outputs <- function(fit, new) if (method == "regression") list(numeric = predict(fit, new), details = venn_abers_details(fit, new)) else lapply(c("log", "brier"), function(loss) list(prob = predict(fit, new, loss = loss), class = predict(fit, new, type = "class", loss = loss), details = venn_abers_details(fit, new, loss = loss)))'
      )
      run(c(prefix,
        'set.seed(840)',
        'd <- data.frame(x = rnorm(360), z = rnorm(360))',
        'd$truth <- if (k == 1L) 2 + d$x + cos(1:360) else factor(rep(c("third", "first", "second")[seq_len(k)], length.out = 360), levels = c("third", "first", "second")[seq_len(k)])',
        'd$x[c(1, 141, 241, 341)] <- NA_real_',
        'spec <- if (k == 1L) parsnip::linear_reg(engine = "lm") else if (k == 2L) parsnip::logistic_reg(engine = "glm") else parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE)',
        'base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ ., d[1:120, ]) |> recipes::step_impute_median(x) |> recipes::step_normalize(recipes::all_numeric_predictors())) |> workflows::add_model(spec)',
        'if (post) base <- workflows::add_tailor(base, if (k == 1L) tailor::tailor() |> tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE) else tailor::tailor() |> tailor::adjust_probability_calibration(method = if (k == 2L) "logistic" else "multinomial", smooth = (k == 3L)))',
        'frozen <- if (post) generics::fit(base, d[1:120, ], data_calibration = d[121:220, ]) else generics::fit(base, d[1:120, ])',
        guards,
        'fit <- if (k == 1L) calibrate_workflow(frozen, d[221:320, ], "truth", epsilon = .1) else calibrate_workflow(frozen, d[221:320, ], "truth", classification = method, loss = "brier", precision = 2)',
        'stopifnot(identical(fit$models[[1L]], frozen), all(c("R", "vennabers", "recipes", "workflows") %in% names(fit$versions)))',
        'if (post) stopifnot(all(c("tailor", "probably") %in% names(fit$versions)))',
        'stopifnot(!any(vapply(fit[setdiff(names(fit), "models")], is.data.frame, logical(1))), is.null(fit$folds))',
        'new <- d[321:360, c("x", "z")]', outputs,
        paste0('saveRDS(list(fit = fit, new = new, expected = outputs(fit, new), empty = outputs(fit, new[0, ]), single = outputs(fit, new[1, , drop = FALSE]), parameters = fit$parameters, classes = fit$class_levels), ', deparse(path), ')')
      ))
      run(c(prefix, guards, outputs,
        'for (name in c("venn_abers", "venn_abers_regression")) trace(name, where = asNamespace("vennabers"), tracer = quote(stop("Unexpected calibration learning")), print = FALSE)',
        paste0('saved <- readRDS(', deparse(path), ')'),
        'stopifnot(identical(saved$expected, outputs(saved$fit, saved$new)), identical(saved$empty, outputs(saved$fit, saved$new[0, ])), identical(saved$single, outputs(saved$fit, saved$new[1, , drop = FALSE])), identical(saved$parameters, saved$fit$parameters), identical(saved$classes, saved$fit$class_levels))',
        'if (k == 1L) stopifnot(saved$fit$parameters$requested$epsilon == .1, saved$fit$parameters$effective$m == 5, saved$fit$parameters$effective$epsilon == 10/101) else stopifnot(saved$fit$loss == "brier", saved$fit$parameters$effective$precision == 2, identical(saved$fit$class_levels, c("third", "first", "second")[seq_len(k)]))'
      ))
      unlink(path)
    }
  }
})

  testthat::expect_true(TRUE)
})

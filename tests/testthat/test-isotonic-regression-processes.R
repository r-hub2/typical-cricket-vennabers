testthat::test_that("isotonic rejects incompatible methods and undersized P before fitting", {
  skip_if_no_integration()
  d <- data.frame(x = 1:24, truth = 1:24)
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  va <- function(data) rsample::make_splits(list(analysis = 1:16, assessment = 17:24), data)
  post <- list(method = "isotonic_numeric", strategy = "holdout", splits = function(data)
    rsample::make_splits(list(analysis = 1:15, assessment = 16L), data))
  spec <- calibrated_workflow(base, "truth", splits = va, postprocessing = post)
  testthat::expect_error(generics::fit(spec, d), "isotonic_numeric B/P plan:.*at least 2")
  classified <- workflows::workflow() |> workflows::add_formula(group ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  testthat::expect_error(generics::fit(calibrated_workflow(classified, "group",
    postprocessing = list(method = "isotonic_numeric", strategy = "holdout", prop = .3)),
    transform(d, group = factor(rep(c("a", "b"), 12)))), "isotonic_numeric.*model classes require")
  if (requireNamespace("tailor", quietly = TRUE) && requireNamespace("probably", quietly = TRUE)) {
    with_post <- workflows::add_tailor(base, tailor::tailor() |>
      tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE))
    testthat::expect_error(calibrated_workflow(with_post, "truth",
      postprocessing = post), "isotonic_numeric requires a base workflow without a postprocessor")
  }
})

testthat::test_that("isotonic B, P and C learning roles remain separate", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:60), truth = 10 + 3 * sin(1:60) + cos(1:60))
  base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, d) |>
    recipes::step_normalize(x)) |> workflows::add_model(parsnip::linear_reg(engine = "lm"))
  va <- function(data) rsample::make_splits(list(analysis = 1:40, assessment = 41:60), data)
  post <- list(method = "isotonic_numeric", strategy = "holdout", splits = function(data)
    rsample::make_splits(list(analysis = 1:25, assessment = 26:40), data))
  spec <- calibrated_workflow(base, "truth", splits = va, postprocessing = post)
  original <- workflows::extract_fit_engine(generics::fit(spec, d))
  base_fit <- function(state) state$models[[1L]]$workflow
  for (role in c("B", "P", "C")) {
    changed <- d
    changed$truth[original$postprocessing[[1L]][[role]]] <- changed$truth[original$postprocessing[[1L]][[role]]] + 3
    current <- workflows::extract_fit_engine(generics::fit(spec, changed))
    if (role == "B") {
      testthat::expect_false(identical(coef(workflows::extract_fit_engine(base_fit(current))),
                                       coef(workflows::extract_fit_engine(base_fit(original)))))
    } else {
      testthat::expect_identical(coef(workflows::extract_fit_engine(base_fit(current))),
                                 coef(workflows::extract_fit_engine(base_fit(original))))
      testthat::expect_identical(workflows::extract_recipe(base_fit(current))$steps,
                                 workflows::extract_recipe(base_fit(original))$steps)
    }
    if (role == "P") testthat::expect_false(identical(current$models[[1L]]$isotonic,
                                                     original$models[[1L]]$isotonic)) else if (role == "C") {
      testthat::expect_identical(current$models[[1L]]$isotonic, original$models[[1L]]$isotonic)
      testthat::expect_false(identical(current$calibrator, original$calibrator))
    }
  }
  changed <- d; changed$x[1:25] <- changed$x[1:25] + 5
  current <- workflows::extract_fit_engine(generics::fit(spec, changed))
  testthat::expect_false(identical(workflows::extract_recipe(base_fit(current))$steps,
                                   workflows::extract_recipe(base_fit(original))$steps))
  bad <- d; bad$truth[30L] <- Inf
  testthat::expect_error(generics::fit(spec, bad), "regression truth")
})

testthat::test_that("isotonic installed example and RDS predict in fresh R without learning", {
  skip_if_no_integration("tune", "yardstick")
  example <- new.env()
  sys.source(system.file("examples", "isotonic-regression-holdout.R", package = "vennabers"), example)
  testthat::expect_length(example$cases, 2L)
  snapshot <- function(fit, new, schema) lapply(c(0L, 1L, 7L), function(n) {
    batch <- new[seq_len(n), , drop = FALSE]
    native <- list(prediction = predict(fit, batch), details = venn_abers_details(fit, batch))
    list(native = native, decoded = lapply(native, venn_abers_decode, schema = schema))
  })
  cases <- lapply(example$cases, function(case) {
    case$expected <- snapshot(case$fit, case$new, case$schema)
    case
  })
  calls <- new.env(); calls$n <- 0L
  original <- getS3method("predict", "workflow", envir = asNamespace("workflows"))
  counted <- function(object, ...) {
    calls$n <- calls$n + 1L
    original(object, ...)
  }
  registerS3method("predict", "workflow", counted, envir = asNamespace("stats"))
  on.exit(registerS3method("predict", "workflow", original, envir = asNamespace("stats")), add = TRUE)
  for (case in cases) for (n in c(7L, 1L, 0L)) {
    calls$n <- 0L
    venn_abers_details(case$fit, case$new[seq_len(n), , drop = FALSE])
    state <- workflows::extract_fit_engine(case$fit)
    testthat::expect_identical(calls$n, if (n) length(state$models) else 0L)
  }
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  saveRDS(cases, saved)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "library(vennabers); loadNamespace('workflows'); loadNamespace('recipes')",
    "guard <- function(...) stop('Unexpected learning')",
    "registerS3method('fit', 'workflow', guard, envir=asNamespace('generics'))",
    "trace('prep.recipe', where=asNamespace('recipes'), tracer=quote(stop('Unexpected learning')), print=FALSE)",
    paste0("snapshot <- ", paste(deparse(snapshot), collapse = "\n")),
    paste0("cases <- readRDS(", deparse(saved), ")"),
    "for (case in cases) stopifnot(identical(case$expected, snapshot(case$fit, case$new, case$schema)))"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

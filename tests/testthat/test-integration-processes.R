testthat::test_that("integration-processes", {
# Each process uses an installed package; startup files cannot mask load order.
package_root <- find.package("vennabers")
package_library <- dirname(package_root)
if (file.exists(file.path(package_root, "R", "venn_abers.R"))) {
  package_library <- withr::local_tempdir("vennabers-library-")
  status <- system2(file.path(R.home("bin"), "R"),
                    c("CMD", "INSTALL", paste0("--library=", package_library), package_root))
  if (status != 0L) stop("Could not install vennabers for process tests.", call. = FALSE)
}
run_process <- function(code, libraries = .libPaths()) {
  script <- tempfile(fileext = ".R")
  output <- tempfile(fileext = ".log")
  on.exit(unlink(c(script, output)), add = TRUE)
  prefix <- c(paste0(".libPaths(", paste(deparse(libraries), collapse = ""), ")"),
              paste0("package_library <- ", deparse(package_library)))
  writeLines(c(prefix, code), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = output, stderr = output)
  if (status != 0L) stop(paste(readLines(output), collapse = "\n"))
}
run_process(c(
  'stopifnot(!"parsnip" %in% loadedNamespaces())',
  'loadNamespace("vennabers", lib.loc = package_library)',
  'stopifnot(!"parsnip" %in% loadedNamespaces())',
  'stopifnot(nrow(predict(vennabers::venn_abers(c(FALSE, TRUE), c(.1, .9)), .5)) == 1L)',
  'stopifnot(is.double(predict(vennabers::venn_abers_regression(1:3, 1:3), 2)$prediction))',
  'stopifnot(!"parsnip" %in% loadedNamespaces())'
), libraries = .Library)

for (method in c("inductive", "cross")) {
if (all(vapply(c("parsnip", "workflows", "rsample", "recipes", "rlang", "tibble", "withr", "generics", "tidyselect"),
               requireNamespace, logical(1), quietly = TRUE))) {
  construction <- c(
    paste0("method <- ", deparse(method)),
    'data <- data.frame(x = sin(1:80), truth = factor(rep(c("yes", "no", "no", "yes"), 20)))',
    'base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |> workflows::add_model(parsnip::logistic_reg(engine = "glm"))',
    'model <- vennabers::calibrated_workflow(base, "truth", method = method, seed = 42)',
    'second <- vennabers::calibrated_workflow(base, "truth", method = method, seed = 42)',
    'stopifnot(nrow(parsnip::get_fit("vennabers")) == 4L)',
    'fit <- generics::fit(model, data)',
    'stopifnot(identical(predict(fit, data, type = "prob"), vennabers::venn_abers_details(fit, data)$prediction))'
  )
  for (order in list(c('loadNamespace("vennabers")', 'loadNamespace("parsnip")'),
                    c('loadNamespace("parsnip")', 'loadNamespace("vennabers")'))) {
    run_process(c(order, construction, 'stopifnot(!"package:vennabers" %in% search(), !"package:parsnip" %in% search())'))
  }
  payload <- tempfile(fileext = ".rds")
  run_process(c(construction, paste0('saveRDS(list(fit = fit, data = data, prob = predict(fit, data, type = "prob"), details = vennabers::venn_abers_details(fit, data, loss = "brier")), ', deparse(payload), ')')))
  run_process(c('loadNamespace("vennabers")', 'loadNamespace("workflows")', paste0('saved <- readRDS(', deparse(payload), ')'),
                'stopifnot(identical(saved$prob, predict(saved$fit, saved$data, type = "prob")), identical(saved$details, vennabers::venn_abers_details(saved$fit, saved$data, loss = "brier")))'))
  unlink(payload)
  run_process(c(
    'loadNamespace("parsnip")',
    'parsnip::set_new_model("vennabers")',
    'before <- parsnip::get_fit("vennabers")',
    'error <- tryCatch(loadNamespace("vennabers"), error = identity)',
    'stopifnot(inherits(error, "error"), grepl("foreign", conditionMessage(error)), identical(before, parsnip::get_fit("vennabers")))'
  ))
}

if (all(vapply(c("tune", "yardstick", "future", "rsample", "workflows"),
               requireNamespace, logical(1), quietly = TRUE))) {
  payload <- tempfile(fileext = ".rds")
  run_process(c(
    paste0("method <- ", deparse(method)),
    'loadNamespace("vennabers")',
    'data <- data.frame(x = sin(1:120), truth = factor(rep(c("yes", "no", "no", "yes"), 30)))',
    'base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |> workflows::add_model(parsnip::logistic_reg(engine = "glm"))',
    'outer <- withr::with_seed(101, rsample::vfold_cv(data, v = 3, repeats = 2))',
    'prepared <- withr::with_seed(102, lapply(outer$splits, function(s) if (method == "cross") rsample::vfold_cv(rsample::analysis(s), v = 5) else rsample::initial_split(rsample::analysis(s), .8)))',
    'resolver <- vennabers::local_calibration_splits(prepared, names(data))',
    'model <- vennabers::calibrated_workflow(base, "truth", method = method, splits = resolver, seed = 42, loss = "brier")',
    paste0('saveRDS(list(model = model, outer = outer), ', deparse(payload), ')')
  ))
  run_process(c(
    'loadNamespace("vennabers")', 'loadNamespace("stats")',
    paste0('saved <- readRDS(', deparse(payload), ')'),
    'Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")',
    'future::plan(future::multisession, workers = 2)',
    'control <- tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, event_level = "first", parallel_over = "resamples", pkgs = c("vennabers", "stats"))',
    'run <- function(outer) withr::with_seed(103, tune::fit_resamples(saved$model, outer, metrics = yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss), control = control))',
    'one <- run(saved$outer)',
    'two <- run(saved$outer)',
    'stopifnot(nrow(tune::collect_notes(one)) == 0L, nrow(tune::collect_notes(two)) == 0L, identical(tune::collect_predictions(one), tune::collect_predictions(two)), identical(tune::collect_metrics(one), tune::collect_metrics(two)))',
    'reversed <- rsample::manual_rset(rev(saved$outer$splits), tibble::as_tibble(saved$outer[nrow(saved$outer):1, c("id", "id2")]))',
    'three <- run(reversed)',
    'stopifnot(nrow(tune::collect_notes(three)) == 0L)',
    'normalize <- function(r) { p <- as.data.frame(tune::collect_predictions(r, summarize = FALSE)); p <- p[order(p$id, p$id2, p$.row), ]; rownames(p) <- NULL; p }',
    'stopifnot(identical(normalize(one), normalize(three)))',
    'states <- tune::collect_extracts(one)$.extracts',
    'stopifnot(length(states) == 6L, all(vapply(states, inherits, logical(1), "vennabers_fit")))',
    'future::plan(future::sequential)'
  ))
  unlink(payload)
}
}

if (all(vapply(c("parsnip", "workflows", "rsample", "nnet"),
               requireNamespace, logical(1), quietly = TRUE))) {
  for (method in c("inductive", "cross")) for (classification in c("ovo", "ovr")) for (k in 2:4) {
    payload <- tempfile(fileext = ".rds")
    run_process(c(
      'loadNamespace("vennabers")',
      paste0('classification <- ', deparse(classification)), paste0('k <- ', k),
      paste0('method <- ', deparse(method)),
      'data <- data.frame(x = sin(1:120), truth = factor(rep(seq_len(k), length.out = 120)))',
      # parsnip/nnet returns a malformed probability table for two classes.
      # A binomial base model still exercises the explicit VA matrix contract.
      'spec <- if (k == 2L) parsnip::logistic_reg(engine = "glm") else parsnip::multinom_reg(engine = "nnet")',
      'base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |> workflows::add_model(spec)',
      'fit <- generics::fit(vennabers::calibrated_workflow(base, "truth", classification = classification, method = method, seed = 32, loss = "brier", precision = 2), data)',
      'new <- data[1:4, , drop = FALSE]',
      'outputs <- function(fit, new) lapply(c("log", "brier"), function(loss) list(prob = predict(fit, new, type = "prob", loss = loss), class = predict(fit, new, type = "class", loss = loss), details = vennabers::venn_abers_details(fit, new, loss = loss)))',
      paste0('saveRDS(list(fit = fit, new = new, outputs = outputs, expected = outputs(fit, new), empty = outputs(fit, new[0, ])), ', deparse(payload), ')')
    ))
    run_process(c('loadNamespace("vennabers")', 'loadNamespace("workflows")', 'loadNamespace("nnet")',
      paste0('saved <- readRDS(', deparse(payload), ')'),
      'stopifnot(identical(saved$expected, saved$outputs(saved$fit, saved$new)), identical(saved$empty, saved$outputs(saved$fit, saved$new[0, ])))'))
    unlink(payload)
  }
}

if (all(vapply(c("parsnip", "workflows", "rsample", "nnet", "tune", "yardstick", "future"),
               requireNamespace, logical(1), quietly = TRUE))) {
  for (classification in c("ovo", "ovr")) {
    run_process(c(
      'loadNamespace("vennabers")',
      paste0('classification <- ', deparse(classification)),
      'data <- iris; names(data)[5L] <- "truth"',
      'base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |> workflows::add_model(parsnip::multinom_reg(engine = "nnet"))',
      'outer <- withr::with_seed(331, rsample::vfold_cv(data, v = 3, repeats = 2, strata = truth))',
      'prepared <- withr::with_seed(332, lapply(outer$splits, function(s) rsample::vfold_cv(rsample::analysis(s), v = 5, strata = truth)))',
      'resolver <- vennabers::local_calibration_splits(prepared, names(data))',
      'Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")',
      'future::plan(future::multisession, workers = 2)',
      'control <- tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, parallel_over = "resamples", pkgs = c("vennabers", "nnet"))',
      'for (splits in list(NULL, resolver)) {',
      'model <- vennabers::calibrated_workflow(base, "truth", method = "cross", classification = classification, splits = splits, seed = 33, loss = "brier")',
      'run <- function(outer) withr::with_seed(333, tune::fit_resamples(model, outer, metrics = yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss), control = control))',
      'one <- run(outer); two <- run(outer)',
      'stopifnot(nrow(tune::collect_notes(one)) == 0L, nrow(tune::collect_notes(two)) == 0L, identical(tune::collect_predictions(one), tune::collect_predictions(two)), identical(tune::collect_metrics(one), tune::collect_metrics(two)))',
      'reversed <- rsample::manual_rset(rev(outer$splits), tibble::as_tibble(outer[nrow(outer):1, c("id", "id2")]))',
      'three <- run(reversed)',
      'normalize <- function(r) { p <- as.data.frame(tune::collect_predictions(r, summarize = FALSE)); p <- p[order(p$id, p$id2, p$.row), ]; rownames(p) <- NULL; p }',
      'stopifnot(nrow(tune::collect_notes(three)) == 0L, identical(normalize(one), normalize(three)), nrow(normalize(one)) == 2L * nrow(data))',
      'states <- tune::collect_extracts(one)$.extracts',
      'stopifnot(length(states) == 6L, all(vapply(states, function(s) length(s$models), integer(1)) == 5L))',
      '}',
      'future::plan(future::sequential)'
    ))
  }
}

if (all(vapply(c("parsnip", "workflows", "recipes", "tune", "yardstick", "future", "rsample"),
               requireNamespace, logical(1), quietly = TRUE))) {
  for (method in c("inductive", "cross")) {
    construction <- c(
      paste0("method <- ", deparse(method)),
      'data <- data.frame(x = c(NA, sin(2:80)), truth = 10 + 3*sin(1:80) + cos(1:80))',
      'base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, data) |> recipes::step_impute_median(x) |> recipes::step_normalize(x)) |> workflows::add_model(parsnip::linear_reg(engine = "lm"))',
      'model <- vennabers::calibrated_workflow(base, "truth", method = method, seed = 42, epsilon = .3)',
      'fit <- generics::fit(model, data)',
      'stopifnot(identical(predict(fit, data, type = "numeric"), vennabers::venn_abers_details(fit, data)$prediction))'
    )
    for (order in list(c('loadNamespace("vennabers")', 'loadNamespace("parsnip")'),
                      c('loadNamespace("parsnip")', 'loadNamespace("vennabers")'))) {
      run_process(c(order, construction))
    }
    payload <- tempfile(fileext = ".rds")
    run_process(c(construction,
      paste0('saveRDS(list(fit = fit, data = data, prediction = predict(fit, data), details = vennabers::venn_abers_details(fit, data)), ', deparse(payload), ')')))
    run_process(c('loadNamespace("vennabers")', 'loadNamespace("workflows")',
      paste0('saved <- readRDS(', deparse(payload), ')'),
      'stopifnot(identical(saved$prediction, predict(saved$fit, saved$data)), identical(saved$details, vennabers::venn_abers_details(saved$fit, saved$data)))',
      'for (n in 0:1) stopifnot(identical(predict(saved$fit, saved$data[seq_len(n), ]), vennabers::venn_abers_details(saved$fit, saved$data[seq_len(n), ])$prediction))'))
    unlink(payload)
    run_process(c(construction,
      'Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")',
      'outer <- withr::with_seed(301, rsample::vfold_cv(data, v = 3, repeats = 2))',
      'prepared <- withr::with_seed(302, lapply(outer$splits, function(s) if (method == "cross") rsample::vfold_cv(rsample::analysis(s), v = 5) else rsample::initial_split(rsample::analysis(s), .8)))',
      'resolver <- vennabers::local_calibration_splits(prepared, names(data))',
      'future::plan(future::multisession, workers = 2)',
      'for (splits in list(NULL, resolver)) {',
      'model <- vennabers::calibrated_workflow(base, "truth", method = method, splits = splits, seed = 42, epsilon = .3)',
      'control <- tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, parallel_over = "resamples", pkgs = c("vennabers", "stats"))',
      'run <- function(outer) withr::with_seed(303, tune::fit_resamples(model, outer, metrics = yardstick::metric_set(yardstick::rmse, yardstick::mae), control = control))',
      'one <- run(outer); two <- run(outer)',
      'stopifnot(nrow(tune::collect_notes(one)) == 0L, nrow(tune::collect_notes(two)) == 0L, identical(tune::collect_predictions(one), tune::collect_predictions(two)), identical(tune::collect_metrics(one), tune::collect_metrics(two)), nrow(tune::collect_predictions(one)) == 2L*nrow(data))',
      'reversed <- rsample::manual_rset(rev(outer$splits), tibble::as_tibble(outer[nrow(outer):1, c("id", "id2")]))',
      'three <- run(reversed)',
      'normalize <- function(r) { p <- as.data.frame(tune::collect_predictions(r, summarize = FALSE)); p <- p[order(p$id, p$id2, p$.row), ]; rownames(p) <- NULL; p }',
      'stopifnot(nrow(tune::collect_notes(three)) == 0L, identical(normalize(one), normalize(three)))',
      'states <- tune::collect_extracts(one)$.extracts',
      'stopifnot(length(states) == 6L, all(vapply(states, function(s) length(s$models), integer(1)) == if (method == "cross") 5L else 1L))',
      '}',
      'future::plan(future::sequential)'
    ))
  }
}

  testthat::expect_true(TRUE)
})

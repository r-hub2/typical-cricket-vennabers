# The 16 fixed regression cells, using independently fitted public workflows
# and the min-max characterization of increasing least-squares isotonic fits.
regression_composition_data <- function() data.frame(x = sin(1:96), z = cos((1:96) / 4),
  truth = 10 + 3 * sin(1:96) + cos(1:96) + sin((1:96) / 4))
regression_composition_base <- function(d) workflows::workflow() |>
  workflows::add_recipe(recipes::recipe(truth ~ ., d) |>
    recipes::step_normalize(recipes::all_numeric_predictors())) |>
  workflows::add_model(parsnip::linear_reg(engine = "lm"))
regression_composition_cells <- expand.grid(mode = c("matched", "full"), method = c("inductive", "cross"),
  strategy = c("holdout", "ensemble"), post = c("tailor_linear", "isotonic_numeric"), stringsAsFactors = FALSE)
regression_composition_spec <- function(base, cell, splits = NULL) {
  post <- list(method = cell$post, strategy = cell$strategy)
  if (cell$strategy == "holdout") post$prop <- .3 else post$v <- 3L
  args <- list(base, "truth", method = cell$method, prediction_model = cell$mode,
    seed = 24, m = 2, postprocessing = post)
  if (!is.null(splits)) args$splits <- splits else if (cell$method == "cross") args$v <- 2L
  do.call(calibrated_workflow, args)
}

testthat::test_that("16 numeric compositions preserve public postprocessors, VA layers and actual costs", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes")) testthat::skip_if_not_installed(pkg)
  audit <- new.env(); audit$fits <- audit$posts <- audit$predictions <- audit$tailor <- 0L
  original_fit <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  original_predict <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  original_tailor <- getS3method("predict", "tailor", envir = asNamespace("stats"))
  external <- function(object) !inherits(workflows::extract_spec_parsnip(object), "vennabers")
  registerS3method("fit", "workflow", function(object, ...) {
    if (external(object)) audit$fits <- audit$fits + 1L
    original_fit(object, ...)
  }, envir = asNamespace("generics"))
  registerS3method("predict", "workflow", function(object, ...) {
    if (external(object)) audit$predictions <- audit$predictions + 1L
    original_predict(object, ...)
  }, envir = asNamespace("stats"))
  registerS3method("predict", "tailor", function(object, ...) {
    audit$tailor <- audit$tailor + 1L
    original_tailor(object, ...)
  }, envir = asNamespace("stats"))
  on.exit({
    registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics"))
    registerS3method("predict", "workflow", original_predict, envir = asNamespace("stats"))
    registerS3method("predict", "tailor", original_tailor, envir = asNamespace("stats"))
  }, add = TRUE)
  linear <- probably::cal_estimate_linear; isotonic <- stats::isoreg
  testthat::local_mocked_bindings(cal_estimate_linear = function(...) {
    audit$posts <- audit$posts + 1L; linear(...)
  }, .package = "probably")
  testthat::local_mocked_bindings(isoreg = function(...) {
    audit$posts <- audit$posts + 1L; isotonic(...)
  }, .package = "stats")
  close <- function(actual, expected, s, info) {
    testthat::expect_identical(dim(actual), dim(expected), info = info)
    testthat::expect_true(all(is.finite(actual)) && all(is.finite(expected)), info = info)
    testthat::expect_true(all(abs(actual - expected) <= 1e-10 * s + 1e-10 * abs(expected)), info = info)
  }
  for (s in c(1e-9, 1, 1e9)) {
    d <- regression_composition_data(); d$truth <- s * d$truth
    base <- regression_composition_base(d); new <- d[1:7, c("x", "z")]
    public <- function(p, data, method) {
      if (method == "tailor_linear") {
        ref <- generics::fit(workflows::add_tailor(base, tailor::tailor() |>
          tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE)),
          d[p$B, ], data_calibration = d[p$P, ])
        return(predict(ref, data)$.pred)
      }
      ref <- generics::fit(base, d[p$B, ])
      input <- predict(ref, d[p$P, ])$.pred; knots <- sort(unique(input))
      # Tied observations retain their original multiplicity in each interval.
      values <- vapply(seq_along(knots), function(k) max(vapply(seq_len(k), function(a)
        min(vapply(seq.int(k, length(knots)), function(b)
          mean(d$truth[p$P][input >= knots[a] & input <= knots[b]]), numeric(1))), numeric(1))), numeric(1))
      x <- predict(ref, data)$.pred
      if (length(knots) == 1L) rep(values, length(x)) else stats::approx(knots, values, xout = x, rule = 2)$y
    }
    for (i in seq_len(nrow(regression_composition_cells))) {
      cell <- regression_composition_cells[i, ]; info <- paste(s, paste(cell, collapse = " "))
      full <- cell$mode == "full"; cross <- cell$method == "cross"
      # K VA sources; J B/P pairs in each source (one for holdout).
      K <- if (cross) 2L else 1L; J <- if (cell$strategy == "ensemble") 3L else 1L
      cost <- (K + full) * J; retained <- if (full) J else K * J
      audit$fits <- audit$posts <- 0L
      fit <- generics::fit(regression_composition_spec(base, cell), d)
      state <- workflows::extract_fit_engine(fit)
      testthat::expect_identical(c(audit$fits, audit$posts), rep(cost, 2L), info = info)
      testthat::expect_output(print(state), paste("Base fits:", cost))
      testthat::expect_output(print(state), paste("Models per nonempty prediction:", retained))
      plans <- function(id) if (J > 1L) state$postprocessing[[id]]$inner_folds else list(state$postprocessing[[id]])
      leaves <- function(id) if (J > 1L) state$models[[id]]$models else list(state$models[[id]])
      sources <- names(state$calibration)
      testthat::expect_named(state$models, if (full) "Full1" else sources)
      for (id in names(state$models)) {
        ps <- plans(id); models <- leaves(id)
        testthat::expect_length(models, J)
        reference <- lapply(ps, public, data = new, method = cell$post)
        for (j in seq_len(J)) {
          base_reference <- generics::fit(base, d[ps[[j]]$B, ])
          for (batch in list(d[ps[[j]]$P, ], new))
            close(predict(models[[j]]$workflow, batch)$.pred, predict(base_reference, batch)$.pred, s, info)
          testthat::expect_identical(ps[[j]]$method, cell$post)
          testthat::expect_identical(ps[[j]]$rank, if (cell$post == "tailor_linear") 2L else NULL)
          actual <- predict(models[[j]], new, type = "numeric")$.pred
          close(actual, reference[[j]], s, info)
          if (cell$post == "isotonic_numeric") close(actual, reference[[j]], 1, info)
        }
        close(predict(state$models[[id]], new)$.pred, rowMeans(do.call(cbind, reference)), s, info)
      }
      for (id in sources) {
        ps <- plans(id); C <- ps[[1L]]$C
        expected <- rowMeans(do.call(cbind, lapply(ps, public, data = d[C, ], method = cell$post)))
        close(state$calibration[[id]]$prediction, expected, s, info)
      }
      # Freeze exact composition inputs at the VA layer: last-bit mean changes
      # can move exact isotonic ties independently of the mean-layer tolerance.
      held <- unlist(lapply(state$calibration, `[[`, "prediction"), use.names = FALSE)
      truth <- unlist(lapply(state$calibration, `[[`, "truth"), use.names = FALSE)
      ids <- rep(sources, vapply(state$calibration, function(record) length(record$row), 1L))
      core <- if (cross) cross_venn_abers_regression(truth, held, ids, m = 2) else
        venn_abers_regression(truth, held, m = 2)
      input <- vapply(sources, function(id) predict(state$models[[state$sources[[id]]]], new)$.pred, numeric(nrow(new)))
      expected <- if (cross) predict(core, input, type = "details") else predict(core, input[, 1L])
      actual <- venn_abers_details(fit, new)
      if (cross) {
        close(actual$prediction$.pred, expected$prediction$prediction, s, info)
        for (column in c("lower", "upper")) close(actual$aggregate[[column]], expected$prediction[[column]], s, info)
        for (id in sources) for (column in c("prediction", "lower", "upper"))
          close(actual$details[actual$details$fold == id, column][[1L]], expected$folds[[id]][[column]], s, info)
      } else close(as.matrix(actual$details[c("prediction", "lower", "upper")]), as.matrix(expected), s, info)
      testthat::expect_true(all(actual$details$postprocessing_input == "numeric"))
      testthat::expect_null(state$class_levels)
      audit$fits <- audit$posts <- 0L
      for (n in c(0L, 1L, 7L)) {
        audit$predictions <- audit$tailor <- 0L
        small <- venn_abers_details(fit, new[seq_len(n), ])
        testthat::expect_identical(audit$predictions, if (n) retained else 0L, info = info)
        testthat::expect_identical(audit$tailor, if (n && cell$post == "tailor_linear") retained else 0L, info = info)
        testthat::expect_named(small$prediction, ".pred")
        testthat::expect_identical(nrow(small$prediction), n)
        testthat::expect_identical(names(small$details), names(actual$details))
        close(small$prediction$.pred, actual$prediction$.pred[seq_len(n)], s, info)
      }
      testthat::expect_identical(c(audit$fits, audit$posts), c(0L, 0L))
      for (type in c("conf_int", "pred_int")) testthat::expect_error(predict(fit, new, type = type), "type")
    }
  }
})

testthat::test_that("full numeric compositions relearn F from C while excluding C from its source learning", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes")) testthat::skip_if_not_installed(pkg)
  d <- regression_composition_data(); base <- regression_composition_base(d)
  held <- list(73:96, 1:36, 37:72)
  learned <- function(state) {
    model <- state$models$Full1
    models <- if (inherits(model, "vennabers_ensemble")) model$models else list(model)
    lapply(models, function(m) {
      numeric_isotonic <- inherits(m, "vennabers_isotonic_workflow")
      list(engine = coef(workflows::extract_fit_engine(if (numeric_isotonic) m$workflow else m)),
        post = if (numeric_isotonic) m$isotonic else coef(workflows::extract_tailor(m)$adjustments[[1]]$results$fit$estimates[[1]]$estimate))
    })
  }
  cells <- subset(regression_composition_cells, mode == "full")
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]; info <- paste(cell, collapse = " ")
    va <- function(data) {
      splits <- lapply(held, function(C) rsample::make_splits(list(analysis = setdiff(1:96, C), assessment = C), data))
      if (cell$method == "cross") rsample::manual_rset(splits, c("one", "two", "three")) else splits[[1L]]
    }
    spec <- regression_composition_spec(base, cell, splits = va)
    fit <- function(data) workflows::extract_fit_engine(generics::fit(spec, data))
    original <- fit(d); source <- names(original$calibration)[1L]
    for (column in c("x", "truth")) {
      changed <- d
      changed[held[[1L]], column] <- if (column == "x") changed[held[[1L]], column] * 2 + .3 else changed[held[[1L]], column] + 3
      after <- fit(changed)
      testthat::expect_identical(after$postprocessing[[source]], original$postprocessing[[source]], info = info)
      testthat::expect_false(identical(learned(after), learned(original)), info = info)
      testthat::expect_identical(identical(after$calibration[[source]]$prediction, original$calibration[[source]]$prediction),
        column == "truth", info = info)
      if (column == "truth") testthat::expect_false(identical(after$calibrator, original$calibrator), info = info)
      if (cell$method == "cross") for (other in names(original$calibration)[-1L])
        testthat::expect_false(identical(after$calibration[[other]]$prediction, original$calibration[[other]]$prediction), info = info)
    }
  }
})

testthat::test_that("isotonic holdouts keep predictor and truth mutations within B P and C roles", {
  for (pkg in c("workflows", "rsample", "recipes")) testthat::skip_if_not_installed(pkg)
  d <- regression_composition_data(); base <- regression_composition_base(d)
  for (method in c("inductive", "cross")) {
    va <- function(data) {
      splits <- lapply(list(1:32, 33:64, 65:96), function(C)
        rsample::make_splits(list(analysis = setdiff(1:96, C), assessment = C), data))
      if (method == "cross") rsample::manual_rset(splits, c("one", "two", "three")) else splits[[1L]]
    }
    spec <- regression_composition_spec(base, list(method = method, mode = "matched", strategy = "holdout", post = "isotonic_numeric"), va)
    fit <- function(data) workflows::extract_fit_engine(generics::fit(spec, data))
    original <- fit(d); model <- original$models[[1L]]; plan <- original$postprocessing[[1L]]
    for (role in c("B", "P", "C")) for (column in c("x", "truth")) {
      changed <- d; rows <- plan[[role]]; info <- paste(method, role, column)
      changed[rows, column] <- if (column == "x") changed[rows, column]^2 + .5 else changed[rows, column] + 3
      after <- fit(changed); other <- after$models[[1L]]
      testthat::expect_identical(identical(coef(workflows::extract_fit_engine(model$workflow)),
        coef(workflows::extract_fit_engine(other$workflow))), role != "B", info = info)
      if (role != "B") testthat::expect_identical(workflows::extract_recipe(model$workflow)$steps,
        workflows::extract_recipe(other$workflow)$steps, info = info)
      if (role == "P") testthat::expect_false(identical(model$isotonic, other$isotonic), info = info)
      if (role == "C") {
        testthat::expect_identical(model$isotonic, other$isotonic, info = info)
        testthat::expect_identical(identical(original$calibration[[1L]]$prediction, after$calibration[[1L]]$prediction),
          column == "truth", info = info)
        if (column == "truth") testthat::expect_false(identical(original$calibrator, after$calibrator), info = info)
      }
    }
  }
})

testthat::test_that("final numeric postprocessor failures abort after successful sources without fallback", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes")) testthat::skip_if_not_installed(pkg)
  d <- regression_composition_data(); base <- regression_composition_base(d)
  linear <- probably::cal_estimate_linear; isotonic <- stats::isoreg
  calls <- 0L; fail <- 0L
  testthat::local_mocked_bindings(cal_estimate_linear = function(...) {
    calls <<- calls + 1L
    if (calls == fail) stop("final numeric postprocessor failure")
    linear(...)
  }, .package = "probably")
  testthat::local_mocked_bindings(isoreg = function(...) {
    calls <<- calls + 1L
    if (calls == fail) stop("final numeric postprocessor failure")
    isotonic(...)
  }, .package = "stats")
  cells <- subset(regression_composition_cells, mode == "full")
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]; calls <- 0L
    K <- if (cell$method == "cross") 2L else 1L; J <- if (cell$strategy == "ensemble") 3L else 1L
    fail <- K * J + 1L
    spec <- regression_composition_spec(base, cell)
    if (cell$post == "tailor_linear") {
      testthat::expect_warning(testthat::expect_error(generics::fit(spec, d), "Full1 / .*no learned linear state"),
        "Full1 / .*final numeric postprocessor failure")
    } else testthat::expect_error(generics::fit(spec, d), "Full1 / .*final numeric postprocessor failure")
    testthat::expect_identical(calls, fail)
  }
})

testthat::test_that("16 numeric compositions restore saved outer states and reproduce on two workers", {
  for (pkg in c("workflows", "rsample", "tailor", "probably", "recipes", "tune", "yardstick")) testthat::skip_if_not_installed(pkg)
  testthat::skip_if_not(worker_replay_available())
  d <- regression_composition_data(); base <- regression_composition_base(d)
  new <- d[1:7, c("x", "z")]
  snapshot <- function(fit, new, schema) lapply(c(0L, 1L, 7L), function(n) {
    batch <- new[seq_len(n), , drop = FALSE]
    list(prediction = predict(fit, batch), details = vennabers::venn_abers_details(fit, batch))
  })
  environment(snapshot) <- globalenv()
  metrics <- "yardstick::metric_set(yardstick::rmse, yardstick::mae)"
  packages <- c("vennabers", "stats")
  outer <- withr::with_seed(31, rsample::vfold_cv(d, v = 2))
  control <- tune::control_resamples(save_pred = TRUE, extract = identity, parallel_over = "resamples", pkgs = packages)
  records <- lapply(seq_len(nrow(regression_composition_cells)), function(i) {
    model <- regression_composition_spec(base, regression_composition_cells[i, ])
    evaluation <- withr::with_seed(43, tune::fit_resamples(model, outer,
      metrics = eval(parse(text = metrics)), control = control))
    testthat::expect_identical(nrow(tune::collect_notes(evaluation)), 0L)
    list(snapshot = snapshot, replay = list(model = model, outer = outer, evaluation = evaluation, new = new),
      expected = lapply(tune::collect_extracts(evaluation)$.extracts, snapshot, new = new, schema = NULL))
  })
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  saveRDS(records, saved)
  guarded <- c(workflows = "fit.workflow", recipes = "prep.recipe", tailor = "fit.tailor", stats = "isoreg")
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "library(vennabers)",
    paste0("records <- readRDS(", deparse(saved), ")"),
    sprintf("trace('%s', where = asNamespace('%s'), tracer = quote(stop('Unexpected learning')), print = FALSE)", guarded, names(guarded)),
    "for (s in records) stopifnot(identical(s$expected, lapply(tune::collect_extracts(s$replay$evaluation)$.extracts, s$snapshot, new = s$replay$new, schema = NULL)))",
    sprintf("untrace('%s', where = asNamespace('%s'))", guarded, names(guarded)),
    "for (s in records) {",
    worker_replay(metrics, packages, character()),
    "}",
    "cat('replayed', length(records), 'numeric compositions\\n')"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  output <- readLines(log)
  testthat::expect_equal(status, 0L, info = paste(output, collapse = "\n"))
  testthat::expect_true("replayed 16 numeric compositions" %in% output)
})

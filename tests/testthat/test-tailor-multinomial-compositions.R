# Shared Tailor multinomial compositions from three classes: the 16 fixed cells
# (OvO/OvR x holdout/ensemble x Inductive/Cross x matched/full) layered against
# independently fitted public Tailor workflows (smooth = FALSE) and the
# unchanged manual VA core, and the fixed probably/nnet iteration limit.
multinomial_composition_data <- function() withr::with_seed(123, {
  truth <- factor(rep(c("z", "a", "b"), length.out = 360), levels = c("z", "a", "b"))
  data.frame(x = rnorm(360, sd = 2) + c(-1, 1, -1)[truth], z = rnorm(360, sd = 2) + c(-1, -1, 1)[truth],
    w = rnorm(360, sd = 2) + c(1, -1, -1)[truth], truth = truth)
})
multinomial_composition_base <- function(d) workflows::workflow() |>
  workflows::add_recipe(recipes::recipe(truth ~ ., d) |> recipes::step_normalize(recipes::all_numeric_predictors())) |>
  workflows::add_model(parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE, maxit = 1000))
multinomial_canonical_tailor <- function() tailor::tailor() |>
  tailor::adjust_probability_calibration(method = "multinomial", smooth = FALSE)
multinomial_composition_cells <- expand.grid(mode = c("matched", "full"), method = c("inductive", "cross"),
  strategy = c("holdout", "ensemble"), classification = c("ovo", "ovr"), stringsAsFactors = FALSE)
# Automatic VA folds (Cross v = 2) unless explicit VA splits are supplied.
multinomial_composition_spec <- function(base, cell, splits = NULL) {
  post <- list(method = "tailor_multinomial", strategy = cell$strategy)
  if (cell$strategy == "holdout") post$prop <- .3 else post$v <- 3L
  args <- list(base, "truth", method = cell$method, classification = cell$classification,
    prediction_model = cell$mode, seed = 24, postprocessing = post)
  if (!is.null(splits)) args$splits <- splits else if (cell$method == "cross") args$v <- 2L
  do.call(calibrated_workflow, args)
}

testthat::test_that("16 shared multinomial compositions equal public Tailor workflows, VA cores and cost formulas", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes", "nnet")) testthat::skip_if_not_installed(pkg)
  d <- multinomial_composition_data(); base <- multinomial_composition_base(d); classes <- levels(d$truth)
  predictors <- c("x", "z", "w")
  canonical <- multinomial_canonical_tailor()
  public <- function(p) generics::fit(workflows::add_tailor(base, canonical), d[p$B, ], data_calibration = d[p$P, ])
  probability <- function(model, data) {
    p <- as.matrix(predict(model, data, type = "prob")[paste0(".pred_", classes)])
    colnames(p) <- classes
    p
  }
  audit <- new.env(); audit$fits <- audit$post <- audit$predictions <- audit$tailor <- 0L
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
  multinomial <- probably::cal_estimate_multinomial
  testthat::local_mocked_bindings(cal_estimate_multinomial = function(...) {
    audit$post <- audit$post + 1L
    multinomial(...)
  }, .package = "probably")
  for (i in seq_len(nrow(multinomial_composition_cells))) {
    cell <- multinomial_composition_cells[i, ]; info <- paste(cell, collapse = " ")
    full <- cell$mode == "full"; cross <- cell$method == "cross"
    # K VA source folds, J inner B/P pairs per composition.
    K <- if (cross) 2L else 1L; J <- if (cell$strategy == "ensemble") 3L else 1L
    audit$fits <- audit$post <- 0L
    fit <- generics::fit(multinomial_composition_spec(base, cell), d)
    state <- workflows::extract_fit_engine(fit)
    # Engine and Tailor fits cost J, 2J, KJ or (K+1)J and retain J, J, KJ or J pairs.
    cost <- (K + full) * J; retained <- if (full) J else K * J
    testthat::expect_identical(audit$fits, cost, info = info)
    testthat::expect_identical(audit$post, cost, info = info)
    testthat::expect_output(print(state), paste("Base fits:", cost))
    testthat::expect_output(print(state), paste("Models per nonempty prediction:", retained))
    testthat::expect_identical(state$class_levels, classes)
    plans <- function(id) if (J > 1L) state$postprocessing[[id]]$inner_folds else list(state$postprocessing[[id]])
    leaves <- function(id) if (J > 1L) state$models[[id]]$models else list(state$models[[id]])
    sources <- names(state$calibration)
    testthat::expect_named(state$models, if (full) "Full1" else sources)
    # Every retained leaf equals its own public B/P Tailor workflow; the inner
    # mean weights complete class-aligned vectors equally.
    composed <- list()
    for (id in names(state$models)) {
      ps <- plans(id); models <- leaves(id)
      testthat::expect_length(models, J)
      for (j in seq_len(J)) {
        p <- ps[[j]]
        testthat::expect_identical(p$method, "tailor_multinomial")
        testthat::expect_identical(p$rank, 3L)
        testthat::expect_true(p$converged)
        testthat::expect_lte(max(abs(probability(public(p), d) - probability(models[[j]], d))), 1e-10)
      }
      composed[[id]] <- probability(state$models[[id]], d)
      testthat::expect_lte(max(abs(composed[[id]] - Reduce(`+`, lapply(models, probability, data = d)) / J)), 1e-10)
    }
    # Held C values come from the source compositions only, also in full mode.
    for (id in sources) {
      ps <- plans(id); C <- ps[[1L]]$C
      expected <- Reduce(`+`, lapply(ps, function(p) probability(public(p), d[C, ]))) / J
      testthat::expect_lte(max(abs(expected - state$calibration[[id]]$probability)), 1e-10)
    }
    # Unchanged manual VA core on the stored held C values and composed inputs.
    truth <- do.call(c, lapply(state$calibration, `[[`, "truth"))
    held <- do.call(rbind, lapply(state$calibration, `[[`, "probability"))
    colnames(held) <- classes
    ids <- rep(sources, vapply(state$calibration, function(record) length(record$row), 1L))
    core <- if (cross) cross_venn_abers(truth, held, ids, method = cell$classification) else
      venn_abers(truth, held, method = cell$classification)
    input <- if (!cross) composed[[1L]] else if (full) stats::setNames(rep(composed["Full1"], K), sources) else composed[sources]
    for (loss in c("log", "brier")) {
      actual <- venn_abers_details(fit, d[predictors], loss = loss)
      testthat::expect_lte(max(abs(as.matrix(actual$prediction) - predict(core, input, loss = loss))), 1e-10)
      testthat::expect_identical(predict(fit, d, type = "class", loss = loss)$.pred_class,
        predict(core, input, type = "class", loss = loss))
      testthat::expect_true(all(actual$details$postprocessing_method == "tailor_multinomial"))
    }
    # Prediction learns nothing and calls each retained workflow and Tailor once per nonempty batch.
    audit$fits <- audit$post <- 0L
    for (n in c(0L, 1L, 7L)) {
      audit$predictions <- audit$tailor <- 0L
      small <- venn_abers_details(fit, d[seq_len(n), predictors])
      testthat::expect_identical(audit$predictions, if (n) retained else 0L, info = info)
      testthat::expect_identical(audit$tailor, if (n) retained else 0L, info = info)
      testthat::expect_named(small$prediction, paste0(".pred_", classes))
      testthat::expect_identical(nrow(small$prediction), n)
      testthat::expect_true(all(is.finite(as.matrix(small$prediction))))
    }
    testthat::expect_identical(c(audit$fits, audit$post), c(0L, 0L))
  }
})

testthat::test_that("full multinomial compositions learn F from C but keep that source's Q fit", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes", "nnet")) testthat::skip_if_not_installed(pkg)
  d <- multinomial_composition_data(); base <- multinomial_composition_base(d)
  # The first source holds C = 271:360 in both VA methods; under Cross these rows
  # belong to the Q of the other two sources.
  held <- list(271:360, 1:135, 136:270)
  va <- list(inductive = function(data) rsample::make_splits(list(analysis = 1:270, assessment = held[[1L]]), data),
    cross = function(data) rsample::manual_rset(lapply(held, function(rows)
      rsample::make_splits(list(analysis = setdiff(1:360, rows), assessment = rows), data)), c("one", "two", "three")))
  learned <- function(state, id = "Full1") {
    p <- state$postprocessing[[id]]
    lapply(if (is.null(p$inner_folds)) list(p) else p$inner_folds, `[[`, "coefficients")
  }
  C <- held[[1L]]
  for (method in names(va)) for (classification in c("ovo", "ovr")) for (strategy in c("holdout", "ensemble")) {
    cell <- list(mode = "full", method = method, strategy = strategy, classification = classification)
    spec <- multinomial_composition_spec(base, cell, splits = va[[method]])
    fit <- function(data) workflows::extract_fit_engine(generics::fit(spec, data))
    original <- fit(d); source <- names(original$calibration)[1L]
    for (column in c("x", "truth")) {
      changed <- d; info <- paste(method, classification, strategy, column)
      if (column == "x") changed$x[C] <- changed$x[C] * 2 + .3 else
        changed$truth[C] <- levels(d$truth)[as.integer(d$truth[C]) %% 3L + 1L]
      after <- fit(changed)
      testthat::expect_identical(after$postprocessing[[source]], original$postprocessing[[source]], info = info)
      testthat::expect_false(identical(learned(after), learned(original)), info = info)
      if (method == "cross") for (other in names(original$calibration)[-1L]) {
        testthat::expect_false(identical(learned(after, other), learned(original, other)), info = info)
      }
      same_held <- identical(after$calibration[[source]]$probability, original$calibration[[source]]$probability)
      testthat::expect_identical(same_held, column == "truth", info = info)
      if (column == "truth") testthat::expect_false(identical(after$calibrator, original$calibrator), info = info)
    }
  }
})

testthat::test_that("16 shared multinomial compositions restore from RDS without learning and replay on outer workers", {
  for (pkg in c("workflows", "rsample", "tailor", "probably", "recipes", "nnet")) testthat::skip_if_not_installed(pkg)
  testthat::skip_if_not(worker_replay_available())
  d <- multinomial_composition_data(); base <- multinomial_composition_base(d)
  new <- d[1:7, c("x", "z", "w")]
  snapshot <- function(fit, new, schema) lapply(c("log", "brier"), function(loss) lapply(c(0L, 1L, 7L), function(n) {
    x <- new[seq_len(n), , drop = FALSE]
    list(prob = predict(fit, x, type = "prob", loss = loss), class = predict(fit, x, type = "class", loss = loss),
      details = vennabers::venn_abers_details(fit, x, loss = loss))
  }))
  environment(snapshot) <- globalenv()
  metrics <- "yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)"
  packages <- c("vennabers", "stats", "nnet")
  outer <- withr::with_seed(31, rsample::vfold_cv(d, v = 2))
  control <- tune::control_resamples(save_pred = TRUE, extract = identity, parallel_over = "resamples", pkgs = packages)
  records <- lapply(seq_len(nrow(multinomial_composition_cells)), function(i) {
    model <- multinomial_composition_spec(base, multinomial_composition_cells[i, ])
    evaluation <- withr::with_seed(43, tune::fit_resamples(model, outer, metrics = eval(parse(text = metrics)), control = control))
    testthat::expect_identical(nrow(tune::collect_notes(evaluation)), 0L)
    list(snapshot = snapshot, replay = list(model = model, outer = outer, evaluation = evaluation, new = new),
      expected = lapply(tune::collect_extracts(evaluation)$.extracts, snapshot, new = new, schema = NULL))
  })
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  saveRDS(records, saved)
  # Guard every restored extract once, then untrace once: the worker replay
  # inside the loop fits by design, so its helper gets no untrace targets.
  guarded <- c(workflows = "fit.workflow", recipes = "prep.recipe", tailor = "fit.tailor")
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
    "cat('replayed', length(records), 'compositions\\n')"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  output <- readLines(log)
  testthat::expect_equal(status, 0L, info = paste(output, collapse = "\n"))
  testthat::expect_true("replayed 16 compositions" %in% output)
})

testthat::test_that("the fixed Tailor multinomial iteration limit fails weakly informative P without fallback", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "nnet")) testthat::skip_if_not_installed(pkg)
  # Weak signal: B-trained probabilities on P stay within about 0.26-0.42.
  weak <- function(seed, n) withr::with_seed(seed, {
    truth <- factor(rep(c("a", "b", "c"), length.out = n))
    data.frame(x = rnorm(n) + .02 * c(-1, 0, 1)[truth], z = rnorm(n) + .02 * c(1, -1, 0)[truth], truth = truth)
  })
  d <- rbind(weak(8, 600), weak(9, 120)); B <- 1:420; P <- 421:600
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x + z) |>
    workflows::add_model(parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
  # The public Tailor/probably path itself calls nnet::multinom without an
  # iteration setting and stops at its default of 100 iterations; an
  # unrestricted fit of the same input converges to clearly different values.
  reference <- generics::fit(workflows::add_tailor(base, multinomial_canonical_tailor()), d[B, ], data_calibration = d[P, ])
  stopped <- workflows::extract_tailor(reference)$adjustments[[1L]]$results$fit$estimates[[1L]]$estimate
  testthat::expect_identical(stopped$convergence, 1L)
  testthat::expect_null(stopped$call$maxit)
  engine <- workflows::extract_fit_parsnip(reference)
  input <- predict(engine, d[P, ], type = "prob")
  testthat::expect_true(all(as.matrix(input) > .25 & as.matrix(input) < .43))
  input$truth <- d$truth[P]
  limit <- nnet::multinom(truth ~ .pred_a + .pred_b, input, trace = FALSE, maxit = 100L)
  testthat::expect_identical(limit$convergence, 1L)
  testthat::expect_identical(coef(limit), coef(stopped))
  optimum <- nnet::multinom(truth ~ .pred_a + .pred_b, input, trace = FALSE, maxit = 10000L)
  testthat::expect_identical(optimum$convergence, 0L)
  testthat::expect_gt(max(abs(predict(optimum, input, type = "probs") - predict(stopped, input, type = "probs"))), .01)
  calls <- 0L; multinomial <- probably::cal_estimate_multinomial
  testthat::local_mocked_bindings(cal_estimate_multinomial = function(...) {
    calls <<- calls + 1L
    multinomial(...)
  }, .package = "probably")
  va <- function(data) rsample::make_splits(list(analysis = 1:600, assessment = 601:720), data)
  pairs <- list(holdout = function(data) rsample::make_splits(list(analysis = B, assessment = P), data),
    ensemble = function(data) rsample::manual_rset(lapply(list(1:210, 211:420, P), function(p)
      rsample::make_splits(list(analysis = setdiff(1:600, p), assessment = p), data)), c("first", "second", "weak")))
  for (classification in c("ovo", "ovr")) for (strategy in names(pairs)) {
    calls <- 0L
    spec <- calibrated_workflow(base, "truth", classification = classification, splits = va,
      postprocessing = list(method = "tailor_multinomial", strategy = strategy, splits = pairs[[strategy]]))
    label <- if (strategy == "holdout") "Holdout1" else "weak"
    # No retry, no switch to none and no partial ensemble: converged pairs are discarded.
    testthat::expect_error(generics::fit(spec, d), paste0("Inductive1 / ", label,
      " tailor_multinomial \\(B=420, P=180\\).*P multinomial fit did not converge"))
    testthat::expect_identical(calls, if (strategy == "holdout") 1L else 3L)
  }
})

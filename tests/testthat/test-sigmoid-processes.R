testthat::test_that("Sigmoid holdouts produce complete native outer results", {
  for (pkg in c("tune", "yardstick", "recipes")) testthat::skip_if_not_installed(pkg)
  source(system.file("examples", "sigmoid-holdout.R", package = "vennabers"), local = new.env())
  testthat::expect_true(TRUE)
})

testthat::test_that("Sigmoid schemas and full predictive state survive fresh R without learning", {
  for (pkg in c("workflows", "recipes")) testthat::skip_if_not_installed(pkg)
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
  for (input in c("probability", "glm_link")) for (method in c("inductive", "cross")) {
    fit <- generics::fit(calibrated_workflow(base, outcome, schema = schema, method = method, seed = 24,
      postprocessing = list(method = "sigmoid", sigmoid_input = input, strategy = "holdout", prop = .3)), d)
    new <- d[setdiff(names(d), outcome)]
    details <- venn_abers_details(fit, new)
    decoded <- venn_abers_decode(details, schema)
    testthat::expect_identical(colnames(decoded$probability), levels(original[[schema$outcome]]))
    testthat::expect_identical(decoded$postprocessing[[1]]$classes, schema$class_levels)
    testthat::expect_identical(decoded$postprocessing[[1]]$columns, names(original))
    testthat::expect_identical(details$postprocessing[[1]]$schema, schema)
    testthat::expect_identical(decoded$postprocessing[[1]]$sigmoid_class, schema$class_levels[2])
    records[[paste(input, method)]] <- list(fit = fit, new = new, expected = outputs(fit, new))
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
    "trace('optim', where=asNamespace('stats'), tracer=quote(stop('Unexpected learning')), print=FALSE)",
    "trace('glm.fit', where=asNamespace('stats'), tracer=quote(stop('Unexpected learning')), print=FALSE)",
    paste0("outputs <- ", paste(deparse(outputs), collapse = "\n")),
    paste0("records <- readRDS(", deparse(rds), ")"),
    "for (r in records) stopifnot(identical(r$expected, outputs(r$fit, r$new)))"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

testthat::test_that("GLM-link Sigmoid completes all 16 fixed compositions natively, layer by layer and from RDS", {
  for (pkg in c("workflows", "recipes", "rsample", "tune", "yardstick")) testthat::skip_if_not_installed(pkg)
  checker <- new.env()
  sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), checker)
  cl <- factor(rep(c("a", "b", "c"), length.out = 480))
  raw <- withr::with_seed(123, data.frame(x = rnorm(480, sd = 2) + c(-1, 1, -1)[cl],
    z = rnorm(480, sd = 2) + c(-1, -1, 1)[cl], truth = cl))
  metrics <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
  # Independent odds coupling of each class pair's binary VA result.
  couple <- function(q, classes) {
    pairs <- utils::combn(classes, 2, simplify = FALSE)
    p <- vapply(classes, function(cl) {
      relevant <- which(vapply(pairs, function(pair) cl %in% pair, logical(1)))
      odds <- lapply(relevant, function(j) {
        prob <- if (cl == pairs[[j]][2L]) q[, j] else 1 - q[, j]
        (1 - prob) / prob
      })
      1 / (1 + Reduce(`+`, odds))
    }, numeric(nrow(q)))
    p / rowSums(p)
  }
  fit_method <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  predict_method <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  audit <- new.env(); audit$fits <- audit$predictions <- 0L
  base_call <- function(object) !inherits(workflows::extract_spec_parsnip(object), "vennabers")
  registerS3method("fit", "workflow", function(object, ...) {
    if (base_call(object)) audit$fits <- audit$fits + 1L
    fit_method(object, ...)
  }, envir = asNamespace("generics"))
  registerS3method("predict", "workflow", function(object, ...) {
    if (base_call(object)) audit$predictions <- audit$predictions + 1L
    predict_method(object, ...)
  }, envir = asNamespace("stats"))
  on.exit({
    registerS3method("fit", "workflow", fit_method, envir = asNamespace("generics"))
    registerS3method("predict", "workflow", predict_method, envir = asNamespace("stats"))
  }, add = TRUE)
  optimizer <- stats::optim; audit$optim <- 0L
  testthat::local_mocked_bindings(optim = function(...) { audit$optim <- audit$optim + 1L; optimizer(...) }, .package = "stats")
  records <- list()
  for (path in c("binary", "pairwise")) {
    d <- if (path == "binary") droplevels(raw[raw$truth != "c", ]) else raw
    rownames(d) <- NULL; n <- nrow(d); train <- seq_len(n * 3 / 4); test <- setdiff(seq_len(n), train)
    base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x + z, d) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(parsnip::logistic_reg(engine = "glm"))
    outer <- withr::with_seed(42, rsample::vfold_cv(d, v = 2, repeats = 2))
    split <- rsample::make_splits(list(analysis = train, assessment = test), d)
    for (strategy in c("holdout", "ensemble")) for (method in c("inductive", "cross")) for (mode in c("matched", "full")) {
      key <- paste(path, strategy, method, mode)
      post <- list(method = "sigmoid", sigmoid_input = "glm_link", strategy = strategy)
      if (strategy == "holdout") post$prop <- .3 else post$v <- 3L
      args <- list(base, "truth", method = method, prediction_model = mode, seed = 24, postprocessing = post)
      if (method == "cross") args$v <- 5L
      if (path == "pairwise") { args$classification <- "ovo"; args$model_training <- "pairwise" }
      model <- do.call(calibrated_workflow, args)
      result <- tune::fit_resamples(model, outer, metrics = metrics,
        control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, allow_par = FALSE))
      states <- checker$check_outer_identities(result, outer, metrics = 2L)$extracts
      for (h in seq_along(states)) {
        F_rows <- as.integer(outer$splits[[h]], "analysis"); A_rows <- as.integer(outer$splits[[h]], "assessment")
        for (part in if (path == "pairwise") states[[h]]$pair_models else list(states[[h]])) {
          rows <- if (path == "pairwise") part$rows else seq_along(F_rows)
          for (record in part$postprocessing) for (p in if (strategy == "ensemble") record$inner_folds else list(record)) {
            testthat::expect_false(any(F_rows[rows[c(p$B, p$P, p$C)]] %in% A_rows), info = key)
            testthat::expect_length(intersect(p$B, p$P), 0L)
            testthat::expect_length(intersect(p$Q, p$C), 0L)
            testthat::expect_identical(p$input, "glm_link")
          }
        }
      }
      last <- tune::last_fit(model, split, metrics = metrics)
      native <- tune::collect_predictions(last)
      testthat::expect_equal(nrow(tune::collect_notes(last)), 0L)
      testthat::expect_identical(native$.row, test)
      testthat::expect_true(all(is.finite(tune::collect_metrics(last)$.estimate)))
      # Costs per domain with J inner pairs and K=5 VA folds: J, 2J, KJ or (K+1)J.
      audit$fits <- audit$optim <- 0L
      fitted <- generics::fit(model, d[train, ])
      columns <- paste0(".pred_", levels(d$truth))
      testthat::expect_identical(as.matrix(native[columns]), as.matrix(predict(fitted, d[test, ], type = "prob")[columns]))
      J <- if (strategy == "ensemble") 3L else 1L
      K <- if (method == "cross") 5L else 1L
      domains <- if (path == "pairwise") 3L else 1L
      cost <- domains * J * (K + (mode == "full"))
      retained <- domains * J * (if (mode == "full") 1L else K)
      if (key == "pairwise ensemble cross full") testthat::expect_identical(c(cost, retained), c(54L, 9L))
      testthat::expect_equal(c(audit$fits, audit$optim), c(cost, cost), info = key)
      state <- workflows::extract_fit_engine(fitted)
      testthat::expect_output(print(state), paste0("Base fits: ", cost, "\\b"))
      testthat::expect_output(print(state), paste0("Models per nonempty prediction: ", retained, "\\b"))
      new <- d[test[1:7], c("x", "z")]
      for (size in c(0L, 1L, 7L)) {
        audit$predictions <- 0L
        predict(fitted, new[seq_len(size), ], type = "prob")
        testthat::expect_equal(audit$predictions, if (size) retained else 0L, info = key)
      }
      testthat::expect_equal(audit$fits, cost)
      # Leaf from its own B engine link, equal inner mean, then the unchanged VA core.
      q <- list(log = list(), brier = list())
      for (part in if (path == "pairwise") state$pair_models else list(state)) {
        event <- paste0(".pred_", part$event_class)
        composed <- lapply(part$models, function(m) {
          leaves <- if (strategy == "ensemble") m$models else list(m)
          values <- vapply(leaves, function(leaf) {
            link <- stats::predict(workflows::extract_fit_engine(leaf$workflow),
              recipes::bake(workflows::extract_recipe(leaf$workflow), new), type = "link")
            rebuilt <- stats::plogis(-(leaf$coefficients[[1L]] * link / leaf$scale + leaf$coefficients[[2L]]))
            actual <- predict(leaf, new)
            testthat::expect_lte(max(abs(actual[[paste0(".pred_", leaf$classes[2L])]] - rebuilt)), 1e-7)
            actual[[event]]
          }, numeric(nrow(new)))
          mean_value <- rowMeans(matrix(values, nrow(new)))
          actual <- predict(m, new, type = "prob")[[event]]
          testthat::expect_lte(max(abs(actual - mean_value)), 1e-10)
          actual
        })
        calibration <- part$calibration
        ids <- names(calibration)
        # Matched sources still exist: their predictions on C are the held values.
        if (mode == "matched") for (id in ids) {
          domain <- d[train, ][if (path == "pairwise") part$rows else seq_along(train), ]
          C <- part$postprocessing[[id]]$C
          testthat::expect_identical(calibration[[id]]$row, C)
          testthat::expect_lte(max(abs(predict(part$models[[id]], domain[C, ], type = "prob")[[event]] -
            calibration[[id]]$probability)), 1e-10)
        }
        input <- do.call(cbind, lapply(ids, function(id) composed[[if (mode == "full") "Full1" else id]]))
        truth <- unlist(lapply(calibration, `[[`, "truth"), use.names = FALSE)
        held <- unlist(lapply(calibration, `[[`, "probability"), use.names = FALSE)
        core <- if (method == "cross") cross_venn_abers(truth, held, rep(ids, lengths(lapply(calibration, `[[`, "row")))) else
          venn_abers(truth, held)
        for (loss in names(q)) {
          expected <- predict(core, if (method == "cross") `colnames<-`(input, ids) else input[, 1L], loss = loss)$probability
          if (path == "pairwise") q[[loss]][[part$pair_id]] <- expected else
            testthat::expect_lte(max(abs(predict(fitted, new, type = "prob", loss = loss)[[event]] - expected)), 1e-10)
        }
      }
      if (path == "pairwise") for (loss in names(q)) testthat::expect_lte(max(abs(unname(as.matrix(
        predict(fitted, new, type = "prob", loss = loss))) - couple(do.call(cbind, q[[loss]]), levels(d$truth)))), 1e-10)
      if (mode == "full") {
        matched <- workflows::extract_fit_engine(records[[sub("full$", "matched", key)]]$fit)
        left <- if (path == "pairwise") matched$pair_models else list(matched)
        right <- if (path == "pairwise") state$pair_models else list(state)
        for (i in seq_along(right)) {
          testthat::expect_identical(right[[i]]$calibrator, left[[i]]$calibrator)
          testthat::expect_identical(right[[i]]$calibration, left[[i]]$calibration)
          testthat::expect_identical(right[[i]]$folds, left[[i]]$folds)
          testthat::expect_identical(right[[i]]$postprocessing[names(left[[i]]$postprocessing)], left[[i]]$postprocessing)
          testthat::expect_named(right[[i]]$models, "Full1")
        }
      }
      records[[key]] <- list(fit = fitted, new = new)
      if (key == "pairwise ensemble cross full") replay <- list(model = model, outer = outer, new = new)
    }
  }
  snapshot <- function(fit, new, schema = NULL) lapply(c(0L, 1L, 7L), function(n) {
    batch <- new[seq_len(n), , drop = FALSE]
    list(prob = predict(fit, batch, type = "prob"), class = predict(fit, batch, type = "class"),
      details = venn_abers_details(fit, batch))
  })
  for (key in names(records)) records[[key]]$expected <- snapshot(records[[key]]$fit, records[[key]]$new)
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  workers <- worker_replay_available()
  saveRDS(list(records = records, snapshot = snapshot, replay = if (workers) replay), saved)
  guarded <- c(workflows = "fit.workflow", stats = "optim", stats = "glm.fit", recipes = "prep.recipe",
    rsample = "initial_split", rsample = "vfold_cv")
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "library(vennabers); loadNamespace('workflows')",
    sprintf("trace('%s', where = asNamespace('%s'), tracer = quote(stop('Unexpected learning or planning')), print = FALSE)",
      guarded, names(guarded)),
    paste0("s <- readRDS(", deparse(saved), ")"),
    "for (r in s$records) stopifnot(identical(r$expected, s$snapshot(r$fit, r$new)))",
    # The most composed GLM-link cell must also rebuild identically on two outer workers.
    if (workers) worker_replay("yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)",
      c("vennabers", "stats"), guarded)
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

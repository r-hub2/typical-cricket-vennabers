testthat::test_that("all 80 fixed compositions preserve calibration and survive native evaluation and RDS", {
  for (pkg in c("tailor", "probably", "recipes", "tune", "yardstick", "nnet")) testthat::skip_if_not_installed(pkg)
  example <- new.env()
  sys.source(system.file("examples", "full-postprocessing.R", package = "vennabers"), example)
  testthat::expect_length(example$cases, 80L)
  # Start with multinomial Tailor in the fresh process, before another method
  # could load its optional S3 registrations incidentally.
  cases <- example$cases
  cases <- cases[order(!grepl("tailor_multinomial.*full$", names(cases)))]
  components <- function(state) if (identical(state$parameters$requested$model_training, "pairwise")) state$pair_models else list(state)
  # Classification: 1e-10 absolute. Regression: 1e-10*s + 1e-10*|expected| with
  # a positive scale s fixed from held truth before any comparison.
  close_within <- function(actual, expected, s = NULL, absolute = 1e-10) {
    testthat::expect_identical(dim(actual), dim(expected))
    testthat::expect_true(all(is.finite(actual)))
    tolerance <- if (is.null(s)) absolute else 1e-10 * s + 1e-10 * abs(expected)
    testthat::expect_true(all(abs(actual - expected) <= tolerance))
  }
  # Postprocessor layer from each retained leaf's own B engine output.
  rebuild_leaf <- function(m, new, columns, regression) {
    if (inherits(m, "vennabers_sigmoid")) {
      x <- predict(m$workflow, new, type = "prob")[[paste0(".pred_", m$classes[2L])]]
      p <- 1 / (1 + exp(m$coefficients[[1L]] * x / m$scale + m$coefficients[[2L]]))
      return(matrix(c(1 - p, p), ncol = 2L, dimnames = list(NULL, paste0(".pred_", m$classes)))[, columns, drop = FALSE])
    }
    if (inherits(m, "vennabers_isotonic_workflow")) {
      x <- predict(m$workflow, new)$.pred; knots <- m$isotonic
      y <- if (length(knots$x) == 1L) rep(knots$y, length(x)) else stats::approx(knots$x, knots$y, x, rule = 2)$y
      return(matrix(y, dimnames = list(NULL, columns)))
    }
    engine <- predict(workflows::extract_fit_parsnip(m), recipes::bake(workflows::extract_recipe(m), new),
      type = if (regression) "numeric" else "prob")
    as.matrix(predict(workflows::extract_tailor(m), engine)[columns])
  }
  # Independent odds coupling, after each class pair's binary VA result.
  couple <- function(q, classes) {
    pairs <- utils::combn(classes, 2, simplify = FALSE)
    p <- vapply(classes, function(cl) {
      js <- which(vapply(pairs, function(pair) cl %in% pair, logical(1)))
      odds <- lapply(js, function(j) {
        prob <- if (cl == pairs[[j]][2]) q[, j] else 1 - q[, j]
        (1 - prob) / prob
      })
      1 / (1 + Reduce(`+`, odds))
    }, numeric(nrow(q)))
    p / rowSums(p)
  }
  for (name in names(cases)[grepl("full$", names(cases))]) {
    r <- cases[[name]]; state <- workflows::extract_fit_engine(r$fit)
    matched <- workflows::extract_fit_engine(cases[[sub("full$", "matched", name)]]$fit)
    settings <- state$parameters$requested
    regression <- settings$mode == "regression"; cross <- settings$method == "cross"
    pairwise <- identical(settings$model_training, "pairwise")
    new <- r$new[1:7, , drop = FALSE]
    actual <- venn_abers_details(r$fit, new)
    left <- components(matched); right <- components(state)
    probabilities <- list(log = list(), brier = list())
    for (i in seq_along(right)) {
      part <- right[[i]]; old <- left[[i]]
      testthat::expect_identical(part$calibrator, old$calibrator)
      testthat::expect_identical(part$calibration, old$calibration)
      testthat::expect_identical(part$folds, old$folds)
      testthat::expect_identical(part$postprocessing[names(old$postprocessing)], old$postprocessing)
      testthat::expect_named(part$models, "Full1")
      diagnostic <- part$postprocessing$Full1; model <- part$models$Full1
      models <- if (settings$postprocessing$strategy == "ensemble") model$models else list(model)
      plans <- if (settings$postprocessing$strategy == "ensemble") diagnostic$inner_folds else list(diagnostic)
      testthat::expect_equal(length(models), diagnostic$observed_base_fits)
      testthat::expect_equal(length(models), diagnostic$observed_postprocessor_fits)
      for (j in seq_along(plans)) {
        p <- plans[[j]]
        testthat::expect_setequal(c(p$B, p$P), p$Q)
        testthat::expect_length(intersect(p$B, p$P), 0L)
        testthat::expect_length(p$C, 0L)
        testthat::expect_identical(p$schema, r$schema)
        testthat::expect_identical(p$VA_fold, "Full1")
        testthat::expect_identical(p$method, settings$postprocessing$method)
      }
      # Leaf prediction, equal-weight composition, then independent public VA
      # construction from the original held C values (never F-on-C values).
      calibration <- part$calibration
      truth <- do.call(c, lapply(calibration, `[[`, "truth"))
      s <- if (regression) stats::sd(truth)
      columns <- if (regression) ".pred" else paste0(".pred_", part$class_levels)
      leaf <- lapply(models, function(m) as.matrix(predict(m, new, type = if (regression) "numeric" else "prob")[columns]))
      for (j in seq_along(models)) {
        rebuilt <- rebuild_leaf(models[[j]], new, columns, regression)
        if (settings$postprocessing$method == "sigmoid") close_within(leaf[[j]], rebuilt, absolute = 1e-7) else
          close_within(leaf[[j]], rebuilt, if (settings$postprocessing$method == "isotonic_numeric") 1 else s)
      }
      average <- Reduce(`+`, leaf) / length(leaf)
      composed <- as.matrix(predict(model, new, type = if (regression) "numeric" else "prob")[columns])
      close_within(composed, average, s)
      # Use the exact composed values for VA: rounding of inner means can move
      # tied scores in a step function, independently of mean-layer tolerance.
      ids <- rep(names(calibration), vapply(calibration, function(c) length(c$row), 1L))
      if (regression) {
        held <- unlist(lapply(calibration, `[[`, "prediction"), use.names = FALSE)
        core <- if (cross) cross_venn_abers_regression(truth, held, ids) else venn_abers_regression(truth, held)
        expected <- predict(core, if (cross) matrix(composed[, 1], nrow(new), length(calibration),
          dimnames = list(NULL, names(calibration))) else composed[, 1])
        close_within(actual$prediction$.pred, expected$prediction, s)
      } else {
        binary <- pairwise || settings$classification == "binary"
        held <- if (binary) unlist(lapply(calibration, `[[`, "probability"), use.names = FALSE) else
          do.call(rbind, lapply(calibration, `[[`, "probability"))
        args <- list(truth = truth, probability = held)
        if (!binary) args$method <- settings$classification
        if (cross) args$fold <- ids
        core <- do.call(if (cross) cross_venn_abers else venn_abers, args)
        input <- if (binary) composed[, paste0(".pred_", part$event_class)] else {
          colnames(composed) <- part$class_levels
          composed
        }
        if (cross) input <- if (binary) matrix(input, nrow(new), length(calibration), dimnames = list(NULL, names(calibration))) else
          stats::setNames(rep(list(input), length(calibration)), names(calibration))
        for (loss in names(probabilities)) {
          expected <- predict(core, input, loss = loss)
          if (binary) expected <- expected$probability
          if (pairwise) probabilities[[loss]][[i]] <- expected else {
            prediction <- as.matrix(predict(r$fit, new, type = "prob", loss = loss))
            if (binary) prediction <- prediction[, paste0(".pred_", part$event_class)]
            close_within(unname(prediction), unname(expected))
          }
        }
      }
    }
    if (pairwise) for (loss in names(probabilities)) close_within(unname(as.matrix(predict(r$fit, new, type = "prob", loss = loss))),
      unname(couple(do.call(cbind, probabilities[[loss]]), state$class_levels)))
  }
  snapshot <- function(fit, new, schema) lapply(c(0L, 1L, 7L), function(n) {
    batch <- new[seq_len(n), , drop = FALSE]
    details <- venn_abers_details(fit, batch)
    list(prediction = predict(fit, batch, type = if (schema$mode == "regression") "numeric" else "prob"),
      class = if (schema$mode != "regression") predict(fit, batch, type = "class"),
      details = details, decoded = venn_abers_decode(details, schema))
  })
  records <- lapply(cases, function(r) {
    r$expected <- snapshot(r$fit, r$new, r$schema)
    r
  })
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  saveRDS(records, saved)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "library(vennabers); loadNamespace('workflows')",
    "guard <- function(...) stop('Unexpected learning or planning')",
    "registerS3method('fit', 'workflow', guard, envir=asNamespace('generics'))",
    "registerS3method('fit', 'tailor', guard, envir=asNamespace('generics'))",
    "for (f in c('optim', 'glm.fit', 'isoreg')) trace(f, where=asNamespace('stats'), tracer=quote(guard()), print=FALSE)",
    "trace('prep.recipe', where=asNamespace('recipes'), tracer=quote(guard()), print=FALSE)",
    "for (f in c('initial_split', 'vfold_cv')) trace(f, where=asNamespace('rsample'), tracer=quote(guard()), print=FALSE)",
    paste0("snapshot <- ", paste(deparse(snapshot), collapse = "\n")),
    paste0("records <- readRDS(", deparse(saved), ")"),
    "for (r in records) stopifnot(identical(r$expected, snapshot(r$fit, r$new, r$schema)))"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

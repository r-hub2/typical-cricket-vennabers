testthat::test_that("numeric ensembles feed mean complete pair predictions into IVAR", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:60), truth = 10 + 3 * sin(1:60) + cos(1:60))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  fit <- generics::fit(calibrated_workflow(base, "truth", seed = 42,
    postprocessing = list(method = "isotonic_numeric", strategy = "ensemble", v = 3L)), d)
  state <- workflows::extract_fit_engine(fit)
  ensemble <- state$models[[1L]]
  pairs <- vapply(ensemble$models, function(model) predict(model, d["x"])$.pred, numeric(nrow(d)))
  averaged <- predict(ensemble, d["x"], type = "numeric")$.pred
  testthat::expect_equal(averaged, rowMeans(pairs), tolerance = 1e-10)
  core <- venn_abers_regression(state$calibration[[1]]$truth, state$calibration[[1]]$prediction)
  # Freeze identical VA inputs: last-bit mean differences can change exact ties.
  testthat::expect_equal(predict(fit, d["x"])$.pred, predict(core, averaged)$prediction, tolerance = 1e-10)
})

ensemble_regression_data <- function() data.frame(x = sin(1:96), z = cos((1:96) / 4),
  truth = 10 + 3 * sin(1:96) + cos(1:96) + sin((1:96) / 4))
ensemble_regression_base <- function(d) workflows::workflow() |>
  workflows::add_recipe(recipes::recipe(truth ~ ., d) |> recipes::step_normalize(recipes::all_numeric_predictors())) |>
  workflows::add_model(parsnip::linear_reg(engine = "lm"))
ensemble_regression_plan <- function(data) {
  n <- nrow(data)
  # Pure local plan: different Q sizes give J=2 or J=3, with unequal P sizes.
  held <- if (n > 60L) list(1:12, 13:36, 37:n) else list(1:18, 19:n)
  rsample::manual_rset(lapply(held, function(p) rsample::make_splits(
    list(analysis = setdiff(seq_len(n), p), assessment = p), data)), paste0("inner", seq_along(held)))
}

testthat::test_that("IVAR and CVAR average numeric pair outputs before VA across scales and moderation", {
  skip_if_no_integration("tailor", "probably")
  close <- function(actual, expected, s) {
    testthat::expect_identical(dim(actual), dim(expected))
    testthat::expect_true(all(is.finite(actual)) && all(is.finite(expected)))
    testthat::expect_true(all(abs(actual - expected) <= 1e-10 * s + 1e-10 * abs(expected)))
  }
  for (s in c(1e-9, 1, 1e9)) for (method in c("inductive", "cross")) for (post_method in c("tailor_linear", "isotonic_numeric")) {
    d <- ensemble_regression_data(); d$truth <- d$truth * s
    base <- ensemble_regression_base(d)
    held_rows <- list(1:16, 17:48, 49:96)
    splits <- lapply(held_rows, function(C) rsample::make_splits(list(analysis = setdiff(1:96, C), assessment = C), d))
    va <- if (method == "cross") rsample::manual_rset(splits, c("one", "two", "three")) else splits[[1L]]
    for (option in list(list(m = 1), list(m = 7), list(epsilon = .875))) {
      fit <- generics::fit(do.call(calibrated_workflow, c(list(workflow = base, outcome = "truth",
        method = method, splits = va, seed = 42,
        postprocessing = list(method = post_method, strategy = "ensemble", splits = ensemble_regression_plan)), option)), d)
      state <- workflows::extract_fit_engine(fit)
      held <- new <- truth <- ids <- list()
      for (id in names(state$models)) {
        ensemble <- state$models[[id]]; plans <- state$postprocessing[[id]]$inner_folds
        testthat::expect_identical(names(ensemble$models), names(plans))
        testthat::expect_identical(state$postprocessing[[id]]$observed_base_fits, length(plans))
        testthat::expect_identical(state$postprocessing[[id]]$observed_postprocessor_fits, length(plans))
        testthat::expect_identical(sort(unlist(lapply(plans, `[[`, "P"), use.names = FALSE)), sort(plans[[1]]$Q))
        calibration <- predictions <- list()
        for (j in names(plans)) {
          p <- plans[[j]]
          testthat::expect_identical(p$inner_fold, j)
          testthat::expect_identical(p$input, "numeric")
          testthat::expect_identical(p$method, post_method)
          testthat::expect_setequal(p$B, setdiff(p$Q, p$P))
          testthat::expect_length(intersect(p$B, p$P), 0L)
          testthat::expect_length(intersect(p$Q, p$C), 0L)
          if (post_method == "tailor_linear") {
            ref <- generics::fit(workflows::add_tailor(base, tailor::tailor() |>
              tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE)), d[p$B, ], data_calibration = d[p$P, ])
            reference <- function(data) predict(ref, data)$.pred
            testthat::expect_identical(p$rank, 2L)
          } else {
            ref <- generics::fit(base, d[p$B, ])
            input <- predict(ref, d[p$P, ])$.pred
            # Independent min-max characterization of increasing least squares.
            # Each tied input's observations keep their original multiplicity.
            knots <- sort(unique(input))
            values <- vapply(seq_along(knots), function(k) max(vapply(seq_len(k), function(a)
              min(vapply(seq.int(k, length(knots)), function(b)
                mean(d$truth[p$P][input >= knots[a] & input <= knots[b]]), numeric(1))), numeric(1))), numeric(1))
            reference <- function(data) stats::approx(knots, values, xout = predict(ref, data)$.pred, rule = 2)$y
            testthat::expect_null(p$rank)
          }
          calibration[[j]] <- reference(d[p$C, ])
          predictions[[j]] <- reference(d)
          actual <- predict(ensemble$models[[j]], d, type = "numeric")$.pred
          if (post_method == "isotonic_numeric") testthat::expect_true(all(
            abs(actual - predictions[[j]]) <= 1e-10 + 1e-10 * abs(predictions[[j]])))
          close(actual, predictions[[j]], s)
        }
        held[[id]] <- rowMeans(do.call(cbind, calibration))
        new[[id]] <- rowMeans(do.call(cbind, predictions))
        sizes <- vapply(plans, `[[`, 1L, "P_size")
        weighted <- drop(do.call(cbind, predictions) %*% (sizes / sum(sizes)))
        testthat::expect_gt(max(abs(new[[id]] - weighted)) / s, 1e-5)
        close(predict(ensemble, d)$.pred, new[[id]], s)
        close(state$calibration[[id]]$prediction, held[[id]], s)
        truth[[id]] <- d$truth[plans[[1]]$C]; ids[[id]] <- rep(id, length(truth[[id]]))
      }
      if (method == "cross") testthat::expect_identical(lengths(lapply(state$postprocessing, `[[`, "inner_folds")), c(one = 3L, two = 3L, three = 2L))
      core <- if (method == "cross") do.call(cross_venn_abers_regression,
        c(list(truth = unlist(truth), prediction = unlist(held), fold = unlist(ids)), option)) else
        do.call(venn_abers_regression, c(list(truth = truth[[1]], prediction = held[[1]]), option))
      expected <- if (method == "cross") predict(core, do.call(cbind, new), type = "details") else predict(core, new[[1]])
      actual <- venn_abers_details(fit, d[c("x", "z")])
      testthat::expect_identical(actual$prediction, predict(fit, d[c("x", "z")]))
      testthat::expect_null(state$class_levels)
      testthat::expect_true(all(actual$details$postprocessing_input == "numeric"))
      for (id in names(state$models)) {
        rows <- actual$details$fold == id
        plans <- state$postprocessing[[id]]$inner_folds
        testthat::expect_true(all(vapply(actual$details$inner_fold[rows], identical, logical(1), names(plans))))
        for (size in c("B_size", "P_size")) testthat::expect_true(all(vapply(actual$details[[size]][rows],
          identical, logical(1), vapply(plans, `[[`, 1L, size))))
      }
      if (method == "cross") {
        close(actual$prediction$.pred, expected$prediction$prediction, s)
        for (column in c("lower", "upper")) close(actual$aggregate[[column]], expected$prediction[[column]], s)
        for (id in names(expected$folds)) for (column in c("prediction", "lower", "upper"))
          close(actual$details[actual$details$fold == id, column][[1L]], expected$folds[[id]][[column]], s)
      } else close(as.matrix(actual$details[c("prediction", "lower", "upper")]), as.matrix(expected), s)
      for (n in 0:1) {
        small <- venn_abers_details(fit, d[seq_len(n), c("x", "z")])
        testthat::expect_identical(names(small$details), names(actual$details))
        testthat::expect_identical(nrow(small$prediction), n)
        if (n) close(small$prediction$.pred, actual$prediction$.pred[1L], s)
      }
    }
  }
})

testthat::test_that("numeric ensembles count complete fits and predictions without relearning", {
  skip_if_no_integration("tailor", "probably")
  d <- ensemble_regression_data(); base <- ensemble_regression_base(d)
  fit_method <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  predict_method <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  tailor_method <- getS3method("predict", "tailor", envir = asNamespace("stats"))
  isotonic_method <- getS3method("predict", "vennabers_isotonic_workflow")
  audit <- new.env(); audit$fits <- audit$posts <- audit$predictions <- audit$adjustments <- 0L
  audit$broken <- FALSE; audit$draws <- numeric()
  counted_fit <- function(object, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) {
      audit$fits <- audit$fits + 1L; audit$draws <- c(audit$draws, runif(1))
    }
    fit_method(object, ...)
  }
  counted_predict <- function(object, ...) {
    p <- predict_method(object, ...)
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) {
      audit$predictions <- audit$predictions + 1L
      if (audit$broken && audit$predictions == 2L) p$.pred[1L] <- Inf
    }
    p
  }
  registerS3method("fit", "workflow", counted_fit, envir = asNamespace("generics"))
  registerS3method("predict", "workflow", counted_predict, envir = asNamespace("stats"))
  registerS3method("predict", "tailor", function(object, ...) {
    audit$adjustments <- audit$adjustments + 1L; tailor_method(object, ...)
  }, envir = asNamespace("stats"))
  counted_isotonic <- function(object, ...) {
    audit$adjustments <- audit$adjustments + 1L; isotonic_method(object, ...)
  }
  testthat::local_mocked_bindings(predict.vennabers_isotonic_workflow = counted_isotonic, .package = "vennabers")
  registerS3method("predict", "vennabers_isotonic_workflow", counted_isotonic, envir = asNamespace("stats"))
  on.exit({
    registerS3method("fit", "workflow", fit_method, envir = asNamespace("generics"))
    registerS3method("predict", "workflow", predict_method, envir = asNamespace("stats"))
    registerS3method("predict", "tailor", tailor_method, envir = asNamespace("stats"))
    registerS3method("predict", "vennabers_isotonic_workflow", isotonic_method, envir = asNamespace("stats"))
  }, add = TRUE)
  linear <- probably::cal_estimate_linear; isotonic <- stats::isoreg
  testthat::local_mocked_bindings(cal_estimate_linear = function(...) {
    audit$posts <- audit$posts + 1L; linear(...)
  }, .package = "probably")
  testthat::local_mocked_bindings(isoreg = function(...) {
    audit$posts <- audit$posts + 1L; isotonic(...)
  }, .package = "stats")
  for (method in c("inductive", "cross")) for (post_method in c("tailor_linear", "isotonic_numeric")) for (J in 2:3) {
    spec <- calibrated_workflow(base, "truth", method = method, seed = 24,
      postprocessing = list(method = post_method, strategy = "ensemble", v = J))
    expected_draws <- withr::with_seed(24, {
      outer <- if (method == "cross") rsample::vfold_cv(d, v = 5)$splits else list(rsample::initial_split(d, .8))
      for (q in outer) rsample::vfold_cv(rsample::analysis(q), v = J)
      runif(length(outer) * J)
    })
    audit$fits <- audit$posts <- 0L; audit$draws <- numeric()
    set.seed(17); before <- .Random.seed
    fit <- generics::fit(spec, d)
    count <- if (method == "cross") 5L * J else J
    testthat::expect_identical(.Random.seed, before)
    testthat::expect_identical(audit$draws, expected_draws)
    testthat::expect_identical(audit$fits, count)
    testthat::expect_identical(audit$posts, count)
    state <- workflows::extract_fit_engine(fit)
    testthat::expect_output(print(state), paste0("Models per nonempty prediction: ", count))
    expected <- predict(fit, d[c("x", "z")])
    for (n in c(0L, 1L, nrow(d))) {
      audit$predictions <- audit$adjustments <- 0L
      result <- venn_abers_details(fit, d[seq_len(n), c("x", "z")])$prediction
      testthat::expect_identical(audit$predictions, if (n) count else 0L)
      testthat::expect_identical(audit$adjustments, if (n) count else 0L)
      testthat::expect_equal(result, expected[seq_len(n), ], tolerance = 1e-10)
      testthat::expect_identical(audit$fits, count)
      testthat::expect_identical(audit$posts, count)
    }
    ensemble <- state$models[[1]]
    testthat::expect_identical(predict(ensemble, d[FALSE, ])$.pred, numeric())
    testthat::expect_error(predict(ensemble, d, type = "prob"), "type")
    for (remove in c(TRUE, FALSE)) {
      missing <- ensemble
      if (remove) missing$models[[2]] <- NULL else missing$models[2] <- list(NULL)
      testthat::expect_error(predict(missing, d), "incomplete inner ensemble")
    }
    audit$broken <- TRUE; audit$predictions <- 0L
    testthat::expect_error(predict(fit, d), "Fold2.*finite")
    audit$broken <- FALSE
  }
})

testthat::test_that("each numeric pair retains method-specific size and rank rules", {
  skip_if_no_integration("tailor", "probably")
  d <- data.frame(x = rep(1, 16), truth = c(2, 8, 4, 6, 1:12))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  va <- function(data) rsample::make_splits(list(analysis = 1:4, assessment = 5:16), data)
  plan <- function(data) rsample::manual_rset(lapply(list(1:2, 3:4), function(p)
    rsample::make_splits(list(analysis = setdiff(1:4, p), assessment = p), data)), c("first", "second"))
  post <- list(method = "isotonic_numeric", strategy = "ensemble", splits = plan)
  fit <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d))
  model <- workflows::extract_fit_engine(fit)$models[[1]]
  testthat::expect_equal(suppressWarnings(predict(model, data.frame(x = c(-100, 0, 100))))$.pred, rep(5, 3), tolerance = 1e-10)
  flat <- d; flat$truth[1:4] <- 12
  constant <- suppressWarnings(generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), flat))
  testthat::expect_equal(suppressWarnings(predict(workflows::extract_fit_engine(constant)$models[[1]], d))$.pred, rep(12, 16), tolerance = 1e-10)
  post$method <- "tailor_linear"
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d), "first tailor_linear B/P plan.*at least 3")
  va <- function(data) rsample::make_splits(list(analysis = 1:6, assessment = 7:16), data)
  post$splits <- function(data) rsample::manual_rset(lapply(list(1:3, 4:6), function(p)
    rsample::make_splits(list(analysis = setdiff(1:6, p), assessment = p), data)), c("first", "second"))
  testthat::expect_error(suppressWarnings(generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d)), "first tailor_linear.*full design rank 2")
  d$x <- 1:16
  varying <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  testthat::expect_s3_class(generics::fit(calibrated_workflow(varying, "truth", splits = va, postprocessing = post), d), "workflow")
})

testthat::test_that("numeric ensembles reject invalid local plans before learning and abort failed pairs", {
  skip_if_no_integration("tailor", "probably")
  d <- ensemble_regression_data(); base <- ensemble_regression_base(d)
  va <- function(data) rsample::make_splits(list(analysis = 1:72, assessment = 73:96), data)
  post <- list(method = "tailor_linear", strategy = "ensemble", splits = ensemble_regression_plan)
  fit_method <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  fits <- 0L
  counted <- function(object, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) fits <<- fits + 1L
    fit_method(object, ...)
  }
  registerS3method("fit", "workflow", counted, envir = asNamespace("generics"))
  on.exit(registerS3method("fit", "workflow", fit_method, envir = asNamespace("generics")), add = TRUE)
  for (method in c("tailor_linear", "isotonic_numeric")) {
    check <- function(provider, pattern = "B/P plan") {
      bad <- post; bad$method <- method; bad$splits <- provider
      testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = bad), d), pattern)
      testthat::expect_identical(fits, 0L)
    }
    check(function(data) rsample::vfold_cv(data, v = 3, repeats = 2))
    check(function(data) rsample::bootstraps(data, times = 3))
    check(function(data) { p <- ensemble_regression_plan(data); p$splits[[3]] <- p$splits[[2]]; p })
    check(function(data) ensemble_regression_plan(data)[1:2, ])
    check(function(data) { p <- ensemble_regression_plan(data); p$id[2] <- p$id[1]; p })
    for (mutate in list(function(x) { x$x[1] <- 99; x }, function(x) x[nrow(x):1, ],
        function(x) { x$x <- as.integer(x$x); x }, function(x) x[rev(names(x))]))
      check(function(data) ensemble_regression_plan(mutate(data)), "binding")
    check(function(data) ensemble_regression_plan(d), "binding")
    for (indices in list(list(analysis = 12:72, assessment = 1:12),
        list(analysis = c(13:72, 13L), assessment = 1:12), list(analysis = 14:72, assessment = 1:12))) {
      check(function(data) { p <- ensemble_regression_plan(data); p$splits[[3]]$in_id <- indices$analysis
        p$splits[[3]]$out_id <- indices$assessment; p })
    }
    check(function(data) rsample::manual_rset(lapply(list(1L, 2:72), function(p)
      rsample::make_splits(list(analysis = setdiff(1:72, p), assessment = p), data)), c("tiny", "rest")), "tiny.*at least")
    bad <- post; bad$method <- method
    bad$splits <- function(data) if (nrow(data) == 48L) NULL else ensemble_regression_plan(data)
    outer <- rsample::manual_rset(lapply(list(1:16, 17:48, 49:96), function(c)
      rsample::make_splits(list(analysis = setdiff(1:96, c), assessment = c), d)), c("one", "two", "three"))
    testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", method = "cross", splits = outer, postprocessing = bad), d), "three.*B/P plan")
    testthat::expect_identical(fits, 0L)
  }
  original <- probably::cal_estimate_linear
  for (fault in c("none", "missing", "nonconverged", "nonfinite", "warning")) local({
    calls <- 0L
    replacement <- function(...) {
      calls <<- calls + 1L
      if (calls == 2L && fault == "none") return(probably::cal_estimate_none(...))
      result <- original(...)
      if (calls == 2L) {
        if (fault == "missing") result$estimates <- list()
        if (fault == "nonconverged") result$estimates[[1]]$estimate$converged <- FALSE
        if (fault == "nonfinite") result$estimates[[1]]$estimate$coefficients[] <- Inf
        if (fault == "warning") warning("instrumented successful warning")
      }
      result
    }
    testthat::local_mocked_bindings(cal_estimate_linear = replacement, .package = "probably")
    spec <- calibrated_workflow(base, "truth", splits = va, postprocessing = post)
    if (fault == "warning") {
      testthat::expect_warning(valid <- generics::fit(spec, d), "inner2 tailor_linear.*successful warning")
      testthat::expect_match(workflows::extract_fit_engine(valid)$postprocessing[[1]]$inner_folds[[2]]$warnings, "successful warning")
    } else {
      testthat::expect_error(suppressWarnings(generics::fit(spec, d)), "Inductive1 / inner2 tailor_linear.*P")
      testthat::expect_identical(calls, 2L)
    }
  })
  original_iso <- stats::isoreg; calls <- 0L
  testthat::local_mocked_bindings(isoreg = function(...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("instrumented isotonic failure")
    original_iso(...)
  }, .package = "stats")
  post$method <- "isotonic_numeric"
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d), "Inductive1 / inner2 isotonic_numeric.*instrumented isotonic failure")
  testthat::expect_identical(calls, 2L)
})

testthat::test_that("every numeric inner B P and C boundary excludes local held-out learning", {
  skip_if_no_integration("tailor", "probably")
  d <- ensemble_regression_data(); base <- ensemble_regression_base(d)
  for (method in c("inductive", "cross")) for (post_method in c("tailor_linear", "isotonic_numeric")) {
    va <- function(data) {
      splits <- lapply(list(1:32, 33:64, 65:96), function(c)
        rsample::make_splits(list(analysis = setdiff(1:96, c), assessment = c), data))
      if (method == "cross") rsample::manual_rset(splits, c("one", "two", "three")) else splits[[1]]
    }
    spec <- calibrated_workflow(base, "truth", method = method, splits = va, seed = 24,
      postprocessing = list(method = post_method, strategy = "ensemble", splits = ensemble_regression_plan))
    fit <- function(data) workflows::extract_fit_engine(generics::fit(spec, data))
    original <- fit(d)
    base_fit <- function(model) if (post_method == "isotonic_numeric") model$workflow else model
    post_fit <- function(model) if (post_method == "isotonic_numeric") model$isotonic else
      coef(workflows::extract_tailor(model)$adjustments[[1]]$results$fit$estimates[[1]]$estimate)
    for (j in names(original$models[[1]]$models)) {
      p <- original$postprocessing[[1]]$inner_folds[[j]]; model <- original$models[[1]]$models[[j]]
      for (role in c("B", "P", "C")) for (column in c("x", "truth")) {
        changed <- d
        changed[p[[role]], column] <- if (column == "x") changed[p[[role]], column]^2 + .5 else changed[p[[role]], column] + .5
        after <- fit(changed); other <- after$models[[1]]$models[[j]]
        if (role == "B") {
          testthat::expect_false(identical(coef(workflows::extract_fit_engine(base_fit(model))), coef(workflows::extract_fit_engine(base_fit(other)))))
          if (column == "x") testthat::expect_false(identical(workflows::extract_recipe(base_fit(model))$steps, workflows::extract_recipe(base_fit(other))$steps))
        } else {
          testthat::expect_identical(coef(workflows::extract_fit_engine(base_fit(model))), coef(workflows::extract_fit_engine(base_fit(other))))
          testthat::expect_identical(workflows::extract_recipe(base_fit(model))$steps, workflows::extract_recipe(base_fit(other))$steps)
        }
        if (role == "P") testthat::expect_false(identical(post_fit(model), post_fit(other)))
        if (role == "C") {
          testthat::expect_identical(post_fit(model), post_fit(other))
          testthat::expect_identical(original$postprocessing[[1]], after$postprocessing[[1]])
          testthat::expect_false(identical(original$calibration[[1]], after$calibration[[1]]))
        }
      }
    }
  }
})

testthat::test_that("numeric local plans reproduce across call order and accept static matching data", {
  skip_if_no_integration("tailor", "probably")
  d <- ensemble_regression_data(); base <- ensemble_regression_base(d)
  va <- function(data) rsample::make_splits(list(analysis = seq_len(nrow(data) - 16L),
    assessment = seq.int(nrow(data) - 15L, nrow(data))), data)
  for (method in c("tailor_linear", "isotonic_numeric")) {
    post <- list(method = method, strategy = "ensemble", splits = ensemble_regression_plan)
    spec <- calibrated_workflow(base, "truth", splits = va, seed = 42, postprocessing = post)
    inputs <- list(d, d[1:72, ])
    forward <- lapply(inputs, function(data) generics::fit(spec, data))
    backward <- rev(lapply(rev(inputs), function(data) generics::fit(spec, data)))
    for (i in seq_along(inputs)) {
      testthat::expect_identical(predict(forward[[i]], d), predict(backward[[i]], d))
      testthat::expect_identical(workflows::extract_fit_engine(forward[[i]])$postprocessing,
        workflows::extract_fit_engine(backward[[i]])$postprocessing)
      static <- post; static$splits <- ensemble_regression_plan(rsample::analysis(va(inputs[[i]])))
      ref <- generics::fit(calibrated_workflow(base, "truth", splits = va(inputs[[i]]), seed = 42, postprocessing = static), inputs[[i]])
      testthat::expect_identical(predict(ref, d), predict(forward[[i]], d))
    }
  }
})

testthat::test_that("numeric ensemble means preserve finite extremes and cancellation", {
  skip_if_no_integration("tailor", "probably")
  d <- ensemble_regression_data()
  fit <- generics::fit(calibrated_workflow(ensemble_regression_base(d), "truth", seed = 42,
    postprocessing = list(method = "tailor_linear", strategy = "ensemble", v = 3L)), d)
  ensemble <- workflows::extract_fit_engine(fit)$models[[1]]
  # Isolate averaging with fixed outputs at the complete owned prediction seam.
  column <- 0L
  values <- rbind(c(1e308, 1e308, 1e308), c(1e16, 1, -1e16), c(-1e308, -1e308, -1e308))
  supplied <- function(object, new_data, ...) {
    column <<- column + 1L
    tibble::tibble(.pred = values[, column])
  }
  testthat::local_mocked_bindings(predict.vennabers_tailor_workflow = supplied, .package = "vennabers")
  actual <- predict(ensemble, d[1:3, ])$.pred
  testthat::expect_true(all(is.finite(actual)))
  testthat::expect_equal(actual[c(1, 3)], c(1e308, -1e308), tolerance = 1e-15)
  testthat::expect_equal(actual[2], 1/3, tolerance = 1e-15)
})

testthat::test_that("linear Tailor preserves public workflow and IVAR/CVAR parity across scales", {
  for (pkg in c("workflows", "tailor", "probably", "rsample")) testthat::skip_if_not_installed(pkg)
  # Predeclared positive scale s is the unit multiplier, independent of results.
  close <- function(actual, expected, s) {
    testthat::expect_identical(dim(actual), dim(expected))
    testthat::expect_true(all(is.finite(actual)) && all(is.finite(expected)))
    testthat::expect_true(all(abs(actual - expected) <= 1e-10 * s + 1e-10 * abs(expected)))
  }
  calls <- new.env(); calls$n <- 0L
  original <- getS3method("predict", "workflow", envir = asNamespace("workflows"))
  counted <- function(object, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) calls$n <- calls$n + 1L
    original(object, ...)
  }
  registerS3method("predict", "workflow", counted, envir = asNamespace("stats"))
  on.exit(registerS3method("predict", "workflow", original, envir = asNamespace("stats")), add = TRUE)
  for (s in c(1e-9, 1, 1e9)) for (shift in c(0, 100)) for (method in c("inductive", "cross")) {
    d <- data.frame(x = sin(1:120), truth = s * (shift + 3 * sin(1:120) + cos(1:120)))
    base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
      workflows::add_model(parsnip::linear_reg(engine = "lm"))
    held_rows <- list(1:15, 16:60, 61:120)
    splits <- lapply(held_rows, function(C) rsample::make_splits(list(analysis = setdiff(1:120, C), assessment = C), d))
    plan <- if (method == "cross") rsample::manual_rset(splits, c("a", "b", "c")) else splits[[1L]]
    # C's smallest size is 15: m=7 and epsilon=.875 are its upper valid boundary.
    for (option in list(list(m = 1), list(m = 7), list(epsilon = .875))) {
      fit <- generics::fit(do.call(calibrated_workflow, c(list(workflow = base, outcome = "truth",
        method = method, splits = plan, seed = 42,
        postprocessing = list(method = "tailor_linear", strategy = "holdout", prop = .3)), option)), d)
      state <- workflows::extract_fit_engine(fit)
      held <- new <- truth <- ids <- list()
      for (id in names(state$models)) {
        p <- state$postprocessing[[id]]
        ref <- generics::fit(workflows::add_tailor(base, tailor::tailor() |>
          tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE)),
          d[p$B, ], data_calibration = d[p$P, ])
        held[[id]] <- predict(ref, d[p$C, ])$.pred
        new[[id]] <- predict(ref, d)$.pred
        truth[[id]] <- d$truth[p$C]; ids[[id]] <- rep(id, length(p$C))
        close(held[[id]], state$calibration[[id]]$prediction, s)
        close(new[[id]], predict(state$models[[id]], d)$.pred, s)
        testthat::expect_identical(p$rank, 2L)
        testthat::expect_identical(p$input, "numeric")
        testthat::expect_setequal(c(p$B, p$P), p$Q)
        testthat::expect_length(intersect(p$B, p$P), 0L)
        testthat::expect_length(intersect(p$Q, p$C), 0L)
      }
      core <- if (method == "cross") do.call(cross_venn_abers_regression,
        c(list(truth = unlist(truth), prediction = unlist(held), fold = unlist(ids)), option)) else
        do.call(venn_abers_regression, c(list(truth = truth[[1L]], prediction = held[[1L]]), option))
      expected <- if (method == "cross") predict(core, do.call(cbind, new), type = "details") else predict(core, new[[1L]])
      calls$n <- 0L
      actual <- venn_abers_details(fit, d["x"])
      testthat::expect_identical(calls$n, length(state$models))
      testthat::expect_true(all(actual$details$postprocessing_input == "numeric"))
      testthat::expect_null(state$class_levels)
      testthat::expect_identical(actual$prediction, predict(fit, d["x"]))
      if (method == "cross") {
        close(actual$prediction$.pred, expected$prediction$prediction, s)
        close(as.matrix(actual$aggregate[c("lower", "upper")]), as.matrix(expected$prediction[c("lower", "upper")]), s)
        for (id in names(expected$folds)) close(as.matrix(actual$details[actual$details$fold == id,
          c("prediction", "lower", "upper")]), as.matrix(expected$folds[[id]]), s)
        # Independent equal-weight means for points and both bounds, despite unequal C sizes.
        for (column in c("prediction", "lower", "upper")) {
          mean <- rowMeans(vapply(expected$folds, `[[`, numeric(nrow(d)), column))
          value <- if (column == "prediction") actual$prediction$.pred else actual$aggregate[[column]]
          close(value, mean, s)
        }
      } else close(as.matrix(actual$details[c("prediction", "lower", "upper")]), as.matrix(expected), s)
      # Rows of d are calibration points, so their scores are knots, where BLAS
      # rounding that depends on the row count (ATLAS) moves them across the step.
      query <- data.frame(x = .5)
      batch <- venn_abers_details(fit, rbind(d["x"], query))$prediction$.pred[nrow(d) + 1L]
      for (n in 0:1) {
        calls$n <- 0L
        small <- venn_abers_details(fit, query[seq_len(n), , drop = FALSE])
        testthat::expect_identical(calls$n, if (n) length(state$models) else 0L)
        testthat::expect_identical(names(small$details), names(actual$details))
        testthat::expect_identical(nrow(small$prediction), n)
        if (n) close(small$prediction$.pred, batch, s)
      }
      for (type in c("conf_int", "pred_int")) testthat::expect_error(predict(fit, d, type = type), "type")
    }
  }
})

testthat::test_that("numeric zero/one truth stays regression and automatic P plans are unstratified", {
  for (pkg in c("workflows", "tailor", "probably", "rsample")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:120), truth = rep(c(0, 1), 60))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  for (method in c("inductive", "cross")) {
    spec <- calibrated_workflow(base, "truth", method = method, seed = 42,
      postprocessing = list(method = "tailor_linear", strategy = "holdout", prop = .3))
    set.seed(987); before <- .Random.seed
    fit <- generics::fit(spec, d)
    testthat::expect_identical(.Random.seed, before)
    state <- workflows::extract_fit_engine(fit)
    expected <- withr::with_seed(42, {
      va <- if (method == "cross") rsample::vfold_cv(d, v = 5)$splits else list(rsample::initial_split(d, .8))
      lapply(va, function(s) rsample::initial_split(rsample::analysis(s), .7))
    })
    for (i in seq_along(expected)) {
      testthat::expect_identical(state$postprocessing[[i]]$analysis, as.integer(expected[[i]], "analysis"))
      testthat::expect_identical(state$postprocessing[[i]]$assessment, as.integer(expected[[i]], "assessment"))
    }
    testthat::expect_null(state$class_levels)
    testthat::expect_type(predict(fit, d)$.pred, "double")
  }
})

testthat::test_that("linear holdouts separate early plan errors from later failed fits", {
  for (pkg in c("workflows", "tailor", "probably", "rsample")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:120), truth = 10 + 3 * sin(1:120) + cos(1:120))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  va <- function(data) rsample::make_splits(list(analysis = 1:90, assessment = 91:120), data)
  post <- list(method = "tailor_linear", strategy = "holdout",
    splits = function(data) rsample::make_splits(list(analysis = 1:60, assessment = 61:90), data))
  spec <- calibrated_workflow(base, "truth", splits = va, postprocessing = post)
  fit <- generics::fit(spec, d)
  canonical <- tailor::tailor() |> tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE)
  adopted <- generics::fit(calibrated_workflow(workflows::add_tailor(base, canonical), "truth",
    splits = va, postprocessing = post), d)
  testthat::expect_identical(predict(adopted, d), predict(fit, d))
  duplicate <- canonical; duplicate$adjustments <- rep(duplicate$adjustments, 2L)
  for (bad in list(tailor::tailor() |> tailor::adjust_numeric_calibration(method = "linear", smooth = TRUE),
      tailor::tailor() |> tailor::adjust_numeric_calibration(method = "none"),
      duplicate,
      workflows::extract_tailor(workflows::extract_fit_engine(fit)$models[[1L]]))) {
    testthat::expect_error(calibrated_workflow(workflows::add_tailor(base, bad), "truth", postprocessing = post), "[Tt]ailor")
  }
  constant <- d; constant$x <- 1
  testthat::expect_error(suppressWarnings(generics::fit(spec, constant)), "tailor_linear.*full design rank 2")
  original <- probably::cal_estimate_linear
  for (fault in c("none", "missing", "nonconverged", "nonfinite", "warning")) local({
    replacement <- function(...) {
      if (fault == "none") return(probably::cal_estimate_none(...))
      result <- original(...)
      if (fault == "missing") result$estimates <- list()
      if (fault == "nonconverged") result$estimates[[1L]]$estimate$converged <- FALSE
      if (fault == "nonfinite") result$estimates[[1L]]$estimate$coefficients[] <- Inf
      if (fault == "warning") warning("instrumented successful warning")
      result
    }
    testthat::local_mocked_bindings(cal_estimate_linear = replacement, .package = "probably")
    if (fault == "warning") {
      testthat::expect_warning(valid <- generics::fit(spec, d), "Holdout1 tailor_linear.*instrumented successful warning")
      testthat::expect_match(workflows::extract_fit_engine(valid)$postprocessing[[1L]]$warnings, "successful warning")
    } else testthat::expect_error(suppressWarnings(generics::fit(spec, d)), "Holdout1 tailor_linear")
  })
  original_predict <- getS3method("predict", "model_fit", envir = asNamespace("parsnip"))
  for (value in c(Inf, NA_real_, NaN)) local({
    corrupt <- function(object, new_data, type = NULL, opts = list(), ...) {
      p <- original_predict(object, new_data, type = type, opts = opts, ...)
      if (".pred" %in% names(p) && !inherits(object$spec, "vennabers")) p$.pred[1L] <- value
      p
    }
    testthat::local_mocked_bindings(predict.model_fit = corrupt, .package = "parsnip")
    registerS3method("predict", "model_fit", corrupt, envir = asNamespace("stats"))
    on.exit(registerS3method("predict", "model_fit", original_predict, envir = asNamespace("stats")), add = TRUE)
    testthat::expect_error(generics::fit(spec, d), "Tailor numeric input")
    testthat::expect_error(predict(fit, d), "Tailor numeric input")
  })
  # Invalid plans must fail before the external workflow training boundary.
  calls <- new.env(); calls$n <- 0L
  original_fit <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  guard <- function(object, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) calls$n <- calls$n + 1L
    original_fit(object, ...)
  }
  registerS3method("fit", "workflow", guard, envir = asNamespace("generics"))
  on.exit(registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics")), add = TRUE)
  check <- function(provider, pattern = "B/P plan") {
    bad <- post; bad$splits <- provider
    testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = bad), d), pattern)
    testthat::expect_identical(calls$n, 0L)
  }
  check(function(data) NULL)
  for (P in list(integer(), 90L, 89:90)) check(function(data)
    rsample::make_splits(list(analysis = setdiff(1:90, P), assessment = P), data))
  for (positions in list(list(analysis = 1:60, assessment = 60:90), list(analysis = 1:60, assessment = 62:90))) {
    check(function(data) rsample::make_splits(positions, data))
  }
  for (mutate in list(function(x) { x$x[1L] <- 999; x }, function(x) x[nrow(x):1, ],
      function(x) { x$x <- as.integer(x$x); x }, function(x) x[rev(names(x))])) {
    check(function(data) rsample::make_splits(list(analysis = 1:60, assessment = 61:90), mutate(data)), "binding")
  }
  bad_post <- post; bad_post$method <- "tailor_logistic"
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", postprocessing = bad_post), d), "incompatible")
  for (option in list(list(m = 16), list(epsilon = .01), list(epsilon = .999))) {
    testthat::expect_error(generics::fit(do.call(calibrated_workflow,
      c(list(workflow = base, outcome = "truth", splits = va, postprocessing = post), option)), d), "calibration size")
  }
  for (value in c(Inf, NA_real_, NaN)) {
    bad <- d; bad$truth[80L] <- value
    testthat::expect_error(generics::fit(spec, bad), "regression truth")
  }
  testthat::expect_identical(calls$n, 0L)
  # Even a bad last B/P plan in Cross prevents the first B fit.
  cross_post <- post
  cross_post$splits <- function(data) {
    if (!1 %in% as.integer(rownames(data))) return(NULL)
    rsample::initial_split(data, .7)
  }
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", method = "cross", seed = 42,
    postprocessing = cross_post), d), "B/P plan")
  testthat::expect_identical(calls$n, 0L)
  testthat::expect_error(suppressWarnings(generics::fit(spec, constant)), "full design rank 2")
  testthat::expect_identical(calls$n, 1L)
  # Exactly three varying P rows suffice for the two-coefficient fit.
  post$splits <- function(data) rsample::make_splits(list(analysis = 1:87, assessment = 88:90), data)
  testthat::expect_s3_class(generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d), "workflow")
})

testthat::test_that("B/P/C mutations obey the regression learning roles", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:120), truth = 10 + 3 * sin(1:120) + cos(1:120))
  base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, d) |> recipes::step_normalize(x)) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  for (method in c("inductive", "cross")) {
    original <- workflows::extract_fit_engine(generics::fit(calibrated_workflow(base, "truth", method = method, seed = 42,
      postprocessing = list(method = "tailor_linear", strategy = "holdout", prop = .3)), d))
    va <- function(data) {
      splits <- lapply(original$folds, function(p) rsample::make_splits(p[c("analysis", "assessment")], data))
      if (method == "cross") rsample::manual_rset(splits, names(splits)) else splits[[1L]]
    }
    post <- list(method = "tailor_linear", strategy = "holdout", splits = function(data) {
      i <- which(vapply(original$postprocessing, function(p) identical(rownames(data), rownames(d)[p$Q]), logical(1)))
      rsample::make_splits(original$postprocessing[[i]][c("analysis", "assessment")], data)
    })
    for (role in c("B", "P", "C")) for (column in c("x", "truth")) {
      changed <- d
      rows <- original$postprocessing[[1L]][[role]]
      changed[rows, column] <- changed[rows, column] + .5
      other <- workflows::extract_fit_engine(generics::fit(calibrated_workflow(base, "truth", method = method,
        splits = va, postprocessing = post), changed))
      engine <- function(state) stats::coef(workflows::extract_fit_engine(state$models[[1L]]))
      recipe <- function(state) workflows::extract_recipe(state$models[[1L]])$steps
      if (role != "B") {
        testthat::expect_identical(engine(other), engine(original))
        testthat::expect_identical(recipe(other), recipe(original))
      } else if (column == "x") testthat::expect_false(identical(recipe(other), recipe(original))) else
        testthat::expect_false(identical(engine(other), engine(original)))
      if (role == "C") testthat::expect_identical(other$postprocessing[[1L]]$coefficients, original$postprocessing[[1L]]$coefficients)
      if (role == "P") testthat::expect_false(identical(other$postprocessing[[1L]]$coefficients, original$postprocessing[[1L]]$coefficients))
    }
  }
})

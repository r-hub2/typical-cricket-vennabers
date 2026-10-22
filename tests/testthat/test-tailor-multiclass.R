testthat::test_that("shared Tailor holdouts preserve public workflows and multiclass VA", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "nnet")) testthat::skip_if_not_installed(pkg)
  calls <- new.env(); calls$workflow <- calls$tailor <- 0L
  original_workflow <- getS3method("predict", "workflow", envir = asNamespace("workflows"))
  original_tailor <- getS3method("predict", "tailor", envir = asNamespace("tailor"))
  count_workflow <- function(object, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) calls$workflow <- calls$workflow + 1L
    original_workflow(object, ...)
  }
  count_tailor <- function(object, ...) {
    calls$tailor <- calls$tailor + 1L
    original_tailor(object, ...)
  }
  registerS3method("predict", "workflow", count_workflow, envir = asNamespace("stats"))
  registerS3method("predict", "tailor", count_tailor, envir = asNamespace("stats"))
  on.exit({
    registerS3method("predict", "workflow", original_workflow, envir = asNamespace("stats"))
    registerS3method("predict", "tailor", original_tailor, envir = asNamespace("stats"))
  }, add = TRUE)
  for (k in 2:4) {
    classes <- c("z", "a", "b", "c")[seq_len(k)]
    d <- withr::with_seed(123, {
      truth <- factor(rep(classes, length.out = 240), levels = classes)
      data.frame(x = rnorm(240, sd = 2) + c(-1, 1, -1, 1)[truth],
        z = rnorm(240, sd = 2) + c(-1, -1, 1, 1)[truth],
        w = rnorm(240, sd = 2) + c(1, -1, -1, 1)[truth], truth = truth)
    })
    base <- workflows::workflow() |> workflows::add_formula(truth ~ x + z + w) |>
      workflows::add_model(if (k == 2L) parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
    method <- if (k == 2L) "logistic" else "multinomial"
    post <- list(method = paste0("tailor_", method), strategy = "holdout", prop = .3)
    probability <- function(model, data) {
      p <- as.matrix(predict(model, data, type = "prob")[paste0(".pred_", classes)])
      colnames(p) <- classes
      p
    }
    for (va in c("inductive", "cross")) for (classification in c("ovo", "ovr")) for (precision in list(NULL, 0L, 2L)) {
      fit <- generics::fit(calibrated_workflow(base, "truth", method = va, classification = classification,
        precision = precision, seed = 24, postprocessing = post), d)
      state <- workflows::extract_fit_engine(fit)
      held <- new <- truth <- ids <- list()
      testthat::expect_identical(state$class_levels, classes)
      testthat::expect_null(state$event_class)
      for (id in names(state$models)) {
        plan <- state$postprocessing[[id]]
        ref <- generics::fit(workflows::add_tailor(base, tailor::tailor() |>
          tailor::adjust_probability_calibration(method = method, smooth = FALSE)),
          d[plan$B, ], data_calibration = d[plan$P, ])
        held[[id]] <- probability(ref, d[plan$C, ])
        new[[id]] <- probability(ref, d)
        truth[[id]] <- d$truth[plan$C]
        ids[[id]] <- rep(id, length(plan$C))
        testthat::expect_lte(max(abs(held[[id]] - state$calibration[[id]]$probability)), 1e-10)
        testthat::expect_lte(max(abs(new[[id]] - probability(state$models[[id]], d))), 1e-10)
        testthat::expect_identical(plan$method, post$method)
        testthat::expect_identical(plan$rank, k)
        testthat::expect_setequal(c(plan$B, plan$P), plan$Q)
        testthat::expect_length(intersect(plan$B, plan$P), 0L)
        testthat::expect_length(intersect(plan$Q, plan$C), 0L)
      }
      core <- if (va == "cross") cross_venn_abers(do.call(c, truth), do.call(rbind, held), unlist(ids),
        method = classification, precision = precision) else
        venn_abers(truth[[1]], held[[1]], method = classification, precision = precision)
      for (loss in c("log", "brier")) {
        input <- if (va == "cross") new else new[[1]]
        expected <- predict(core, input, loss = loss)
        actual <- venn_abers_details(fit, d[c("x", "z", "w")], loss = loss)
        testthat::expect_lte(max(abs(as.matrix(actual$prediction) - expected)), 1e-10)
        testthat::expect_identical(actual$prediction, predict(fit, d[c("x", "z", "w")], type = "prob", loss = loss))
        testthat::expect_identical(predict(fit, d, type = "class", loss = loss)$.pred_class,
          predict(core, input, type = "class", loss = loss))
        for (n in 0:1) {
          calls$workflow <- calls$tailor <- 0L
          small <- venn_abers_details(fit, d[seq_len(n), c("x", "z", "w")], loss = loss)
          expected_calls <- if (n) length(state$models) else 0L
          testthat::expect_identical(calls$workflow, expected_calls)
          testthat::expect_identical(calls$tailor, expected_calls)
          testthat::expect_identical(nrow(small$prediction), n)
          if (n) testthat::expect_lte(max(abs(as.matrix(small$prediction) - expected[1L, , drop = FALSE])), 1e-10)
          testthat::expect_identical(names(small$details), names(actual$details))
        }
      }
    }
  }
})

testthat::test_that("multiclass Tailor rejects incompatible methods, plans and failed fits", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "nnet")) testthat::skip_if_not_installed(pkg)
  d <- withr::with_seed(123, {
    truth <- factor(rep(c("z", "a", "b"), 80), levels = c("z", "a", "b"))
    data.frame(x = rnorm(240, sd = 2) + c(-1, 1, -1, 1)[truth],
        z = rnorm(240, sd = 2) + c(-1, -1, 1, 1)[truth],
        w = rnorm(240, sd = 2) + c(1, -1, -1, 1)[truth], truth = truth)
  })
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x + z + w) |>
    workflows::add_model(parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
  va <- function(data) rsample::make_splits(list(analysis = 1:180, assessment = 181:240), data)
  post <- list(method = "tailor_multinomial", strategy = "holdout",
    splits = function(data) rsample::make_splits(list(analysis = 1:120, assessment = 121:180), data))
  spec <- calibrated_workflow(base, "truth", classification = "ovo", splits = va, postprocessing = post)
  fit <- generics::fit(spec, d)
  canonical <- tailor::tailor() |> tailor::adjust_probability_calibration(method = "multinomial", smooth = FALSE)
  adopted <- generics::fit(calibrated_workflow(workflows::add_tailor(base, canonical), "truth",
    classification = "ovo", splits = va, postprocessing = post), d)
  testthat::expect_identical(predict(fit, d, type = "prob"), predict(adopted, d, type = "prob"))
  testthat::expect_length(workflows::extract_tailor(workflows::extract_fit_engine(adopted)$models[[1]])$adjustments, 1L)
  for (bad in list(tailor::tailor() |> tailor::adjust_probability_calibration(method = "logistic", smooth = FALSE),
                  tailor::tailor() |> tailor::adjust_probability_calibration(method = "multinomial", smooth = TRUE),
                  canonical |> tailor::adjust_probability_threshold(.5),
                  workflows::extract_tailor(workflows::extract_fit_engine(fit)$models[[1]]))) {
    testthat::expect_error(calibrated_workflow(workflows::add_tailor(base, bad), "truth",
      classification = "ovr", postprocessing = post), "[Tt]ailor")
  }
  constant <- d; constant$x <- constant$z <- constant$w <- 1
  testthat::expect_error(suppressWarnings(generics::fit(spec, constant)), "tailor_multinomial.*full design rank 3")
  original <- probably::cal_estimate_multinomial
  for (fault in c("none", "missing", "nonconverged", "nonfinite", "warning")) local({
    replacement <- function(...) {
      if (fault == "none") return(probably::cal_estimate_none(...))
      result <- original(...)
      if (fault == "missing") result$estimates <- list()
      if (fault == "nonconverged") result$estimates[[1]]$estimate$convergence <- 1L
      if (fault == "nonfinite") result$estimates[[1]]$estimate$wts[] <- NA_real_
      if (fault == "warning") warning("instrumented successful warning")
      result
    }
    testthat::local_mocked_bindings(cal_estimate_multinomial = replacement, .package = "probably")
    if (fault == "warning") {
      testthat::expect_warning(valid <- generics::fit(spec, d), "tailor_multinomial.*instrumented successful warning")
      testthat::expect_match(workflows::extract_fit_engine(valid)$postprocessing[[1]]$warnings, "instrumented successful warning")
    } else testthat::expect_error(suppressWarnings(generics::fit(spec, d)), "Holdout1 tailor_multinomial")
  })
  original_predict <- getS3method("predict", "model_fit", envir = asNamespace("parsnip"))
  for (fault in c("missing", "duplicate", "nonfinite", "range", "sum", "tolerance", "reverse")) local({
    corrupt <- function(object, new_data, type = NULL, opts = list(), ...) {
      p <- original_predict(object, new_data, type = type, opts = opts, ...)
      if (identical(type, "prob") && !inherits(object$spec, "vennabers")) {
        if (fault == "missing") p <- p[-1L]
        if (fault == "duplicate") names(p)[2L] <- names(p)[1L]
        if (fault == "nonfinite") p[1L, 1L] <- Inf
        if (fault == "range") p[1L, 1L] <- -.1
        if (fault == "sum") p <- p * (1 + 2e-8)
        if (fault == "tolerance") p <- p * (1 + 5e-9)
        if (fault == "reverse") p <- p[rev(seq_along(p))]
      }
      p
    }
    testthat::local_mocked_bindings(predict.model_fit = corrupt, .package = "parsnip")
    registerS3method("predict", "model_fit", corrupt, envir = asNamespace("stats"))
    on.exit(registerS3method("predict", "model_fit", original_predict, envir = asNamespace("stats")), add = TRUE)
    if (fault %in% c("tolerance", "reverse")) {
      valid <- generics::fit(spec, d)
      testthat::expect_identical(names(predict(valid, d, type = "prob")), paste0(".pred_", levels(d$truth)))
    } else {
      testthat::expect_error(generics::fit(spec, d), "probabilit|column|Column")
      testthat::expect_error(predict(fit, d, type = "prob"), "probabilit|column|Column")
    }
  })
  # The external workflow boundary must not be reached by invalid local plans.
  original_fit <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  guard <- function(object, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) stop("UNEXPECTED LEARNING")
    original_fit(object, ...)
  }
  registerS3method("fit", "workflow", guard, envir = asNamespace("generics"))
  on.exit(registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics")), add = TRUE)
  for (classification in c("ovo", "ovr")) {
    for (k in 2:4) {
      bad_data <- d
      bad_data$truth <- factor(rep(letters[seq_len(k)], length.out = nrow(d)))
      bad_post <- list(method = if (k == 2L) "tailor_multinomial" else "tailor_logistic", strategy = "holdout", prop = .3)
      testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", classification = classification,
        splits = va, postprocessing = bad_post), bad_data), "require explicitly selecting tailor_")
    }
    for (P in list(integer(), 180L, 178:180, which(d$truth[1:180] != "b"), which(d$truth[1:180] == "b"))) {
      bad_post <- post
      bad_post$splits <- function(data) rsample::make_splits(list(analysis = setdiff(1:180, P), assessment = P), data)
      testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", classification = classification,
        splits = va, postprocessing = bad_post), d), "Holdout1 tailor_multinomial B/P plan")
    }
  }
})

testthat::test_that("multiclass P and C mutations preserve the appropriate learned state", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes", "nnet")) testthat::skip_if_not_installed(pkg)
  for (k in 2:4) {
    d <- withr::with_seed(123, {
      truth <- factor(rep(letters[seq_len(k)], length.out = 240))
      data.frame(x = rnorm(240, sd = 2) + c(-1, 1, -1, 1)[truth],
        z = rnorm(240, sd = 2) + c(-1, -1, 1, 1)[truth],
        w = rnorm(240, sd = 2) + c(1, -1, -1, 1)[truth], truth = truth)
    })
    base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ ., d) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(if (k == 2L) parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
    post <- list(method = if (k == 2L) "tailor_logistic" else "tailor_multinomial", strategy = "holdout", prop = .4)
    for (method in c("inductive", "cross")) for (classification in c("ovo", "ovr")) {
      spec <- calibrated_workflow(base, "truth", method = method, classification = classification, seed = 24, postprocessing = post)
      fit <- generics::fit(spec, d)
      state <- workflows::extract_fit_engine(fit)
      # Replay the stored plans so mutations test data roles, not random splits.
      va_provider <- function(data) {
        parts <- lapply(state$folds, function(f) rsample::make_splits(f[c("analysis", "assessment")], data))
        if (method == "cross") rsample::manual_rset(parts, names(parts)) else parts[[1L]]
      }
      post$splits <- function(data) {
        id <- which(vapply(state$postprocessing, function(p) identical(rownames(data), rownames(d)[p$Q]), logical(1)))
        rsample::make_splits(state$postprocessing[[id]][c("analysis", "assessment")], data)
      }
      post$prop <- NULL
      replay <- calibrated_workflow(base, "truth", method = method, classification = classification,
        splits = va_provider, seed = 24, postprocessing = post)
      p <- state$postprocessing[[1L]]
      for (role in c("P", "C")) for (column in c("x", "truth")) {
        changed <- d
        if (column == "x") changed$x[p[[role]]] <- changed$x[p[[role]]] + .1 else {
          rows <- p[[role]][match(levels(d$truth)[1:2], as.character(d$truth[p[[role]]]))]
          changed$truth[rows] <- rev(changed$truth[rows])
        }
        other <- workflows::extract_fit_engine(generics::fit(replay, changed))
        testthat::expect_identical(workflows::extract_recipe(state$models[[1]])$steps,
          workflows::extract_recipe(other$models[[1]])$steps)
        testthat::expect_identical(stats::coef(workflows::extract_fit_engine(state$models[[1]])),
          stats::coef(workflows::extract_fit_engine(other$models[[1]])))
        if (role == "P") testthat::expect_false(identical(p$coefficients, other$postprocessing[[1]]$coefficients)) else
          testthat::expect_identical(p$coefficients, other$postprocessing[[1]]$coefficients)
      }
      post$splits <- NULL; post$prop <- .4
    }
  }
})

testthat::test_that("missing Cross C classes follow VA rules independently of B and P", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "nnet")) testthat::skip_if_not_installed(pkg)
  d <- withr::with_seed(123, {
    truth <- factor(rep(c("z", "a", "b"), 80), levels = c("z", "a", "b"))
    data.frame(x = rnorm(240, sd = 2) + c(-1, 1, -1)[truth],
      z = rnorm(240, sd = 2) + c(-1, -1, 1)[truth],
      w = rnorm(240, sd = 2) + c(1, -1, -1)[truth], truth = truth)
  })
  base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |>
    workflows::add_model(parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
  for (classification in c("ovo", "ovr")) {
    first <- which(d$truth %in% if (classification == "ovo") c("z", "a") else "z")[1:12]
    remaining <- setdiff(1:240, first)
    held <- c(list(first), split(remaining, cut(seq_along(remaining), 2L, labels = FALSE)))
    plans <- rsample::manual_rset(lapply(held, function(h) rsample::make_splits(
      list(analysis = setdiff(1:240, h), assessment = h), d)), c("First", "Second", "Third"))
    fit <- generics::fit(calibrated_workflow(base, "truth", method = "cross", classification = classification,
      splits = plans, seed = 24, postprocessing = list(method = "tailor_multinomial", strategy = "holdout", prop = .4)), d)
    state <- workflows::extract_fit_engine(fit)
    testthat::expect_lt(length(unique(state$calibration$First$truth)), 3L)
    for (p in state$postprocessing) {
      testthat::expect_setequal(as.character(d$truth[p$B]), levels(d$truth))
      testthat::expect_setequal(as.character(d$truth[p$P]), levels(d$truth))
    }
    testthat::expect_true(all(is.finite(as.matrix(predict(fit, d, type = "prob")))))
  }
})

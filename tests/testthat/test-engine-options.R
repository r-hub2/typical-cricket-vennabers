# Fixed engine options from parsnip::set_engine() reach every local base fit,
# which learns only from its own rows. No VA-level fit-data interface exists.

options_data <- function(n = 240L, classes = c("a", "b")) {
  truth <- factor(rep(classes, length.out = n), levels = classes)
  data.frame(x = sin(1:n) + as.integer(truth) / 2, z = cos(1:n * 1.7), truth = truth)
}

options_base <- function() {
  workflows::workflow() |>
    workflows::add_formula(truth ~ .) |>
    workflows::add_model(parsnip::decision_tree(mode = "classification") |>
      parsnip::set_engine("rpart", xval = 0L, maxsurrogate = 0L))
}

# Fits with rpart recorded as original rows (`x` is unique) and engine options.
recorded_fit <- function(d, ...) {
  log <- new.env()
  log$fits <- list()
  original <- rpart::rpart
  testthat::local_mocked_bindings(rpart = function(formula, data, ...) {
    log$fits[[length(log$fits) + 1L]] <- list(rows = sort(match(data$x, d$x)), options = list(...))
    original(formula, data, ...)
  }, .package = "rpart")
  fit <- generics::fit(calibrated_workflow(options_base(), "truth", seed = 42, ...), d)
  list(state = workflows::extract_fit_engine(fit), fits = log$fits)
}

expect_fixed_options <- function(recorded, rows, info) {
  testthat::expect_identical(lapply(recorded$fits, `[[`, "rows"), unname(lapply(rows, sort)), info = info)
  for (f in recorded$fits) {
    testthat::expect_identical(f$options, list(xval = 0L, maxsurrogate = 0L), info = info)
  }
}

testthat::test_that("set_engine() options reach every source, pair, full and B fit on its local rows", {
  skip_if_no_integration("rpart")
  d <- options_data()
  for (method in c("inductive", "cross")) for (mode in c("matched", "full")) {
    info <- paste(method, mode)
    recorded <- do.call(recorded_fit, c(list(d, method = method, prediction_model = mode),
                                        if (method == "cross") list(v = 2L)))
    sources <- lapply(recorded$state$folds, `[[`, "analysis")
    expect_fixed_options(recorded, c(sources, if (mode == "full") list(seq_len(nrow(d)))), info)
    # Stored models carry the option.
    for (model in recorded$state$models) {
      testthat::expect_identical(workflows::extract_fit_engine(model)$control$xval, 0L, info = info)
    }
  }

  # Every pair's sources first, then each pair's Full1 on all its rows.
  multi <- options_data(classes = c("a", "b", "c"))
  recorded <- recorded_fit(multi, classification = "ovo", model_training = "pairwise",
                           method = "cross", v = 2L, prediction_model = "full")
  pairs <- unname(recorded$state$pair_models)
  sources <- unlist(lapply(pairs, function(pair) lapply(pair$folds, function(f) pair$rows[f$analysis])), recursive = FALSE)
  expect_fixed_options(recorded, c(sources, lapply(pairs, `[[`, "rows")), "pairwise")
  testthat::expect_length(recorded$fits, 9L)

  # Postprocessing fits the engine on each B only; P and C never reach it.
  for (strategy in c("holdout", "ensemble")) {
    post <- c(list(method = "sigmoid", sigmoid_input = "probability", strategy = strategy),
              if (strategy == "holdout") list(prop = .3) else list(v = 2L))
    recorded <- recorded_fit(d, method = "cross", v = 2L, prediction_model = "full", postprocessing = post)
    plans <- unlist(lapply(recorded$state$postprocessing, function(p) if (is.null(p$inner_folds)) list(p) else p$inner_folds),
                    recursive = FALSE)
    expect_fixed_options(recorded, unname(lapply(plans, `[[`, "B")), strategy)
    testthat::expect_length(recorded$fits, if (strategy == "holdout") 3L else 6L)
  }
})

testthat::test_that("fit options stay in set_engine(): outer dots are rejected", {
  skip_if_no_integration("rpart")
  d <- options_data()
  testthat::expect_error(calibrated_workflow(options_base(), "truth", xval = 0L),
                         "`...` must be empty; unsupported argument: xval.", fixed = TRUE)
  testthat::expect_error(generics::fit(calibrated_workflow(options_base(), "truth"), d, xval = 0L),
                         "must be empty")
})

testthat::test_that("native XGBoost early stopping translates to local internal validation", {
  skip_if_no_integration()
  # Translation needs no installed xgboost; the real fit is checked by
  # tools/xgboost-early-stopping.R.
  spec <- parsnip::boost_tree(mode = "classification", trees = 100L, stop_iter = 10L) |>
    parsnip::set_engine("xgboost", validation = 0.2)
  fit <- parsnip::translate(spec)$method$fit
  testthat::expect_identical(fit$func, c(pkg = "parsnip", fun = "xgb_train"))
  testthat::expect_identical(lapply(fit$args[c("nrounds", "early_stop", "validation")], rlang::eval_tidy),
                             list(nrounds = 100L, early_stop = 10L, validation = 0.2))
})

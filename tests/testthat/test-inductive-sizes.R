testthat::test_that("independent Inductive sizes resolve sklearn rounding on the local domain", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:100), truth = factor(rep(c("yes", "no", "no", "yes"), 25)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  sizes <- function(...) {
    state <- workflows::extract_fit_engine(generics::fit(calibrated_workflow(base, "truth", seed = 3, ...), d))
    fold <- state$folds$Inductive1
    c(length(fold$analysis), length(fold$assessment), length(fold$unassigned))
  }
  # Training floors, calibration ceils, whole numbers are exact, a missing size is the complement.
  testthat::expect_identical(sizes(training_size = .555, calibration_size = .071), c(55L, 8L, 37L))
  testthat::expect_identical(sizes(training_size = .7), c(70L, 30L, 0L))
  testthat::expect_identical(sizes(calibration_size = .155), c(84L, 16L, 0L))
  testthat::expect_identical(sizes(calibration_size = 15), c(85L, 15L, 0L))
  testthat::expect_identical(sizes(training_size = .5, calibration_size = 10), c(50L, 10L, 40L))
  testthat::expect_identical(sizes(training_size = 10, calibration_size = 1, strata = FALSE), c(10L, 1L, 89L))

  # Exact totals under stratification, where rsample rounds each class down.
  three <- data.frame(x = sin(1:21), truth = factor(rep(c("a", "b", "c"), 7)))
  multi <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  rsample_rows <- length(as.integer(withr::with_seed(1, rsample::initial_split(three, prop = .8, strata = truth, pool = 0)), data = "analysis"))
  testthat::expect_identical(rsample_rows, 15L)
  for (strata in list(NULL, FALSE)) {
    state <- workflows::extract_fit_engine(generics::fit(calibrated_workflow(multi, "truth", classification = "ovr",
      seed = 1, training_size = .8, strata = strata), three))
    fold <- state$folds$Inductive1
    testthat::expect_identical(c(length(fold$analysis), length(fold$assessment)), c(16L, 5L))
    if (is.null(strata)) {
      # Largest remainders, ties in level order: 16 = 6 + 5 + 5, calibration takes the rest.
      testthat::expect_identical(as.vector(table(three$truth[fold$analysis])), c(6L, 5L, 5L))
    }
  }
  state <- workflows::extract_fit_engine(generics::fit(calibrated_workflow(multi, "truth", classification = "ovr",
    seed = 1, training_size = 9, calibration_size = 7), three))
  fold <- state$folds$Inductive1
  testthat::expect_identical(as.vector(table(three$truth[fold$analysis])), c(3L, 3L, 3L))
  testthat::expect_identical(as.vector(table(three$truth[fold$assessment])), c(3L, 2L, 2L))
  testthat::expect_identical(sort(unlist(fold[c("analysis", "assessment", "unassigned")], use.names = FALSE)), 1:21)

  # The same seed gives the same draw; a different seed a different one.
  plan <- function(seed) workflows::extract_fit_engine(generics::fit(calibrated_workflow(base, "truth", seed = seed,
    training_size = 50, calibration_size = 20), d))$folds$Inductive1
  testthat::expect_identical(plan(9), plan(9))
  testthat::expect_false(identical(plan(9)$analysis, plan(10)$analysis))
})

testthat::test_that("the 100 -> 50/20/30 scenario keeps U out of matched and inside Full1", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:100) + rep(c(1, 0, 0, 1), 25), truth = factor(rep(c("yes", "no", "no", "yes"), 25)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (mode in c("matched", "full")) {
    fit <- generics::fit(calibrated_workflow(base, "truth", seed = 42, training_size = 50, calibration_size = 20,
      prediction_model = mode), d)
    state <- workflows::extract_fit_engine(fit)
    fold <- state$folds$Inductive1
    testthat::expect_identical(lengths(fold[c("analysis", "assessment", "unassigned", "training_size", "calibration_size")]),
      c(analysis = 50L, assessment = 20L, unassigned = 30L, training_size = 1L, calibration_size = 1L))
    testthat::expect_identical(c(fold$training_size, fold$calibration_size), c(50L, 20L))
    testthat::expect_identical(sort(c(fold$analysis, fold$assessment, fold$unassigned)), 1:100)
    # Stratified: 25 of each class in training, 10 in calibration.
    testthat::expect_identical(as.vector(table(d$truth[fold$analysis])), c(25L, 25L))
    testthat::expect_identical(as.vector(table(d$truth[fold$assessment])), c(10L, 10L))
    event <- state$event_class
    source <- generics::fit(base, d[fold$analysis, ])
    core <- venn_abers(d$truth[fold$assessment] == event,
      predict(source, d[fold$assessment, ], type = "prob")[[paste0(".pred_", event)]])
    testthat::expect_identical(state$calibration$Inductive1$row, fold$assessment)
    predictor <- if (mode == "full") generics::fit(base, d) else source
    expected <- predict(core, predict(predictor, d, type = "prob")[[paste0(".pred_", event)]])$probability
    testthat::expect_lte(max(abs(predict(fit, d, type = "prob")[[paste0(".pred_", event)]] - expected)), 1e-10)
    testthat::expect_identical(unname(state$sources), if (mode == "full") "Full1" else "Inductive1")
    testthat::expect_identical(state$parameters$requested[c("training_size", "calibration_size")],
      list(training_size = 50, calibration_size = 20))
    testthat::expect_identical(state$parameters$effective[c("training_size", "calibration_size", "unassigned_size")],
      list(training_size = 50L, calibration_size = 20L, unassigned_size = 30L))
    testthat::expect_null(state$parameters$effective$calibration_prop)
    printed <- utils::capture.output(print(state))
    testthat::expect_true("Inductive rows: training 50 | calibration 20 | unassigned 30 " %in% printed)
  }
})

testthat::test_that("regression and postprocessing learn only on training; Full1 uses all of F", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:100), y = cos(1:100) + sin(1:100))
  base <- workflows::workflow() |> workflows::add_formula(y ~ x) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  fit <- generics::fit(calibrated_workflow(base, "y", seed = 5, training_size = .5, calibration_size = .2), d)
  state <- workflows::extract_fit_engine(fit)
  fold <- state$folds$Inductive1
  testthat::expect_identical(lengths(fold[c("analysis", "assessment", "unassigned")]),
    c(analysis = 50L, assessment = 20L, unassigned = 30L))
  source <- generics::fit(base, d[fold$analysis, ])
  core <- venn_abers_regression(d$y[fold$assessment], predict(source, d[fold$assessment, ])$.pred)
  expected <- predict(core, predict(source, d)$.pred)$prediction
  testthat::expect_lte(max(abs(predict(fit, d)$.pred - expected)), 1e-10)

  b <- data.frame(x = sin(1:100) + rep(c(1, 0, 0, 1), 25), truth = factor(rep(c("yes", "no", "no", "yes"), 25)))
  logistic <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  post <- list(method = "sigmoid", strategy = "holdout", prop = .25, sigmoid_input = "probability")
  for (mode in c("matched", "full")) {
    state <- workflows::extract_fit_engine(generics::fit(calibrated_workflow(logistic, "truth", seed = 6,
      training_size = 50, calibration_size = 20, prediction_model = mode, postprocessing = post), b))
    fold <- state$folds$Inductive1
    held <- state$postprocessing$Inductive1
    # B/P is drawn inside training only; neither C nor U takes part.
    testthat::expect_identical(held$Q, fold$analysis)
    testthat::expect_identical(sort(c(held$B, held$P)), fold$analysis)
    if (mode == "full") {
      testthat::expect_identical(state$postprocessing$Full1$Q, 1:100)
      testthat::expect_identical(state$postprocessing$Full1$B_size + state$postprocessing$Full1$P_size, 100L)
    }
  }
})

testthat::test_that("default calls, saved fits and older fits without U keep their contract", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:100), truth = factor(rep(c("yes", "no", "no", "yes"), 25)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  # Both sizes NULL: the calibration_prop draw and the stored shape are unchanged.
  for (prop in c(.2, .35)) {
    spec <- if (prop == .2) calibrated_workflow(base, "truth", seed = 42) else
      calibrated_workflow(base, "truth", seed = 42, calibration_prop = prop)
    state <- workflows::extract_fit_engine(generics::fit(spec, d))
    split <- withr::with_seed(42, rsample::initial_split(d, prop = 1 - prop, strata = truth, pool = 0))
    testthat::expect_identical(state$folds$Inductive1[c("analysis", "assessment", "training_size", "calibration_size")],
      list(analysis = as.integer(split, data = "analysis"), assessment = as.integer(split, data = "assessment"),
           training_size = length(as.integer(split, data = "analysis")),
           calibration_size = length(as.integer(split, data = "assessment"))))
    testthat::expect_false(any(c("training_size", "calibration_size") %in% names(state$parameters$requested)))
    testthat::expect_false(any(c("training_size", "calibration_size", "unassigned_size") %in% names(state$parameters$effective)))
    testthat::expect_identical(state$parameters$effective$calibration_prop, prop)
    testthat::expect_false(any(grepl("Inductive rows", utils::capture.output(print(state)))))
  }

  fit <- generics::fit(calibrated_workflow(base, "truth", seed = 42, training_size = 50, calibration_size = 20,
    prediction_model = "full"), d)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  testthat::expect_identical(predict(readRDS(path), d, type = "prob"), predict(fit, d, type = "prob"))
  # A fit saved before U existed has no field: U is empty and nothing else changes.
  older <- fit
  older$fit$fit$fit$folds$Inductive1$unassigned <- NULL
  older$fit$fit$fit$parameters$effective[c("training_size", "calibration_size", "unassigned_size")] <- NULL
  testthat::expect_identical(predict(older, d, type = "prob"), predict(fit, d, type = "prob"))
  testthat::expect_identical(venn_abers_details(older, d), venn_abers_details(fit, d))
  testthat::expect_false(any(grepl("Inductive rows", utils::capture.output(print(workflows::extract_fit_engine(older))))))
})

testthat::test_that("invalid sizes and conflicts fail with the domain and both sizes", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:100), truth = factor(rep(c("yes", "no", "no", "yes"), 25)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (size in list(0, -1, 1.5, 1 + 1e-9, Inf, NA_real_, TRUE, "10", numeric(), c(.1, .2), matrix(.5))) {
    testthat::expect_error(calibrated_workflow(base, "truth", training_size = size), "`training_size` must be NULL")
    testthat::expect_error(calibrated_workflow(base, "truth", calibration_size = size), "`calibration_size` must be NULL")
  }
  testthat::expect_error(calibrated_workflow(base, "truth", training_size = quote(tune())), "cannot be tuned")
  testthat::expect_error(calibrated_workflow(base, "truth", calibration_size = quote(tune("c"))), "cannot be tuned")
  testthat::expect_error(calibrated_workflow(base, "truth", training_size = .7, calibration_size = .4), "sum to more than 1")
  testthat::expect_error(calibrated_workflow(base, "truth", method = "cross", training_size = .5), "method = 'inductive'")
  testthat::expect_error(calibrated_workflow(base, "truth", calibration_prop = .2, calibration_size = .2), "not both")
  testthat::expect_error(calibrated_workflow(base, "truth", training_size = .5,
    splits = function(data) rsample::initial_split(data)), "partition all rows")
  testthat::expect_error(calibrated_workflow(base, "truth", training_size = .5, training = 1), "training")

  fails <- function(text, ...) testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", seed = 1, ...), d), text, fixed = TRUE)
  fails("`training_size` = 150 (150 rows) request 150 rows, but only the 100 local rows are available.", training_size = 150)
  fails("`training_size` = 50 (50 rows) and `calibration_size` = 60 (60 rows) request 110 rows, but only the 100 local rows are available.",
    training_size = 50, calibration_size = 60)
  fails("`training_size` = 0.5 (50 rows) and `calibration_size` = 51 (51 rows) request 101 rows", training_size = .5, calibration_size = 51)
  fails("`training_size` = 0.001 (0 rows) leaves 0 training rows of the 100 local rows", training_size = .001)
  fails("`training_size` = 100 (100 rows) leaves 0 calibration rows of the 100 local rows", training_size = 100)
  fails("`calibration_size` = 100 (100 rows) leaves 0 training rows", calibration_size = 100)

  # After the size check, the existing checks run unchanged: every class must reach C.
  three <- data.frame(x = sin(1:30), truth = factor(rep(c("a", "b", "c"), 10)))
  multi <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::multinom_reg(engine = "nnet"))
  testthat::expect_error(generics::fit(calibrated_workflow(multi, "truth", classification = "ovr", seed = 1,
    training_size = 20, calibration_size = 2), three), "class")
})

testthat::test_that("pairwise sizes resolve on each local pair with U per pair and Full1 on all pair rows", {
  skip_if_no_integration()
  # Unequal classes give pairs of 70, 60 and 50 rows.
  d <- withr::with_seed(11, data.frame(x = stats::rnorm(90) + rep(c(0, .8, 1.6), c(40, 30, 20)),
    truth = factor(rep(c("a", "b", "c"), c(40, 30, 20)))))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  pairwise <- function(...) calibrated_workflow(base, "truth", classification = "ovo",
    model_training = "pairwise", seed = 7, ...)
  pairs <- list(c("a", "b"), c("a", "c"), c("b", "c"))
  for (mode in c("matched", "full")) {
    state <- workflows::extract_fit_engine(generics::fit(pairwise(training_size = .5, calibration_size = 10,
      prediction_model = mode), d))
    expected <- list(training_size = c(pair_1 = 35L, pair_2 = 30L, pair_3 = 25L),
      calibration_size = c(pair_1 = 10L, pair_2 = 10L, pair_3 = 10L),
      unassigned_size = c(pair_1 = 25L, pair_2 = 20L, pair_3 = 15L))
    testthat::expect_identical(state$parameters$effective[names(expected)], expected)
    testthat::expect_null(state$folds)
    for (i in seq_along(pairs)) {
      component <- state$pair_models[[i]]
      fold <- component$folds$Inductive1
      local <- d$truth[component$rows]
      testthat::expect_identical(component$rows, which(d$truth %in% pairs[[i]]))
      testthat::expect_identical(lengths(fold[c("analysis", "assessment", "unassigned")]),
        c(analysis = expected$training_size[[i]], assessment = 10L, unassigned = expected$unassigned_size[[i]]))
      # Positions are local to the pair, which U completes.
      testthat::expect_identical(sort(c(fold$analysis, fold$assessment, fold$unassigned)), seq_along(local))
      testthat::expect_identical(component$calibration$Inductive1$row, fold$assessment)
      # Exact stratified totals per pair, as in the shared domain.
      testthat::expect_identical(sum(table(droplevels(local[fold$analysis]))), expected$training_size[[i]])
      testthat::expect_true(all(table(droplevels(local[fold$assessment])) > 0L))
      learned <- function(id) length(workflows::extract_fit_engine(component$models[[id]])$y)
      if (mode == "full") {
        # Full1 learns on all pair rows, including the pair's C and U.
        testthat::expect_identical(names(component$models), "Full1")
        testthat::expect_identical(learned("Full1"), length(local))
      } else testthat::expect_identical(learned("Inductive1"), expected$training_size[[i]])
    }
    printed <- utils::capture.output(print(state))
    testthat::expect_true(all(c("pair_1 Inductive rows: training 35 | calibration 10 | unassigned 25 ",
      "pair_3 Inductive rows: training 25 | calibration 10 | unassigned 15 ") %in% printed))
  }

  # Absolute sizes are strict per pair: only the 50-row pair [b, c] is too small.
  testthat::expect_error(generics::fit(pairwise(training_size = 30, calibration_size = 25), d),
    "Pair [b, c]: `training_size` = 30 (30 rows) and `calibration_size` = 25 (25 rows) request 55 rows, but only the 50 local pair rows are available.",
    fixed = TRUE)
  # Explicit pair plans, providers and prepared plans still partition all rows.
  pair_plans <- lapply(stats::setNames(pairs, paste0("pair_", 1:3)), function(pair) {
    local <- d[d$truth %in% pair, ]
    local$truth <- factor(local$truth, levels = pair)
    rsample::initial_split(local)
  })
  testthat::expect_error(pairwise(training_size = .5, splits = pair_plans), "partition all rows")
  testthat::expect_error(pairwise(calibration_size = 10, splits = function(local) rsample::initial_split(local)),
    "partition all rows")
  prepared <- local_calibration_splits(list(rsample::initial_split(d)), names(d))
  testthat::expect_error(calibrated_workflow(base, "truth", classification = "ovr", training_size = .5,
    splits = prepared), "partition all rows")

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  fit <- generics::fit(pairwise(training_size = 20, calibration_size = 10, prediction_model = "full"), d)
  saveRDS(fit, path)
  testthat::expect_identical(predict(readRDS(path), d, type = "prob"), predict(fit, d, type = "prob"))
  testthat::expect_identical(venn_abers_details(readRDS(path), d), venn_abers_details(fit, d))
})

# Records every automatic size resolution: the domain's row count and the
# resolved training, calibration and unassigned counts.
record_size_resolutions <- function(env = parent.frame()) {
  resolve <- vennabers:::.inductive_size_positions
  log <- new.env()
  log$calls <- list()
  testthat::local_mocked_bindings(.inductive_size_positions = function(domain, settings, local) {
    plan <- resolve(domain, settings, local)
    log$calls[[length(log$calls) + 1L]] <- c(n = nrow(local), lengths(plan))
    plan
  }, .package = "vennabers", .env = env)
  log
}

# sklearn rounding of the requested proportions on each recorded domain.
expect_size_rule <- function(calls, training, calibration) {
  for (call in calls) {
    n <- call[["n"]]
    testthat::expect_identical(unname(call[c("analysis", "assessment", "unassigned")]),
      as.integer(c(floor(n * training), ceiling(n * calibration), n - floor(n * training) - ceiling(n * calibration))))
  }
}

testthat::test_that("Grid searches with Venn-Abers candidates resolve proportions on their own rows", {
  skip_if_no_integration("tune", "yardstick", "dials")
  d <- withr::with_seed(12, data.frame(x = stats::rnorm(120), truth = factor(rep(c("yes", "no"), 60))))
  d$x <- d$x + (d$truth == "yes")
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x, d) |> recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  folds <- function(data) rsample::vfold_cv(data, v = 2, strata = truth, pool = 0)
  model <- list(method = "grid", grid = data.frame(degree = 1:2), resamples = folds)

  # Joint model search inside the 60 training rows: candidates resolve on 30-row resample analyses.
  log <- record_size_resolutions()
  fit <- generics::fit(calibrated_workflow(base, "truth", seed = 4, training_size = .5, calibration_size = .25,
    tuning = list(strategy = "joint", model = model)), d)
  ns <- vapply(log$calls, `[[`, 0, "n")
  testthat::expect_identical(ns, c(120, rep(30, 4)))
  expect_size_rule(log$calls, .5, .25)
  state <- workflows::extract_fit_engine(fit)
  testthat::expect_identical(lengths(state$folds$Inductive1[c("analysis", "assessment", "unassigned")]),
    c(analysis = 60L, assessment = 30L, unassigned = 30L))
  testthat::expect_identical(venn_abers_tuning(fit)$fit_roles, state$folds)

  # Pairwise joint grid: the cached pair candidates resolve on each resample's pair rows.
  three <- withr::with_seed(13, data.frame(x = stats::rnorm(120), truth = factor(rep(c("a", "b", "c"), 40))))
  three$x <- three$x + as.integer(three$truth)
  pair_base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x, three) |> recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  log <- record_size_resolutions()
  fit <- generics::fit(calibrated_workflow(pair_base, "truth", classification = "ovo", model_training = "pairwise",
    seed = 5, training_size = .5, calibration_size = .25, prediction_model = "full",
    tuning = list(strategy = "joint", model = list(method = "grid", grid = data.frame(degree = 1:2)),
      joint_resamples = folds)), three)
  ns <- vapply(log$calls, `[[`, 0, "n")
  # Each 60-row resample analysis has 40-row pairs; the final fit has 80-row pairs.
  testthat::expect_identical(ns, c(rep(40, 6), rep(80, 3)))
  expect_size_rule(log$calls, .5, .25)
  state <- workflows::extract_fit_engine(fit)
  testthat::expect_identical(state$parameters$effective$unassigned_size, c(pair_1 = 20L, pair_2 = 20L, pair_3 = 20L))
  testthat::expect_identical(names(venn_abers_tuning(fit)$joint$selected$pairs), paste0("pair_", 1:3))

  # Row counts cannot follow candidates onto resamples; two-stage candidates have no Venn-Abers step.
  calibration <- list(method = "grid", grid = data.frame(va_loss = c("log", "brier")), resamples = folds)
  rejected <- "supply `training_size` and `calibration_size` as proportions, not row counts"
  testthat::expect_error(calibrated_workflow(base, "truth", training_size = 60,
    tuning = list(strategy = "joint", model = model)), rejected, fixed = TRUE)
  testthat::expect_error(calibrated_workflow(base, "truth", training_size = .5, calibration_size = 20,
    tuning = list(strategy = "two_stage", model = model, calibration = calibration, selection_prop = .5),
    loss = tune::tune("va_loss")), rejected, fixed = TRUE)
  testthat::expect_error(calibrated_workflow(pair_base, "truth", classification = "ovo", model_training = "pairwise",
    calibration_size = 10, tuning = list(strategy = "joint", model = list(method = "grid", grid = data.frame(degree = 1:2)),
      joint_resamples = folds)), rejected, fixed = TRUE)
  log <- record_size_resolutions()
  fit <- generics::fit(calibrated_workflow(base, "truth", seed = 4, training_size = 50, calibration_size = 20,
    tuning = list(strategy = "two_stage", model = model)), d)
  testthat::expect_identical(log$calls, list(c(n = 120L, analysis = 50L, assessment = 20L, unassigned = 50L)))
  testthat::expect_identical(lengths(venn_abers_tuning(fit)$fit_roles$Inductive1[c("analysis", "assessment", "unassigned")]),
    c(analysis = 50L, assessment = 20L, unassigned = 50L))

  # Pairwise two-stage searches keep absolute row counts per pair.
  fit <- generics::fit(calibrated_workflow(pair_base, "truth", classification = "ovo", model_training = "pairwise",
    seed = 5, training_size = 40, calibration_size = 20, tuning = list(strategy = "two_stage",
      model = list(method = "grid", grid = data.frame(degree = 1:2), resamples = folds))), three)
  testthat::expect_identical(workflows::extract_fit_engine(fit)$parameters$effective$unassigned_size,
    c(pair_1 = 20L, pair_2 = 20L, pair_3 = 20L))
})

testthat::test_that("Bayes calibration selection resolves proportions on S candidates and on F", {
  skip_if_no_integration("tune", "yardstick", "dials")
  d <- withr::with_seed(14, data.frame(x = stats::rnorm(120), truth = factor(rep(c("yes", "no"), 60))))
  d$x <- d$x + (d$truth == "yes")
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  folds <- function(data) rsample::vfold_cv(data, v = 2, strata = truth, pool = 0)
  calibration <- list(method = "bayes", initial = 2L, iter = 1L, resamples = folds,
    param_info = dials::parameters(dials::new_qual_param(type = "character", values = c("log", "brier"),
      label = c(va_loss = "Loss"))),
    control = tune::control_bayes(seed = 5, save_pred = TRUE, no_improve = 2L))
  log <- record_size_resolutions()
  fit <- generics::fit(calibrated_workflow(base, "truth", seed = 8, loss = tune::tune("va_loss"),
    training_size = .6, calibration_size = .2, prediction_model = "full",
    tuning = list(strategy = "joint", selection_prop = .5, calibration = calibration)), d)
  ns <- vapply(log$calls, `[[`, 0, "n")
  # S candidates resolve on 30-row resample analyses, the final fit on the 60 rows of F.
  testthat::expect_true(length(ns) > 2L)
  testthat::expect_true(all(utils::head(ns, -1L) == 30))
  testthat::expect_identical(utils::tail(ns, 1L), 60)
  expect_size_rule(log$calls, .6, .2)
  selection <- venn_abers_tuning(fit)
  state <- workflows::extract_fit_engine(fit)
  fold <- state$folds$Inductive1
  testthat::expect_length(selection$calibration$final_rows, 60L)
  testthat::expect_identical(lengths(fold[c("analysis", "assessment", "unassigned")]),
    c(analysis = 36L, assessment = 12L, unassigned = 12L))
  testthat::expect_identical(sort(c(fold$analysis, fold$assessment, fold$unassigned)), 1:60)
  testthat::expect_identical(state$parameters$effective[c("training_size", "calibration_size", "unassigned_size")],
    list(training_size = 36L, calibration_size = 12L, unassigned_size = 12L))
  testthat::expect_identical(names(state$models), "Full1")

  # A proportion that leaves F empty names F; the S candidates had rows to spare.
  ten <- function(data) rsample::make_splits(list(analysis = 1:10, assessment = 11:nrow(data)), data)
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", seed = 8, loss = tune::tune("va_loss"),
    training_size = .08, tuning = list(strategy = "joint", selection_split = ten,
      calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")), resamples = folds))), d),
    "`training_size` = 0.08 (0 rows) leaves 0 training rows of the 10 final-fit F rows", fixed = TRUE)
})

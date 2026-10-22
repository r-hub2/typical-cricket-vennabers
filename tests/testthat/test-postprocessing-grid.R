# Records every non-proxy base workflow fit with its B rows, P rows and
# whether a Tailor postprocessor was attached.
record_postprocessing_fits <- function(env = parent.frame()) {
  audit <- new.env()
  audit$fits <- list()
  # Observe B/P at the owned composition boundary, then count the public B fit.
  fit_tailor <- .fit_tailor_workflow
  testthat::local_mocked_bindings(.fit_tailor_workflow = function(workflow, data, data_calibration, ...) {
    audit$calibration <- data_calibration$row_id
    on.exit(audit$calibration <- NULL)
    fit_tailor(workflow, data, data_calibration, ...)
  }, .package = "vennabers", .env = env)
  fit_method <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  registerS3method("fit", "workflow", function(object, data, ..., data_calibration = NULL) {
    if (!inherits(workflows::extract_spec_parsnip(object), c("vennabers", "vennabers_grid"))) {
      post <- tryCatch(workflows::extract_postprocessor(object), error = function(e) NULL)
      audit$fits[[length(audit$fits) + 1L]] <- list(rows = data$row_id,
        calibration = if (!is.null(audit$calibration)) audit$calibration else data_calibration$row_id,
        tailor = !is.null(audit$calibration) || !is.null(post))
    }
    fit_method(object, data, ..., data_calibration = data_calibration)
  }, envir = asNamespace("generics"))
  withr::defer(registerS3method("fit", "workflow", fit_method, envir = asNamespace("generics")), envir = env)
  audit
}

grid_binary_data <- function(n = 240L) {
  withr::with_seed(11, {
    x <- rnorm(n, sd = 1.5)
    truth <- factor(ifelse(stats::runif(n) < stats::plogis(x - .4 * x^2 + .5), "yes", "no"), levels = c("no", "yes"))
    data.frame(row_id = seq_len(n), x = x, truth = truth)
  })
}

grid_tuned_base <- function(data, engine = parsnip::logistic_reg(engine = "glm")) {
  recipe <- recipes::recipe(truth ~ ., data = data) |>
    recipes::update_role(row_id, new_role = "id") |>
    recipes::step_poly(x, degree = tune::tune("degree"))
  workflows::workflow() |> workflows::add_recipe(recipe) |> workflows::add_model(engine)
}

# Positional plans keep every role fixed when target values are mutated.
positional_va <- function(data) {
  n <- nrow(data)
  rsample::make_splits(list(analysis = seq_len(n - n %/% 4L), assessment = seq.int(n - n %/% 4L + 1L, n)), data)
}
positional_holdout <- function(data) {
  n <- nrow(data)
  rsample::make_splits(list(analysis = seq_len(n - n %/% 3L), assessment = seq.int(n - n %/% 3L + 1L, n)), data)
}
positional_ensemble <- function(data) {
  n <- nrow(data)
  half <- list(seq_len(n %/% 2L), seq.int(n %/% 2L + 1L, n))
  rsample::manual_rset(lapply(half, function(P)
    rsample::make_splits(list(analysis = setdiff(seq_len(n), P), assessment = P), data)), c("first", "second"))
}
unstratified_selection <- function(d) rsample::vfold_cv(d, v = 2)
stratified_selection <- function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)

flip_truth <- function(data, rows) {
  data$truth[rows] <- ifelse(data$truth[rows] == "yes", "no", "yes")
  data
}

testthat::test_that("local grid search stays inside B, leaves P and C out and keeps nested candidates fixed", {
  for (pkg in c("tune", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- grid_binary_data()
  base <- grid_tuned_base(d)
  post <- list(method = "tailor_logistic", strategy = "holdout", splits = positional_holdout)
  for (strategy in c("two_stage", "joint")) {
    config <- list(strategy = strategy, model = list(method = "grid",
      grid = data.frame(degree = c(1L, 2L)), resamples = unstratified_selection))
    spec <- calibrated_workflow(base, "truth", splits = positional_va, seed = 5,
      tuning = config, postprocessing = post)
    audit <- record_postprocessing_fits()
    fitted <- generics::fit(spec, d)
    state <- workflows::extract_fit_engine(fitted)
    selection <- venn_abers_tuning(fitted)
    diagnostic <- state$postprocessing$Inductive1
    local <- selection$models$Inductive1$Holdout1
    testthat::expect_named(selection$models$Inductive1, "Holdout1")
    testthat::expect_identical(local$B, diagnostic$B)
    testthat::expect_identical(local$fold_id, "Inductive1")
    testthat::expect_identical(diagnostic$selected, local$selected)
    testthat::expect_identical(sort(unlist(lapply(local$resamples, `[[`, "assessment"), use.names = FALSE)),
      seq_along(diagnostic$B))
    fits <- audit$fits
    final <- fits[[length(fits)]]
    testthat::expect_identical(final$rows, diagnostic$B)
    testthat::expect_identical(final$calibration, diagnostic$P)
    search <- fits[-length(fits)]
    touched <- unlist(lapply(search, function(f) c(f$rows, f$calibration)), use.names = FALSE)
    testthat::expect_true(all(touched %in% diagnostic$B))
    testthat::expect_false(any(touched %in% c(diagnostic$P, diagnostic$C)))
    # Two stages: 2 candidates x 2 resamples of the raw base plus one final fit.
    # Joint: each candidate composition has one fixed inner pair, no further search.
    testthat::expect_length(fits, 5L)
    testthat::expect_identical(vapply(search, `[[`, TRUE, "tailor"), rep(strategy == "joint", 4L))
    testthat::expect_identical(selection$cost$observed_base_fits, 5L)
    testthat::expect_identical(local$cost$observed_base_fits, 4L)
    # A plan function leaves the nested joint VA fold count unknown in advance.
    testthat::expect_identical(local$budget$gross_base_fit_estimate, if (strategy == "joint") NA_real_ else 4)
    testthat::expect_true(all(vapply(unlist(local$decisions, recursive = FALSE), is.null, TRUE)))
    testthat::expect_identical(local$candidates$status, rep("eligible", 2L))
    # The unchanged VA core reproduces the output from the held C predictions.
    calibration <- state$calibration$Inductive1
    core <- venn_abers(calibration$truth, calibration$probability)
    p <- predict(state$models$Inductive1, d, type = "prob")$.pred_no
    testthat::expect_lte(max(abs(predict(core, p, loss = state$loss)$probability -
      predict(fitted, d, type = "prob")$.pred_no)), 1e-10)
    engine_coef <- function(s) stats::coef(workflows::extract_fit_engine(s$models$Inductive1))
    tailor_coef <- function(s) s$postprocessing$Inductive1$coefficients
    for (role in c("P", "C")) {
      changed <- generics::fit(spec, flip_truth(d, diagnostic[[role]]))
      other <- workflows::extract_fit_engine(changed)
      again <- venn_abers_tuning(changed)$models$Inductive1$Holdout1
      info <- paste(strategy, role)
      testthat::expect_identical(again[c("chosen", "scores", "resamples", "selected")],
        local[c("chosen", "scores", "resamples", "selected")], info = info)
      testthat::expect_identical(engine_coef(other), engine_coef(state), info = info)
      testthat::expect_identical(identical(tailor_coef(other), tailor_coef(state)), role == "C", info = info)
      # Held C predictions change only through the P-trained postprocessor.
      testthat::expect_identical(identical(other$calibration$Inductive1$probability,
        state$calibration$Inductive1$probability), role == "C", info = info)
    }
  }
})

testthat::test_that("every inner B, including each F pair, runs its own search", {
  testthat::skip_if_not_installed("tune")
  d <- grid_binary_data()
  base <- grid_tuned_base(d)
  post <- list(method = "sigmoid", sigmoid_input = "probability", strategy = "ensemble", splits = positional_ensemble)
  config <- list(strategy = "two_stage", model = list(method = "grid",
    grid = data.frame(degree = c(1L, 2L)), resamples = unstratified_selection))
  spec <- calibrated_workflow(base, "truth", splits = positional_va, seed = 5,
    prediction_model = "full", tuning = config, postprocessing = post)
  audit <- record_postprocessing_fits()
  fitted <- generics::fit(spec, d)
  state <- workflows::extract_fit_engine(fitted)
  selection <- venn_abers_tuning(fitted)
  testthat::expect_named(selection$models, c("Inductive1", "Full1"))
  for (path in names(selection$models)) {
    testthat::expect_named(selection$models[[path]], c("first", "second"))
    for (inner in c("first", "second")) {
      local <- selection$models[[path]][[inner]]
      plan <- state$postprocessing[[path]]$inner_folds[[inner]]
      testthat::expect_identical(local$B, plan$B, info = paste(path, inner))
      testthat::expect_identical(plan$selected, local$selected, info = paste(path, inner))
      testthat::expect_identical(local$inner_fold, inner)
    }
  }
  # Two paths x two inner pairs x (2 candidates x 2 resamples + 1 final fit).
  testthat::expect_length(audit$fits, 20L)
  testthat::expect_identical(selection$cost$observed_base_fits, 20L)
  # F excludes nothing but its own P: F searches may see held C rows.
  full <- state$postprocessing$Full1$inner_folds$first
  testthat::expect_identical(full$Q, seq_len(nrow(d)))
  # Mutating P_first leaves the first search alone; the second B contains P_first.
  source_first <- state$postprocessing$Inductive1$inner_folds$first
  changed <- venn_abers_tuning(generics::fit(spec, flip_truth(d, source_first$P)))
  testthat::expect_identical(changed$models$Inductive1$first[c("chosen", "scores")],
    selection$models$Inductive1$first[c("chosen", "scores")])
  testthat::expect_false(identical(changed$models$Inductive1$second$scores,
    selection$models$Inductive1$second$scores))
})

testthat::test_that("postprocessing grid covers Cross, full and both strategies with complete cost records", {
  for (pkg in c("tune", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- grid_binary_data(400L)
  base <- grid_tuned_base(d)
  expected <- c(two_stage.holdout.matched = 10L, two_stage.holdout.full = 15L,
    two_stage.ensemble.matched = 20L, two_stage.ensemble.full = 30L,
    joint.holdout.matched = 18L, joint.holdout.full = 39L)
  for (strategy in c("two_stage", "joint")) for (kind in c("holdout", "ensemble")) for (model in c("matched", "full")) {
    post <- if (kind == "holdout") list(method = "tailor_logistic", strategy = "holdout", prop = .3) else
      list(method = "sigmoid", sigmoid_input = "probability", strategy = "ensemble", v = 2)
    config <- list(strategy = strategy, model = list(method = "grid",
      grid = data.frame(degree = c(1L, 2L)), resamples = stratified_selection))
    fitted <- generics::fit(calibrated_workflow(base, "truth", method = "cross", v = 2, seed = 3,
      prediction_model = model, tuning = config, postprocessing = post), d)
    state <- workflows::extract_fit_engine(fitted)
    selection <- venn_abers_tuning(fitted)
    name <- paste(strategy, kind, model, sep = ".")
    inner <- if (kind == "holdout") "Holdout1" else c("Fold1", "Fold2")
    paths <- c("Fold1", "Fold2", if (model == "full") "Full1")
    testthat::expect_named(selection$models, paths)
    for (path in paths) testthat::expect_named(selection$models[[path]], inner, info = name)
    budget <- selection$models$Fold1[[1L]]$budget
    testthat::expect_identical(budget$inner_pairs, if (strategy == "joint" && kind == "ensemble") 2L else 1L, info = name)
    searches <- unlist(lapply(selection$models, function(path) lapply(path, function(x) x$cost$observed_base_fits)))
    testthat::expect_identical(selection$cost$observed_base_fits, sum(searches) + length(searches), info = name)
    # Joint ensemble candidates may fail inside a nested composition and are then excluded.
    if (!(strategy == "joint" && kind == "ensemble")) {
      testthat::expect_identical(selection$cost$observed_base_fits, expected[[name]], info = name)
    }
    testthat::expect_true(all(is.finite(as.matrix(predict(fitted, d, type = "prob")))), info = name)
    details <- venn_abers_details(fitted, d[1:3, ])
    testthat::expect_true(all(details$details$model_source %in% names(state$models)), info = name)
  }
})

grid_multiclass_data <- function(n = 360L, k = 3L) {
  withr::with_seed(21, {
    truth <- factor(rep(letters[seq_len(k)], length.out = n))
    data.frame(row_id = seq_len(n), x = rnorm(n, sd = 2) + c(-1, 1, -1, 1)[truth],
      z = rnorm(n, sd = 2) + c(-1, -1, 1, 1)[truth], w = rnorm(n, sd = 2) + c(1, -1, -1, 1)[truth], truth = truth)
  })
}

testthat::test_that("pairwise joint grid with postprocessing selects in every local pair B", {
  for (pkg in c("tune", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- grid_multiclass_data(900L)[c("row_id", "x", "truth")]
  base <- grid_tuned_base(d)
  post <- list(method = "tailor_logistic", strategy = "ensemble", v = 2)
  model <- list(method = "grid", grid = data.frame(degree = c(1L, 2L)), resamples = stratified_selection)
  args <- list(base, "truth", classification = "ovo", model_training = "pairwise",
    seed = 9, prediction_model = "full", postprocessing = post,
    tuning = list(strategy = "joint", model = model, joint_resamples = stratified_selection))
  fitted <- generics::fit(do.call(calibrated_workflow, args), d)
  state <- workflows::extract_fit_engine(fitted)
  selection <- venn_abers_tuning(fitted)
  # No global Cartesian pair grid: one complete assessment of the local rule.
  testthat::expect_identical(selection$joint$candidates$id, "joint_0001")
  testthat::expect_identical(selection$joint$candidates$status, "eligible")
  testthat::expect_named(selection$models, paste0("pair_", 1:3))
  for (pair in names(selection$models)) {
    component <- state$pair_models[[pair]]
    testthat::expect_named(selection$models[[pair]], c("Inductive1", "Full1"))
    for (path in names(selection$models[[pair]])) for (inner in c("Fold1", "Fold2")) {
      local <- selection$models[[pair]][[path]][[inner]]
      info <- paste(pair, path, inner)
      testthat::expect_identical(local$pair_id, pair, info = info)
      testthat::expect_identical(local$fold_id, path, info = info)
      testthat::expect_identical(local$B, component$postprocessing[[path]]$inner_folds[[inner]]$B, info = info)
      # Selection resamples cover exactly this pair's B rows, never other classes.
      testthat::expect_identical(sort(unlist(lapply(local$resamples, `[[`, "assessment"), use.names = FALSE)),
        seq_along(local$B), info = info)
      testthat::expect_true(all(d$truth[component$rows[local$B]] %in% component$class_levels), info = info)
      testthat::expect_true(local$chosen %in% local$candidates$id, info = info)
    }
  }
  testthat::expect_true(all(is.finite(as.matrix(predict(fitted, d, type = "prob")))))
  # Without local resamples the pairwise grid cannot fall back to a global choice.
  args$tuning$model$resamples <- NULL
  testthat::expect_error(do.call(calibrated_workflow, args), "resamples")
  args$tuning <- list(strategy = "two_stage", model = model)
  two <- venn_abers_tuning(generics::fit(do.call(calibrated_workflow, args), d))
  testthat::expect_null(two$joint)
  testthat::expect_named(two$models$pair_2$Full1, c("Fold1", "Fold2"))
  testthat::expect_true(all(two$models$pair_2$Full1$Fold1$candidates$status == "eligible"))
})

testthat::test_that("shared multiclass calibration grid keeps S out of every final B and P", {
  for (pkg in c("tune", "tailor", "probably", "nnet")) testthat::skip_if_not_installed(pkg)
  d <- grid_multiclass_data(2400L)
  base <- grid_tuned_base(d, parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE, maxit = 1000))
  for (strategy in c("joint", "two_stage")) {
    config <- list(strategy = strategy, selection_prop = .5,
      model = list(method = "grid", grid = data.frame(degree = c(1L, 2L)), resamples = stratified_selection),
      calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")), resamples = stratified_selection))
    fitted <- generics::fit(calibrated_workflow(base, "truth", classification = "ovr", seed = 4,
      loss = tune::tune("va_loss"), tuning = config,
      postprocessing = list(method = "tailor_multinomial", strategy = "holdout", prop = .3)), d)
    state <- workflows::extract_fit_engine(fitted)
    selection <- venn_abers_tuning(fitted)
    roles <- selection$calibration
    testthat::expect_length(intersect(roles$selection_rows, roles$final_rows), 0L)
    plan <- state$postprocessing$Inductive1
    global <- function(rows) roles$final_rows[rows]
    testthat::expect_length(intersect(global(c(plan$B, plan$P, plan$C)), roles$selection_rows), 0L)
    testthat::expect_identical(selection$models$Inductive1$Holdout1$B, plan$B)
    testthat::expect_identical(roles$candidates$status, c("eligible", "eligible"))
    # Candidate compositions each ran their own local B search.
    decisions <- roles$decisions[[1L]][[1L]]
    testthat::expect_named(decisions$Inductive1, "Holdout1")
    testthat::expect_true(decisions$Inductive1$Holdout1$chosen %in% c("model_001", "model_002"))
    if (strategy == "two_stage") for (id in names(roles$resamples)) {
      choices <- lapply(roles$decisions, function(candidate) candidate[[id]]$Inductive1$Holdout1$chosen)
      testthat::expect_identical(choices[[1L]], choices[[2L]])
    }
    testthat::expect_identical(state$loss, roles$selected[[1L]])
  }
})

testthat::test_that("regression grids choose inside B for isotonic and linear Tailor", {
  for (pkg in c("tune", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- withr::with_seed(8, {
    x <- rnorm(240)
    data.frame(row_id = seq_len(240), x = x, truth = 3 + x - .5 * x^2 + rnorm(240, sd = .5))
  })
  base <- grid_tuned_base(d, parsnip::linear_reg(engine = "lm"))
  for (method in c("isotonic_numeric", "tailor_linear")) for (strategy in c("two_stage", "joint")) {
    config <- list(strategy = strategy, model = list(method = "grid",
      grid = data.frame(degree = c(1L, 2L)), resamples = unstratified_selection))
    spec <- calibrated_workflow(base, "truth", method = "cross", v = 2, seed = 2, prediction_model = "full",
      tuning = config, postprocessing = list(method = method, strategy = "holdout", prop = .3))
    fitted <- generics::fit(spec, d)
    state <- workflows::extract_fit_engine(fitted)
    selection <- venn_abers_tuning(fitted)
    info <- paste(method, strategy)
    for (path in c("Fold1", "Fold2", "Full1")) {
      plan <- state$postprocessing[[path]]
      local <- selection$models[[path]]$Holdout1
      testthat::expect_identical(local$B, plan$B, info = info)
      testthat::expect_identical(sort(unlist(lapply(local$resamples, `[[`, "assessment"), use.names = FALSE)),
        seq_along(plan$B), info = info)
      testthat::expect_identical(plan$selected, local$selected, info = info)
      testthat::expect_true(local$chosen %in% local$candidates$id[local$candidates$status == "eligible"], info = info)
    }
    testthat::expect_true(all(is.finite(predict(fitted, d)$.pred)), info = info)
    testthat::expect_identical(venn_abers_details(fitted, d[1:2, ])$details$inner_fold, rep("Holdout1", 4L), info = info)
  }
})

testthat::test_that("tuned postprocessing validates static plans and has no final fallback", {
  for (pkg in c("tune", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- grid_binary_data()
  base <- grid_tuned_base(d)
  model <- list(method = "grid", grid = data.frame(degree = c(1L, 2L)), resamples = unstratified_selection)
  static <- list(method = "tailor_logistic", strategy = "holdout",
    splits = rsample::initial_split(d[seq_len(180), ], prop = .7))
  testthat::expect_error(calibrated_workflow(base, "truth", postprocessing = static,
    tuning = list(strategy = "joint", model = model)), "provider function")
  testthat::expect_error(calibrated_workflow(base, "truth", postprocessing = static, loss = tune::tune("va_loss"),
    tuning = list(strategy = "two_stage", model = model, selection_prop = .2,
      calibration = list(method = "grid", grid = data.frame(va_loss = "log"), resamples = unstratified_selection))),
    "provider function")
  testthat::expect_s3_class(calibrated_workflow(base, "truth", postprocessing = static,
    tuning = list(strategy = "two_stage", model = model)), "workflow")
  # Every candidate succeeds on the smaller selection data; the selected final
  # B fit fails and the error surfaces without trying the other candidate.
  sizes <- integer()
  fit_method <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  registerS3method("fit", "workflow", function(object, data, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), c("vennabers", "vennabers_grid"))) {
      sizes <<- c(sizes, nrow(data))
      if (nrow(data) > 120L) stop("planned final B failure")
    }
    fit_method(object, data, ...)
  }, envir = asNamespace("generics"))
  on.exit(registerS3method("fit", "workflow", fit_method, envir = asNamespace("generics")), add = TRUE)
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", seed = 1,
    tuning = list(strategy = "two_stage", model = model),
    postprocessing = list(method = "sigmoid", sigmoid_input = "probability", strategy = "holdout", prop = .25)), d),
    "Holdout1 sigmoid.*planned final B failure")
  testthat::expect_identical(sum(sizes > 120L), 1L)
})

testthat::test_that("each B search has its own RNG stream, so another B's data cannot shift it", {
  for (pkg in c("tune", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- grid_binary_data(360L)
  base <- grid_tuned_base(d)
  # Three positional VA folds: P of Fold3 (end of the second block) lies in B of Fold1.
  va <- function(data) {
    blocks <- split(seq_len(nrow(data)), rep(1:3, each = nrow(data) / 3))
    rsample::manual_rset(lapply(blocks, function(C)
      rsample::make_splits(list(analysis = setdiff(seq_len(nrow(data)), C), assessment = C), data)), paste0("Fold", 1:3))
  }
  config <- list(strategy = "two_stage", model = list(method = "grid",
    grid = data.frame(degree = c(1L, 2L)), resamples = stratified_selection))
  spec <- calibrated_workflow(base, "truth", method = "cross", splits = va, seed = 5, tuning = config,
    postprocessing = list(method = "tailor_logistic", strategy = "holdout", splits = positional_holdout))
  state <- workflows::extract_fit_engine(generics::fit(spec, d))
  plan <- state$postprocessing
  testthat::expect_true(any(plan$Fold3$P %in% plan$Fold1$B))
  selection <- state$tuning$models
  seeds <- vapply(selection, function(x) x$Holdout1$search_seed, 1L)
  testthat::expect_false(anyDuplicated(seeds) > 0L)
  changed <- workflows::extract_fit_engine(generics::fit(spec, flip_truth(d, plan$Fold3$P)))
  testthat::expect_identical(changed$tuning$models$Fold3$Holdout1[c("resamples", "scores", "chosen", "search_seed")],
    selection$Fold3$Holdout1[c("resamples", "scores", "chosen", "search_seed")])
  testthat::expect_identical(stats::coef(workflows::extract_fit_engine(changed$models$Fold3)),
    stats::coef(workflows::extract_fit_engine(state$models$Fold3)))
})

testthat::test_that("pairwise joint calibration selection and GLM-link Sigmoid work with local searches", {
  for (pkg in c("tune", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- grid_multiclass_data(1200L)[c("row_id", "x", "truth")]
  base <- grid_tuned_base(d)
  config <- list(strategy = "joint", selection_prop = .5,
    model = list(method = "grid", grid = data.frame(degree = c(1L, 2L)), resamples = stratified_selection),
    calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")), resamples = stratified_selection))
  fitted <- generics::fit(calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
    seed = 12, loss = tune::tune("va_loss"), splits = function(x) rsample::initial_split(x, prop = .8, strata = truth, pool = 0),
    tuning = config, postprocessing = list(method = "sigmoid", sigmoid_input = "glm_link", strategy = "holdout", prop = .3)), d)
  selection <- venn_abers_tuning(fitted)
  state <- workflows::extract_fit_engine(fitted)
  roles <- selection$calibration
  testthat::expect_identical(roles$candidates$status, c("eligible", "eligible"))
  for (pair in names(state$pair_models)) {
    component <- state$pair_models[[pair]]
    local <- selection$models[[pair]]$Inductive1$Holdout1
    testthat::expect_identical(local$B, component$postprocessing$Inductive1$B)
    testthat::expect_identical(component$postprocessing$Inductive1$input, "glm_link")
    global <- roles$final_rows[component$rows[c(local$B, component$postprocessing$Inductive1$P)]]
    testthat::expect_length(intersect(global, roles$selection_rows), 0L)
  }
  testthat::expect_true(all(is.finite(as.matrix(predict(fitted, d, type = "prob")))))
})

testthat::test_that("postprocessed Grid searches exclude S and held P/C predictors with learned recipes", {
  for (pkg in c("tune", "recipes", "rpart")) testthat::skip_if_not_installed(pkg)
  d <- withr::with_seed(84, data.frame(x = rnorm(480), truth = rnorm(480)))
  base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, d) |>
    recipes::step_normalize(x) |> recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::decision_tree(tree_depth = tune::tune("depth"), min_n = 10) |>
      parsnip::set_mode("regression") |> parsnip::set_engine("rpart"))
  keys <- c("chosen", "selected", "candidates", "scores", "resamples", "search_seed", "budget", "cost")
  for (strategy in c("joint", "two_stage")) {
    spec <- calibrated_workflow(base, "truth", method = "cross", splits = positional_ensemble,
      prediction_model = "full", seed = 31, epsilon = tune::tune("moderation"),
      postprocessing = list(method = "isotonic_numeric", strategy = "ensemble", splits = positional_ensemble),
      tuning = list(strategy = strategy,
        selection_split = function(d) rsample::make_splits(list(analysis = 1:240, assessment = 241:480), d),
        model = list(method = "grid", grid = data.frame(degree = 1:2, depth = 1:2), resamples = unstratified_selection),
        calibration = list(method = "grid", grid = data.frame(moderation = .3), resamples = unstratified_selection)))
    fit <- generics::fit(spec, d)
    state <- workflows::extract_fit_engine(fit)
    original <- venn_abers_tuning(fit)
    testthat::expect_named(original$models, c("first", "second", "Full1"))
    changed <- d
    changed$x[241:480] <- changed$x[241:480] + 100
    changed$truth[241:480] <- changed$truth[241:480] - 100
    other <- generics::fit(spec, changed)
    for (n in c(0L, 1L, 7L)) {
      testthat::expect_identical(predict(other, d[seq_len(n), ]), predict(fit, d[seq_len(n), ]), info = paste(strategy, n))
      testthat::expect_identical(venn_abers_details(other, d[seq_len(n), ]), venn_abers_details(fit, d[seq_len(n), ]),
        info = paste(strategy, n))
    }
    testthat::expect_false(identical(venn_abers_tuning(other)$calibration$scores, original$calibration$scores), info = strategy)
    for (path in names(original$models)) for (inner in c("first", "second")) {
      local <- original$models[[path]][[inner]]
      plan <- state$postprocessing[[path]]$inner_folds[[inner]]
      testthat::expect_identical(venn_abers_tuning(other)$models[[path]][[inner]][keys], local[keys], info = paste(strategy, path, inner))
      testthat::expect_equal(as.data.frame(local$native$splits[[1L]]$data), d[plan$B, ], ignore_attr = TRUE)
      testthat::expect_identical(plan$selected, local$selected)
      testthat::expect_true(all(vapply(unlist(local$decisions, recursive = FALSE), is.null, TRUE)))
    }
    for (inner in c("first", "second")) testthat::expect_equal(
      rlang::eval_tidy(workflows::extract_spec_parsnip(state$models$Full1$models[[inner]]$workflow)$args$min_n), 10)
    # In a Cross/full ensemble, P_first belongs to another B; its own search
    # must remain unchanged while those sibling and full searches see the data.
    for (role in c("P", "C")) {
      changed <- d
      rows <- state$postprocessing$first$inner_folds$first[[role]]
      changed$x[rows] <- changed$x[rows] + 20
      changed$truth[rows] <- changed$truth[rows] + 30
      other <- generics::fit(spec, changed)
      again <- venn_abers_tuning(other)
      testthat::expect_identical(workflows::extract_fit_engine(other)$postprocessing$first$inner_folds$first,
        state$postprocessing$first$inner_folds$first, info = paste(strategy, role))
      testthat::expect_identical(again$models$first$first[keys], original$models$first$first[keys], info = paste(strategy, role))
      testthat::expect_identical(tune::collect_predictions(again$models$first$first$native),
        tune::collect_predictions(original$models$first$first$native), info = paste(strategy, role))
      full_inner <- if (role == "P") "first" else "second"
      testthat::expect_false(identical(again$models$Full1[[full_inner]]$scores,
        original$models$Full1[[full_inner]]$scores), info = paste(strategy, role))
    }
  }
})

testthat::test_that("Grid excludes incomplete inner ensembles, isolates reports and keeps stable ties", {
  for (pkg in c("tune", "recipes")) testthat::skip_if_not_installed(pkg)
  d <- withr::with_seed(8, data.frame(x = rnorm(240), truth = rnorm(240)))
  base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, d) |>
    recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(parsnip::linear_reg(engine = "lm"))
  # Degree 100 cannot be fitted in any inner B; degree 2 fits but its public
  # predictions drop a row. Either way the invalid candidate comes first, so
  # the exact tie behind it must still go to the first eligible Grid ID.
  original_predict <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  registerS3method("predict", "workflow", function(object, new_data, ...) {
    result <- original_predict(object, new_data, ...)
    if (!inherits(workflows::extract_spec_parsnip(object), c("vennabers", "vennabers_grid")) &&
        identical(workflows::extract_preprocessor(object)$steps[[1L]]$degree, 2L)) result <- result[-1L, , drop = FALSE]
    result
  }, envir = asNamespace("stats"))
  on.exit(registerS3method("predict", "workflow", original_predict, envir = asNamespace("stats")), add = TRUE)
  report <- yardstick::new_numeric_metric(function(...) stop("report unavailable"), "minimize")
  for (strategy in c("joint", "two_stage")) for (invalid in c(100L, 2L)) {
    fit <- generics::fit(calibrated_workflow(base, "truth", splits = positional_va, seed = 5,
      postprocessing = list(method = "isotonic_numeric", strategy = "ensemble", splits = positional_ensemble),
      tuning = list(strategy = strategy, model = list(method = "grid",
        grid = data.frame(degree = c(invalid, 1L, 1L)), resamples = unstratified_selection,
        report_metrics = list(broken = report)))), d)
    choices <- venn_abers_tuning(fit)$models$Inductive1
    testthat::expect_named(choices, c("first", "second"))
    for (local in choices) {
      testthat::expect_identical(local$candidates$status, c("excluded", "eligible", "eligible"))
      testthat::expect_identical(local$chosen, "model_002")
      testthat::expect_identical(local$candidates$estimate[2L], local$candidates$estimate[3L])
      testthat::expect_true(is.na(local$candidates$estimate[1L]))
      testthat::expect_match(local$candidates$reason[1L], "omit, duplicate, or mismatch")
      # Joint candidates run the whole inner ensemble; two_stage the raw workflow.
      notes <- paste(tune::collect_notes(local$native)$note, collapse = " ")
      if (strategy == "joint") testthat::expect_match(notes, "first isotonic_numeric") else testthat::expect_true(nzchar(notes))
      testthat::expect_match(local$report_metrics[[2L]][[1L]]$broken, "report unavailable")
      testthat::expect_identical(local$cost$excluded_candidates, 1L)
    }
    testthat::expect_true(all(is.finite(predict(fit, d[1:7, ])$.pred)))
  }
})

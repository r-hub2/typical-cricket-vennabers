testthat::test_that("holdout and ensemble fits consume named inner pairs without shape detection", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:80), truth = factor(rep(c("a", "b"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  va <- function(data) rsample::make_splits(list(analysis = 1:60, assessment = 61:80), data)
  inner <- function(data) {
    P <- list(seq_len(nrow(data) %/% 3L), seq.int(nrow(data) %/% 3L + 1L, nrow(data)))
    rsample::manual_rset(lapply(P, function(rows) rsample::make_splits(
      list(analysis = setdiff(seq_len(nrow(data)), rows), assessment = rows), data)), c("B", "P"))
  }
  fit_pair <- .fit_postprocessing_pair
  seen <- list()
  testthat::local_mocked_bindings(.fit_postprocessing_pair = function(workflow, data, plan, fold,
      domain, settings, id, inner_id, ...) {
    members <- if (settings$postprocessing$strategy == "holdout") list(Holdout1 = fold$postprocessing) else
      fold$postprocessing
    testthat::expect_identical(plan, members[[inner_id]])
    seen[[length(seen) + 1L]] <<- paste(id, inner_id, sep = "/")
    fit_pair(workflow, data, plan, fold, domain, settings, id, inner_id, ...)
  })
  for (strategy in c("holdout", "ensemble")) for (prediction_model in c("matched", "full")) {
    seen <- list()
    post <- list(method = "sigmoid", sigmoid_input = "glm_link", strategy = strategy,
      splits = if (strategy == "ensemble") inner else function(data) inner(data)$splits[[1L]])
    fit <- generics::fit(calibrated_workflow(base, "truth", splits = va, seed = 42,
      prediction_model = prediction_model, postprocessing = post), d)
    state <- workflows::extract_fit_engine(fit)
    expected <- if (strategy == "holdout") "Holdout1" else c("B", "P")
    plan <- state$folds$Inductive1$postprocessing
    diagnostics <- state$postprocessing$Inductive1
    # As in v0.9.0, a holdout plan is the bare B/P pair; an ensemble plan names its pairs.
    if (strategy == "holdout") {
      testthat::expect_named(plan, c("Q", "B", "P", "analysis", "assessment"))
      testthat::expect_identical(plan, diagnostics[names(plan)])
    } else {
      testthat::expect_named(plan, expected)
      testthat::expect_identical(plan, lapply(diagnostics$inner_folds, function(p) p[names(plan[[1L]])]))
    }
    sources <- if (prediction_model == "full") c("Inductive1", "Full1") else "Inductive1"
    testthat::expect_identical(unlist(seen), unlist(lapply(sources, function(id) paste(id, expected, sep = "/"))))
    source <- if (prediction_model == "full") "Full1" else "Inductive1"
    record <- state$postprocessing[[source]]
    # The v0.9.0 record format: an ensemble record holds its pairs only under
    # inner_folds, and details derive the per-pair list-columns from them.
    # Q has 60 rows and F 80; each inner plan holds out a first third as P.
    sizes <- if (prediction_model == "full") c(54L, 26L) else c(40L, 20L)
    if (strategy == "ensemble") {
      testthat::expect_false(any(c("inner_fold", "B_size", "P_size") %in% names(record)))
      testthat::expect_identical(vapply(record$inner_folds, `[[`, 1L, "B_size"), c(B = sizes[[1L]], P = sizes[[2L]]))
      testthat::expect_identical(vapply(record$inner_folds, `[[`, 1L, "P_size"), c(B = sizes[[2L]], P = sizes[[1L]]))
    } else {
      testthat::expect_identical(c(record$B_size, record$P_size), sizes)
    }
    for (n in c(0L, 1L, 7L)) {
      details <- venn_abers_details(fit, d[seq_len(n), "x", drop = FALSE])
      testthat::expect_identical(details$prediction, predict(fit, d[seq_len(n), "x", drop = FALSE], type = "prob"))
      expected_fields <- if (strategy == "holdout") list(inner_fold = rep("Holdout1", n),
        B_size = stats::setNames(rep(record$B_size, n), rep(source, n)),
        P_size = stats::setNames(rep(record$P_size, n), rep(source, n))) else list(
        inner_fold = rep(list(c("B", "P")), n),
        B_size = rep(list(vapply(record$inner_folds, `[[`, 1L, "B_size")), n),
        P_size = rep(list(vapply(record$inner_folds, `[[`, 1L, "P_size")), n))
      for (field in names(expected_fields)) {
        testthat::expect_identical(details$details[[field]], expected_fields[[field]])
      }
    }
  }
})

testthat::test_that("binary ensembles average complete pairs before the unchanged VA core", {
  skip_if_no_integration("tailor", "probably")
  d <- data.frame(x = sin(1:160), truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  # Unequal P sizes distinguish equal pair weights from row-count weights.
  inner <- function(data) {
    n <- nrow(data)
    rows <- as.vector(t(matrix(order(data$truth), ncol = 2)))
    held <- list(rows[1:20], rows[21:50], rows[51:n])
    rsample::manual_rset(lapply(held, function(p) rsample::make_splits(
      list(analysis = setdiff(seq_len(n), p), assessment = p), data)), c("small", "medium", "large"))
  }
  for (method in c("inductive", "cross")) for (kind in c("tailor_logistic", "probability", "glm_link")) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = "ensemble", splits = inner)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    fit <- generics::fit(calibrated_workflow(base, "truth", method = method, seed = 24, postprocessing = post), d)
    state <- workflows::extract_fit_engine(fit)
    held <- new <- truth <- ids <- list()
    for (id in names(state$models)) {
      ensemble <- state$models[[id]]
      diagnostics <- state$postprocessing[[id]]
      testthat::expect_length(ensemble$models, 3)
      testthat::expect_identical(names(ensemble$models), c("small", "medium", "large"))
      predictions <- calibration <- list()
      for (j in names(ensemble$models)) {
        p <- diagnostics$inner_folds[[j]]
        testthat::expect_setequal(p$B, setdiff(p$Q, p$P))
        testthat::expect_length(intersect(p$B, p$P), 0)
        testthat::expect_length(intersect(p$Q, p$C), 0)
        # Refit through the public holdout seam on the same B/P/C rows.
        refpost <- post; refpost$strategy <- "holdout"
        refpost$splits <- rsample::make_splits(p[c("analysis", "assessment")], d[p$Q, ])
        va <- rsample::make_splits(list(analysis = p$Q, assessment = p$C), d)
        ref <- generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = refpost), d)
        complete <- workflows::extract_fit_engine(ref)$models[[1]]
        predictions[[j]] <- predict(complete, d["x"], type = "prob")$.pred_a
        calibration[[j]] <- predict(complete, d[p$C, "x", drop = FALSE], type = "prob")$.pred_a
      }
      new[[id]] <- unname(rowMeans(do.call(cbind, predictions)))
      held[[id]] <- unname(rowMeans(do.call(cbind, calibration)))
      testthat::expect_equal(predict(ensemble, d["x"], type = "prob")$.pred_a, new[[id]], tolerance = 1e-10)
      testthat::expect_equal(state$calibration[[id]]$probability, held[[id]], tolerance = 1e-10)
      testthat::expect_identical(diagnostics$observed_base_fits, 3L)
      testthat::expect_identical(diagnostics$observed_postprocessor_fits, 3L)
      truth[[id]] <- d$truth[diagnostics$C] == "a"
      ids[[id]] <- rep(id, length(diagnostics$C))
    }
    core <- if (method == "cross") cross_venn_abers(unlist(truth), unlist(held), unlist(ids)) else
      venn_abers(truth[[1]], held[[1]])
    for (loss in c("log", "brier")) {
      expected <- predict(core, if (method == "cross") do.call(cbind, new) else new[[1]], loss = loss)
      testthat::expect_equal(predict(fit, d["x"], type = "prob", loss = loss)$.pred_a, expected$probability, tolerance = 1e-10)
    }
  }
})

testthat::test_that("ensemble details preserve every pair and support both events and rounding", {
  skip_if_no_integration("tailor", "probably")
  d <- data.frame(x = sin(1:160), truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (kind in c("tailor_logistic", "probability", "glm_link")) for (method in c("inductive", "cross"))
    for (event in c("first", "second")) for (precision in list(NULL, 2L)) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = "ensemble", v = 3L)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    fit <- generics::fit(calibrated_workflow(base, "truth", method = method, event_level = event,
      precision = precision, seed = 24, postprocessing = post), d)
    state <- workflows::extract_fit_engine(fit)
    column <- paste0(".pred_", state$event_class)
    held <- truth <- ids <- new <- list()
    for (id in names(state$models)) {
      p <- state$postprocessing[[id]]
      testthat::expect_identical(sort(unlist(lapply(p$inner_folds, `[[`, "P"), use.names = FALSE)), sort(p$Q))
      held[[id]] <- predict(state$models[[id]], d[p$C, ], type = "prob")[[column]]
      new[[id]] <- predict(state$models[[id]], d, type = "prob")[[column]]
      truth[[id]] <- d$truth[p$C] == state$event_class
      ids[[id]] <- rep(id, length(p$C))
    }
    core <- if (method == "cross") cross_venn_abers(unlist(truth), unlist(held), unlist(ids), precision = precision) else
      venn_abers(truth[[1]], held[[1]], precision = precision)
    for (loss in c("log", "brier")) {
      expected <- predict(core, if (method == "cross") do.call(cbind, new) else new[[1]], loss = loss)
      details <- venn_abers_details(fit, d["x"], loss = loss)
      testthat::expect_lte(max(abs(details$prediction[[column]] - expected$probability)), 1e-10)
      testthat::expect_identical(details$prediction, predict(fit, d["x"], type = "prob", loss = loss))
      testthat::expect_identical(details$postprocessing, state$postprocessing)
      testthat::expect_identical(details$details$inner_fold[[1]], names(state$postprocessing[[1]]$inner_folds))
      testthat::expect_equal(details$details$B_size[[1]], vapply(state$postprocessing[[1]]$inner_folds, `[[`, 1L, "B_size"))
      testthat::expect_identical(as.character(predict(fit, d["x"], type = "class", loss = loss)$.pred_class),
        ifelse(expected$probability >= .5, state$event_class, setdiff(levels(d$truth), state$event_class)))
      for (n in 0:1) {
        small <- venn_abers_details(fit, d[seq_len(n), "x", drop = FALSE], loss = loss)
        testthat::expect_equal(nrow(small$prediction), n)
        testthat::expect_identical(names(small$details), names(details$details))
      }
    }
    count <- if (method == "cross") 15L else 3L
    testthat::expect_output(print(state), paste0("Base fits: ", count))
    testthat::expect_output(print(state), paste0("Models per nonempty prediction: ", count))
  }
})

testthat::test_that("ensemble requests and local partitions fail before any learning", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:80), truth = factor(rep(c("a", "b"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  post <- list(method = "sigmoid", sigmoid_input = "probability", strategy = "ensemble", v = 3)
  for (v in list(NULL, TRUE, NA_real_, Inf, numeric(), c(2, 3), 1, 2.5, matrix(3), .Machine$integer.max + 1)) {
    bad <- post; bad["v"] <- list(v)
    testthat::expect_error(calibrated_workflow(base, "truth", postprocessing = bad), "postprocessing.*v")
  }
  for (bad in list(post[-4], c(post, prop = .3), c(post, list(splits = function(data) NULL)),
                  c(post, v = 3), c(post, unknown = 1))) {
    testthat::expect_error(calibrated_workflow(base, "truth", postprocessing = bad), "[Pp]ostprocessing")
  }
  va <- function(data) rsample::make_splits(list(analysis = 1:60, assessment = 61:80), data)
  plan <- function(data) {
    rsample::manual_rset(lapply(list(1:20, 21:40, 41:60), function(p)
      rsample::make_splits(list(analysis = setdiff(1:60, p), assessment = p), data)), c("one", "two", "three"))
  }
  original <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  guard <- function(object, data, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) stop("UNEXPECTED LEARNING")
    original(object, data, ...)
  }
  registerS3method("fit", "workflow", guard, envir = asNamespace("generics"))
  on.exit(registerS3method("fit", "workflow", original, envir = asNamespace("generics")), add = TRUE)
  check <- function(provider, pattern = "B/P plan") {
    request <- post; request$v <- NULL; request$splits <- provider
    testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = request), d), pattern)
  }
  check(function(data) list(analysis = 1:40, assessment = 41:60))
  check(function(data) rsample::vfold_cv(data, v = 3, repeats = 2))
  check(function(data) rsample::bootstraps(data, times = 3))
  check(function(data) { p <- plan(data); p$splits[[3]] <- p$splits[[2]]; p })
  check(function(data) { p <- plan(data); p$splits <- p$splits[1:2]; p[1:2, ] })
  check(function(data) { p <- plan(data); p$id[2] <- p$id[1]; p })
  for (mutate in list(function(x) { x$x[1] <- 99; x }, function(x) x[nrow(x):1, ],
                     function(x) { x$x <- as.integer(x$x); x }, function(x) x[rev(names(x))],
                     function(x) { levels(x$truth) <- rev(levels(x$truth)); x })) {
    check(function(data) plan(mutate(data)), "binding")
  }
  check(function(data) plan(d), "binding")
  for (indices in list(list(analysis = c(21:60, 21L), assessment = 1:20),
                       list(analysis = 20:60, assessment = 1:20),
                       list(analysis = 22:60, assessment = 1:20),
                       list(analysis = 21:60, assessment = c(1:19, 81L)),
                       list(analysis = seq(2L, 60L, 2L), assessment = seq(1L, 59L, 2L)))) {
    check(function(data) { p <- plan(data); p$splits[[3]]$in_id <- indices$analysis
      p$splits[[3]]$out_id <- indices$assessment; p })
  }
  # Valid partitions with class-deficient B/P and too-small P are rejected too.
  check(function(data) rsample::manual_rset(lapply(list(seq(1L, 59L, 2L), seq(2L, 60L, 2L)), function(p)
    rsample::make_splits(list(analysis = setdiff(1:60, p), assessment = p), data)), c("odd", "even")), "must contain all")
  request <- post; request$v <- 100L
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = request), d), "not reduced")
  # The final local Q must be validated before the first pair can train.
  request <- post; request$v <- NULL
  request$splits <- function(data) if (data$x[1] == d$x[1]) NULL else rsample::vfold_cv(data, v = 2, strata = truth)
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", method = "cross", seed = 24, postprocessing = request), d), "B/P plan")
})

testthat::test_that("inner fits consume one stable RNG stream and preserve all required pairs", {
  skip_if_no_integration("tailor", "probably")
  d <- data.frame(x = sin(1:160), truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  fit_method <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  predict_method <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  audit <- new.env(); audit$draws <- numeric(); audit$predictions <- 0L
  instrumented_fit <- function(object, data, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) audit$draws <- c(audit$draws, runif(1))
    fit_method(object, data, ...)
  }
  instrumented_predict <- function(object, new_data, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) audit$predictions <- audit$predictions + 1L
    predict_method(object, new_data, ...)
  }
  registerS3method("fit", "workflow", instrumented_fit, envir = asNamespace("generics"))
  registerS3method("predict", "workflow", instrumented_predict, envir = asNamespace("stats"))
  on.exit({
    registerS3method("fit", "workflow", fit_method, envir = asNamespace("generics"))
    registerS3method("predict", "workflow", predict_method, envir = asNamespace("stats"))
  }, add = TRUE)
  original_optim <- stats::optim
  original_logistic <- probably::cal_estimate_logistic
  audit$post <- 0L
  testthat::local_mocked_bindings(optim = function(...) { audit$post <- audit$post + 1L; original_optim(...) }, .package = "stats")
  testthat::local_mocked_bindings(cal_estimate_logistic = function(...) {
    audit$post <- audit$post + 1L; original_logistic(...)
  }, .package = "probably")
  for (method in c("inductive", "cross")) for (kind in c("tailor_logistic", "probability", "glm_link")) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = "ensemble", v = 3L)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    spec <- calibrated_workflow(base, "truth", method = method, seed = 24, postprocessing = post)
    # Independently consume just the automatic plans, followed by one draw per fit.
    expected <- withr::with_seed(24, {
      outer <- if (method == "cross") rsample::vfold_cv(d, v = 5, strata = truth, pool = 0)$splits else
        list(rsample::initial_split(d, .8, strata = truth, pool = 0))
      for (q in outer) rsample::vfold_cv(rsample::analysis(q), v = 3, strata = truth, pool = 0)
      runif(length(outer) * 3)
    })
    set.seed(17); before <- .Random.seed
    audit$draws <- numeric(); audit$post <- 0L
    fit <- generics::fit(spec, d)
    testthat::expect_identical(.Random.seed, before)
    testthat::expect_identical(audit$draws, expected)
    testthat::expect_equal(audit$post, length(expected))
    state <- workflows::extract_fit_engine(fit)
    testthat::expect_equal(sum(vapply(state$models, function(x) length(x$models), 1L)), length(expected))
    for (n in c(0L, 1L, 160L)) {
      audit$predictions <- 0L; post_before <- audit$post; fit_before <- audit$draws
      predict(fit, d[seq_len(n), "x", drop = FALSE], type = "prob")
      testthat::expect_equal(audit$predictions, if (n) length(expected) else 0L)
      testthat::expect_identical(audit$post, post_before)
      testthat::expect_identical(audit$draws, fit_before)
    }
    # NULL seed uses exactly the same caller stream, without a restart per pair.
    set.seed(24); audit$draws <- numeric()
    unseeded <- generics::fit(calibrated_workflow(base, "truth", method = method, postprocessing = post), d)
    after <- .Random.seed
    testthat::expect_identical(audit$draws, expected)
    set.seed(24)
    outer <- if (method == "cross") rsample::vfold_cv(d, v = 5, strata = truth, pool = 0)$splits else
      list(rsample::initial_split(d, .8, strata = truth, pool = 0))
    for (q in outer) rsample::vfold_cv(rsample::analysis(q), v = 3, strata = truth, pool = 0)
    runif(length(expected))
    testthat::expect_identical(after, .Random.seed)
    testthat::expect_identical(predict(fit, d, type = "prob"), predict(unseeded, d, type = "prob"))
  }
})

testthat::test_that("pure local plans allow different J and replay in reverse call order", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:160), row_id = 1:160, truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  outer <- rsample::manual_rset(lapply(list(1:40, 41:90, 91:160), function(c)
    rsample::make_splits(list(analysis = setdiff(1:160, c), assessment = c), d)), c("one", "two", "three"))
  provider <- function(data) withr::with_seed(37,
    rsample::vfold_cv(data, v = if (nrow(data) > 100) 3L else 2L, strata = truth, pool = 0))
  prepared <- lapply(outer$splits, function(s) provider(rsample::analysis(s)))
  resolver <- local_calibration_splits(prepared, names(d))
  resolver <- unserialize(serialize(resolver, NULL))
  for (i in rev(seq_along(prepared))) {
    testthat::expect_identical(unname(lapply(resolver(rsample::analysis(outer$splits[[i]]))$splits, as.integer, "assessment")),
      lapply(prepared[[i]]$splits, as.integer, "assessment"))
  }
  fit <- function(splits, provider) generics::fit(calibrated_workflow(base, "truth", method = "cross", splits = splits,
    postprocessing = list(method = "sigmoid", sigmoid_input = "probability", strategy = "ensemble", splits = provider)), d)
  first <- fit(outer, provider)
  replay <- fit(outer[3:1, ], resolver)
  original <- workflows::extract_fit_engine(first); reverse <- workflows::extract_fit_engine(replay)
  testthat::expect_equal(unname(vapply(original$models, function(x) length(x$models), 1L)), c(3L, 3L, 2L))
  for (id in names(original$models)) {
    testthat::expect_identical(original$postprocessing[[id]], reverse$postprocessing[[id]])
    testthat::expect_identical(predict(original$models[[id]], d, type = "prob"), predict(reverse$models[[id]], d, type = "prob"))
  }
  testthat::expect_equal(predict(first, d, type = "prob"), predict(replay, d, type = "prob"), tolerance = 1e-10)
})

testthat::test_that("every local pair excludes its P from B learning and its C from all learning", {
  skip_if_no_integration("tailor", "probably")
  d <- data.frame(x = sin(1:160), truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, d) |> recipes::step_normalize(x)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  inner <- function(data) rsample::manual_rset(lapply(split(seq_len(nrow(data)), (seq_len(nrow(data)) - 1L) %/% 40L), function(p)
    rsample::make_splits(list(analysis = setdiff(seq_len(nrow(data)), p), assessment = p), data)), c("one", "two", "three"))
  va <- function(data) rsample::make_splits(list(analysis = 1:120, assessment = 121:160), data)
  for (kind in c("tailor_logistic", "probability", "glm_link")) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = "ensemble", splits = inner)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    spec <- calibrated_workflow(base, "truth", splits = va, postprocessing = post)
    fit <- function(data) workflows::extract_fit_engine(generics::fit(spec, data))
    original <- fit(d)
    one_class_c <- d; one_class_c$truth[121:160] <- "a"
    testthat::expect_true(all(fit(one_class_c)$calibration[[1]]$truth))
    base_model <- function(model) if (kind == "tailor_logistic") model else model$workflow
    for (j in names(original$models[[1]]$models)) {
      p <- original$postprocessing[[1]]$inner_folds[[j]]
      model <- base_model(original$models[[1]]$models[[j]])
      for (role in c("B", "P", "C")) for (variable in c("x", "truth")) {
        changed <- d; rows <- p[[role]]
        if (variable == "x") changed$x[rows] <- changed$x[rows] * 2 + .1 else
          changed$truth[rows] <- ifelse(changed$truth[rows] == "a", "b", "a")
        after <- fit(changed)
        other <- base_model(after$models[[1]]$models[[j]])
        if (role != "B") {
          testthat::expect_identical(coef(workflows::extract_fit_engine(model)), coef(workflows::extract_fit_engine(other)))
          testthat::expect_identical(workflows::extract_recipe(model)$steps, workflows::extract_recipe(other)$steps)
        } else {
          testthat::expect_false(identical(coef(workflows::extract_fit_engine(model)), coef(workflows::extract_fit_engine(other))))
        }
        if (role == "C") testthat::expect_identical(original$postprocessing, after$postprocessing)
        if (role == "P") testthat::expect_false(identical(p$coefficients, after$postprocessing[[1]]$inner_folds[[j]]$coefficients))
      }
    }
  }
})

testthat::test_that("one failed or missing inner pair fails the entire ensemble with context", {
  skip_if_no_integration("tailor", "probably")
  d <- data.frame(x = sin(1:160), truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  original_fit <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  original_predict <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  for (kind in c("tailor_logistic", "probability", "glm_link")) local({
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = "ensemble", v = 3L)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    spec <- calibrated_workflow(base, "truth", seed = 24, postprocessing = post)
    fit <- generics::fit(spec, d)
    ensemble <- workflows::extract_fit_engine(fit)$models[[1]]
    missing <- ensemble; missing$models[[2]] <- NULL
    testthat::expect_error(predict(missing, d, type = "prob"), "Inductive1.*incomplete")
    missing <- ensemble; missing$models[2] <- list(NULL)
    testthat::expect_error(predict(missing, d, type = "prob"), "Inductive1.*incomplete")
    calls <- 0L
    failed_fit <- function(object, data, ...) {
      if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) {
        calls <<- calls + 1L
        if (calls == 2L) stop("injected B failure")
      }
      original_fit(object, data, ...)
    }
    registerS3method("fit", "workflow", failed_fit, envir = asNamespace("generics"))
    on.exit(registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics")), add = TRUE)
    testthat::expect_error(generics::fit(spec, d), "Inductive1 / Fold2.*injected B failure")
    registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics"))
    calls <- 0L
    failed_predict <- function(object, new_data, ...) {
      if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) {
        calls <<- calls + 1L
        if (calls == 2L) stop("injected prediction failure")
      }
      original_predict(object, new_data, ...)
    }
    registerS3method("predict", "workflow", failed_predict, envir = asNamespace("stats"))
    on.exit(registerS3method("predict", "workflow", original_predict, envir = asNamespace("stats")), add = TRUE)
    testthat::expect_error(predict(fit, d, type = "prob"), "Inductive1.*Fold2.*injected prediction failure")
    calls <- 0L
    testthat::expect_error(generics::fit(spec, d), "Inductive1.*Fold2.*injected prediction failure")
  })
})

testthat::test_that("ensemble pairs align named outputs and reject failed postprocessors", {
  skip_if_no_integration("tailor", "probably")
  d <- data.frame(x = sin(1:160), truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  original <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  for (kind in c("tailor_logistic", "probability", "glm_link")) local({
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = "ensemble", v = 3L)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    spec <- calibrated_workflow(base, "truth", seed = 24, postprocessing = post)
    fit <- generics::fit(spec, d)
    expected <- predict(fit, d, type = "prob")
    calls <- 0L; broken <- FALSE
    reordered <- function(object, new_data, type = NULL, ...) {
      p <- original(object, new_data, type = type, ...)
      if (inherits(workflows::extract_spec_parsnip(object), "vennabers")) return(p)
      calls <<- calls + 1L
      if (calls %% 2L == 0L) {
        if (broken) { if (type == "raw") p[1] <- Inf else p[[1]][1] <- NaN }
        if (type == "prob") p <- p[rev(names(p))]
      }
      p
    }
    registerS3method("predict", "workflow", reordered, envir = asNamespace("stats"))
    on.exit(registerS3method("predict", "workflow", original, envir = asNamespace("stats")), add = TRUE)
    testthat::expect_identical(predict(fit, d, type = "prob"), expected)
    broken <- TRUE; calls <- 0L
    testthat::expect_error(predict(fit, d, type = "prob"), "Inductive1.*Fold2")
    registerS3method("predict", "workflow", original, envir = asNamespace("stats"))
    calls <- 0L
    if (kind == "tailor_logistic") {
      estimate <- probably::cal_estimate_logistic
      testthat::local_mocked_bindings(cal_estimate_logistic = function(..., smooth = FALSE) {
        calls <<- calls + 1L
        if (calls == 2L) probably::cal_estimate_none(...) else estimate(..., smooth = smooth)
      }, .package = "probably")
    } else {
      optimizer <- stats::optim
      testthat::local_mocked_bindings(optim = function(...) {
        calls <<- calls + 1L
        result <- optimizer(...)
        if (calls == 2L) result$convergence <- 1L
        result
      }, .package = "stats")
    }
    testthat::expect_error(generics::fit(spec, d), "Inductive1 / Fold2.*P")
  })
})

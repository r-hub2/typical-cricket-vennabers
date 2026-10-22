testthat::test_that("Tailor holdouts reproduce complete public workflows and the VA core", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(seq_len(160)), truth = factor(rep(c("yes", "no", "no", "yes"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  va <- function(data) rsample::make_splits(list(analysis = 1:120, assessment = 121:160), data)
  post <- list(method = "tailor_logistic", strategy = "holdout",
    splits = function(data) rsample::make_splits(list(analysis = 1:80, assessment = 81:120), data))
  fit <- generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d)
  reference <- generics::fit(workflows::add_tailor(base, tailor::tailor() |>
    tailor::adjust_probability_calibration(method = "logistic", smooth = FALSE)),
    d[1:80, ], data_calibration = d[81:120, ])
  held <- predict(reference, d[121:160, ], type = "prob")$.pred_no
  core <- venn_abers(d$truth[121:160] == "no", held)
  expected <- predict(core, predict(reference, d, type = "prob")$.pred_no)
  testthat::expect_equal(predict(fit, d["x"], type = "prob")$.pred_no, expected$probability, tolerance = 1e-10)
  testthat::expect_equal(workflows::extract_fit_engine(fit)$calibration$Inductive1$probability, unname(held), tolerance = 1e-10)
})

testthat::test_that("Tailor cannot hide invalid base probabilities during fit or predict", {
  for (pkg in c("workflows", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:160), truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  spec <- calibrated_workflow(base, "truth", seed = 24,
    postprocessing = list(method = "tailor_logistic", strategy = "holdout", prop = .3))
  fit <- generics::fit(spec, d)
  original <- getS3method("predict", "model_fit", envir = asNamespace("parsnip"))
  corrupt <- function(object, new_data, type = NULL, opts = list(), ...) {
    p <- original(object, new_data, type = type, opts = opts, ...)
    if (identical(type, "prob") && !inherits(object$spec, "vennabers")) p$.pred_b <- p$.pred_b + .2
    p
  }
  testthat::local_mocked_bindings(predict.model_fit = corrupt, .package = "parsnip")
  registerS3method("predict", "model_fit", corrupt, envir = asNamespace("stats"))
  on.exit(registerS3method("predict", "model_fit", original, envir = asNamespace("stats")), add = TRUE)
  testthat::expect_error(generics::fit(spec, d), "Tailor.*probabilit")
  testthat::expect_error(predict(fit, d["x"], type = "prob"), "Tailor.*probabilit")
})

testthat::test_that("only the requested, learned and converged Tailor method is accepted", {
  for (pkg in c("workflows", "tailor", "probably", "rsample")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:160), truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  post <- list(method = "tailor_logistic", strategy = "holdout", prop = .3)
  canonical <- tailor::tailor() |> tailor::adjust_probability_calibration(method = "logistic", smooth = FALSE)
  configured <- workflows::add_tailor(base, canonical)
  fit <- generics::fit(calibrated_workflow(configured, "truth", seed = 24, postprocessing = post), d)
  auto <- generics::fit(calibrated_workflow(base, "truth", seed = 24, postprocessing = post), d)
  testthat::expect_identical(predict(fit, d, type = "prob"), predict(auto, d, type = "prob"))
  testthat::expect_length(workflows::extract_tailor(workflows::extract_fit_engine(fit)$models[[1]])$adjustments, 1)
  testthat::expect_error(calibrated_workflow(configured, "truth"), "postprocessor")
  for (bad in list(tailor::tailor() |> tailor::adjust_probability_threshold(.5),
                  tailor::tailor() |> tailor::adjust_probability_calibration(method = "logistic", smooth = TRUE),
                  canonical |> tailor::adjust_probability_threshold(.5),
                  tailor::tailor() |> tailor::adjust_probability_calibration(method = "logistic", smooth = FALSE, extra = 1))) {
    testthat::expect_error(calibrated_workflow(workflows::add_tailor(base, bad), "truth", postprocessing = post), "[Tt]ailor")
  }
  learned <- workflows::extract_tailor(workflows::extract_fit_engine(fit)$models[[1]])
  testthat::expect_error(calibrated_workflow(workflows::add_tailor(base, learned), "truth", postprocessing = post), "train|fit")
  constant <- d; constant$x <- 1
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", seed = 24, postprocessing = post), constant), "P requires full design rank")
  if (requireNamespace("nnet", quietly = TRUE)) {
    failed_engine <- workflows::update_model(base,
      parsnip::mlp(hidden_units = 2, epochs = 1, mode = "classification") |>
        parsnip::set_engine("nnet", trace = FALSE))
    testthat::expect_error(generics::fit(calibrated_workflow(failed_engine, "truth", seed = 24, postprocessing = post), d),
      "B engine did not converge")
  }
  # Fault injection at the external probably boundary: finite predictions alone
  # must not hide absent calibration, a method switch, or failed convergence.
  original <- probably::cal_estimate_logistic
  for (fault in c("none", "missing", "nonconverged", "warning")) local({
    replacement <- function(...) {
      if (fault == "none") return(probably::cal_estimate_none(...))
      result <- original(...)
      if (fault == "missing") result$estimates <- list()
      if (fault == "nonconverged") result$estimates[[1]]$estimate$converged <- FALSE
      if (fault == "warning") warning("instrumented successful warning")
      result
    }
    testthat::local_mocked_bindings(cal_estimate_logistic = replacement, .package = "probably")
    if (fault == "warning") {
      testthat::expect_warning(valid <- generics::fit(calibrated_workflow(base, "truth", seed = 24, postprocessing = post), d),
        "Inductive1 / Holdout1 tailor_logistic.*instrumented successful warning")
      testthat::expect_match(workflows::extract_fit_engine(valid)$postprocessing[[1]]$warnings, "instrumented successful warning")
    } else {
      testthat::expect_error(suppressWarnings(generics::fit(calibrated_workflow(base, "truth", seed = 24, postprocessing = post), d)),
        "Inductive1 / Holdout1 tailor_logistic")
    }
  })
})

testthat::test_that("both VA methods preserve Tailor predictions, options, plans and RNG", {
  for (pkg in c("workflows", "tailor", "probably", "rsample")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:160), truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  post <- list(method = "tailor_logistic", strategy = "holdout", prop = .3)
  for (method in c("inductive", "cross")) for (event in c("first", "second")) for (precision in list(NULL, 2L)) {
    set.seed(2026); rng <- .Random.seed
    fit <- generics::fit(calibrated_workflow(base, "truth", method = method, seed = 24,
      event_level = event, precision = precision, postprocessing = post), d)
    testthat::expect_identical(.Random.seed, rng)
    state <- workflows::extract_fit_engine(fit)
    held <- truth <- new <- ids <- list()
    for (id in names(state$models)) {
      plan <- state$postprocessing[[id]]
      testthat::expect_setequal(c(plan$B, plan$P), plan$Q)
      testthat::expect_length(intersect(plan$B, plan$P), 0L)
      testthat::expect_length(intersect(plan$Q, plan$C), 0L)
      ref <- generics::fit(workflows::add_tailor(base, tailor::tailor() |>
        tailor::adjust_probability_calibration(method = "logistic", smooth = FALSE)),
        d[plan$B, ], data_calibration = d[plan$P, ])
      column <- paste0(".pred_", state$event_class)
      held[[id]] <- unname(predict(ref, d[plan$C, ], type = "prob")[[column]])
      new[[id]] <- unname(predict(ref, d, type = "prob")[[column]])
      truth[[id]] <- d$truth[plan$C] == state$event_class
      ids[[id]] <- rep(id, length(plan$C))
      testthat::expect_lte(max(abs(held[[id]] - state$calibration[[id]]$probability)), 1e-10)
      testthat::expect_identical(plan$observed_base_fits, 1L)
      testthat::expect_identical(plan$observed_postprocessor_fits, 1L)
    }
    core <- if (method == "cross") cross_venn_abers(unlist(truth), unlist(held), unlist(ids), precision = precision) else
      venn_abers(truth[[1]], held[[1]], precision = precision)
    for (loss in c("log", "brier")) {
      expected <- predict(core, if (method == "cross") do.call(cbind, new) else new[[1]], loss = loss)
      actual <- venn_abers_details(fit, d["x"], loss = loss)
      testthat::expect_lte(max(abs(actual$prediction[[column]] - expected$probability)), 1e-10)
      testthat::expect_identical(actual$prediction, predict(fit, d["x"], type = "prob", loss = loss))
      testthat::expect_identical(as.character(predict(fit, d["x"], type = "class", loss = loss)$.pred_class),
        ifelse(expected$probability >= .5, state$event_class, setdiff(levels(d$truth), state$event_class)))
      for (n in 0:1) {
        small <- venn_abers_details(fit, d[seq_len(n), "x", drop = FALSE], loss = loss)
        testthat::expect_equal(nrow(small$prediction), n)
        testthat::expect_identical(names(small$details), names(actual$details))
      }
    }
    # Replay the actual automatic plans, preserving Q-local coordinates.
    va <- lapply(state$folds, function(f) rsample::make_splits(f[c("analysis", "assessment")], d))
    va <- if (method == "cross") rsample::manual_rset(va, names(va)) else va[[1]]
    inner <- lapply(state$postprocessing, function(p) rsample::make_splits(p[c("analysis", "assessment")], d[p$Q, ]))
    post$splits <- local_calibration_splits(inner, names(d)); post$prop <- NULL
    replay <- generics::fit(calibrated_workflow(base, "truth", method = method, splits = va,
      event_level = event, precision = precision, postprocessing = post), d)
    testthat::expect_identical(predict(fit, d, type = "prob"), predict(replay, d, type = "prob"))
    post <- list(method = "tailor_logistic", strategy = "holdout", prop = .3)
  }
  set.seed(12)
  va <- rsample::initial_split(d, .8, strata = truth, pool = 0)
  q <- d[as.integer(va, "analysis"), ]
  inner <- rsample::initial_split(q, .7, strata = truth, pool = 0)
  expected_rng <- .Random.seed
  set.seed(12)
  fit <- generics::fit(calibrated_workflow(base, "truth", postprocessing = post), d)
  testthat::expect_identical(.Random.seed, expected_rng)
  p <- workflows::extract_fit_engine(fit)$postprocessing$Inductive1
  testthat::expect_identical(p$assessment, as.integer(inner, "assessment"))
})

testthat::test_that("P, C and prediction batches respect complete workflow learning boundaries", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:160), row_id = 1:160, truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, d) |>
    recipes::step_normalize(x)) |> workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  va <- function(data) rsample::make_splits(list(analysis = 1:120, assessment = 121:160), data)
  post <- list(method = "tailor_logistic", strategy = "holdout",
    splits = function(data) rsample::make_splits(list(analysis = 1:80, assessment = 81:120), data))
  spec <- calibrated_workflow(base, "truth", splits = va, postprocessing = post)
  fit <- generics::fit(spec, d)
  state <- workflows::extract_fit_engine(fit)
  changed <- d
  changed$truth[81:120] <- ifelse(changed$truth[81:120] == "a", "b", "a")
  changed_fit <- generics::fit(spec, changed)
  other <- workflows::extract_fit_engine(changed_fit)
  original_model <- state$models[[1]]; changed_model <- other$models[[1]]
  testthat::expect_equal(stats::coef(workflows::extract_fit_engine(original_model)), stats::coef(workflows::extract_fit_engine(changed_model)))
  testthat::expect_identical(workflows::extract_recipe(original_model)$steps, workflows::extract_recipe(changed_model)$steps)
  testthat::expect_false(identical(state$postprocessing[[1]]$coefficients, other$postprocessing[[1]]$coefficients))
  changed <- d; changed$x[81:120] <- changed$x[81:120] * 1000
  other <- workflows::extract_fit_engine(generics::fit(spec, changed))
  testthat::expect_equal(stats::coef(workflows::extract_fit_engine(other$models[[1]])), stats::coef(workflows::extract_fit_engine(original_model)))
  testthat::expect_identical(workflows::extract_recipe(other$models[[1]])$steps, workflows::extract_recipe(original_model)$steps)
  changed <- d; changed$x[121:160] <- 1e3; changed$truth[121:160] <- "a"
  other <- workflows::extract_fit_engine(generics::fit(spec, changed))
  testthat::expect_equal(stats::coef(workflows::extract_fit_engine(other$models[[1]])), stats::coef(workflows::extract_fit_engine(original_model)))
  testthat::expect_identical(other$postprocessing[[1]]$coefficients, state$postprocessing[[1]]$coefficients)
  testthat::expect_identical(workflows::extract_recipe(other$models[[1]])$steps, workflows::extract_recipe(original_model)$steps)
  testthat::expect_true(all(other$calibration[[1]]$truth)) # One-class C remains valid.
  cross <- generics::fit(calibrated_workflow(base, "truth", method = "cross", seed = 24,
    postprocessing = list(method = "tailor_logistic", strategy = "holdout", prop = .3)), d)
  # Count only complete source workflow and Tailor predictions, at public boundaries.
  audit <- new.env(); audit$workflow <- audit$tailor <- 0L
  original_predict <- getS3method("predict", "workflow")
  original_tailor <- getS3method("predict", "tailor")
  original_fit <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  trace("predict.workflow", where = asNamespace("workflows"), tracer = bquote({
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) assign("workflow", get("workflow", .(audit)) + 1L, .(audit))
  }), print = FALSE)
  trace("predict.tailor", where = asNamespace("tailor"), tracer = bquote(assign("tailor", get("tailor", .(audit)) + 1L, .(audit))), print = FALSE)
  trace("fit.workflow", where = asNamespace("workflows"), tracer = quote(stop("Unexpected learning")), print = FALSE)
  registerS3method("predict", "workflow", get("predict.workflow", asNamespace("workflows")), envir = asNamespace("stats"))
  registerS3method("predict", "tailor", get("predict.tailor", asNamespace("tailor")), envir = asNamespace("stats"))
  registerS3method("fit", "workflow", get("fit.workflow", asNamespace("workflows")), envir = asNamespace("generics"))
  on.exit({
    untrace("predict.workflow", where = asNamespace("workflows"))
    untrace("predict.tailor", where = asNamespace("tailor"))
    untrace("fit.workflow", where = asNamespace("workflows"))
    registerS3method("predict", "workflow", original_predict, envir = asNamespace("stats"))
    registerS3method("predict", "tailor", original_tailor, envir = asNamespace("stats"))
    registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics"))
  }, add = TRUE)
  for (fitted in list(fit, cross)) for (n in c(0L, 1L, 160L)) {
    audit$workflow <- audit$tailor <- 0L
    predict(fitted, d[seq_len(n), c("x", "row_id")], type = "prob")
    expected <- if (n) length(workflows::extract_fit_engine(fitted)$models) else 0L
    testthat::expect_equal(audit$workflow, expected)
    testthat::expect_equal(audit$tailor, expected)
  }
})

testthat::test_that("postprocessing rejects invalid requests and all local plans before learning", {
  for (pkg in c("workflows", "tailor", "probably", "rsample")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:80), truth = factor(rep(c("a", "b"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  post <- list(method = "tailor_logistic", strategy = "holdout", prop = .25)
  for (bad in list(list(), c(post, prop = .4), c(post, unknown = 1), post[-1], post[-2],
                  c(post, list(sigmoid_input = NULL)), c(post, list(v = NULL)), c(post, list(splits = NULL)))) {
    testthat::expect_error(calibrated_workflow(base, "truth", postprocessing = bad), "[Pp]ostprocessing")
  }
  for (size in list(NULL, TRUE, NA_real_, Inf, numeric(), c(.1, .2), 0, 1, matrix(.2))) {
    bad <- post; bad["prop"] <- list(size)
    testthat::expect_error(calibrated_workflow(base, "truth", postprocessing = bad), "prop")
  }
  testthat::expect_error(calibrated_workflow(base, "truth", tuning = list(), postprocessing = post), "tuning")
  # A boundary guard ensures schema/plan rejection precedes any recipe or engine fit.
  original_fit <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  trace("fit.workflow", where = asNamespace("workflows"), tracer = quote({
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) stop("UNEXPECTED LEARNING")
  }), print = FALSE)
  registerS3method("fit", "workflow", get("fit.workflow", asNamespace("workflows")), envir = asNamespace("generics"))
  on.exit({
    untrace("fit.workflow", where = asNamespace("workflows"))
    registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics"))
  }, add = TRUE)
  va <- function(data) rsample::make_splits(list(analysis = 1:60, assessment = 61:80), data)
  check <- function(provider, pattern) {
    testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", splits = va,
      postprocessing = list(method = "tailor_logistic", strategy = "holdout", splits = provider)), d), pattern)
  }
  check(function(data) NULL, "B/P plan")
  for (mutate in list(function(x) { x$x[1] <- 99; x }, function(x) x[nrow(x):1, ],
                     function(x) { x$x <- as.integer(x$x); x },
                     function(x) { levels(x$truth) <- rev(levels(x$truth)); x })) {
    check(function(data) rsample::make_splits(list(analysis = 1:40, assessment = 41:60), mutate(data)), "binding")
  }
  for (positions in list(list(analysis = 1:40, assessment = 40:60), list(analysis = 1:40, assessment = 42:60),
                         list(analysis = 1:59, assessment = 60L), list(analysis = 1:58, assessment = 59:60),
                         list(analysis = seq(1, 59, 2), assessment = seq(2, 60, 2)))) {
    check(function(data) rsample::make_splits(positions, data), "B/P plan")
  }
  # A late bad Cross plan also fails before the first source fit.
  provider <- function(data) {
    if (data$x[1] == d$x[1]) return(NULL)
    rsample::initial_split(data, .7, strata = truth, pool = 0)
  }
  testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", method = "cross", seed = 42,
    postprocessing = list(method = "tailor_logistic", strategy = "holdout", splits = provider)), d), "B/P plan")
})

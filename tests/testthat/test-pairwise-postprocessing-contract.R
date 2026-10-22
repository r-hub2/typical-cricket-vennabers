testthat::test_that("pair-local plans reject invalid bindings and class occupancy before learning", {
  for (pkg in c("workflows", "rsample", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:180), truth = factor(rep(c("a", "b", "c"), 60)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  va <- function(data) rsample::make_splits(list(analysis = 1:80, assessment = 81:120), data)
  valid <- function(data) rsample::manual_rset(lapply(list(1:40, 41:80), function(p)
    rsample::make_splits(list(analysis = setdiff(1:80, p), assessment = p), data)), c("one", "two"))
  original <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  guard <- function(object, data, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) stop("UNEXPECTED LEARNING")
    original(object, data, ...)
  }
  registerS3method("fit", "workflow", guard, envir = asNamespace("generics"))
  on.exit(registerS3method("fit", "workflow", original, envir = asNamespace("generics")), add = TRUE)
  check <- function(provider, pattern, strategy = "ensemble", kind = "sigmoid") {
    post <- list(method = kind, strategy = strategy, splits = provider)
    if (kind == "sigmoid") post$sigmoid_input <- "probability"
    testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", classification = "ovo",
      model_training = "pairwise", splits = va, postprocessing = post), d), pattern)
  }
  # Fail in the last class pair after earlier plans have been validated.
  for (mutation in list(function(x) { x$x[1] <- 99; x }, function(x) x[80:1, ],
    function(x) { x$truth <- factor(x$truth, levels = rev(levels(x$truth))); x },
    function(x) x[rev(names(x))])) {
    check(function(data) {
      if (identical(levels(data$truth), c("b", "c"))) data <- mutation(data)
      valid(data)
    }, "pair_3 Inductive1 / ensemble.*binding")
  }
  for (bad in list(function(data) valid(d),
    function(data) list(analysis = 1:40, assessment = 41:80),
    function(data) rsample::vfold_cv(data, v = 2, repeats = 2),
    function(data) { p <- valid(data); p$splits[[2]] <- p$splits[[1]]; p })) {
    check(bad, "pair_1 Inductive1 / ensemble.*B/P plan")
  }
  check(function(data) {
    p <- valid(data)
    b <- which(data$truth == levels(data$truth)[1])
    p$splits[[2]] <- rsample::make_splits(list(analysis = b, assessment = setdiff(1:80, b)), data)
    p$splits[[1]] <- rsample::make_splits(list(analysis = setdiff(1:80, b), assessment = b), data)
    p
  }, "pair_1 Inductive1 / one.*B .*all model classes")
  check(function(data) rsample::make_splits(list(analysis = 3:80, assessment = 1:2), data),
    "pair_1 Inductive1 / Holdout1.*P has 2 rows.*at least 3", "holdout", "tailor_logistic")
  check(function(data) rsample::make_splits(list(analysis = 2:80, assessment = 1L), data),
    "pair_1 Inductive1 / Holdout1.*P .*all model classes", "holdout")
  check(valid, "pair_1 tailor_multinomial B/P plan.*tailor_logistic", kind = "tailor_multinomial")
})

testthat::test_that("pairwise complete paths preserve RNG and actual fit and prediction costs", {
  for (pkg in c("workflows", "rsample", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:240), truth = factor(rep(c("a", "b", "c", "a"), 60)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  original_fit <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  original_predict <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  draws <- numeric(); calls <- post_fits <- 0L
  fit_audit <- function(object, data, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) {
      stopifnot(nlevels(data$truth) == 2L, all(table(data$truth) > 0))
      draws <<- c(draws, runif(1))
    }
    original_fit(object, data, ...)
  }
  predict_audit <- function(object, new_data, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) calls <<- calls + 1L
    original_predict(object, new_data, ...)
  }
  registerS3method("fit", "workflow", fit_audit, envir = asNamespace("generics"))
  registerS3method("predict", "workflow", predict_audit, envir = asNamespace("stats"))
  on.exit({
    registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics"))
    registerS3method("predict", "workflow", original_predict, envir = asNamespace("stats"))
  }, add = TRUE)
  optimizer <- stats::optim; logistic <- probably::cal_estimate_logistic
  testthat::local_mocked_bindings(optim = function(...) {
    post_fits <<- post_fits + 1L; optimizer(...)
  }, .package = "stats")
  testthat::local_mocked_bindings(cal_estimate_logistic = function(...) {
    post_fits <<- post_fits + 1L; logistic(...)
  }, .package = "probably")
  for (method in c("inductive", "cross")) for (strategy in c("holdout", "ensemble"))
    for (kind in c("tailor_logistic", "probability", "glm_link")) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = strategy)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    if (strategy == "holdout") post$prop <- .3 else post$v <- 2L
    args <- list(base, "truth", classification = "ovo", model_training = "pairwise",
      method = method, seed = 24, postprocessing = post)
    if (method == "cross") args$v <- 2L
    # Existing VA pair plans use isolated seeds; B/P plans and fits then consume
    # one enclosing stream in pair, VA-fold, inner-fold order.
    expected <- withr::with_seed(24, {
      qs <- list()
      for (pair in combn(levels(d$truth), 2, simplify = FALSE)) {
        local <- d[d$truth %in% pair, ]; local$truth <- factor(local$truth, levels = pair)
        outer <- withr::with_seed(24, if (method == "cross")
          rsample::vfold_cv(local, v = 2, strata = truth, pool = 0)$splits else
          list(rsample::initial_split(local, .8, strata = truth, pool = 0)))
        qs <- c(qs, lapply(outer, rsample::analysis))
      }
      for (q in qs) if (strategy == "holdout") rsample::initial_split(q, .7, strata = truth, pool = 0) else
        rsample::vfold_cv(q, v = 2, strata = truth, pool = 0)
      runif(length(qs) * if (strategy == "ensemble") 2L else 1L)
    })
    draws <- numeric(); post_fits <- 0L
    set.seed(19); before <- .Random.seed
    fit <- generics::fit(do.call(calibrated_workflow, args), d)
    testthat::expect_identical(.Random.seed, before)
    testthat::expect_identical(draws, expected)
    testthat::expect_equal(post_fits, length(expected))
    state <- workflows::extract_fit_engine(fit)
    cost <- sum(vapply(state$pair_models, function(p)
      sum(vapply(p$postprocessing, `[[`, 1L, "observed_base_fits")), 1L))
    testthat::expect_equal(cost, length(expected))
    for (n in c(0L, 1L, 7L)) {
      calls <- 0L; learned <- draws; processors <- post_fits
      details <- venn_abers_details(fit, d[seq_len(n), "x", drop = FALSE])
      testthat::expect_equal(calls, if (n) cost else 0L)
      testthat::expect_identical(draws, learned)
      testthat::expect_identical(post_fits, processors)
      testthat::expect_equal(nrow(details$prediction), n)
      testthat::expect_true(all(details$details$postprocessing_method == post$method))
    }
    args$seed <- NULL
    # Without a seed, pair VA splits consume the caller stream as before.
    set.seed(24); first <- generics::fit(do.call(calibrated_workflow, args), d)
    rng <- .Random.seed
    set.seed(24); second <- generics::fit(do.call(calibrated_workflow, args), d)
    testthat::expect_identical(.Random.seed, rng)
    testthat::expect_identical(predict(first, d, type = "prob"), predict(second, d, type = "prob"))
  }
})

testthat::test_that("pair-local B P C boundaries survive role and foreign-class mutations", {
  for (pkg in c("workflows", "recipes", "rsample", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:180), truth = factor(rep(c("a", "b", "c"), 60)))
  base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, d) |> recipes::step_normalize(x)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  va <- function(data) rsample::make_splits(list(analysis = 1:80, assessment = 81:120), data)
  inner <- function(data) rsample::manual_rset(lapply(list(1:40, 41:80), function(p)
    rsample::make_splits(list(analysis = setdiff(1:80, p), assessment = p), data)), c("one", "two"))
  for (kind in c("tailor_logistic", "probability", "glm_link")) for (strategy in c("holdout", "ensemble")) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = strategy,
      splits = if (strategy == "ensemble") inner else function(data) inner(data)$splits[[1]])
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    spec <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise", splits = va, postprocessing = post)
    fitted <- generics::fit(spec, d)
    fit <- function(data) workflows::extract_fit_engine(generics::fit(spec, data))$pair_models$pair_3
    original <- workflows::extract_fit_engine(fitted)$pair_models$pair_3
    p <- original$postprocessing[[1]]
    members <- if (strategy == "ensemble") p$inner_folds else list(Holdout1 = p)
    leaf <- function(component, j) {
      model <- component$models[[1]]
      if (strategy == "ensemble") model <- model$models[[j]]
      if (kind != "tailor_logistic") model <- model$workflow
      model
    }
    for (j in names(members)) for (role in c("B", "P", "C")) for (variable in c("x", "truth")) {
      changed <- d; rows <- original$rows[members[[j]][[role]]]
      if (variable == "x") changed$x[rows] <- changed$x[rows] * 2 + .3 else
        changed$truth[rows] <- ifelse(changed$truth[rows] == "b", "c", "b")
      after <- fit(changed)
      before_model <- leaf(original, j); after_model <- leaf(after, j)
      before_coef <- coef(workflows::extract_fit_engine(before_model))
      after_coef <- coef(workflows::extract_fit_engine(after_model))
      if (role == "B") testthat::expect_false(identical(before_coef, after_coef)) else {
        testthat::expect_identical(before_coef, after_coef)
        testthat::expect_identical(workflows::extract_recipe(before_model)$steps, workflows::extract_recipe(after_model)$steps)
      }
      after_p <- after$postprocessing[[1]]
      if (strategy == "ensemble") after_p <- after_p$inner_folds[[j]]
      if (role == "P") testthat::expect_false(identical(members[[j]]$coefficients, after_p$coefficients))
      if (role == "C") testthat::expect_identical(original$postprocessing, after$postprocessing)
    }
    foreign <- d; foreign$x[foreign$truth == "a"] <- 1000
    after <- fit(foreign)
    testthat::expect_identical(original$postprocessing, after$postprocessing)
    testthat::expect_identical(original$calibration, after$calibration)
    for (j in names(members)) testthat::expect_identical(
      coef(workflows::extract_fit_engine(leaf(original, j))), coef(workflows::extract_fit_engine(leaf(after, j))))
  }
})

testthat::test_that("incomplete or failed class pairs never yield a global prediction", {
  for (pkg in c("workflows", "rsample", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:240), truth = factor(rep(c("a", "b", "c", "a"), 60)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  original_fit <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  original_predict <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  on.exit({
    registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics"))
    registerS3method("predict", "workflow", original_predict, envir = asNamespace("stats"))
  }, add = TRUE)
  for (kind in c("tailor_logistic", "probability", "glm_link")) for (strategy in c("holdout", "ensemble")) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = strategy)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    if (strategy == "holdout") post$prop <- .3 else post$v <- 2L
    spec <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
      method = "cross", v = 2, seed = 24, postprocessing = post)
    fit <- generics::fit(spec, d)
    original <- workflows::extract_fit_parsnip(fit)
    bad <- original; bad$fit$pair_models[[2]] <- NULL
    testthat::expect_error(predict(bad, d, type = "prob"), "Pair/class mapping")
    bad <- original; bad$fit$pair_models[[2]]$models[1] <- list(NULL)
    testthat::expect_error(predict(bad, d, type = "prob"), "pair_2 Fold1.*identity")
    bad <- original
    if (strategy == "ensemble") {
      bad$fit$pair_models[[2]]$models[[1]]$models[[2]] <- bad$fit$pair_models[[1]]$models[[1]]$models[[2]]
      testthat::expect_error(predict(bad, d, type = "prob"), "pair_2 Fold1 / Fold2.*identity")
      bad <- original; bad$fit$pair_models[[2]]$models[[1]]$models[[2]] <- NULL
      testthat::expect_error(predict(bad, d, type = "prob"), "pair_2 Fold1.*incomplete")
    } else {
      bad$fit$pair_models[[2]]$models[[1]] <- bad$fit$pair_models[[1]]$models[[1]]
      testthat::expect_error(predict(bad, d, type = "prob"), "pair_2 Fold1.*identity")
    }
    calls <- 0L
    target <- if (strategy == "ensemble") 6L else 3L
    inner <- if (strategy == "ensemble") "Fold2" else "Holdout1"
    failed <- function(object, data, ...) {
      if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) {
        calls <<- calls + 1L
        if (calls == target) stop("injected B failure")
      }
      original_fit(object, data, ...)
    }
    registerS3method("fit", "workflow", failed, envir = asNamespace("generics"))
    testthat::expect_error(generics::fit(spec, d), paste0("pair_2 Fold1 / ", inner, ".*injected B failure"))
    registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics"))
    calls <- 0L; broken <- FALSE
    reordered <- function(object, new_data, type = NULL, ...) {
      p <- original_predict(object, new_data, type = type, ...)
      if (inherits(workflows::extract_spec_parsnip(object), "vennabers")) return(p)
      calls <<- calls + 1L
      if (calls == target && broken) stop("injected prediction failure")
      if (type == "prob") p <- p[rev(names(p))]
      p
    }
    expected <- predict(fit, d, type = "prob")
    registerS3method("predict", "workflow", reordered, envir = asNamespace("stats"))
    testthat::expect_identical(predict(fit, d, type = "prob"), expected)
    calls <- 0L; broken <- TRUE
    testthat::expect_error(predict(fit, d, type = "prob"), paste0("pair_2 Fold1.*", if (strategy == "ensemble") inner else "", ".*injected prediction failure"))
    registerS3method("predict", "workflow", original_predict, envir = asNamespace("stats"))
  }
})

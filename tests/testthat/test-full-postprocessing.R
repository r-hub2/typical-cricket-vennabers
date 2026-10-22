testthat::test_that("full postprocessing preserves held calibration and retains only its new ensemble", {
  skip_if_no_integration("tailor", "probably")
  d <- data.frame(x = sin(1:180), truth = factor(rep(c("a", "b", "c"), 60)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  args <- list(base, "truth", classification = "ovo", model_training = "pairwise",
    method = "cross", v = 5, seed = 42,
    postprocessing = list(method = "tailor_logistic", strategy = "ensemble", v = 3))
  fit_method <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  predict_method <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  base_fits <- post_fits <- predictions <- 0L
  draws <- numeric()
  registerS3method("fit", "workflow", function(object, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) {
      base_fits <<- base_fits + 1L; draws <<- c(draws, runif(1))
    }
    fit_method(object, ...)
  }, envir = asNamespace("generics"))
  registerS3method("predict", "workflow", function(object, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) predictions <<- predictions + 1L
    predict_method(object, ...)
  }, envir = asNamespace("stats"))
  on.exit({
    registerS3method("fit", "workflow", fit_method, envir = asNamespace("generics"))
    registerS3method("predict", "workflow", predict_method, envir = asNamespace("stats"))
  }, add = TRUE)
  logistic <- probably::cal_estimate_logistic
  testthat::local_mocked_bindings(cal_estimate_logistic = function(...) {
    post_fits <<- post_fits + 1L; logistic(...)
  }, .package = "probably")
  matched <- generics::fit(do.call(calibrated_workflow, args), d)
  source_draws <- draws
  base_fits <- post_fits <- 0L; draws <- numeric()
  args$prediction_model <- "full"
  set.seed(71); rng <- .Random.seed
  full <- generics::fit(do.call(calibrated_workflow, args), d)
  testthat::expect_identical(.Random.seed, rng)
  testthat::expect_identical(draws[seq_along(source_draws)], source_draws)
  testthat::expect_equal(base_fits, 54L)
  testthat::expect_equal(post_fits, 54L)
  before <- workflows::extract_fit_engine(matched)
  after <- workflows::extract_fit_engine(full)
  for (id in names(after$pair_models)) {
    left <- before$pair_models[[id]]; right <- after$pair_models[[id]]
    testthat::expect_identical(right$calibration, left$calibration)
    testthat::expect_identical(right$calibrator, left$calibrator)
    testthat::expect_identical(right$folds, left$folds)
    testthat::expect_identical(right$postprocessing[names(left$postprocessing)], left$postprocessing)
    testthat::expect_named(right$models, "Full1")
    testthat::expect_length(right$models$Full1$models, 3L)
    testthat::expect_true(all(right$sources == "Full1"))
    testthat::expect_identical(right$postprocessing$Full1$Q, seq_along(right$rows))
    testthat::expect_length(right$postprocessing$Full1$C, 0L)
  }
  cost <- sum(vapply(after$pair_models, function(p)
    sum(vapply(p$postprocessing, `[[`, 1L, "observed_base_fits")), 1L))
  testthat::expect_equal(cost, 54L)
  testthat::expect_output(print(after), "Base fits: 54")
  testthat::expect_output(print(after), "Models per nonempty prediction: 9")
  for (n in c(0L, 1L, 7L)) {
    predictions <- 0L
    venn_abers_details(full, d[seq_len(n), "x", drop = FALSE])
    testthat::expect_equal(predictions, if (n) 9L else 0L)
    testthat::expect_equal(base_fits, 54L)
    testthat::expect_equal(post_fits, 54L)
  }
  details <- venn_abers_details(full, d[1:7, "x", drop = FALSE])
  testthat::expect_identical(details$prediction, predict(full, d[1:7, ], type = "prob"))
  for (id in names(after$pair_models)) {
    rows <- details$details$pair_id == id
    testthat::expect_true(all(details$details$model_source[rows] == paste0(id, "/Full1")))
    sizes <- vapply(after$pair_models[[id]]$postprocessing$Full1$inner_folds, `[[`, 1L, "B_size")
    testthat::expect_true(all(vapply(details$details$B_size[rows], identical, logical(1), sizes)))
  }
})

testthat::test_that("full domains use their own local inner counts and fail without fallback", {
  skip_if_no_integration("tailor", "probably")
  d <- data.frame(x = sin(1:240), truth = factor(rep(c("a", "b", "c", "a"), 60)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  va <- function(data) {
    n <- nrow(data); h <- seq_len(n %/% 3L)
    rsample::manual_rset(lapply(list(h, setdiff(seq_len(n), h)), function(C)
      rsample::make_splits(list(analysis = setdiff(seq_len(n), C), assessment = C), data)), c("small", "large"))
  }
  inner <- function(data) withr::with_seed(9, rsample::vfold_cv(data,
    v = if (nrow(data) > 100L) 3L else 2L, strata = truth, pool = 0))
  for (kind in c("tailor_logistic", "probability", "glm_link")) for (strategy in c("holdout", "ensemble")) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = strategy,
      splits = if (strategy == "ensemble") inner else function(data) {
        split <- inner(data)$splits[[1]]
        rsample::make_splits(list(analysis = as.integer(split, "analysis"), assessment = as.integer(split, "assessment")), data)
      })
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    spec <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
      method = "cross", splits = va, seed = 24, prediction_model = "full", postprocessing = post)
    fit <- generics::fit(spec, d); state <- workflows::extract_fit_engine(fit)
    counts <- lapply(state$pair_models, function(p) vapply(p$postprocessing, `[[`, 1L, "observed_base_fits"))
    expected <- if (strategy == "ensemble") list(c(3L, 2L, 3L), c(3L, 2L, 3L), c(2L, 2L, 3L)) else rep(list(rep(1L, 3L)), 3L)
    testthat::expect_identical(lapply(unname(counts), unname), expected)
    testthat::expect_output(print(state), paste("Base fits:", sum(unlist(expected))))
    details <- venn_abers_details(fit, d[1:7, "x", drop = FALSE])
    for (id in names(state$pair_models)) {
      p <- state$pair_models[[id]]$postprocessing$Full1
      sizes <- if (strategy == "ensemble") vapply(p$inner_folds, `[[`, 1L, "P_size") else p$P_size
      rows <- details$details$pair_id == id
      testthat::expect_true(all(vapply(as.list(details$details$P_size[rows]), identical, logical(1), sizes)))
    }
    # Public workflow fit boundary: every source succeeds, then the first F
    # engine or learned postprocessor fails. Neither may return source models.
    original <- getS3method("fit", "workflow", envir = asNamespace("generics"))
    seen <- 0L; fail <- sum(vapply(counts, function(x) sum(x[names(x) != "Full1"]), 1L)) + 1L
    registerS3method("fit", "workflow", function(object, ...) {
      if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) {
        seen <<- seen + 1L
        if (seen == fail) stop("final engine failure")
      }
      original(object, ...)
    }, envir = asNamespace("generics"))
    tryCatch(testthat::expect_error(generics::fit(spec, d), "pair_1 Full1 / .*final engine failure"),
      finally = registerS3method("fit", "workflow", original, envir = asNamespace("generics")))
    testthat::expect_equal(seen, fail)
    bad <- workflows::extract_fit_parsnip(fit)
    if (strategy == "ensemble") bad$fit$pair_models[[1]]$models$Full1$models[[2]] <- NULL else
      bad$fit$pair_models[[1]]$models$Full1 <- NULL
    testthat::expect_error(predict(bad, d, type = "prob"), "pair_1.*(incomplete|identity)")
    testthat::expect_error(predict(fit, d, type = "prob", prediction_model = "matched"), "Unknown|must be empty|unused|does not accept")
  }
})

testthat::test_that("C can change F but never its original Q fits, while F B and P keep their roles", {
  skip_if_no_integration("tailor", "probably")
  va <- function(data) rsample::make_splits(list(analysis = 1:108, assessment = 109:144), data)
  inner <- function(data) {
    n <- nrow(data); ends <- c(0L, n %/% 6L, n %/% 2L, n)
    rsample::manual_rset(lapply(1:3, function(j) {
      P <- seq.int(ends[j] + 1L, ends[j + 1L])
      rsample::make_splits(list(analysis = setdiff(seq_len(n), P), assessment = P), data)
    }), c("one", "two", "three"))
  }
  for (kind in c("tailor_logistic", "probability", "glm_link", "tailor_linear", "isotonic_numeric"))
    for (strategy in c("holdout", "ensemble")) {
    numeric <- kind %in% c("tailor_linear", "isotonic_numeric")
    d <- data.frame(x = sin(1:144), truth = if (numeric) 10 + sin(1:144) + cos(1:144) else factor(rep(c("a", "b"), 72)))
    base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x, d) |> recipes::step_normalize(x)) |>
      workflows::add_model(if (numeric) parsnip::linear_reg(engine = "lm") else parsnip::logistic_reg(engine = "glm"))
    post <- list(method = if (kind %in% c("probability", "glm_link")) "sigmoid" else kind, strategy = strategy,
      splits = if (strategy == "ensemble") inner else function(data) inner(data)$splits[[3]])
    if (post$method == "sigmoid") post$sigmoid_input <- kind
    spec <- calibrated_workflow(base, "truth", splits = va, seed = 42, prediction_model = "full", postprocessing = post)
    fit <- function(data) workflows::extract_fit_engine(generics::fit(spec, data))
    original <- fit(d)
    member <- function(s) {
      model <- s$models$Full1
      if (strategy == "ensemble") model <- model$models[[1]]
      model
    }
    workflow <- function(s) {
      model <- member(s)
      if (inherits(model, "workflow")) model else model$workflow
    }
    learned <- function(s) {
      m <- member(s)
      diagnostic <- s$postprocessing$Full1
      if (strategy == "ensemble") diagnostic <- diagnostic$inner_folds[[1]]
      post_state <- if (inherits(m, "vennabers_isotonic_workflow")) m$isotonic else diagnostic$coefficients
      list(engine = coef(workflows::extract_fit_engine(workflow(s))), postprocessor = post_state)
    }
    for (column in c("x", "truth")) {
      changed <- d
      if (column == "x") changed$x[109:144] <- changed$x[109:144] * 2 + .3 else
        if (numeric) changed$truth[109:144] <- changed$truth[109:144] + 3 else
          changed$truth[109:144] <- ifelse(changed$truth[109:144] == "a", "b", "a")
      after <- fit(changed)
      testthat::expect_identical(after$postprocessing$Inductive1, original$postprocessing$Inductive1)
      testthat::expect_false(identical(learned(after), learned(original)), info = paste(kind, strategy, column))
      if (column == "truth") {
        column_name <- if (numeric) "prediction" else "probability"
        testthat::expect_identical(after$calibration[[1]][[column_name]], original$calibration[[1]][[column_name]])
      }
    }
    diagnostic <- original$postprocessing$Full1
    if (strategy == "ensemble") diagnostic <- diagnostic$inner_folds[[1]]
    for (role in c("B", "P")) for (column in c("x", "truth")) {
      changed <- d; rows <- diagnostic[[role]]
      # Not affine: step_normalize, estimated on B, would absorb an affine change of all B rows.
      if (column == "x") changed$x[rows] <- exp(changed$x[rows]) else
        if (numeric) changed$truth[rows] <- changed$truth[rows] + 3 else
          changed$truth[rows] <- ifelse(changed$truth[rows] == "a", "b", "a")
      after <- fit(changed)
      same <- identical(coef(workflows::extract_fit_engine(workflow(after))), coef(workflows::extract_fit_engine(workflow(original))))
      testthat::expect_identical(same, role == "P")
      if (role == "P") testthat::expect_identical(workflows::extract_recipe(workflow(after))$steps,
        workflows::extract_recipe(workflow(original))$steps)
      testthat::expect_false(identical(learned(after), learned(original)), info = paste(kind, strategy, column))
    }
  }
})

testthat::test_that("failed final postprocessors abort after successful source fits", {
  skip_if_no_integration("tailor", "probably")
  d <- data.frame(x = sin(1:120), truth = factor(rep(c("a", "b"), 60)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  logistic <- probably::cal_estimate_logistic; optimizer <- stats::optim
  seen <- 0L
  testthat::local_mocked_bindings(cal_estimate_logistic = function(...) {
    seen <<- seen + 1L
    if (seen == 4L) stop("final postprocessor failure")
    logistic(...)
  }, .package = "probably")
  testthat::local_mocked_bindings(optim = function(...) {
    seen <<- seen + 1L
    if (seen == 4L) stop("final postprocessor failure")
    optimizer(...)
  }, .package = "stats")
  for (kind in c("tailor_logistic", "sigmoid")) {
    seen <- 0L
    post <- list(method = kind, strategy = "ensemble", v = 3)
    if (kind == "sigmoid") post$sigmoid_input <- "glm_link"
    spec <- calibrated_workflow(base, "truth", seed = 42, prediction_model = "full", postprocessing = post)
    if (kind == "tailor_logistic") {
      testthat::expect_warning(testthat::expect_error(generics::fit(spec, d), "Full1 / Fold1.*no learned logistic state"),
        "Full1 / Fold1.*final postprocessor failure")
    } else testthat::expect_error(generics::fit(spec, d), "Full1 / Fold1.*final postprocessor failure")
    testthat::expect_equal(seen, 4L)
  }
})

testthat::test_that("a caller fold named Full1 keeps its source diagnostics", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:180), truth = factor(rep(c("a", "b", "c"), 60)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  va <- function(data) {
    n <- nrow(data); C <- seq_len(n %/% 2L)
    rsample::manual_rset(lapply(list(C, setdiff(seq_len(n), C)), function(rows)
      rsample::make_splits(list(analysis = setdiff(seq_len(n), rows), assessment = rows), data)), c("Full1", "Full1.1"))
  }
  for (pairwise in c(FALSE, TRUE)) {
    local <- if (pairwise) d else droplevels(d[d$truth != "c", ])
    args <- list(base, "truth", method = "cross", splits = va, seed = 42,
      postprocessing = list(method = "sigmoid", sigmoid_input = "glm_link", strategy = "ensemble", v = 2L))
    if (pairwise) { args$classification <- "ovo"; args$model_training <- "pairwise" }
    matched <- generics::fit(do.call(calibrated_workflow, args), local)
    args$prediction_model <- "full"
    full <- generics::fit(do.call(calibrated_workflow, args), local)
    before <- workflows::extract_fit_engine(matched); after <- workflows::extract_fit_engine(full)
    left <- if (pairwise) before$pair_models else list(before)
    right <- if (pairwise) after$pair_models else list(after)
    for (i in seq_along(right)) {
      testthat::expect_identical(right[[i]]$postprocessing[names(left[[i]]$postprocessing)], left[[i]]$postprocessing)
      testthat::expect_length(right[[i]]$postprocessing, 3L)
      testthat::expect_identical(right[[i]]$calibration, left[[i]]$calibration)
      testthat::expect_false(names(right[[i]]$models) %in% names(left[[i]]$postprocessing))
    }
    details <- venn_abers_details(full, local[1:7, "x", drop = FALSE])
    testthat::expect_identical(details$prediction, predict(full, local[1:7, ], type = "prob"))
    for (i in seq_along(right)) {
      part <- right[[i]]
      id <- names(part$models)
      label <- if (pairwise) paste(part$pair_id, id, sep = "/") else id
      testthat::expect_true(label %in% details$details$model_source)
    }
  }
})

testthat::test_that("pairwise and shared multiclass F compositions keep B/P roles and predict from RDS without plans", {
  skip_if_no_integration("tailor", "probably", "nnet")
  cl <- factor(rep(c("a", "b", "c"), 80))
  d <- withr::with_seed(5, data.frame(x = rnorm(240, sd = 2) + as.integer(cl), z = rnorm(240, sd = 2) - as.integer(cl), truth = cl))
  # Positional plans keep every role fixed under mutation; a fresh process sets
  # the flag so any plan call during prediction fails.
  va <- function(data) {
    if (nzchar(Sys.getenv("VENNABERS_PLAN_GUARD"))) stop("plan function called")
    n <- nrow(data)
    rsample::make_splits(list(analysis = seq_len(n - n %/% 4L), assessment = seq.int(n - n %/% 4L + 1L, n)), data)
  }
  inner <- function(data) {
    if (nzchar(Sys.getenv("VENNABERS_PLAN_GUARD"))) stop("plan function called")
    n <- nrow(data); half <- list(seq_len(n %/% 2L), seq.int(n %/% 2L + 1L, n))
    rsample::manual_rset(lapply(half, function(P)
      rsample::make_splits(list(analysis = setdiff(seq_len(n), P), assessment = P), data)), c("first", "second"))
  }
  cases <- list(
    pairwise_tailor = list(training = "pairwise", classification = "ovo", method = "tailor_logistic"),
    pairwise_sigmoid = list(training = "pairwise", classification = "ovo", method = "sigmoid"),
    ovo3 = list(training = "shared", classification = "ovo", method = "tailor_multinomial"),
    ovr3 = list(training = "shared", classification = "ovr", method = "tailor_multinomial"))
  fits <- list(); new <- d[1:7, c("x", "z")]
  for (name in names(cases)) {
    case <- cases[[name]]; pairwise <- case$training == "pairwise"
    engine <- if (pairwise) parsnip::logistic_reg(engine = "glm") else
      parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE, maxit = 1000)
    base <- workflows::workflow() |> workflows::add_recipe(recipes::recipe(truth ~ x + z, d) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |> workflows::add_model(engine)
    post <- list(method = case$method, strategy = "ensemble", splits = inner)
    if (case$method == "sigmoid") post$sigmoid_input <- "probability"
    args <- list(base, "truth", classification = case$classification, splits = va, seed = 42,
      prediction_model = "full", postprocessing = post)
    if (pairwise) args$model_training <- "pairwise"
    spec <- do.call(calibrated_workflow, args)
    fits[[name]] <- generics::fit(spec, d)
    state <- function(data) {
      s <- workflows::extract_fit_engine(generics::fit(spec, data))
      if (pairwise) s$pair_models else list(s)
    }
    original <- state(d)
    for (i in seq_along(original)) {
      part <- original[[i]]; rows <- if (pairwise) part$rows else seq_len(nrow(d))
      member <- part$models$Full1$models[[1L]]
      engine_of <- function(m) coef(workflows::extract_fit_engine(if (inherits(m, "workflow")) m else m$workflow))
      plan <- part$postprocessing$Full1$inner_folds[[1L]]
      for (role in c("B", "P")) for (column in c("x", "truth")) {
        changed <- d; target <- rows[plan[[role]]]
        if (column == "x") changed$x[target] <- changed$x[target] * 2 + .3 else {
          # Swap within the domain's own classes so pair memberships stay intact.
          levels <- part$class_levels[1:2]
          changed$truth[target] <- ifelse(changed$truth[target] == levels[1L], levels[2L],
            ifelse(changed$truth[target] == levels[2L], levels[1L], as.character(changed$truth[target])))
        }
        other <- state(changed)[[i]]$models$Full1$models[[1L]]
        info <- paste(name, i, role, column)
        testthat::expect_identical(identical(engine_of(other), engine_of(member)), role == "P", info = info)
        testthat::expect_false(identical(predict(other, new, type = "prob"), predict(member, new, type = "prob")), info = info)
      }
    }
  }
  saved <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile()
  on.exit(unlink(c(saved, script, log)), add = TRUE)
  saveRDS(list(fits = fits, new = new, expected = lapply(fits, function(fit)
    list(prob = predict(fit, new, type = "prob"), details = venn_abers_details(fit, new)))), saved)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "library(vennabers); loadNamespace('workflows'); Sys.setenv(VENNABERS_PLAN_GUARD = '1')",
    "guard <- function(...) stop('Unexpected learning')",
    "registerS3method('fit', 'workflow', guard, envir=asNamespace('generics'))",
    "registerS3method('fit', 'tailor', guard, envir=asNamespace('generics'))",
    paste0("saved <- readRDS(", deparse(saved), ")"),
    "for (name in names(saved$fits)) {",
    "  fit <- saved$fits[[name]]",
    "  stopifnot(identical(predict(fit, saved$new, type = 'prob'), saved$expected[[name]]$prob))",
    "  stopifnot(identical(venn_abers_details(fit, saved$new), saved$expected[[name]]$details))",
    "}"), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  testthat::expect_equal(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

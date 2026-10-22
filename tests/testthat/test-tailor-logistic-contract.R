# Tailor logistic contract gaps of the 32 fixed compositions: shared two-class
# OvO/OvR keeps multiclass rules, local plans fail early, learned-state
# failures fail late without retry, none or omitted pairs, and data roles hold.
two_class_data <- function() withr::with_seed(123, {
  truth <- factor(rep(c("z", "a"), length.out = 240), levels = c("z", "a"))
  data.frame(x = rnorm(240, sd = 2) + c(-1, 1)[truth], z = rnorm(240, sd = 2) + c(1, -1)[truth],
    truth = truth)
})
two_class_base <- function(d, engine = parsnip::logistic_reg(engine = "glm")) workflows::workflow() |>
  workflows::add_recipe(recipes::recipe(truth ~ ., d) |> recipes::step_normalize(recipes::all_numeric_predictors())) |>
  workflows::add_model(engine)
two_class_va <- function(data) rsample::make_splits(list(analysis = 1:180, assessment = 181:240), data)
two_class_post <- function(strategy) list(method = "tailor_logistic", strategy = strategy,
  splits = if (strategy == "holdout") function(data) rsample::make_splits(list(analysis = 1:120, assessment = 121:180), data) else
    function(data) rsample::manual_rset(lapply(list(1:60, 61:120, 121:180), function(p)
      rsample::make_splits(list(analysis = setdiff(1:180, p), assessment = p), data)), c("one", "two", "three")))

testthat::test_that("shared two-class Tailor logistic keeps multiclass class and tie rules", {
  for (pkg in c("workflows", "tailor", "probably", "rsample")) testthat::skip_if_not_installed(pkg)
  # Varied B/P rows train both learners; every C row and the new row share x = 0,
  # so held scores tie and balanced C truth gives exactly one half. Without an
  # intercept the base probability at x = 0 is exactly 1/2, which keeps the
  # Tailor score of one row equal to that of a batch on every BLAS: OpenBLAS
  # AVX-512 kernels round leftover rows differently from 4-row blocks.
  d <- data.frame(x = c(sin(1:120) * 2, rep(0, 40)), truth = factor(rep(c("z", "a"), 80), levels = c("z", "a")))
  d$x[1:120] <- d$x[1:120] + ifelse(d$truth[1:120] == "a", .5, -.5)
  base <- workflows::workflow() |> workflows::add_variables(outcomes = truth, predictors = x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"), formula = truth ~ 0 + x)
  va <- function(data) rsample::make_splits(list(analysis = 1:120, assessment = 121:160), data)
  for (strategy in c("holdout", "ensemble")) {
    post <- list(method = "tailor_logistic", strategy = strategy, splits = if (strategy == "holdout")
      function(data) rsample::make_splits(list(analysis = 1:80, assessment = 81:120), data) else
      function(data) rsample::manual_rset(lapply(list(1:40, 41:80, 81:120), function(p)
        rsample::make_splits(list(analysis = setdiff(1:120, p), assessment = p), data)), c("one", "two", "three")))
    expected <- c(first = "z", second = "a", ovo = "z", ovr = "z")
    for (path in names(expected)) {
      args <- list(base, "truth", splits = va, postprocessing = post)
      if (path %in% c("first", "second")) args$event_level <- path else args$classification <- path
      fit <- generics::fit(do.call(calibrated_workflow, args), d)
      state <- workflows::extract_fit_engine(fit)
      shared <- path %in% c("ovo", "ovr")
      info <- paste(strategy, path)
      testthat::expect_identical(state$parameters$requested$classification, if (shared) path else "binary", info = info)
      testthat::expect_identical(is.null(state$event_class), shared, info = info)
      for (loss in c("log", "brier")) {
        testthat::expect_identical(unname(unlist(predict(fit, data.frame(x = 0), type = "prob", loss = loss))), c(.5, .5), info = info)
        testthat::expect_identical(predict(fit, data.frame(x = 0), type = "class", loss = loss)$.pred_class,
          factor(expected[[path]], levels = c("z", "a")), info = info)
        details <- venn_abers_details(fit, data.frame(x = 0), loss = loss)$details
        testthat::expect_identical("method" %in% names(details), shared, info = info)
        if (shared) testthat::expect_true(all(details$method == path), info = info)
      }
    }
  }
})

testthat::test_that("shared two-class Tailor logistic plans fail before any learning", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes")) testthat::skip_if_not_installed(pkg)
  d <- two_class_data(); base <- two_class_base(d)
  original <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  registerS3method("fit", "workflow", function(object, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) stop("UNEXPECTED LEARNING")
    original(object, ...)
  }, envir = asNamespace("generics"))
  on.exit(registerS3method("fit", "workflow", original, envir = asNamespace("generics")), add = TRUE)
  bad <- list(
    list(P = 179:180, pattern = "P has 2 rows; requires at least 3 \\(more than design rank 2\\)"),
    list(P = which(d$truth[1:180] == "a"), pattern = "B \\(90 rows\\) must contain all model classes: z, a"),
    list(P = which(d$truth[1:180] == "z")[1:10], pattern = "P \\(10 rows\\) must contain all model classes: z, a"))
  for (classification in c("ovo", "ovr")) for (strategy in c("holdout", "ensemble")) for (case in bad) {
    post <- list(method = "tailor_logistic", strategy = strategy, splits = if (strategy == "holdout")
      function(data) rsample::make_splits(list(analysis = setdiff(1:180, case$P), assessment = case$P), data) else
      function(data) rsample::manual_rset(lapply(list(case$P, setdiff(1:180, case$P)), function(p)
        rsample::make_splits(list(analysis = setdiff(1:180, p), assessment = p), data)), c("bad", "rest")))
    inner <- if (strategy == "holdout") "Holdout1" else "bad"
    testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", classification = classification,
      splits = two_class_va, postprocessing = post), d),
      paste0("Inductive1 / ", inner, " tailor_logistic B/P plan: ", case$pattern))
  }
})

testthat::test_that("shared two-class Tailor logistic refuses failed fits without retry, none or omitted pairs", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes")) testthat::skip_if_not_installed(pkg)
  d <- two_class_data(); base <- two_class_base(d)
  original <- probably::cal_estimate_logistic
  for (classification in c("ovo", "ovr")) for (strategy in c("holdout", "ensemble")) {
    spec <- calibrated_workflow(base, "truth", classification = classification, splits = two_class_va,
      postprocessing = two_class_post(strategy))
    # Fail the last learned P: an ensemble must not skip it or return two pairs.
    target <- if (strategy == "holdout") 1L else 3L
    inner <- if (strategy == "holdout") "Holdout1" else "three"
    for (fault in c("none", "missing", "rank", "nonconverged", "nonfinite")) local({
      calls <- 0L
      testthat::local_mocked_bindings(cal_estimate_logistic = function(..., smooth = FALSE) {
        calls <<- calls + 1L
        if (calls == target && fault == "none") return(probably::cal_estimate_none(...))
        result <- original(..., smooth = smooth)
        if (calls == target) {
          if (fault == "missing") result$estimates <- list()
          if (fault == "rank") result$estimates[[1]]$estimate$rank <- 1L
          if (fault == "nonconverged") result$estimates[[1]]$estimate$converged <- FALSE
          if (fault == "nonfinite") result$estimates[[1]]$estimate$coefficients[] <- NA_real_
        }
        result
      }, .package = "probably")
      # An empty estimate list already breaks Tailor's own fit; the rest reach the P check.
      # A falsified rank also makes Tailor's own prediction warn before that check.
      refused <- function() testthat::expect_error(generics::fit(spec, d), paste0("Inductive1 / ", inner, " tailor_logistic",
        if (fault != "missing") ".*: P "), info = paste(classification, strategy, fault))
      if (fault == "rank") testthat::expect_warning(refused(), "prediction from rank-deficient fit") else refused()
      testthat::expect_identical(calls, target, info = paste(classification, strategy, fault))
    })
  }
  # Real learned-state failures: constant P input, a separable tiny P, an unconverged B engine.
  constant <- d; constant$x <- constant$z <- 1
  formula_base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  # Overlapping B, a ten-row P separated by a narrow gap, then C.
  separable_data <- withr::with_seed(1, data.frame(x = c(rnorm(60), seq(-1, -.025, length.out = 5),
    seq(.025, 1, length.out = 5), rnorm(20)), truth = factor(c(rep(c("z", "a"), 30), rep(c("z", "a"), each = 5),
    rep(c("z", "a"), 10)), levels = c("z", "a"))))
  tiny <- list(method = "tailor_logistic", strategy = "holdout",
    splits = function(data) rsample::make_splits(list(analysis = 1:60, assessment = 61:70), data))
  tiny_va <- function(data) rsample::make_splits(list(analysis = 1:70, assessment = 71:90), data)
  unconverged_base <- two_class_base(d, parsnip::logistic_reg(engine = "glm") |> parsnip::set_engine("glm", control = list(maxit = 1)))
  for (classification in c("ovo", "ovr")) local({
    testthat::expect_error(generics::fit(calibrated_workflow(formula_base, "truth", classification = classification,
      splits = two_class_va, postprocessing = two_class_post("holdout")), constant),
      "Holdout1 tailor_logistic.*P requires full design rank 2")
    calls <- 0L
    testthat::local_mocked_bindings(cal_estimate_logistic = function(...) { calls <<- calls + 1L; original(...) }, .package = "probably")
    # Both glm warnings keep their context; the structural convergence check refuses the fit.
    testthat::expect_warning(testthat::expect_warning(testthat::expect_error(generics::fit(calibrated_workflow(formula_base, "truth",
      classification = classification, splits = tiny_va, postprocessing = tiny), separable_data),
      "Inductive1 / Holdout1 tailor_logistic \\(B=60, P=10\\): P logistic fit did not converge to finite coefficients"),
      "Holdout1 tailor_logistic.*algorithm did not converge"), "Holdout1 tailor_logistic.*fitted probabilities numerically 0 or 1")
    testthat::expect_identical(calls, 1L)
    testthat::expect_warning(testthat::expect_error(generics::fit(calibrated_workflow(unconverged_base, "truth",
      classification = classification, splits = two_class_va, postprocessing = two_class_post("ensemble")), d),
      "Inductive1 / one tailor_logistic.*B engine did not converge"), "Inductive1 / one tailor_logistic.*algorithm did not converge")
  })
})

testthat::test_that("Tailor's own fallback to no calibration is refused on every logistic path", {
  for (pkg in c("workflows", "tailor", "probably", "rsample")) testthat::skip_if_not_installed(pkg)
  # Without a schema, probably cannot build its formula for the class "a b";
  # Tailor then warns and stores no calibration, which must never be used.
  d <- data.frame(x = sin(1:240), truth = factor(rep(c("a b", "c", "c", "a b"), 60), levels = c("a b", "c")))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (path in c("binary", "ovo", "ovr", "pairwise")) for (strategy in c("holdout", "ensemble")) {
    data <- d
    if (path == "pairwise") data$truth <- factor(rep(c("a b", "c", "e", "a b"), 60), levels = c("a b", "c", "e"))
    post <- list(method = "tailor_logistic", strategy = strategy)
    if (strategy == "holdout") post$prop <- .3 else post$v <- 2L
    args <- list(base, "truth", seed = 24, postprocessing = post)
    if (path != "binary") args$classification <- if (path == "pairwise") "ovo" else path
    if (path == "pairwise") args$model_training <- "pairwise"
    context <- paste0(if (path == "pairwise") "pair_1 ", "Inductive1 / ", if (strategy == "holdout") "Holdout1" else "Fold1",
      " tailor_logistic")
    testthat::expect_warning(testthat::expect_error(generics::fit(do.call(calibrated_workflow, args), data),
      paste0(context, ".*P has no learned logistic state")), paste0(context, ".*No calibration is applied"))
  }
})

testthat::test_that("shared two-class Tailor logistic B, P and C mutations keep their learning roles", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes")) testthat::skip_if_not_installed(pkg)
  d <- two_class_data(); base <- two_class_base(d)
  for (classification in c("ovo", "ovr")) for (strategy in c("holdout", "ensemble")) {
    spec <- calibrated_workflow(base, "truth", classification = classification, splits = two_class_va,
      postprocessing = two_class_post(strategy))
    state <- function(data) workflows::extract_fit_engine(generics::fit(spec, data))
    original <- state(d)
    record <- original$postprocessing$Inductive1
    members <- if (strategy == "holdout") list(Holdout1 = record) else record$inner_folds
    leaf <- function(s, j) if (strategy == "holdout") s$models$Inductive1 else s$models$Inductive1$models[[j]]
    pair <- function(s, j) if (strategy == "holdout") s$postprocessing$Inductive1 else s$postprocessing$Inductive1$inner_folds[[j]]
    for (j in names(members)) for (role in c("B", "P", "C")) for (column in c("x", "truth")) {
      changed <- d; rows <- members[[j]][[role]]
      if (column == "x") changed$x[rows] <- changed$x[rows] * 2 + .3 else
        changed$truth[rows] <- ifelse(changed$truth[rows] == "z", "a", "z")
      after <- state(changed)
      info <- paste(classification, strategy, j, role, column)
      same_engine <- identical(coef(workflows::extract_fit_engine(leaf(original, j))), coef(workflows::extract_fit_engine(leaf(after, j))))
      testthat::expect_identical(same_engine, role != "B", info = info)
      if (role != "B") testthat::expect_identical(workflows::extract_recipe(leaf(original, j))$steps,
        workflows::extract_recipe(leaf(after, j))$steps, info = info)
      testthat::expect_identical(identical(pair(original, j)$coefficients, pair(after, j)$coefficients), role == "C", info = info)
      if (role == "C") {
        testthat::expect_identical(after$postprocessing, original$postprocessing, info = info)
        testthat::expect_false(identical(after$calibration, original$calibration), info = info)
      }
    }
  }
})

testthat::test_that("shared two-class full compositions keep F roles and the source C values", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes")) testthat::skip_if_not_installed(pkg)
  d <- two_class_data(); base <- two_class_base(d)
  # Positional plans relative to each local Q serve the source Q and F alike.
  relative <- function(strategy) list(method = "tailor_logistic", strategy = strategy, splits = function(data) {
    n <- nrow(data); thirds <- split(seq_len(n), rep(1:3, each = n / 3))
    if (strategy == "holdout") return(rsample::make_splits(list(analysis = unlist(thirds[1:2]), assessment = thirds[[3]]), data))
    rsample::manual_rset(lapply(thirds, function(p) rsample::make_splits(list(analysis = setdiff(seq_len(n), p),
      assessment = p), data)), c("one", "two", "three"))
  })
  for (classification in c("ovo", "ovr")) for (strategy in c("holdout", "ensemble")) {
    spec <- calibrated_workflow(base, "truth", classification = classification, splits = two_class_va,
      prediction_model = "full", postprocessing = relative(strategy))
    state <- function(data) workflows::extract_fit_engine(generics::fit(spec, data))
    original <- state(d)
    final <- original$postprocessing$Full1
    member <- function(s) if (strategy == "holdout") s$models$Full1 else s$models$Full1$models[[1L]]
    plan <- function(s) if (strategy == "holdout") s$postprocessing$Full1 else s$postprocessing$Full1$inner_folds[[1L]]
    engine <- function(s) coef(workflows::extract_fit_engine(member(s)))
    testthat::expect_named(original$models, "Full1")
    testthat::expect_identical(plan(original)$Q, seq_len(nrow(d)))
    testthat::expect_length(plan(original)$C, 0L)
    # F's own B trains the engine, F's own P only the Tailor.
    for (role in c("B", "P")) for (column in c("x", "truth")) {
      changed <- d; rows <- plan(original)[[role]]
      if (column == "x") changed$x[rows] <- changed$x[rows] * 2 + .3 else
        changed$truth[rows] <- ifelse(changed$truth[rows] == "z", "a", "z")
      after <- state(changed); info <- paste(classification, strategy, role, column)
      testthat::expect_identical(identical(engine(after), engine(original)), role == "P", info = info)
      if (role == "P") testthat::expect_identical(workflows::extract_recipe(member(after))$steps,
        workflows::extract_recipe(member(original))$steps, info = info)
      testthat::expect_false(identical(plan(after)$coefficients, plan(original)$coefficients), info = info)
    }
    # Source C belongs to F: it changes the final composition, never the source fit or held values.
    for (column in c("x", "truth")) {
      changed <- d
      if (column == "x") changed$x[181:240] <- changed$x[181:240] * 2 + .3 else
        changed$truth[181:240] <- ifelse(changed$truth[181:240] == "z", "a", "z")
      after <- state(changed); info <- paste(classification, strategy, "C", column)
      testthat::expect_identical(after$postprocessing$Inductive1, original$postprocessing$Inductive1, info = info)
      testthat::expect_false(identical(list(engine(after), plan(after)$coefficients),
        list(engine(original), plan(original)$coefficients)), info = info)
      if (column == "truth") testthat::expect_identical(after$calibration$Inductive1$probability,
        original$calibration$Inductive1$probability, info = info)
    }
  }
})

testthat::test_that("pairwise Tailor logistic refuses a failed class-pair postprocessor without a global result", {
  for (pkg in c("workflows", "tailor", "probably", "rsample")) testthat::skip_if_not_installed(pkg)
  d <- data.frame(x = sin(1:240), truth = factor(rep(c("a", "b", "c", "a"), 60)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  original <- probably::cal_estimate_logistic
  for (strategy in c("holdout", "ensemble")) for (fault in c("none", "nonconverged")) local({
    post <- list(method = "tailor_logistic", strategy = strategy)
    if (strategy == "holdout") post$prop <- .3 else post$v <- 2L
    spec <- calibrated_workflow(base, "truth", classification = "ovo", model_training = "pairwise",
      method = "cross", v = 2, seed = 24, postprocessing = post)
    # Class pairs are fitted in pair/VA-fold/inner-fold order; fail pair_2 Fold2.
    target <- if (strategy == "holdout") 4L else 8L
    inner <- if (strategy == "holdout") "Holdout1" else "Fold2"
    calls <- 0L
    testthat::local_mocked_bindings(cal_estimate_logistic = function(..., smooth = FALSE) {
      calls <<- calls + 1L
      if (calls == target && fault == "none") return(probably::cal_estimate_none(...))
      result <- original(..., smooth = smooth)
      if (calls == target) result$estimates[[1]]$estimate$converged <- FALSE
      result
    }, .package = "probably")
    testthat::expect_error(generics::fit(spec, d), paste0("pair_2 Fold2 / ", inner, " tailor_logistic.*: P "),
      info = paste(strategy, fault))
    testthat::expect_identical(calls, target, info = paste(strategy, fault))
  })
})

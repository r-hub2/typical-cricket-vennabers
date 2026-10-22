testthat::test_that("binary sigmoid holdouts use the complete recipe and both input types", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:160), truth = factor(rep(c("a", "b", "b", "a"), 40)))
  base <- workflows::workflow() |> workflows::add_recipe(
    recipes::recipe(truth ~ x, d) |> recipes::step_normalize(x)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  va <- function(data) rsample::make_splits(list(analysis = 1:120, assessment = 121:160), data)
  for (input in c("probability", "glm_link")) {
    post <- list(method = "sigmoid", strategy = "holdout", sigmoid_input = input,
      splits = function(data) rsample::make_splits(list(analysis = 1:80, assessment = 81:120), data))
    fit <- generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d)
    state <- workflows::extract_fit_engine(fit)
    reference <- generics::fit(base, d[1:80, ])
    score <- function(data) if (input == "probability") predict(reference, data, type = "prob")$.pred_b else
      as.double(predict(reference, data, type = "raw", opts = list(type = "link")))
    y <- d$truth[81:120] == "b"
    target <- ifelse(y, (sum(y) + 1) / (sum(y) + 2), 1 / (sum(!y) + 2))
    mapping <- glm(target ~ score(d[81:120, ]), family = quasibinomial(), control = glm.control(epsilon = 1e-12))
    mapped <- function(data) plogis(coef(mapping)[1] + coef(mapping)[2] * score(data))
    model <- state$models[[1]]
    testthat::expect_lte(max(abs(predict(model, d, type = "prob")$.pred_b - mapped(d))), 1e-7)
    held <- predict(model, d[121:160, ], type = "prob")$.pred_a
    core <- venn_abers(d$truth[121:160] == "a", held)
    expected <- predict(core, predict(model, d, type = "prob")$.pred_a)$probability
    testthat::expect_lte(max(abs(predict(fit, d["x"], type = "prob")$.pred_a - expected)), 1e-10)
    testthat::expect_identical(state$postprocessing[[1]]$input, input)
  }
})

testthat::test_that("identical sigmoid inputs match sklearn, including rare and constant P", {
  skip_if_no_integration()
  reference <- read.csv(testthat::test_path("extdata", "sigmoid-reference.csv"))
  original <- getS3method("predict", "workflow", envir = asNamespace("workflows"))
  # Fix the public base prediction boundary to the exact Python reference inputs.
  controlled <- function(object, new_data, type = NULL, ...) {
    if (inherits(workflows::extract_spec_parsnip(object), "vennabers")) return(original(object, new_data, type, ...))
    if (type == "raw") return(new_data$score)
    data.frame(.pred_b = new_data$score, .pred_a = 1 - new_data$score)
  }
  registerS3method("predict", "workflow", controlled, envir = asNamespace("stats"))
  on.exit(registerS3method("predict", "workflow", original, envir = asNamespace("stats")), add = TRUE)
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (input in unique(reference$input)) for (case in unique(reference$case)) {
    ref <- reference[reference$input == input & reference$case == case, ]
    p <- ref[ref$role == "P", ]; new <- ref[ref$role == "prediction", ]
    d <- data.frame(x = sin(seq_len(40 + nrow(p))), score = c(rep(.5, 40), p$score),
      truth = factor(c(rep(c("a", "b", "b", "a"), 10), ifelse(p$truth, "b", "a"))))
    va <- function(data) rsample::make_splits(list(analysis = c(1:20, 41:nrow(data)), assessment = 21:40), data)
    post <- list(method = "sigmoid", strategy = "holdout", sigmoid_input = input,
      splits = function(data) rsample::make_splits(list(analysis = 1:20, assessment = 21:nrow(data)), data))
    fit <- generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d)
    model <- workflows::extract_fit_engine(fit)$models[[1]]
    for (r in list(p, new)) {
      actual <- predict(model, data.frame(x = 0, score = r$score), type = "prob")$.pred_b
      testthat::expect_lte(max(abs(actual - r$probability)), 1e-7, label = paste(input, case))
    }
  }
})

testthat::test_that("sigmoid accepts a stationary fit whose line search stops at machine precision", {
  skip_if_no_integration()
  # Pairwise P inputs from the v0.9 tuning acceptance: with sklearn's ftol,
  # L-BFGS-B ends in ABNORMAL_TERMINATION_IN_LNSRCH at the optimum.
  p <- read.csv(testthat::test_path("extdata", "sigmoid-line-search.csv"))
  original <- getS3method("predict", "workflow", envir = asNamespace("workflows"))
  controlled <- function(object, new_data, type = NULL, ...) {
    if (inherits(workflows::extract_spec_parsnip(object), "vennabers")) return(original(object, new_data, type, ...))
    data.frame(.pred_b = new_data$score, .pred_a = 1 - new_data$score)
  }
  registerS3method("predict", "workflow", controlled, envir = asNamespace("stats"))
  on.exit(registerS3method("predict", "workflow", original, envir = asNamespace("stats")), add = TRUE)
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  d <- data.frame(x = sin(seq_len(40 + nrow(p))), score = c(rep(.5, 40), p$score),
    truth = factor(c(rep(c("a", "b", "b", "a"), 10), ifelse(p$truth == 1L, "b", "a"))))
  va <- function(data) rsample::make_splits(list(analysis = c(1:20, 41:nrow(data)), assessment = 21:40), data)
  post <- list(method = "sigmoid", strategy = "holdout", sigmoid_input = "probability",
    splits = function(data) rsample::make_splits(list(analysis = 1:20, assessment = 21:nrow(data)), data))
  fit <- generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d)
  model <- workflows::extract_fit_engine(fit)$models[[1]]
  # Exact optimum of the same smoothed-target logistic objective.
  n1 <- sum(p$truth); n0 <- nrow(p) - n1
  target <- ifelse(p$truth == 1L, (n1 + 1) / (n1 + 2), 1 / (n0 + 2))
  reference <- suppressWarnings(stats::glm(target ~ p$score, family = stats::quasibinomial(),
    control = stats::glm.control(epsilon = 1e-14, maxit = 100)))
  grid <- seq(0, 1, by = .05)
  expected <- stats::plogis(coef(reference)[[1L]] + coef(reference)[[2L]] * grid)
  actual <- predict(model, data.frame(x = 0, score = grid), type = "prob")$.pred_b
  testthat::expect_lte(max(abs(actual - expected)), 1e-7)
})

testthat::test_that("a line-search stop is kept only at a stationary point; other stops fail", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:80), truth = factor(rep(c("a", "b", "b", "a"), 20)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  original_optim <- stats::optim
  # The optimizer's own result, reported as ABNORMAL_TERMINATION_IN_LNSRCH (52).
  fit_after_stop <- function(spec, shift) testthat::with_mocked_bindings(
    workflows::extract_fit_engine(generics::fit(spec, d))$postprocessing[[1]],
    optim = function(...) {
      result <- original_optim(...)
      result$par <- result$par + shift
      result$convergence <- 52L
      result
    }, .package = "stats")
  for (input in c("probability", "glm_link")) {
    spec <- calibrated_workflow(base, "truth", seed = 24,
      postprocessing = list(method = "sigmoid", strategy = "holdout", prop = .3, sigmoid_input = input))
    converged <- workflows::extract_fit_engine(generics::fit(spec, d))$postprocessing[[1]]
    kept <- fit_after_stop(spec, c(0, 0))
    testthat::expect_identical(kept$coefficients, converged$coefficients)
    testthat::expect_true(kept$converged)
    # Away from the optimum, or at a nonfinite point, the same stop is nonconvergence.
    for (shift in list(c(.05, 0), c(0, -.05), c(NaN, 0))) {
      testthat::expect_error(fit_after_stop(spec, shift), "Inductive1 / Holdout1 sigmoid.*P fit did not converge",
        info = paste(input, toString(shift)))
    }
  }
})

testthat::test_that("sigmoid rejects incompatible workflows and invalid local learning inputs", {
  skip_if_no_integration("tailor")
  d <- data.frame(x = sin(1:80), truth = factor(rep(c("a", "b"), 40)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  post <- list(method = "sigmoid", strategy = "holdout", prop = .3, sigmoid_input = "glm_link")
  for (input in list(NULL, NA_character_, "link", c("probability", "glm_link"))) {
    bad <- post; bad["sigmoid_input"] <- list(input)
    testthat::expect_error(calibrated_workflow(base, "truth", postprocessing = bad), "sigmoid_input")
  }
  configured <- workflows::add_tailor(base, tailor::tailor() |> tailor::adjust_probability_threshold(.5))
  testthat::expect_error(calibrated_workflow(configured, "truth", postprocessing = post), "sigmoid B.*without a postprocessor")
  for (family in list(binomial("probit"), quasibinomial(), "gaussian")) {
    bad <- workflows::update_model(base, parsnip::logistic_reg() |> parsnip::set_engine("glm", family = family))
    testthat::expect_error(calibrated_workflow(bad, "truth", postprocessing = post), "sigmoid glm_link B.*binomial")
  }
  bad <- workflows::update_model(base, parsnip::logistic_reg(engine = "glmnet"))
  testthat::expect_error(calibrated_workflow(bad, "truth", postprocessing = post), "sigmoid glm_link B")
  spec <- calibrated_workflow(base, "truth", seed = 24, postprocessing = post)
  fitted <- generics::fit(spec, d)
  original <- getS3method("predict", "workflow", envir = asNamespace("workflows"))
  for (input in c("glm_link", "probability")) for (fault in c("nonfinite", "missing")) local({
    corrupt <- function(object, new_data, type = NULL, ...) {
      p <- original(object, new_data, type = type, ...)
      if (inherits(workflows::extract_spec_parsnip(object), "vennabers")) return(p)
      if (fault == "missing") return(if (type == "raw") p[-1L] else p[-1L, ])
      if (type == "raw") p[1] <- Inf else p[[1]][1] <- NaN
      p
    }
    registerS3method("predict", "workflow", corrupt, envir = asNamespace("stats"))
    on.exit(registerS3method("predict", "workflow", original, envir = asNamespace("stats")), add = TRUE)
    requested <- post; requested$sigmoid_input <- input
    spec <- calibrated_workflow(base, "truth", seed = 24, postprocessing = requested)
    testthat::expect_error(generics::fit(spec, d), "Inductive1 / Holdout1 sigmoid.*P")
    if (input == "glm_link") testthat::expect_error(predict(fitted, d, type = "prob"), "Inductive1 / Holdout1 sigmoid.*prediction")
  })
  # A failed optimizer cannot be accepted merely because it returns coefficients.
  original_optim <- stats::optim
  testthat::local_mocked_bindings(optim = function(...) {
    result <- original_optim(...); result$convergence <- 1L; result
  }, .package = "stats")
  testthat::expect_error(generics::fit(spec, d), "Inductive1 / Holdout1 sigmoid.*P fit did not converge")
})

testthat::test_that("both events and VA methods preserve sigmoid maps, options and data roles", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:160), truth = factor(rep(c("b", "a", "a", "b"), 40), levels = c("b", "a")))
  base <- workflows::workflow() |> workflows::add_recipe(
    recipes::recipe(truth ~ x, d) |> recipes::step_normalize(x)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  for (input in c("probability", "glm_link")) for (method in c("inductive", "cross"))
    for (event in c("first", "second")) for (precision in list(NULL, 2L)) {
    post <- list(method = "sigmoid", strategy = "holdout", prop = .3, sigmoid_input = input)
    spec <- calibrated_workflow(base, "truth", method = method, event_level = event,
      precision = precision, seed = 24, postprocessing = post)
    set.seed(102); rng <- .Random.seed
    fitted <- generics::fit(spec, d)
    testthat::expect_identical(.Random.seed, rng)
    state <- workflows::extract_fit_engine(fitted)
    held <- truth <- new <- ids <- list()
    column <- paste0(".pred_", state$event_class)
    for (id in names(state$models)) {
      model <- state$models[[id]]; plan <- state$postprocessing[[id]]
      held[[id]] <- predict(model, d[plan$C, ], type = "prob")[[column]]
      new[[id]] <- predict(model, d, type = "prob")[[column]]
      truth[[id]] <- d$truth[plan$C] == state$event_class
      ids[[id]] <- rep(id, length(plan$C))
      testthat::expect_lte(max(abs(held[[id]] - state$calibration[[id]]$probability)), 1e-10)
      testthat::expect_length(intersect(plan$B, plan$P), 0)
      testthat::expect_length(intersect(plan$Q, plan$C), 0)
      testthat::expect_setequal(c(plan$B, plan$P), plan$Q)
      testthat::expect_identical(plan$sigmoid_class, levels(d$truth)[2])
    }
    core <- if (method == "cross") cross_venn_abers(unlist(truth), unlist(held), unlist(ids), precision = precision) else
      venn_abers(truth[[1]], held[[1]], precision = precision)
    for (loss in c("log", "brier")) {
      expected <- predict(core, if (method == "cross") do.call(cbind, new) else new[[1]], loss = loss)
      details <- venn_abers_details(fitted, d["x"], loss = loss)
      testthat::expect_lte(max(abs(details$prediction[[column]] - expected$probability)), 1e-10)
      testthat::expect_identical(details$prediction, predict(fitted, d["x"], type = "prob", loss = loss))
      for (n in 0:1) {
        small <- venn_abers_details(fitted, d[seq_len(n), "x", drop = FALSE], loss = loss)
        testthat::expect_equal(nrow(small$prediction), n)
        testthat::expect_identical(names(small$details), names(details$details))
      }
    }
  }
  # Freeze plans so mutations change only the intended learning role.
  va <- function(data) rsample::make_splits(list(analysis = 1:120, assessment = 121:160), data)
  for (input in c("probability", "glm_link")) {
    post <- list(method = "sigmoid", strategy = "holdout", sigmoid_input = input,
      splits = function(data) rsample::make_splits(list(analysis = 1:80, assessment = 81:120), data))
    spec <- calibrated_workflow(base, "truth", splits = va, postprocessing = post)
    state <- workflows::extract_fit_engine(generics::fit(spec, d))
    for (role in c("P", "C")) {
      changed <- d; rows <- state$postprocessing[[1]][[role]]
      changed$x[rows] <- changed$x[rows] * 1000
      changed$truth[rows] <- if (role == "C") "a" else ifelse(changed$truth[rows] == "a", "b", "a")
      after <- workflows::extract_fit_engine(generics::fit(spec, changed))
      testthat::expect_equal(coef(workflows::extract_fit_engine(after$models[[1]]$workflow)),
        coef(workflows::extract_fit_engine(state$models[[1]]$workflow)))
      testthat::expect_identical(workflows::extract_recipe(after$models[[1]]$workflow)$steps,
        workflows::extract_recipe(state$models[[1]]$workflow)$steps)
      if (role == "C") testthat::expect_identical(after$postprocessing[[1]]$coefficients, state$postprocessing[[1]]$coefficients)
      else testthat::expect_false(identical(after$postprocessing[[1]]$coefficients, state$postprocessing[[1]]$coefficients))
    }
  }
})

testthat::test_that("sigmoid learns once on B/P and predicts once per nonempty source batch", {
  skip_if_no_integration()
  d <- data.frame(x = sin(1:80), truth = factor(rep(c("a", "b", "b", "a"), 20)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  original_predict <- getS3method("predict", "workflow", envir = asNamespace("workflows"))
  original_fit <- getS3method("fit", "workflow", envir = asNamespace("workflows"))
  audit <- new.env(); audit$fits <- audit$predictions <- list()
  audit$learning <- TRUE
  fitted_boundary <- function(object, data, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) {
      if (!audit$learning) stop("Unexpected learning")
      audit$fits[[length(audit$fits) + 1L]] <- data$x
    }
    original_fit(object, data, ...)
  }
  prediction_boundary <- function(object, new_data, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) {
      audit$predictions[[length(audit$predictions) + 1L]] <- new_data$x
    }
    original_predict(object, new_data, ...)
  }
  registerS3method("fit", "workflow", fitted_boundary, envir = asNamespace("generics"))
  registerS3method("predict", "workflow", prediction_boundary, envir = asNamespace("stats"))
  on.exit({
    registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics"))
    registerS3method("predict", "workflow", original_predict, envir = asNamespace("stats"))
  }, add = TRUE)
  for (input in c("probability", "glm_link")) for (method in c("inductive", "cross")) {
    audit$learning <- TRUE; audit$fits <- audit$predictions <- list()
    fit <- generics::fit(calibrated_workflow(base, "truth", method = method, seed = 24,
      postprocessing = list(method = "sigmoid", strategy = "holdout", prop = .3, sigmoid_input = input)), d)
    plans <- workflows::extract_fit_engine(fit)$postprocessing
    testthat::expect_equal(length(audit$fits), length(plans))
    testthat::expect_equal(length(audit$predictions), 2 * length(plans))
    for (i in seq_along(plans)) {
      testthat::expect_identical(audit$fits[[i]], d$x[plans[[i]]$B])
      testthat::expect_identical(audit$predictions[[2*i-1]], d$x[plans[[i]]$P])
      testthat::expect_identical(audit$predictions[[2*i]], d$x[plans[[i]]$C])
    }
    audit$learning <- FALSE
    for (n in c(0L, 1L, nrow(d))) {
      audit$predictions <- list()
      predict(fit, d[seq_len(n), "x", drop = FALSE], type = "prob")
      testthat::expect_equal(length(audit$predictions), if (n) length(plans) else 0L)
    }
  }
  audit$learning <- TRUE
  for (positions in list(list(analysis = 1:59, assessment = 60L),
                         list(analysis = c(1:57, 60L), assessment = 58:59),
                         list(analysis = 1:40, assessment = 40:60))) {
    audit$fits <- list()
    va <- function(data) rsample::make_splits(list(analysis = 1:60, assessment = 61:80), data)
    post <- list(method = "sigmoid", strategy = "holdout", sigmoid_input = "probability",
      splits = function(data) rsample::make_splits(positions, data))
    testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", splits = va, postprocessing = post), d),
      "Inductive1 / Holdout1 sigmoid B/P plan")
    testthat::expect_length(audit$fits, 0L)
  }
})

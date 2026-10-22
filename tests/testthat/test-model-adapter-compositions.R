# Adapters in the existing class-pair, postprocessing and search compositions:
# each case matches the native engine, and the callbacks see only the rows of
# their local role.
adapter_factories <- function() {
  example <- new.env(parent = globalenv())
  utils::capture.output(sys.source(system.file("examples", "model-adapter.R", package = "vennabers"), example))
  example
}

# Records every callback fit in `audit$fits`; otherwise the example's GLM adapter.
audited_glm <- function(audit) model_adapter("classification",
  fit = function(x, y) {
    audit$fits[[length(audit$fits) + 1L]] <- list(x = x, y = y)
    list(model = stats::glm(y ~ ., data = cbind(x, y = y), family = stats::binomial()), classes = levels(y))
  },
  predict = function(object, new_data) {
    p <- stats::predict(object$model, new_data, type = "response")
    matrix(c(1 - p, p), ncol = 2L, dimnames = list(NULL, object$classes))
  })

adapter_multiclass <- function() {
  classes <- factor(rep(c("a", "b", "c"), 70))
  d <- data.frame(x = 1.5 * sin(1:210) + as.integer(classes), z = cos(1:210), truth = classes)
  list(train = d[1:180, ], new = d[181:210, c("x", "z")])
}

testthat::test_that("a binary callback carries pairwise OvO on each pair's rows and original levels", {
  skip_if_no_integration()
  e <- adapter_factories()
  d <- adapter_multiclass()
  native <- e$base(parsnip::logistic_reg(engine = "glm"), d$train)
  for (method in c("inductive", "cross")) for (mode in c("matched", "full")) {
    audit <- new.env()
    audit$fits <- list()
    args <- list(outcome = "truth", classification = "ovo", model_training = "pairwise",
                 method = method, prediction_model = mode, seed = 42)
    fit <- generics::fit(do.call(calibrated_workflow, c(list(e$base(audited_glm(audit), d$train)), args)), d$train)
    expected <- generics::fit(do.call(calibrated_workflow, c(list(native), args)), d$train)
    info <- paste(method, mode)
    testthat::expect_equal(venn_abers_details(fit, d$new), venn_abers_details(expected, d$new), info = info)
    testthat::expect_identical(predict(fit, d$new, type = "class"), predict(expected, d$new, type = "class"), info = info)
    # Every pair's source folds first, then each pair's Full1 on all its rows.
    state <- workflows::extract_fit_engine(fit)
    pairs <- unname(state$pair_models)
    sources <- if (method == "cross") 5L else 1L
    full <- if (mode == "full") pairs
    expected_rows <- c(unlist(lapply(pairs, function(pair) vapply(pair$folds, function(f) length(f$analysis), 1L)), use.names = FALSE),
                       vapply(full, function(pair) length(pair$rows), 1L))
    testthat::expect_identical(vapply(audit$fits, function(f) nrow(f$x), 1L), expected_rows, info = info)
    levels_seen <- lapply(audit$fits, function(f) levels(f$y))
    testthat::expect_identical(levels_seen, c(rep(lapply(pairs, `[[`, "class_levels"), each = sources),
                                              lapply(full, `[[`, "class_levels")), info = info)
    for (f in audit$fits) testthat::expect_identical(nlevels(droplevels(f$y)), 2L, info = info)
  }
  # Fit and predict failures name the class pair and its fold.
  only <- function(fit, predict) generics::fit(calibrated_workflow(e$base(model_adapter("classification", fit, predict), d$train),
    "truth", classification = "ovo", model_training = "pairwise", method = "cross", v = 2, seed = 1), d$train)
  half <- function(object, new_data) matrix(.5, nrow(new_data), 2L, dimnames = list(NULL, object))
  testthat::expect_error(only(function(x, y) if (identical(levels(y), c("a", "c"))) stop("pair broke") else levels(y), half),
    "pair_2 Fold1 base workflow training: model_adapter\\(\\) `fit` callback: pair broke")
  testthat::expect_error(only(function(x, y) levels(y),
    function(object, new_data) if (identical(object, c("b", "c"))) stop("no pair") else half(object, new_data)),
    "pair_3 Fold1 base workflow probability prediction.*model_adapter\\(\\) `predict` callback: no pair")
  testthat::expect_error(only(function(x, y) levels(y), function(object, new_data) half(c("a", "b"), new_data)),
    "pair_2 Fold1 base workflow probability prediction.*returned column for unknown class .b.")
})

testthat::test_that("adapter postprocessing keeps its B/P roles and inner ensembles", {
  skip_if_no_integration("tailor", "probably", "nnet")
  e <- adapter_factories()
  binary <- list(train = e$binary[1:160, ], new = e$binary[161:200, c("x", "z")])
  multi <- adapter_multiclass()
  numeric <- list(train = e$regression[1:120, ], new = e$regression[121:150, c("x", "z")])
  glm <- list(adapter = e$glm_adapter, native = parsnip::logistic_reg(engine = "glm"), data = binary)
  multinom <- list(adapter = e$multinom_adapter, native = parsnip::multinom_reg(engine = "nnet"), data = multi,
                   args = list(classification = "ovr"))
  lm <- list(adapter = e$lm_adapter, native = parsnip::linear_reg(engine = "lm"), data = numeric)
  cases <- list(
    list(glm, list(method = "sigmoid", sigmoid_input = "probability"), "inductive"),
    list(glm, list(method = "tailor_logistic"), "cross"),
    list(multinom, list(method = "tailor_multinomial"), "inductive"),
    list(lm, list(method = "tailor_linear"), "cross"),
    list(lm, list(method = "isotonic_numeric"), "inductive"))
  for (case in cases) for (strategy in c("holdout", "ensemble")) {
    engine <- case[[1L]]
    post <- c(case[[2L]], list(strategy = strategy), if (strategy == "holdout") list(prop = .3) else list(v = 2L))
    args <- c(list(outcome = "truth", method = case[[3L]], seed = 42, postprocessing = post), engine$args)
    if (case[[3L]] == "cross") args$v <- 2L
    calls <- new.env()
    calls$rows <- integer()
    counting <- engine$adapter()
    callbacks <- rlang::eval_tidy(counting$eng_args$parameters)
    learn <- callbacks$fit
    callbacks$fit <- function(x, y) {
      calls$rows <- c(calls$rows, nrow(x))
      learn(x, y)
    }
    counting <- parsnip::set_engine(counting, "callbacks", parameters = callbacks)
    d <- engine$data
    fit <- generics::fit(do.call(calibrated_workflow, c(list(e$base(counting, d$train)), args)), d$train)
    expected <- generics::fit(do.call(calibrated_workflow, c(list(e$base(engine$native, d$train)), args)), d$train)
    info <- paste(post$method, strategy, case[[3L]])
    testthat::expect_equal(venn_abers_details(fit, d$new), venn_abers_details(expected, d$new), tolerance = 1e-6, info = info)
    # One callback fit on each B, holdout or inner ensemble member; P and C never reach it.
    state <- workflows::extract_fit_engine(fit)
    plans <- unlist(lapply(state$postprocessing, function(p) if (is.null(p$inner_folds)) list(p) else p$inner_folds), recursive = FALSE)
    testthat::expect_identical(calls$rows, vapply(plans, function(p) length(p$B), 1L, USE.NAMES = FALSE), info = info)
    testthat::expect_identical(sum(vapply(state$postprocessing, `[[`, 1L, "observed_base_fits")), length(plans), info = info)
    testthat::expect_true(all(venn_abers_details(fit, d$new)$details$postprocessing_method == post$method), info = info)
  }
  # Inner ensemble failures keep the source and inner fold.
  second <- local({
    calls <- 0L
    glm_fit <- rlang::eval_tidy(e$glm_adapter()$eng_args$parameters)
    function(x, y) {
      calls <<- calls + 1L
      if (calls == 2L) stop("inner two") else glm_fit$fit(x, y)
    }
  })
  failing <- model_adapter("classification", second, rlang::eval_tidy(e$glm_adapter()$eng_args$parameters)$predict)
  testthat::expect_error(generics::fit(calibrated_workflow(e$base(failing, binary$train), "truth", seed = 1,
    postprocessing = list(method = "sigmoid", sigmoid_input = "probability", strategy = "ensemble", v = 2L)), binary$train),
    "Inductive1 / Fold2 sigmoid \\(B=64, P=64\\): model_adapter\\(\\) `fit` callback: inner two")
  # Task and method rules are unchanged: the adapter fails exactly like its native engine.
  rejects <- function(adapter, native, data, ...) {
    message <- function(spec) tryCatch({
      generics::fit(calibrated_workflow(e$base(spec, data), "truth", ...), data)
      NA_character_
    }, error = conditionMessage)
    found <- message(adapter)
    testthat::expect_false(is.na(found))
    testthat::expect_identical(found, message(native))
  }
  rejects(e$lm_adapter(), lm$native, numeric$train,
          postprocessing = list(method = "tailor_logistic", strategy = "holdout", prop = .3))
  rejects(e$lm_adapter(), lm$native, numeric$train,
          postprocessing = list(method = "sigmoid", sigmoid_input = "probability", strategy = "holdout", prop = .3))
  rejects(e$glm_adapter(), glm$native, binary$train,
          postprocessing = list(method = "tailor_multinomial", strategy = "holdout", prop = .3))
  rejects(e$glm_adapter(), glm$native, binary$train,
          postprocessing = list(method = "isotonic_numeric", strategy = "holdout", prop = .3))
})

testthat::test_that("glm_link stays reserved for the native GLM engine", {
  skip_if_no_integration()
  e <- adapter_factories()
  binary <- e$binary[1:160, ]
  multi <- adapter_multiclass()$train
  # The example adapter calls stats::glm() itself; it is still not the native engine.
  for (strategy in c("holdout", "ensemble")) {
    post <- c(list(method = "sigmoid", sigmoid_input = "glm_link", strategy = strategy),
              if (strategy == "holdout") list(prop = .3) else list(v = 2L))
    testthat::expect_error(calibrated_workflow(e$base(e$glm_adapter(), binary), "truth", postprocessing = post),
                           "sigmoid glm_link B requires a binomial GLM workflow")
    testthat::expect_error(calibrated_workflow(e$base(e$glm_adapter(), multi), "truth", classification = "ovo",
                                               model_training = "pairwise", postprocessing = post),
                           "sigmoid glm_link B requires a binomial GLM workflow")
  }
  testthat::expect_no_error(calibrated_workflow(e$base(parsnip::logistic_reg(engine = "glm"), binary), "truth",
    postprocessing = list(method = "sigmoid", sigmoid_input = "glm_link", strategy = "holdout", prop = .3)))
})

testthat::test_that("grid and Bayes searches tune calibration and recipe parameters, never callbacks", {
  skip_if_no_integration("tune", "yardstick", "dials")
  e <- adapter_factories()
  d <- e$binary[1:160, ]
  new <- e$binary[161:200, c("x", "z")]
  tuned <- function(spec) workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x + z, d) |> recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(spec)
  folds <- function(data) rsample::vfold_cv(data, v = 2)
  same <- function(tuning, info, ...) {
    adapter <- generics::fit(calibrated_workflow(tuned(e$glm_adapter()), "truth", seed = 3, tuning = tuning, ...), d)
    native <- generics::fit(calibrated_workflow(tuned(parsnip::logistic_reg(engine = "glm")), "truth", seed = 3, tuning = tuning, ...), d)
    testthat::expect_equal(venn_abers_details(adapter, new), venn_abers_details(native, new), info = info)
    a <- venn_abers_tuning(adapter)
    n <- venn_abers_tuning(native)
    testthat::expect_identical(lapply(a$models, `[[`, "selected"), lapply(n$models, `[[`, "selected"), info = info)
    testthat::expect_identical(a$calibration$selected, n$calibration$selected, info = info)
  }
  for (strategy in c("joint", "two_stage")) {
    same(list(strategy = strategy, selection_prop = .5,
              model = list(method = "grid", grid = data.frame(degree = 1:2), resamples = folds),
              calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")), resamples = folds)),
         paste("grid", strategy), loss = tune::tune("va_loss"))
  }
  degree <- dials::parameters(dials::new_quant_param(type = "integer", range = c(1L, 3L), inclusive = c(TRUE, TRUE),
                                                     label = c(degree = "Degree")))
  same(list(strategy = "two_stage", model = list(method = "bayes", initial = 3L, iter = 2L, param_info = degree,
              resamples = folds, control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 2L))), "bayes")
  # The adapter declares no parameters: grid columns for it and replacements fail.
  testthat::expect_identical(nrow(generics::tune_args(e$glm_adapter())), 0L)
  testthat::expect_error(calibrated_workflow(tuned(e$glm_adapter()), "truth", tuning = list(strategy = "two_stage",
    model = list(method = "grid", grid = data.frame(degree = 1:2, penalty = c(.1, .2)), resamples = folds))),
    "needs rows with exactly these ID columns: degree")
  testthat::expect_error(calibrated_workflow(e$base(e$glm_adapter(), d), "truth", tuning = list(strategy = "two_stage",
    model = list(method = "grid", grid = data.frame(penalty = c(.1, .2)), resamples = folds))),
    "`tuning\\$model` is required exactly when the base workflow has tune\\(\\) placeholders")
  spec <- e$glm_adapter()
  testthat::expect_identical(stats::update(spec), spec)
  testthat::expect_error(stats::update(e$glm_adapter(), penalty = 1),
    "model_adapter\\(\\) specs have no tunable or replaceable parameters")
})

adapter_glm <- function(audit = NULL, reorder = FALSE) model_adapter("classification",
  fit = function(x, y) {
    if (!is.null(audit)) audit$fits[[length(audit$fits) + 1L]] <- list(x = x, y = y)
    list(model = stats::glm(y ~ ., data = cbind(x, y = y), family = stats::binomial()), classes = levels(y))
  },
  predict = function(object, new_data) {
    p <- stats::predict(object$model, new_data, type = "response")
    p <- matrix(c(1 - p, p), ncol = 2L, dimnames = list(NULL, object$classes))
    if (reorder) p[, 2:1, drop = FALSE] else p
  })

adapter_base <- function(spec, data) workflows::workflow() |>
  workflows::add_recipe(recipes::recipe(truth ~ ., data) |> recipes::step_normalize(recipes::all_numeric_predictors())) |>
  workflows::add_model(spec)

adapter_binary <- function() {
  d <- data.frame(x = sin(1:160) + rep(c(1, 0, 0, 1), 40), z = cos(1:160),
                  truth = factor(rep(c("yes", "no", "no", "yes"), 40), levels = c("yes", "no")))
  list(train = d[1:120, ], new = d[121:160, c("x", "z")])
}

testthat::test_that("model_adapter() builds an unfitted parsnip spec and rejects malformed constructors", {
  skip_if_no_integration()
  fit <- function(x, y) NULL
  predict <- function(object, new_data) NULL
  spec <- model_adapter("regression", fit, predict)
  testthat::expect_s3_class(spec, c("vennabers_adapter", "model_spec"))
  testthat::expect_identical(c(spec$mode, spec$engine), c("regression", "callbacks"))
  callbacks <- rlang::eval_tidy(spec$eng_args$parameters)
  testthat::expect_identical(callbacks[c("fit", "predict")], list(fit = fit, predict = predict))
  testthat::expect_true("vennabers" %in% generics::required_pkgs(spec))
  testthat::expect_identical(nrow(generics::tune_args(spec)), 0L)
  testthat::expect_error(model_adapter("censored regression", fit, predict), "`mode` must be one of: classification, regression")
  testthat::expect_error(model_adapter(c("classification", "regression"), fit, predict), "`mode` must be one of")
  testthat::expect_error(model_adapter("classification", "glm", predict), "`fit` must be a function\\(x, y\\)")
  testthat::expect_error(model_adapter("classification", fit, NULL), "`predict` must be a function\\(object, new_data\\)")
  testthat::expect_error(model_adapter("classification", fit), "Supply `mode`, `fit` and `predict`")
  testthat::expect_error(model_adapter("classification", fit, pred = predict), "Unknown argument: .pred.")
  testthat::expect_error(model_adapter("classification", fit, predict, seed = 1), "unused argument")
})

testthat::test_that("callbacks learn on prepared local rows and calibrate like the native GLM", {
  skip_if_no_integration()
  d <- adapter_binary()
  native <- adapter_base(parsnip::logistic_reg(engine = "glm"), d$train)
  for (method in c("inductive", "cross")) for (mode in c("matched", "full")) {
    audit <- new.env()
    audit$fits <- list()
    args <- list(outcome = "truth", method = method, prediction_model = mode, seed = 42)
    fit <- generics::fit(do.call(calibrated_workflow, c(list(adapter_base(adapter_glm(audit), d$train)), args)), d$train)
    expected <- generics::fit(do.call(calibrated_workflow, c(list(native), args)), d$train)
    info <- paste(method, mode)
    testthat::expect_equal(predict(fit, d$new, type = "prob"), predict(expected, d$new, type = "prob"), info = info)
    testthat::expect_identical(predict(fit, d$new, type = "class"), predict(expected, d$new, type = "class"), info = info)
    testthat::expect_equal(venn_abers_details(fit, d$new), venn_abers_details(expected, d$new), info = info)
    # One fresh callback fit per source fold, then Full1 on all local rows.
    state <- workflows::extract_fit_engine(fit)
    rows <- c(lapply(state$folds, function(f) length(f$analysis)), if (mode == "full") list(Full1 = 120L))
    testthat::expect_identical(lapply(audit$fits, function(f) nrow(f$x)), unname(rows), info = info)
    for (seen in audit$fits) {
      testthat::expect_identical(class(seen$x), "data.frame")
      testthat::expect_identical(names(seen$x), c("x", "z"))
      # The recipe was prepared on exactly these rows before the callback ran.
      testthat::expect_equal(unname(colMeans(seen$x)), c(0, 0), tolerance = 1e-12)
      testthat::expect_identical(levels(seen$y), c("yes", "no"))
    }
  }
  # Columns align by name, not position.
  swapped <- generics::fit(calibrated_workflow(adapter_base(adapter_glm(reorder = TRUE), d$train), "truth", seed = 42), d$train)
  plain <- generics::fit(calibrated_workflow(adapter_base(adapter_glm(), d$train), "truth", seed = 42), d$train)
  testthat::expect_identical(predict(swapped, d$new, type = "prob"), predict(plain, d$new, type = "prob"))
  # A plain adapter workflow predicts like the native model, classes by first maximum.
  raw <- generics::fit(adapter_base(adapter_glm(), d$train), d$train)
  reference <- generics::fit(native, d$train)
  testthat::expect_equal(predict(raw, d$new, type = "prob"), predict(reference, d$new, type = "prob"))
  testthat::expect_identical(predict(raw, d$new), predict(reference, d$new))
  testthat::expect_identical(names(predict(raw, d$new[0, ], type = "prob")), c(".pred_yes", ".pred_no"))
})

testthat::test_that("shared multiclass and regression adapters match their native engines", {
  skip_if_no_integration("nnet")
  classes <- factor(rep(c("a", "b", "c"), 50))
  multi <- data.frame(x = sin(1:150) + as.integer(classes), z = cos(1:150), truth = classes)
  multinom <- model_adapter("classification",
    fit = function(x, y) nnet::multinom(y ~ ., data = cbind(x, y = y), trace = FALSE),
    predict = function(object, new_data) {
      p <- stats::predict(object, new_data, type = "probs")
      matrix(p, nrow(new_data), length(object$lev), dimnames = list(NULL, object$lev))
    })
  native <- adapter_base(parsnip::multinom_reg(engine = "nnet"), multi[1:120, ])
  for (classification in c("ovo", "ovr")) for (method in c("inductive", "cross")) {
    args <- list(outcome = "truth", classification = classification, method = method, seed = 7)
    fit <- generics::fit(do.call(calibrated_workflow, c(list(adapter_base(multinom, multi[1:120, ])), args)), multi[1:120, ])
    expected <- generics::fit(do.call(calibrated_workflow, c(list(native), args)), multi[1:120, ])
    new <- multi[121:150, c("x", "z")]
    testthat::expect_equal(predict(fit, new, type = "prob"), predict(expected, new, type = "prob"),
                           tolerance = 1e-6, info = paste(classification, method))
    testthat::expect_identical(predict(fit, new[1, ], type = "prob")[0, ], predict(fit, new[0, ], type = "prob"))
  }
  regression <- data.frame(x = sin(1:120), z = cos(1:120))
  regression$truth <- 10 + 3 * regression$x + regression$z^2
  lm_adapter <- model_adapter("regression", fit = function(x, y) stats::lm(y ~ ., data = cbind(x, y = y)),
    predict = function(object, new_data) unname(stats::predict(object, new_data)))
  reference <- adapter_base(parsnip::linear_reg(engine = "lm"), regression[1:90, ])
  for (method in c("inductive", "cross")) for (mode in c("matched", "full")) {
    args <- list(outcome = "truth", method = method, prediction_model = mode, m = 2, seed = 5)
    fit <- generics::fit(do.call(calibrated_workflow, c(list(adapter_base(lm_adapter, regression[1:90, ])), args)), regression[1:90, ])
    expected <- generics::fit(do.call(calibrated_workflow, c(list(reference), args)), regression[1:90, ])
    new <- regression[91:120, c("x", "z")]
    testthat::expect_equal(venn_abers_details(fit, new), venn_abers_details(expected, new), info = paste(method, mode))
  }
})

testthat::test_that("prediction contract violations fail concretely with fold and full-fit context", {
  skip_if_no_integration()
  d <- adapter_binary()
  fails <- function(predict, pattern, fit = function(x, y) levels(y), ..., mode = "classification", data = d$train) {
    spec <- model_adapter(mode, fit, predict)
    testthat::expect_error(generics::fit(calibrated_workflow(adapter_base(spec, data), "truth", seed = 1, ...), data), pattern)
  }
  probs <- function(object, n, p = .3) matrix(c(rep(p, n), rep(1 - p, n)), n, 2L, dimnames = list(NULL, object))
  fails(function(object, new_data) probs(object, nrow(new_data))[, 1L],
        "Inductive1 base workflow probability prediction.*must return a numeric matrix with columns named by the classes .yes., .no.; received a numeric of length")
  fails(function(object, new_data) unname(probs(object, nrow(new_data))), "must return a numeric matrix with columns named")
  fails(function(object, new_data) as.data.frame(probs(object, nrow(new_data))), "received a data.frame with dimensions")
  fails(function(object, new_data) probs(object, nrow(new_data))[, 1L, drop = FALSE], "returned no column for class .no.")
  fails(function(object, new_data) probs(c("yes", "yes"), nrow(new_data)), "returned duplicate column .yes.")
  fails(function(object, new_data) probs(c("yes", "No"), nrow(new_data)), "returned column for unknown class .No.")
  fails(function(object, new_data) probs(object, nrow(new_data) + 1L), "returned 25 rows for 24 input rows")
  fails(function(object, new_data) probs(object, nrow(new_data), p = 1.2), "must return finite probabilities in \\[0, 1\\] that sum to 1")
  fails(function(object, new_data) probs(object, nrow(new_data), p = NA), "finite probabilities")
  fails(function(object, new_data) probs(object, nrow(new_data)) * .99, "sum to 1 in every row")
  fails(function(object, new_data) stop("no model"), "Inductive1 base workflow probability prediction.*model_adapter\\(\\) `predict` callback: no model")
  fails(function(object, new_data) probs(object, nrow(new_data)), "Inductive1 base workflow training: model_adapter\\(\\) `fit` callback: bad start",
        fit = function(x, y) stop("bad start"))
  second <- local({
    calls <- 0L
    function(x, y) {
      calls <<- calls + 1L
      if (calls == 2L) stop("second fold") else levels(y)
    }
  })
  fails(function(object, new_data) probs(object, nrow(new_data)), "Fold2 base workflow training: model_adapter\\(\\) `fit` callback: second fold",
        fit = second, method = "cross")
  # Full1 is fitted after every source and keeps its own context.
  fails(function(object, new_data) probs(object, nrow(new_data)), "Full1 base workflow training: model_adapter\\(\\) `fit` callback: all rows",
        fit = function(x, y) if (nrow(x) == 120L) stop("all rows") else levels(y), prediction_model = "full")
  regression <- data.frame(x = sin(1:40), truth = cos(1:40))
  for (case in list(list(function(object, new_data) matrix(0, nrow(new_data), 1L), "received a matrix with dimensions 8 x 1"),
                    list(function(object, new_data) numeric(nrow(new_data) - 1L), "received a numeric of length 7 for 8 input rows"),
                    list(function(object, new_data) c(Inf, numeric(nrow(new_data) - 1L)), "one finite value per input row"),
                    list(function(object, new_data) as.character(seq_len(nrow(new_data))), "received a character"))) {
    fails(case[[1L]], paste0("Inductive1 base workflow numeric prediction: model_adapter\\(\\) `predict` must return a numeric vector.*", case[[2L]]),
          fit = function(x, y) NULL, mode = "regression", data = regression)
  }
})

testthat::test_that("adapters keep the local seed, frozen adoption, outer resamples and RDS", {
  skip_if_no_integration()
  d <- adapter_binary()
  # Callbacks drawing from R's generator follow the local calibrated_workflow() seed.
  jitter <- model_adapter("classification",
    fit = function(x, y) list(shift = stats::runif(1, -.1, .1), classes = levels(y), mean = mean(x$x)),
    predict = function(object, new_data) {
      p <- stats::plogis(new_data$x - object$mean + object$shift)
      matrix(c(p, 1 - p), ncol = 2L, dimnames = list(NULL, object$classes))
    })
  seeded <- function(seed) workflows::extract_fit_engine(generics::fit(calibrated_workflow(adapter_base(jitter, d$train), "truth",
    method = "cross", seed = seed), d$train))
  shift <- function(state) vapply(state$models, function(m) workflows::extract_fit_engine(m)$state$shift, 0)
  testthat::expect_identical(shift(seeded(11)), shift(seeded(11)))
  testthat::expect_false(identical(shift(seeded(11)), shift(seeded(12))))

  # A fitted adapter workflow is an ordinary workflow for frozen adoption.
  native <- adapter_base(parsnip::logistic_reg(engine = "glm"), d$train)
  frozen <- calibrate_workflow(generics::fit(adapter_base(adapter_glm(), d$train[1:80, ]), d$train[1:80, ]), d$train[81:120, ], "truth")
  expected <- calibrate_workflow(generics::fit(native, d$train[1:80, ]), d$train[81:120, ], "truth")
  testthat::expect_equal(predict(frozen, d$new), predict(expected, d$new))

  # Outer resamples refit the whole calibrated adapter workflow per analysis set.
  model <- calibrated_workflow(adapter_base(adapter_glm(), d$train), "truth", seed = 42)
  if (requireNamespace("tune", quietly = TRUE) && requireNamespace("yardstick", quietly = TRUE)) {
    outer <- withr::with_seed(3, rsample::vfold_cv(d$train, v = 2))
    run <- function(workflow) tune::collect_predictions(tune::fit_resamples(workflow, outer,
      metrics = yardstick::metric_set(yardstick::mn_log_loss), control = tune::control_resamples(save_pred = TRUE)))
    testthat::expect_equal(run(model), run(calibrated_workflow(native, "truth", seed = 42)))
  }

  # The saved fit carries its callbacks and fitted states.
  fit <- generics::fit(model, d$train)
  restored <- unserialize(serialize(fit, NULL))
  testthat::expect_identical(predict(restored, d$new, type = "prob"), predict(fit, d$new, type = "prob"))
  testthat::expect_identical(venn_abers_details(restored, d$new), venn_abers_details(fit, d$new))
})

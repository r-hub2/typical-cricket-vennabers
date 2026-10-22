testthat::test_that("Tailor configuration remains canonical when added or adopted", {
  for (pkg in c("workflows", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  for (method in c("logistic", "multinomial", "linear")) {
    numeric <- method == "linear"
    base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
      workflows::add_model(if (numeric) parsnip::linear_reg(engine = "lm") else
        parsnip::logistic_reg(engine = "glm"))
    canonical <- if (numeric) tailor::tailor() |>
      tailor::adjust_numeric_calibration(method = method, smooth = FALSE) else
      tailor::tailor() |> tailor::adjust_probability_calibration(method = method, smooth = FALSE)
    for (workflow in list(base, workflows::add_tailor(base, canonical))) {
      adapted <- .postprocessing_workflow(workflow, list(method = paste0("tailor_", method)))
      testthat::expect_identical(workflows::extract_postprocessor(adapted), canonical)
    }
  }
})

testthat::test_that("the owned boundary preserves fitted Tailor and adapts prediction batches", {
  for (pkg in c("workflows", "tailor", "probably", "nnet")) testthat::skip_if_not_installed(pkg)
  for (method in c("logistic", "multinomial", "linear")) {
    classes <- if (method == "linear") NULL else c("z", "a", "b")[seq_len(if (method == "logistic") 2L else 3L)]
    d <- withr::with_seed(123, {
      x <- matrix(rnorm(720), ncol = 3L)
      data.frame(x = x[, 1L], y = x[, 2L], z = x[, 3L], truth = if (is.null(classes))
        x[, 1L] + x[, 2L]^2 else factor(rep(classes, length.out = 240), levels = classes))
    })
    base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |>
      workflows::add_model(if (method == "linear") parsnip::linear_reg(engine = "lm") else
        if (method == "logistic") parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
    configured <- .postprocessing_workflow(base, list(method = paste0("tailor_", method)))
    domain <- list(outcome = "truth", classes = classes)
    result <- .fit_tailor_workflow(configured, d[1:120, ], d[121:180, ], domain, method, "boundary")
    reference <- generics::fit(configured, d[1:120, ], data_calibration = d[121:180, ])
    learned <- workflows::extract_tailor(result$model)
    testthat::expect_equal(learned, workflows::extract_tailor(reference))
    testthat::expect_identical(workflows::extract_postprocessor(result$model),
      workflows::extract_postprocessor(configured))
    new <- d[181:190, setdiff(names(d), "truth")]
    for (type in if (is.null(classes)) "numeric" else c("prob", "class")) {
      expected <- predict(reference, new, type = type)
      for (n in c(0L, 1L, nrow(new))) {
        actual <- predict(result$model, new[seq_len(n), , drop = FALSE], type = type)
        testthat::expect_identical(names(actual), names(expected))
        testthat::expect_equal(nrow(actual), n)
        if (n) testthat::expect_equal(actual, expected[seq_len(n), , drop = FALSE], tolerance = 1e-10)
      }
    }
    testthat::expect_identical(workflows::extract_tailor(result$model), learned)
  }
  original <- getS3method("predict", "tailor", envir = asNamespace("stats"))
  registerS3method("predict", "tailor", function(object, new_data, ...) {
    output <- original(object, new_data, ...)
    output[c(seq_len(nrow(output)), 1L), , drop = FALSE]
  }, envir = asNamespace("stats"))
  on.exit(registerS3method("predict", "tailor", original, envir = asNamespace("stats")))
  testthat::expect_error(predict(result$model, new, type = "numeric"), "Tailor output requires one prediction per input row")
})

testthat::test_that("prediction-response", {
library(vennabers)

# Reject diagnostic requests while still running the real core prediction method.
predict.response_only <- function(object, ..., type = "response") {
  stopifnot(type == "response")
  NextMethod()
}

packages <- c("parsnip", "workflows", "rsample", "recipes", "withr", "nnet")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  for (mode in c("ovo", "ovr", "binary", "regression")) {
    data <- data.frame(x = iris$Sepal.Length, truth = iris$Species)
    spec <- parsnip::multinom_reg() |> parsnip::set_engine("nnet", trace = FALSE)
    if (mode == "binary") {
      data$truth <- factor(ifelse(data$truth == "setosa", "yes", "no"))
      spec <- parsnip::logistic_reg(engine = "glm")
    } else if (mode == "regression") {
      data$truth <- iris$Petal.Length
      spec <- parsnip::linear_reg(engine = "lm")
    }
    base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
      workflows::add_model(spec)
    for (method in c("inductive", "cross")) {
      options <- if (mode == "regression") list() else list(classification = mode)
      fitted <- generics::fit(do.call(calibrated_workflow,
        c(list(base, "truth", method = method, seed = 42), options)), data)
      for (loss in if (mode == "regression") list(NULL) else c("log", "brier")) {
        for (n in c(0L, 1L, 7L)) {
          new <- data[seq_len(n), , drop = FALSE]
          args <- if (is.null(loss)) list() else list(loss = loss)
          details <- do.call(venn_abers_details, c(list(fitted, new), args))
          type <- if (mode == "regression") "numeric" else "prob"
          stopifnot(identical(details$prediction,
            do.call(predict, c(list(fitted, new, type = type), args))))

          state <- workflows::extract_fit_engine(fitted)
          state$loss <- loss
          class(state$calibrator) <- c("response_only", class(state$calibrator))
          result <- vennabers:::.integration_prediction(state, new["x"])
          stopifnot(identical(names(result), "prediction"),
                    identical(result$prediction, details$prediction))
          if (mode == "regression") {
            stopifnot(identical(vennabers:::predict_vennabers_numeric(state, new["x"]),
                                details$prediction$.pred))
          } else {
            stopifnot(identical(vennabers:::predict_vennabers_prob(state, new["x"]),
                                details$prediction),
                      identical(vennabers:::predict_vennabers_class(state, new["x"]),
                        do.call(predict, c(list(fitted, new, type = "class"), args))$.pred_class))
          }
        }
      }
    }
  }
}

  testthat::expect_true(TRUE)
})

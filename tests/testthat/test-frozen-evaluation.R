testthat::test_that("frozen-evaluation", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

packages <- c("parsnip", "workflows", "hardhat", "rsample", "recipes", "tune", "yardstick", "nnet",
              "rlang", "tibble", "withr", "generics", "tidyselect")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  # Native tidymodels evaluation trains the complete base workflows. Frozen adoption
  # then calibrates each returned workflow on separate calibration rows, which
  # neither the outer analysis nor the outer assessment data contain.
  post_available <- all(vapply(c("tailor", "probably", "mgcv"), requireNamespace, logical(1), quietly = TRUE))
  guards <- list(c("stats", "glm.fit"), c("stats", "lm.fit"), c("nnet", "multinom"),
                 c("workflows", "fit.workflow"), c("recipes", "prep.recipe"))
  if (post_available) guards <- c(guards, list(c("tailor", "fit.tailor"), c("mgcv", "gam")))
  active <- FALSE
  guard <- function(on) if (on != active) {
    for (g in guards) if (on) trace(g[2L], where = asNamespace(g[1L]), tracer = quote(stop("Unexpected learning")),
                                    print = FALSE) else untrace(g[2L], where = asNamespace(g[1L]))
    active <<- on
  }
  on.exit(guard(FALSE), add = TRUE)
  for (method in c("binary", "ovo", "ovr", "regression")) for (post in if (post_available) c(FALSE, TRUE) else FALSE) {
    regression <- method == "regression"
    classes <- if (method == "binary") c("yes", "no") else c("third", "first", "second")
    set.seed(907)
    d <- data.frame(x = rnorm(480), z = rnorm(480))
    score <- exp(cbind(0, d$x, -d$z)[, seq_along(classes), drop = FALSE])
    d$truth <- if (regression) 2 + d$x + cos(1:480) else factor(classes[apply(score / rowSums(score), 1L,
      function(p) sample.int(length(p), 1L, prob = p))], levels = classes)
    outer_data <- d[c(1:240, 361:480), ]
    calibration <- d[241:360, ]
    spec <- if (regression) parsnip::linear_reg(engine = "lm") else if (method == "binary")
      parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE)
    base <- workflows::workflow() |>
      workflows::add_recipe(recipes::recipe(truth ~ ., outer_data) |>
        recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(spec)
    if (post) base <- workflows::add_tailor(base, if (regression) tailor::tailor() |>
      tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE) else tailor::tailor() |>
        tailor::adjust_probability_calibration(method = if (method == "binary") "logistic" else "multinomial",
                                              smooth = method != "binary"))
    metrics <- if (regression) yardstick::metric_set(yardstick::rmse) else
      yardstick::metric_set(yardstick::mn_log_loss, yardstick::brier_class)
    columns <- if (regression) ".pred" else paste0(".pred_", classes)
    type <- if (regression) "numeric" else "prob"
    adopt <- function(workflow, data) {
      if (regression) calibrate_workflow(workflow, data, "truth", epsilon = .1) else
        calibrate_workflow(workflow, data, "truth", classification = method, loss = "brier")
    }
    manual <- function(workflow, data) {
      p <- predict(workflow, data[c("x", "z")], type = type)
      if (regression) return(venn_abers_regression(data$truth, p$.pred, epsilon = .1))
      p <- as.matrix(p[columns])
      colnames(p) <- classes
      if (method == "binary") venn_abers(data$truth == classes[1L], p[, 1L]) else
        venn_abers(data$truth, p, method = method)
    }

    # Repeated outer resamples: every saved workflow is the natively assessed one.
    folds <- withr::with_seed(31, rsample::vfold_cv(d[1:240, ], v = 2, repeats = 2))
    evaluation <- withr::with_seed(32, tune::fit_resamples(base, folds, metrics = metrics,
      control = tune::control_resamples(save_pred = TRUE, extract = function(x) x)))
    native <- tune::collect_predictions(evaluation)
    extracts <- tune::collect_extracts(evaluation)
    assert_conditions(nrow(tune::collect_notes(evaluation)) == 0L, nrow(extracts) == nrow(folds),
              identical(paste(extracts$id, extracts$id2), paste(folds$id, folds$id2)),
              nrow(tune::collect_metrics(evaluation, summarize = FALSE)) == nrow(folds) * if (regression) 1L else 2L)
    guard(TRUE)
    keys <- character()
    for (i in seq_len(nrow(folds))) {
      workflow <- extracts$.extracts[[i]]
      rows <- rsample::complement(folds$splits[[i]])
      assessed <- d[rows, c("x", "z")]
      assessed_native <- native[native$id == folds$id[i] & native$id2 == folds$id2[i], ]
      assessed_native <- assessed_native[order(assessed_native$.row), ]
      adopted <- adopt(workflow, calibration)
      calibrated <- predict(adopted, assessed)
      assert_conditions(identical(assessed_native$.row, sort(rows)),
                max(abs(as.matrix(assessed_native[columns]) - as.matrix(predict(workflow, d[assessed_native$.row, c("x", "z")], type = type)[columns]))) <= 1e-10,
                identical(adopted$models[[1L]], workflow), identical(adopted$calibrator, manual(workflow, calibration)),
                is.null(adopted$folds), identical(names(adopted$models), "Inductive1"),
                nrow(calibrated) == length(rows), all(is.finite(as.matrix(calibrated))))
      keys <- c(keys, paste(folds$id[i], folds$id2[i], rows, extracts$.config[i]))
    }
    stopifnot(!anyDuplicated(keys), setequal(keys, paste(native$id, native$id2, native$.row, native$.config)))
    guard(FALSE)

    # last_fit: the extracted workflow learned only from its analysis rows.
    split <- withr::with_seed(33, rsample::initial_split(outer_data, prop = 1 / 2))
    final <- withr::with_seed(34, tune::last_fit(base, split, metrics = metrics))
    test <- rsample::testing(split)
    changed <- outer_data
    changed$truth[-split$in_id] <- if (regression) -changed$truth[-split$in_id] else rev(changed$truth[-split$in_id])
    changed_split <- withr::with_seed(33, rsample::initial_split(changed, prop = 1 / 2))
    other <- withr::with_seed(34, tune::last_fit(base, changed_split, metrics = metrics))
    stopifnot(identical(changed_split$in_id, split$in_id))
    guard(TRUE)
    workflow <- tune::extract_workflow(final)
    native <- tune::collect_predictions(final)
    adopted <- adopt(workflow, calibration)
    adopted_other <- adopt(tune::extract_workflow(other), calibration)
    assert_conditions(identical(native$.row, setdiff(seq_len(nrow(outer_data)), split$in_id)),
              max(abs(as.matrix(native[columns]) - as.matrix(predict(workflow, test[c("x", "z")], type = type)[columns]))) <= 1e-10,
              identical(adopted$calibrator, manual(workflow, calibration)),
              # Outer test targets reach neither the extracted workflow nor its calibration.
              identical(predict(adopted_other, test), predict(adopted, test)),
              identical(venn_abers_details(adopted_other, test), venn_abers_details(adopted, test)))

    # yardstick scores the adopted probabilities with its own Brier scale:
    # equal for binary, half the summed multiclass score of brier_score().
    if (!regression) {
      p <- predict(adopted, test[c("x", "z")])
      estimate <- if (method == "binary") p$.pred_yes else structure(as.matrix(p[columns]), dimnames = list(NULL, classes))
      score <- if (method == "binary") brier_score(test$truth == "yes", estimate) else brier_score(test$truth, estimate)
      stopifnot(abs(yardstick::brier_class_vec(test$truth, estimate) - score / if (method == "binary") 1 else 2) <= 1e-10)
    }

    # Narrow loss-only choice on separate validation rows: existing predict and
    # yardstick calls, no refit and no change to the adopted state.
    if (!regression) {
      validation <- test[seq_len(nrow(test) %/% 2L), ]
      assessment <- test[-seq_len(nrow(test) %/% 2L), ]
      before <- adopted
      losses <- vapply(c("log", "brier"), function(loss) {
        p <- predict(adopted, validation[c("x", "z")], type = "prob", loss = loss)
        if (method == "binary") yardstick::mn_log_loss_vec(validation$truth, p$.pred_yes) else
          yardstick::mn_log_loss_vec(validation$truth, as.matrix(p[columns]))
      }, numeric(1))
      chosen <- c("log", "brier")[which.min(losses)]
      assert_conditions(all(is.finite(losses)), identical(adopted, before),
                identical(predict(adopted, assessment, loss = chosen),
                          predict(calibrate_workflow(workflow, calibration, "truth", classification = method, loss = chosen),
                                  assessment)))
    }
    guard(FALSE)
  }
})

  testthat::expect_true(TRUE)
})

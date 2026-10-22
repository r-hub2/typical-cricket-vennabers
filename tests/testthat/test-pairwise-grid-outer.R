testthat::test_that("nested pairwise grids score all classes outside testthat mode", {
  testthat::skip_if_not_installed("tune")
  testthat::skip_if_not_installed("rpart")
  # tune enables yardstick's column-selection cache only outside testthat mode.
  withr::local_envvar(TESTTHAT = "false")
  withr::local_seed(804L)
  data <- iris
  names(data)[5L] <- "truth"
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ .) |>
    workflows::add_model(parsnip::decision_tree(tree_depth = tune::tune("depth")) |>
      parsnip::set_mode("classification") |> parsnip::set_engine("rpart"))
  folds <- function(d) rsample::vfold_cv(d, v = 2L, strata = truth, pool = 0)
  config <- list(strategy = "two_stage", selection_prop = .5,
    model = list(method = "grid", grid = data.frame(depth = 1:2), resamples = folds),
    calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")),
      resamples = folds))
  calibrated <- calibrated_workflow(base, "truth", method = "cross", v = 2L,
    prediction_model = "full", classification = "ovo", model_training = "pairwise",
    loss = tune::tune("va_loss"), tuning = config, seed = 806L)
  result <- tune::fit_resamples(calibrated, folds(data),
    metrics = yardstick::metric_set(yardstick::mn_log_loss),
    control = tune::control_resamples(save_pred = TRUE, allow_par = FALSE,
      extract = function(w) venn_abers_tuning(w)$calibration))
  testthat::expect_equal(nrow(tune::collect_notes(result)), 0L)
  metrics <- tune::collect_metrics(result, summarize = FALSE)
  testthat::expect_equal(nrow(metrics), 2L)
  testthat::expect_true(all(is.finite(metrics$.estimate)))
  predictions <- tune::collect_predictions(result)
  testthat::expect_identical(sort(predictions$.row), seq_len(nrow(data)))
  columns <- paste0(".pred_", levels(data$truth))
  for (id in result$id) {
    fold <- predictions[predictions$id == id, ]
    probability <- as.matrix(fold[columns])
    colnames(probability) <- levels(data$truth)
    testthat::expect_equal(metrics$.estimate[metrics$id == id],
      log_loss(fold$truth, probability), tolerance = 1e-10)
  }
  for (extract in result$.extracts) {
    selection <- extract$.extracts[[1L]]
    testthat::expect_identical(selection$candidates$status, rep("eligible", 2L))
    testthat::expect_true(all(is.finite(selection$candidates$estimate)))
    testthat::expect_equal(nrow(tune::collect_notes(selection$native)), 0L)
    native <- tune::collect_metrics(selection$native, summarize = FALSE)
    for (candidate_id in selection$candidates$id) {
      scores <- selection$scores[[candidate_id]]
      candidate <- native[native$candidate_id == candidate_id, ]
      testthat::expect_equal(candidate$.estimate[match(names(scores), candidate$id)],
        unname(scores), tolerance = 1e-10)
    }
  }
})

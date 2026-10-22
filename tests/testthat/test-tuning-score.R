testthat::test_that("Grid IDs and Bayes values share strict candidate scoring", {
  testthat::skip_if_not_installed("yardstick")
  data <- data.frame(truth = c(1, 2, 3, 4))
  plans <- list("Repeat1/Fold1" = list(assessment = 1:2),
                "Repeat1/Fold2" = list(assessment = 3:4))
  candidates <- data.frame(id = c("model_001", "model_002"))
  search <- list(primary_metric = yardstick::rmse, report_metrics = list(
    broken = yardstick::new_numeric_metric(function(...) stop("report failed"), "minimize")))
  counter <- list(n = 7L)
  score <- function(result, keys) .score_tuning_candidates(search, result, data, plans,
    "truth", "first", "binary", candidates, keys, counter, 3L, list(candidates = 2L))
  results <- list()
  for (method in c("grid", "bayes")) {
    keys <- if (method == "grid") data.frame(candidate_id = candidates$id) else
      data.frame(degree = c(1L, 2L), loss = c("log", "brier"))
    native <- data.frame(id = c("Repeat1", "Repeat1"), id2 = c("Fold1", "Fold2"))
    native$.predictions <- lapply(plans, function(plan) {
      rows <- rev(plan$assessment)
      predictions <- do.call(rbind, lapply(1:2, function(i)
        cbind(keys[rep(i, length(rows)), , drop = FALSE],
          data.frame(.row = rows, truth = data$truth[rows], .pred = data$truth[rows] + i))))
      predictions[c(1L, 3L, 2L, 4L), , drop = FALSE]
    })
    native$.extracts <- lapply(plans, function(plan) {
      extracts <- keys[2:1, , drop = FALSE]
      extracts$.extracts <- list(list(decision = "second"), list(decision = "first"))
      extracts
    })
    results[[method]] <- score(native, keys)
    testthat::expect_identical(results[[method]]$candidates$status, rep("eligible", 2L))
    testthat::expect_equal(results[[method]]$candidates$estimate, c(1, 2))
    testthat::expect_identical(results[[method]]$decisions$model_001,
      stats::setNames(list("first", "first"), names(plans)))
    testthat::expect_identical(results[[method]]$report_metrics$model_001[[1L]]$broken,
      "report failed")
    testthat::expect_identical(results[[method]]$cost$observed_base_fits, 4L)

    # Damage only one candidate in one fold; partial scores cannot make it eligible.
    corruptions <- list(
      missing = function(p) p[-1L, ],
      duplicate = function(p) rbind(p, p[1L, ]),
      row = function(p) { p$.row[1L] <- NA_integer_; p },
      truth = function(p) { p$truth[1L] <- -1; p },
      nonfinite = function(p) { p$.pred[1L] <- Inf; p })
    for (corrupt in corruptions) {
      broken <- native
      broken$.predictions[[2L]] <- corrupt(broken$.predictions[[2L]])
      scored <- score(broken, keys)
      testthat::expect_identical(scored$candidates$status, c("excluded", "eligible"))
      testthat::expect_true(is.na(scored$candidates$estimate[1L]))
      testthat::expect_match(scored$candidates$reason[1L], "Repeat1/Fold2")
      testthat::expect_identical(scored$eligible, 2L)
    }
    broken <- native
    broken$id2 <- rev(broken$id2)
    testthat::expect_true(all(grepl("Native resample IDs changed", score(broken, keys)$candidates$reason)))
    broken <- native
    broken$.predictions[2L] <- list(NULL)
    testthat::expect_identical(score(broken, keys)$candidates$status, rep("excluded", 2L))
    search$primary_metric <- yardstick::new_numeric_metric(
      function(...) data.frame(.estimate = NA_real_), "minimize")
    testthat::expect_true(all(grepl("finite estimate", score(native, keys)$candidates$reason)))
    search$primary_metric <- yardstick::rmse
  }
  testthat::expect_identical(results$grid, results$bayes)
})

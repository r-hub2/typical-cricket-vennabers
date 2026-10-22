testthat::test_that("binary-options-parity", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

fixture <- function(name) {
  path <- file.path("extdata", name)
  stopifnot(nzchar(path))
  read.csv(path, stringsAsFactors = FALSE)
}
calibration <- fixture("cross-calibration.csv")
prediction <- fixture("cross-predictions.csv")
fold_reference <- fixture("binary-options-folds.csv")
combined_reference <- fixture("binary-options-combined.csv")
case_ids <- unique(combined_reference$case_id)
stopifnot(length(case_ids) == 25L,
          sum(startsWith(case_ids, "random_")) == 20L,
          setequal(unique(combined_reference$precision), c("none", "0", "2", "-1")))

maximum_error <- 0
for (setting in unique(combined_reference$precision)) {
  precision <- if (setting == "none") NULL else as.numeric(setting)
  for (case_id in case_ids) {
    cal <- calibration[calibration$case_id == case_id, ]
    pred <- prediction[prediction$case_id == case_id, ]
    expected <- combined_reference[combined_reference$case_id == case_id &
                                     combined_reference$precision == setting, ]
    expected_folds <- fold_reference[fold_reference$case_id == case_id &
                                      fold_reference$precision == setting, ]
    identities <- unique(cal$fold)
    input <- matrix(NA_real_, nrow(expected), length(identities),
                    dimnames = list(NULL, identities))
    for (id in identities) input[, id] <- pred$probability[pred$fold == id]
    fit <- cross_venn_abers(cal$truth, cal$probability, cal$fold, precision)
    endpoints <- NULL
    for (loss in c("log", "brier")) {
      actual <- predict(fit, input, loss = loss, type = "details")
      assert_conditions(identical(actual$folds$row, expected_folds$row), identical(actual$folds$fold, expected_folds$fold),
          identical(actual$prediction, predict(fit, input, loss = loss)), identical(predict(fit, input, loss = loss,
              type = "class"), actual$prediction$probability >= 0.5), identical(predict(fit, input[, rev(identities)],
              loss = loss, type = "details"), actual), all(is.finite(as.matrix(actual$folds[3:5]))), all(actual$folds[3:5] >=
              0 & actual$folds[3:5] <= 1), all(actual$folds$lower <= actual$folds$upper), all(is.finite(actual$prediction$probability)),
          all(actual$prediction$probability >= 0 & actual$prediction$probability <= 1))
      if (!is.null(endpoints)) stopifnot(identical(actual$folds[3:4], endpoints))
      endpoints <- actual$folds[3:4]
      maximum_error <- max(maximum_error,
                           abs(actual$prediction$probability - expected[[loss]]),
                           abs(as.matrix(endpoints) - as.matrix(expected_folds[5:6])),
                           abs(actual$folds$probability - expected_folds[[loss]]))
      for (id in identities) {
        a <- cal[cal$fold == id, ]
        b <- expected_folds[expected_folds$fold == id, ]
        inductive <- predict(venn_abers(a$truth, a$probability, precision), input[, id], loss = loss)
        maximum_error <- max(maximum_error,
                             abs(as.matrix(inductive[1:2]) - as.matrix(b[5:6])),
                             abs(inductive$probability - b[[loss]]))
        stopifnot(identical(inductive$probability, actual$folds$probability[actual$folds$fold == id]))
      }
      # Reconstruct both aggregate choices solely from public fold diagnostics.
      pieces <- actual$folds
      reconstructed <- vapply(seq_len(nrow(input)), function(row) {
        f <- pieces[pieces$row == row, ]
        if (loss == "log") plogis(mean(log(f$upper)) - mean(log1p(-f$lower))) else mean(f$probability)
      }, 0)
      maximum_error <- max(maximum_error, abs(reconstructed - expected[[loss]]))
      if (maximum_error > 1e-10) stop("Binary option parity failed: ", case_id, "/", setting, "/", loss)
    }
  }
}
cat("Binary option reference maximum absolute error:", maximum_error, "\n")

  testthat::expect_true(TRUE)
})

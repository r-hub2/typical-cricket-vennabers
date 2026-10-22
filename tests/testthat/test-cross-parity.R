testthat::test_that("cross-parity", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

fixture <- function(name) {
  path <- file.path("extdata", paste0("cross-", name, ".csv"))
  stopifnot(nzchar(path))
  read.csv(path, stringsAsFactors = FALSE)
}
calibration <- fixture("calibration")
prediction <- fixture("predictions")
combined <- fixture("combined")
case_ids <- unique(calibration$case_id)
random_ids <- case_ids[startsWith(case_ids, "random_")]
random_calibration <- calibration[calibration$case_id %in% random_ids, ]
random_folds <- split(random_calibration,
                      paste(random_calibration$case_id, random_calibration$fold))
sizes <- vapply(random_folds, nrow, 0L)
class_counts <- vapply(random_folds, function(x) length(unique(x$truth)), 0L)
fold_counts <- vapply(split(random_calibration$fold, random_calibration$case_id),
                      function(x) length(unique(x)), 0L)
assert_conditions(length(case_ids) == 105L, length(random_ids) == 100L, identical(sort(unique(fold_counts)),
    2:10), identical(sort(unique(sizes)), 1:200), identical(sort(unique(class_counts)), 1:2), any(vapply(random_folds,
    function(x) anyDuplicated(x$probability) > 0, FALSE)), any(vapply(random_folds, function(x) nrow(x) >
    1L && !anyDuplicated(x$probability), FALSE)), identical(as.integer(table(combined$case_id[combined$case_id %in%
    random_ids])), rep(25L, 100)), identical(sort(unique(prediction$case_id)), sort(case_ids)), identical(sort(unique(combined$case_id)),
    sort(case_ids)))

maximum_endpoint_error <- maximum_probability_error <- 0
for (case_id in case_ids) {
  cal <- calibration[calibration$case_id == case_id, ]
  pred <- prediction[prediction$case_id == case_id, ]
  expected <- combined[combined$case_id == case_id, ]
  identities <- unique(cal$fold)
  input <- matrix(NA_real_, nrow(expected), length(identities),
                  dimnames = list(NULL, identities))
  lower <- upper <- input
  for (id in identities) {
    a <- cal[cal$fold == id, ]
    b <- pred[pred$fold == id, ]
    stopifnot(identical(a$row, seq_len(nrow(a))),
              identical(b$row, seq_len(nrow(expected))))
    input[, id] <- b$probability
    # The existing public IVAP boundary supplies the endpoints composed by
    # cross prediction; this does not freeze the cross object's storage layout.
    endpoints <- stats::predict(venn_abers(a$truth, a$probability), b$probability)
    assert_conditions(identical(dim(endpoints), c(nrow(expected), 3L)), identical(names(endpoints), c("lower",
        "upper", "probability")), all(vapply(endpoints, is.double, FALSE)), all(is.finite(as.matrix(endpoints))),
        all(endpoints >= 0 & endpoints <= 1))
    error <- max(abs(as.matrix(endpoints[c("lower", "upper")]) -
                       as.matrix(b[c("lower", "upper")])))
    if (error > 1e-10) stop("Fold endpoint parity failed: ", case_id, "/", id)
    maximum_endpoint_error <- max(maximum_endpoint_error, error)
    lower[, id] <- b$lower
    upper[, id] <- b$upper
  }
  stopifnot(identical(expected$row, seq_len(nrow(expected))))
  fit <- cross_venn_abers(cal$truth, cal$probability, cal$fold)
  actual <- stats::predict(fit, input)
  assert_conditions(is.data.frame(actual), identical(names(actual), "probability"), identical(dim(actual),
      c(nrow(expected), 1L)), is.double(actual$probability), all(is.finite(actual$probability)), all(actual >=
      0 & actual <= 1))
  error <- max(abs(actual$probability - expected$probability))
  if (error > 1e-10) stop("Cross probability parity failed: ", case_id)
  maximum_probability_error <- max(maximum_probability_error, error)

  if (startsWith(case_id, "random_")) {
    stopifnot(nrow(input) == 25L,
              all(vapply(seq_len(ncol(input))[-1], function(k) {
                !identical(input[, 1], input[, k])
              }, FALSE)))
  }
  # Every stored case also checks name matching and observation order.
  stopifnot(identical(stats::predict(fit, input[, rev(identities)]), actual),
            identical(stats::predict(fit, input[nrow(input):1, ])$probability,
                      rev(actual$probability)))
  if (case_id == "crossing_pair") {
    stopifnot(all(1 - exp(rowMeans(log1p(-lower))) > exp(rowMeans(log(upper)))))
  }
  if (case_id == "near_knots") {
    # CSV round-tripping must preserve the exact knot and its neighbours.
    assert_conditions(input[1, 1] < 0.5, input[2, 1] == 0.5, input[3, 1] > 0.5, input[3, 1] - 0.5 == .Machine$double.eps/2,
        cal$probability[2] == input[3, 1])
  }
  if (case_id == "unequal") {
    averaged <- rowMeans(upper / (1 - lower + upper))
    weights <- as.numeric(table(factor(cal$fold, levels = identities)))
    weighted <- stats::plogis(as.vector((log(upper) - log1p(-lower)) %*% weights / sum(weights)))
    pooled <- stats::predict(venn_abers(cal$truth, cal$probability), input[, 1])$probability
    for (alternative in list(averaged, weighted, pooled)) {
      stopifnot(max(abs(alternative - expected$probability)) > 1e-5)
    }
  }
}
cat("Cross reference maximum absolute errors: endpoints =", maximum_endpoint_error,
    "; probability =", maximum_probability_error, "\n")

  testthat::expect_true(TRUE)
})

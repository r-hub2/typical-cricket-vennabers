testthat::test_that("outer identity audit requires exact (id, id2, .row, .config) coverage", {
  for (pkg in c("tune", "yardstick", "rsample", "workflows", "parsnip", "withr", "tibble")) testthat::skip_if_not_installed(pkg)
  checker <- new.env()
  sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), checker)
  d <- withr::with_seed(1, data.frame(x = stats::rnorm(60), truth = factor(rep(c("a", "b"), 30))))
  outer <- withr::with_seed(2, rsample::vfold_cv(d, v = 2, repeats = 2))
  model <- workflows::workflow(truth ~ x, parsnip::logistic_reg())
  result <- withr::with_seed(3, tune::fit_resamples(model, outer, metrics = yardstick::metric_set(yardstick::mn_log_loss),
    control = tune::control_resamples(save_pred = TRUE, extract = function(x) nrow(x$pre$mold$predictors))))
  audited <- checker$check_outer_identities(result, outer, metrics = 1L)
  # Extracts come back aligned with outer rows, whatever order tune reports.
  testthat::expect_identical(unlist(audited$extracts), vapply(outer$splits, function(s) nrow(rsample::analysis(s)), 1L))
  testthat::expect_identical(nrow(audited$predictions), 2L * nrow(d))

  broken <- function(edit) {
    changed <- result
    changed$.predictions[[3L]] <- edit(changed$.predictions[[3L]])
    changed
  }
  expect_rejected <- function(changed, metrics = 1L) testthat::expect_error(checker$check_outer_identities(changed, outer, metrics), "is not TRUE")
  expect_rejected(broken(function(p) p[-1L, ]))
  expect_rejected(broken(function(p) rbind(p, p[1L, ])))
  expect_rejected(broken(function(p) transform(p, .row = .row %% nrow(d) + 1L)))
  expect_rejected(broken(function(p) transform(p, .config = "pre0_mod1_post0")))
  expect_rejected(result, metrics = 2L)
  flat <- result
  flat$id2 <- NULL
  expect_rejected(flat)
  noted <- result
  noted$.notes[[1L]] <- tibble::tibble(location = "x", type = "warning", note = "planted", trace = list(NULL))
  expect_rejected(noted)
  unextracted <- result
  unextracted$.extracts <- NULL
  testthat::expect_error(checker$check_outer_identities(unextracted, outer, 1L))
  missing_fit <- result
  missing_fit$.extracts[[3L]]$.extracts[1L] <- list(NULL)
  expect_rejected(missing_fit)
  wrong_extract_config <- result
  wrong_extract_config$.extracts[[3L]]$.config <- "pre0_mod1_post0"
  expect_rejected(wrong_extract_config)
})

testthat::test_that("Bayes regression workflows replay on outer workers and restore without learning", {
  check_tuning_replay("bayes", "regression")
})

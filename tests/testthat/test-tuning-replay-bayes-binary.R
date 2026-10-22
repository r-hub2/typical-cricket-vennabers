testthat::test_that("Bayes binary workflows replay on outer workers and restore without learning", {
  check_tuning_replay("bayes", "binary")
})

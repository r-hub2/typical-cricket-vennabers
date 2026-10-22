testthat::test_that("Bayes ovr workflows replay on outer workers and restore without learning", {
  check_tuning_replay("bayes", "ovr")
})

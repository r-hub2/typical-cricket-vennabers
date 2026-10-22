testthat::test_that("Bayes ovo workflows replay on outer workers and restore without learning", {
  check_tuning_replay("bayes", "ovo")
})

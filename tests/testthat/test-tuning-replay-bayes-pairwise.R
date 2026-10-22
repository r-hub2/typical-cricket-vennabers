testthat::test_that("Bayes pairwise workflows replay on outer workers and restore without learning", {
  check_tuning_replay("bayes", "pairwise")
})

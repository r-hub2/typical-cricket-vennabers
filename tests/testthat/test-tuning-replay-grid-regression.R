testthat::test_that("Grid regression workflows replay on outer workers and restore without learning", {
  check_tuning_replay("grid", "regression")
})

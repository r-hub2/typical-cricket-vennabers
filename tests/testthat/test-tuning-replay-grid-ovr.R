testthat::test_that("Grid ovr workflows replay on outer workers and restore without learning", {
  check_tuning_replay("grid", "ovr")
})

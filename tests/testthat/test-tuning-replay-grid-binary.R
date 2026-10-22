testthat::test_that("Grid binary workflows replay on outer workers and restore without learning", {
  check_tuning_replay("grid", "binary")
})

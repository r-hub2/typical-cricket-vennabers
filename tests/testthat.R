# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
library(testthat)
library(vennabers)

# Long process/worker and integration matrices remain available for local
# release verification. The maintainer approved the short CRAN profile.
filter <- NULL
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  long <- readLines("cran-long-tests.txt")
  stopifnot(!anyDuplicated(long), all(file.exists(file.path("testthat", long))))
  filter <- paste0("^(", paste(tools::file_path_sans_ext(sub("^test-", "", long)), collapse = "|"), ")$")
  cat("CRAN short profile: excluding", length(long), "long test files listed in cran-long-tests.txt\n")
}
test_check("vennabers", filter = filter, invert = !is.null(filter))

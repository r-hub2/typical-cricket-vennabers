# Source checkouts read man/ directly; R CMD check reads the installed pages.
man <- testthat::test_path("..", "..", "man")
pages <- if (dir.exists(man)) tools::Rd_db(dir = dirname(man)) else tools::Rd_db("vennabers")
tags <- function(rd) vapply(rd, attr, "", "Rd_tag")
section <- function(rd, tag) paste(unlist(rd[tags(rd) == tag]), collapse = "")
exports <- getNamespaceExports("vennabers")

testthat::test_that("every exported function page has examples and a described value", {
  documented <- Filter(function(rd) {
    aliases <- unlist(rd[tags(rd) == "\\alias"])
    "\\usage" %in% tags(rd) && !grepl("internal", section(rd, "\\keyword")) &&
      any(aliases %in% exports)
  }, pages)
  testthat::expect_true(all(c("venn_abers_details.Rd", "venn_abers_tuning.Rd") %in% names(documented)))
  for (page in names(documented)) {
    rd <- documented[[page]]
    testthat::expect_true("\\examples" %in% tags(rd), info = page)
    testthat::expect_true("\\value" %in% tags(rd), info = page)
    testthat::expect_false(grepl("opaque", section(rd, "\\value"), ignore.case = TRUE), info = page)
  }
})

testthat::test_that("calibrator values name their class and point to their methods", {
  for (page in c("venn_abers.Rd", "cross_venn_abers.Rd")) {
    value <- section(pages[[page]], "\\value")
    calibrator <- sub("\\.Rd$", "", page)
    testthat::expect_true(grepl(paste0("class \"", calibrator, "\""), value, fixed = TRUE), info = page)
    for (method in c("predict", "plot", "print")) {
      testthat::expect_true(grepl(paste0(method, ".", calibrator), value, fixed = TRUE), info = paste(page, method))
    }
  }
})

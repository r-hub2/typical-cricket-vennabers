testthat::test_that("the installed parity index names a statement kind and existing evidence for every capability", {
  index <- system.file("PARITY.md", package = "vennabers")
  testthat::expect_true(nzchar(index))
  lines <- readLines(index, encoding = "UTF-8")
  kinds <- c("N", "T", "R", "E", "D", "L")
  # Every kind is defined before the table uses it.
  for (kind in kinds) testthat::expect_true(any(startsWith(lines, paste0("- **", kind, " "))), info = kind)
  start <- which(lines == "## Capability index")
  end <- which(startsWith(lines, "## "))
  end <- min(end[end > start])
  rows <- grep("^\\| ", lines[seq.int(start + 1L, end - 1L)], value = TRUE)
  rows <- rows[-(1:2)]
  testthat::expect_gte(length(rows), 20L)
  cells <- lapply(strsplit(rows, " \\| "), trimws)
  testthat::expect_true(all(lengths(cells) == 3L))
  for (cell in cells) {
    used <- strsplit(cell[[2L]], ", ", fixed = TRUE)[[1L]]
    testthat::expect_true(length(used) > 0L && all(used %in% kinds), info = cell[[1L]])
    testthat::expect_match(cell[[3L]], "`", info = cell[[1L]])
  }
  # Every referenced test, fixture and installed example exists.
  refs <- unique(unlist(regmatches(lines, gregexpr("`[^`]+`", lines))))
  refs <- gsub("`", "", refs, fixed = TRUE)
  tests <- grep("^test-.*[.]R$", refs, value = TRUE)
  fixtures <- sub("^extdata/", "", grep("^extdata/", refs, value = TRUE))
  examples <- sub("^examples/", "", grep("^examples/.*[.]R$", refs, value = TRUE))
  # Lower bounds guard against a truncated index or a broken reference pattern.
  testthat::expect_gte(length(tests), 30L)
  testthat::expect_gte(length(fixtures), 5L)
  testthat::expect_gte(length(examples), 20L)
  for (file in tests) testthat::expect_true(file.exists(testthat::test_path(file)), info = file)
  for (file in fixtures) {
    testthat::expect_true(length(Sys.glob(testthat::test_path("extdata", file))) > 0L, info = file)
  }
  for (file in examples) {
    testthat::expect_true(nzchar(system.file("examples", file, package = "vennabers")), info = file)
  }
})

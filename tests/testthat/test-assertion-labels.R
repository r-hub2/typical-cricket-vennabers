testthat::test_that("assertion-labels", {
walk_calls <- function(x, calls = list()) {
  if (!is.recursive(x)) return(calls)
  if (is.call(x) && identical(as.character(x[[1L]]), "stopifnot")) {
    arguments <- as.list(x)[-1L]
    labels <- names(arguments)
    if (is.null(labels)) labels <- rep("", length(arguments))
    if (length(arguments) >= 5L) calls[[length(calls) + 1L]] <- labels
  }
  children <- as.list(x)
  for (i in seq_along(children)) if (length(children[[i]])) {
    calls <- walk_calls(children[[i]], calls)
  }
  calls
}

root <- Sys.getenv("TESTTHAT_WD")
examples <- file.path(root, "inst", "examples")
if (!dir.exists(examples)) examples <- system.file("examples", package = "vennabers")
files <- c(list.files(".", "^test-.*\\.R$", full.names = TRUE),
           list.files(examples, "\\.R$", full.names = TRUE))
for (file in files) {
  labels <- walk_calls(parse(file))
  if (length(labels)) {
    stop(file, " contains a multi-condition stopifnot().", call. = FALSE)
  }
}

source(if (file.exists(file.path(examples, "assertions.R"))) {
  file.path(examples, "assertions.R")
} else {
  system.file("examples", "assertions.R", package = "vennabers")
})

failure <- tryCatch({
  cases <- list(list(name = "first"))
  for (case in cases) for (precision in list(NULL)) for (iteration in 1L) {
    assert_conditions(TRUE, TRUE, TRUE, TRUE, FALSE)
  }
  NULL
}, error = identity)
if (!all(inherits(failure, "error"),
         grepl("FALSE", conditionMessage(failure), fixed = TRUE),
         grepl("case=cases[[1L]]", conditionMessage(failure), fixed = TRUE),
         grepl("iteration=1", conditionMessage(failure), fixed = TRUE),
         grepl("precision=NULL", conditionMessage(failure), fixed = TRUE))) {
  stop("The assertion diagnostic omitted a failed condition or loop value.", call. = FALSE)
}

  testthat::expect_true(TRUE)
})

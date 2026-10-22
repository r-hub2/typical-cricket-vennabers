testthat::test_that("one-hot class matrices represent existing class decisions", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
error <- function(expr, pattern) {
  result <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result), fixed = TRUE))
}
one_hot <- function(classes, index) {
  result <- matrix(0, length(index), length(classes), dimnames = list(NULL, classes))
  result[cbind(seq_along(index), index)] <- 1
  result
}

# Manual binary: the logical event decision (`>= .5`) maps to FALSE/TRUE columns,
# so an exact 0.5 tie lands in the event column TRUE.
binary <- venn_abers(c(FALSE, TRUE, FALSE, TRUE), c(.2, .2, .8, .8))
probability <- predict(binary, c(.2, .5, .8))$probability
decision <- predict(binary, c(.2, .5, .8), type = "class")
encoded <- venn_abers_one_hot(decision)
assert_conditions(probability[[2L]] == .5, identical(decision, c(FALSE, TRUE, TRUE)),
                  identical(encoded, one_hot(c("FALSE", "TRUE"), c(1L, 2L, 2L))),
                  is.double(encoded), is.null(rownames(encoded)))

# Manual multiclass: exact (1/3, 1/3, 1/3) ties keep the first stored level.
classes <- c("cat", "dog", "bird")
p <- matrix(1 / 3, 6, 3, dimnames = list(NULL, classes))
truth <- factor(rep(classes, 2), levels = classes)
for (method in c("ovo", "ovr")) {
  fit <- venn_abers(truth, p, method = method)
  tied <- predict(fit, p[1:2, ])
  decision <- predict(fit, p[1:2, ], type = "class")
  assert_conditions(all(tied == tied[1L, 1L]),
                    identical(venn_abers_one_hot(decision), one_hot(classes, c(1L, 1L))))
}

# Zero and one rows keep every column; a one-class batch keeps all K columns.
empty <- venn_abers_one_hot(predict(binary, numeric(), type = "class"))
single <- venn_abers_one_hot(factor("dog", levels = classes))
batch <- venn_abers_one_hot(factor(c("dog", "dog"), levels = classes))
assert_conditions(identical(dim(empty), c(0L, 2L)), identical(colnames(empty), c("FALSE", "TRUE")),
                  identical(venn_abers_one_hot(predict(fit, p[0, , drop = FALSE], type = "class")),
                            one_hot(classes, integer())),
                  is.matrix(single), identical(single, one_hot(classes, 2L)),
                  identical(batch, one_hot(classes, c(2L, 2L))),
                  identical(venn_abers_one_hot(factor("only")), one_hot("only", 1L)))

# Missing decisions give missing rows in place; ordered levels only fix the column order.
missing <- venn_abers_one_hot(factor(c("bird", NA, "cat"), levels = classes))
expected <- one_hot(classes, c(3L, 1L, 1L)); expected[2L, ] <- NA_real_
logical_missing <- venn_abers_one_hot(c(NA, TRUE))
ordered <- venn_abers_one_hot(factor(c("high", "low"), levels = c("low", "mid", "high"), ordered = TRUE))
assert_conditions(identical(missing, expected), identical(dim(logical_missing), c(2L, 2L)),
                  all(is.na(logical_missing[1L, ])), identical(logical_missing[2L, ], c("FALSE" = 0, "TRUE" = 1)),
                  identical(ordered, one_hot(c("low", "mid", "high"), c(3L, 1L))))

# Anything other than an existing class decision is rejected.
error(venn_abers_one_hot(factor(character())), "at least one level")
error(venn_abers_one_hot(p), "`x` must be a factor")
error(venn_abers_one_hot(probability), "`x` must be a factor")
error(venn_abers_one_hot(c("cat", "dog")), "`x` must be a factor")
error(venn_abers_one_hot(matrix(TRUE, 2, 2)), "`x` must be a factor")
error(venn_abers_one_hot(list(factor("cat"))), "`x` must be a factor")
error(venn_abers_one_hot(NULL), "`x` must be a factor")
error(venn_abers_one_hot(data.frame(.pred_cat = .4, .pred_dog = .6)), "factor `.pred_class`")
error(venn_abers_one_hot(data.frame(.pred_class = "cat")), "factor `.pred_class`")

packages <- c("parsnip", "workflows", "rsample", "recipes", "tibble", "generics")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  # A null model gives exactly tied workflow probabilities on balanced data.
  spec <- parsnip::null_model(mode = "classification") |> parsnip::set_engine("parsnip")
  d <- data.frame(x = 1:12, truth = factor(rep(c("yes", "no"), 6), levels = c("yes", "no")))
  split <- rsample::make_splits(list(analysis = 1:8, assessment = 9:12), d)
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |> workflows::add_model(spec)
  for (event_level in c("first", "second")) {
    fitted <- generics::fit(calibrated_workflow(base, "truth", splits = split, event_level = event_level), d)
    new <- d[1:2, "x", drop = FALSE]
    prediction <- predict(fitted, new, type = "class")
    event <- if (event_level == "first") 1L else 2L
    assert_conditions(all(as.matrix(predict(fitted, new, type = "prob")) == .5),
                      identical(venn_abers_one_hot(prediction), one_hot(c("yes", "no"), c(event, event))),
                      identical(venn_abers_one_hot(tibble::add_column(prediction, id = 1:2)),
                                venn_abers_one_hot(prediction$.pred_class)))
  }
  dm <- data.frame(x = 1:18, truth = factor(rep(c("a", "b", "c"), 6)))
  splitm <- rsample::make_splits(list(analysis = 1:12, assessment = 13:18), dm)
  basem <- workflows::workflow() |> workflows::add_formula(truth ~ x) |> workflows::add_model(spec)
  for (classification in c("ovo", "ovr")) {
    fitted <- generics::fit(calibrated_workflow(basem, "truth", classification = classification,
                                                splits = splitm), dm)
    new <- dm[1:2, "x", drop = FALSE]
    probability <- as.matrix(predict(fitted, new, type = "prob"))
    assert_conditions(all(probability == probability[1L, 1L]),
                      identical(venn_abers_one_hot(predict(fitted, new, type = "class")),
                                one_hot(c("a", "b", "c"), c(1L, 1L))),
                      identical(venn_abers_one_hot(predict(fitted, dm[0, "x", drop = FALSE], type = "class")),
                                one_hot(c("a", "b", "c"), integer())))
  }

  # Schema-encoded decisions keep their codes; decoding first restores original names.
  schema <- venn_abers_schema(d, "truth")
  development <- venn_abers_encode(d, schema)
  formula <- reformulate(setdiff(names(development), schema$encoded_outcome), schema$encoded_outcome)
  encoded_base <- workflows::workflow() |> workflows::add_formula(formula) |> workflows::add_model(spec)
  model <- calibrated_workflow(encoded_base, schema$encoded_outcome, schema = schema,
                               splits = rsample::make_splits(list(analysis = 1:8, assessment = 9:12), development))
  fitted <- generics::fit(model, development)
  prediction <- predict(fitted, venn_abers_encode(d[1:3, "x", drop = FALSE], schema), type = "class")
  coded <- venn_abers_one_hot(prediction)
  decoded <- venn_abers_one_hot(venn_abers_decode(prediction, schema)$class)
  assert_conditions(identical(colnames(coded), unname(schema$classes)),
                    identical(colnames(decoded), c("yes", "no")),
                    identical(unname(coded), unname(decoded)))
}

  testthat::expect_true(TRUE)
})

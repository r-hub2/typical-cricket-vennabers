testthat::test_that("schema", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
error <- function(expr, pattern = "") {
  result <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
}
# Schema construction/encoding does not require optional integration packages.
d <- data.frame(x = c(.13, .81, .37), z = factor(c("a b", "a-b", "a b")),
                y = factor(c("class", ".pred_yes", "class"), levels = c("class", ".pred_yes")))
names(d) <- c("vaV0001", "ümlaut", ".pred_class")
s <- venn_abers_schema(d, ".pred_class", event_level = "second")
e <- venn_abers_encode(d, s)
assert_conditions(identical(s, venn_abers_schema(d[3:1, ], ".pred_class", event_level = "second")),
          identical(s$encoded_outcome, unname(s$columns[[".pred_class"]])),
          identical(e[[s$columns[["vaV0001"]]]], d[["vaV0001"]]),
          identical(e[[s$columns[["ümlaut"]]]], d[["ümlaut"]]),
          identical(levels(e[[s$encoded_outcome]]), unname(s$classes)),
          identical(s$event_class, ".pred_yes"))
reordered <- d[c(3, 2, 1)]
reordered[[s$outcome]] <- factor(reordered[[s$outcome]], levels = rev(levels(d[[s$outcome]])))
stopifnot(identical(venn_abers_encode(reordered, s), e),
          nrow(venn_abers_encode(d[0, ], s)) == 0L,
          !s$encoded_outcome %in% names(venn_abers_encode(d[1:2], s)))
error(venn_abers_encode(d[1], s), "missing")
error(venn_abers_encode(cbind(d, extra = 1), s), "Unknown")
for (bad_names in list(c("x", "x", "y"), c("", "z", "y"), c(NA, "z", "y"))) {
  bad <- d; names(bad) <- bad_names
  error(venn_abers_schema(bad, "y"), "name")
  error(venn_abers_encode(bad, s), "name")
}
bad <- d; bad[[s$outcome]] <- factor(c("new", "class", "class"))
error(venn_abers_encode(bad, s), "Unknown")
for (field in names(s)) {
  bad <- s; bad[[field]] <- NULL
  error(venn_abers_encode(d, bad), "schema")
}
bad <- s; bad$columns <- bad$columns[-1]
error(venn_abers_encode(d, bad), "schema")
bad <- s; bad$classes[2] <- bad$classes[1]
error(venn_abers_encode(d, bad), "schema")
# Decoding uses exact probability names, retaining dimensions and separate roles.
p <- setNames(data.frame(c(.1, .9, .3), c(.9, .1, .7)), paste0(".pred_", unname(s$classes)))
p$.pred_class <- factor(c(s$classes[2], s$classes[1], s$classes[2]), levels = unname(s$classes))
p[[s$encoded_outcome]] <- e[[s$encoded_outcome]]
p$.row <- c(3L, 1L, 2L); p$id <- c("Repeat1", "Repeat1", "Repeat2")
p$id2 <- c("Fold2", "Fold1", "Fold1"); p$.config <- rep("pre0_mod0_post0", 3)
for (rows in list(1:3, 1L, integer())) {
  native <- p[rows, rev(names(p)), drop = FALSE]
  decoded <- venn_abers_decode(native, s)
  assert_conditions(identical(colnames(decoded$probability), s$class_levels),
            identical(dim(decoded$probability), c(length(rows), 2L)),
            identical(unname(decoded$probability[, 1]), c(.1, .9, .3)[rows]),
            identical(decoded$class, factor(c(".pred_yes", "class", ".pred_yes")[rows], levels = s$class_levels)),
            identical(decoded$truth, d[[s$outcome]][rows]), identical(decoded$outcome, s$outcome),
            identical(decoded$identity, native[c(".config", "id2", "id", ".row")]))
}
only <- venn_abers_decode(p[1:2], s)
stopifnot(is.null(only$class), is.null(only$truth), is.null(only$details), is.null(only$identity))
error(venn_abers_decode(p[-1], s), "probabilit")
bad <- p; names(bad)[2] <- names(bad)[1]
error(venn_abers_decode(bad, s), "name")
bad <- p; names(bad)[1] <- ".pred_unknown"
error(venn_abers_decode(bad, s), "probabilit")
bad <- p; bad$.pred_class <- factor("vaL9999")
error(venn_abers_decode(bad, s), "Unknown")
bad <- p; attr(bad, "vennabers_schema") <- venn_abers_schema(d, s$outcome)
error(venn_abers_decode(bad, s), "different schema")
# A code-looking original is encoded by its own identity, even at another position.
shifted <- d[c(3, 2, 1)]
shifted[[2]] <- as.Date("2020-01-01") + 1:3
ss <- venn_abers_schema(shifted, ".pred_class")
se <- venn_abers_encode(shifted, ss)
stopifnot(ss$columns[["vaV0001"]] != "vaV0001", identical(se[[3]], d[[1]]),
          identical(se[[2]], shifted[[2]]), identical(rownames(se), rownames(shifted)))
for (value in c("", NA_character_)) {
  bad <- d; levels(bad[[s$outcome]]) <- c("class", value)
  error(venn_abers_schema(bad, s$outcome))
}
bad <- d; bad[[s$outcome]][1] <- NA
stopifnot(is.na(venn_abers_encode(bad, s)[[s$encoded_outcome]][1]))

# Multiclass retains the declared universe, including unobserved classes.
for (method in c("ovo", "ovr")) {
  multi <- d
  multi[[s$outcome]] <- factor(multi[[s$outcome]], levels = c(".pred_yes", "class", "vaL0001", "a b"))
  ms <- venn_abers_schema(multi, s$outcome, classification = method)
  me <- venn_abers_encode(multi, ms)
  stopifnot(identical(ms$class_levels, levels(multi[[s$outcome]])),
            identical(levels(me[[ms$encoded_outcome]]), unname(ms$classes)),
            is.null(ms$event_level), is.null(ms$event_class))
  error(venn_abers_schema(multi, s$outcome, classification = method, event_level = "first"), "event_level")
}
# Numeric 0/1 is regression; encoding never turns it into factor levels.
r <- d; r[[s$outcome]] <- c(0L, 1L, 0L)
rs <- venn_abers_schema(r, s$outcome)
re <- venn_abers_encode(r, rs)
assert_conditions(rs$mode == "regression", is.null(rs$class_levels), is.null(rs$classification),
          is.null(rs$event_class), identical(re[[rs$encoded_outcome]], r[[s$outcome]]))
error(venn_abers_schema(r, s$outcome, classification = "binary"), "Regression")
error(venn_abers_schema(r, s$outcome, event_level = "first"), "Regression")
# Decode only existing roles; numbers and diagnostic identities are unchanged.
for (method in c("ovo", "ovr")) {
  ms <- venn_abers_schema(multi, s$outcome, classification = method)
  native <- setNames(data.frame(.1, .2, .3, .4), paste0(".pred_", unname(ms$classes)))
  native$.pred_class <- factor(ms$classes[4], levels = rev(unname(ms$classes)))
  native[[ms$encoded_outcome]] <- factor(ms$classes[2], levels = unname(ms$classes))
  diagnostic <- data.frame(row = 1L, fold = "F1", method = method,
                           event_class = ms$classes[2], lower = .1, upper = .3, probability = .2)
  if (method == "ovo") {
    diagnostic$class_1 <- unname(ms$classes[1]); diagnostic$class_2 <- unname(ms$classes[2])
  } else diagnostic$class <- unname(ms$classes[2])
  for (n in 0:1) {
    decoded <- venn_abers_decode(list(prediction = native[seq_len(n), rev(names(native))],
                                     details = diagnostic[seq_len(n), ]), ms)
    assert_conditions(identical(dim(decoded$probability), c(n, 4L)),
              identical(colnames(decoded$probability), ms$class_levels),
              identical(levels(decoded$class), ms$class_levels),
              identical(levels(decoded$truth), ms$class_levels),
              identical(decoded$details[c("row", "fold", "method", "lower", "upper", "probability")],
                        diagnostic[seq_len(n), c("row", "fold", "method", "lower", "upper", "probability")]))
    for (cl in intersect(c("class_1", "class_2", "class", "event_class"), names(diagnostic))) {
      stopifnot(identical(decoded$details[[cl]], ms$class_levels[match(diagnostic[[cl]][seq_len(n)], ms$classes)]))
    }
  }
  error(venn_abers_decode(native[-1], ms), "probabilit")
  bad <- diagnostic; bad$event_class <- "unknown"
  error(venn_abers_decode(list(prediction = native, details = bad), ms), "Unknown")
}
for (outcome in c(".pred", ".pred_class", ".pred_vaL0001", ".row", ".config", "id", "id2", "a b")) {
  names(r)[3] <- outcome
  rs <- venn_abers_schema(r, outcome)
  native <- data.frame(.pred = c(1.25, -2), id = c("R1", "R1"), id2 = c("F1", "F2"),
                       .row = 1:2, .config = "pre0_mod0_post0")
  native[[rs$encoded_outcome]] <- c(0L, 1L)
  for (n in 0:2) {
    x <- native[seq_len(n), ]
    diagnostic <- data.frame(row = seq_len(n), fold = rep("F1", n),
                             prediction = x$.pred, lower = x$.pred - 1, upper = x$.pred + 1)
    aggregate <- diagnostic[c("row", "lower", "upper")]
    decoded <- venn_abers_decode(list(prediction = x, details = diagnostic, aggregate = aggregate), rs)
    assert_conditions(identical(decoded$prediction, x$.pred), identical(decoded$truth, x[[rs$encoded_outcome]]),
              identical(decoded$outcome, outcome), identical(decoded$details, diagnostic),
              identical(decoded$aggregate, aggregate), is.null(decoded$probability), is.null(decoded$class),
              identical(decoded$identity, x[c("id", "id2", ".row", ".config")]))
  }
  only <- venn_abers_decode(native[".pred"], rs)
  stopifnot(identical(names(only), c("outcome", "prediction")))
  bad <- native; bad$.pred <- "1"
  error(venn_abers_decode(bad, rs), "numeric")
  bad <- native; bad$.pred_class <- factor("vaL0001")
  error(venn_abers_decode(bad, rs), "Unknown")
}
# Every new schema mode rejects incomplete mappings and foreign bindings.
for (schema in list(ms, rs)) {
  original <- if (schema$mode == "regression") r else multi
  encoded <- venn_abers_encode(original, schema)
  reordered <- original[rev(names(original))]
  if (schema$mode != "regression") reordered[[schema$outcome]] <- factor(reordered[[schema$outcome]], levels = rev(schema$class_levels))
  stopifnot(identical(venn_abers_encode(reordered, schema), encoded))
  for (field in names(schema)) {
    broken <- schema; broken[[field]] <- NULL
    error(venn_abers_encode(original, broken), "schema")
  }
  broken <- schema; broken$columns[2] <- broken$columns[1]
  error(venn_abers_encode(original, broken), "schema")
  error(venn_abers_encode(original[-1], schema), "missing")
  error(venn_abers_encode(cbind(original, unknown = 1), schema), "Unknown")
  foreign <- original; names(foreign)[1] <- "foreign predictor"
  foreign_schema <- if (schema$mode == "regression") venn_abers_schema(foreign, schema$outcome) else
    venn_abers_schema(foreign, schema$outcome, classification = schema$classification)
  attr(original, "vennabers_schema") <- foreign_schema
  error(venn_abers_encode(original, schema), "different schema")
}

  testthat::expect_true(TRUE)
})

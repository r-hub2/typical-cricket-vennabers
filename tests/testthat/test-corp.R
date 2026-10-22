testthat::test_that("CORP reliability diagrams and Brier score decompositions", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
error <- function(expr, pattern) {
  result <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result), fixed = TRUE))
}
path <- tempfile(fileext = ".pdf")
grDevices::pdf(path)
on.exit({ grDevices::dev.off(); unlink(path) })
layout <- graphics::par("mfrow")
columns <- c("probability", "estimate", "count", "band_lower", "band_upper",
             "valid_replicates", "band", "level", "replicates", "seed")
components <- c("brier", "mcb", "dsc", "unc")

# PAV pools the descending events 1, 1, 0, 0 into one block, the tie at 0.3
# included, and keeps one row per unique probability in increasing order.
truth <- c(0, 1, 1, 0, 1)
probability <- c(.5, .3, .1, .3, .7)
drawn <- withVisible(corp_plot(truth, probability, band = "none"))
table <- drawn$value
assert_conditions(!drawn$visible, is.data.frame(table), identical(names(table), columns),
                  identical(table$probability, c(.1, .3, .5, .7)),
                  identical(table$estimate, c(.5, .5, .5, 1)),
                  identical(table$count, c(1L, 2L, 1L, 1L)),
                  identical(rownames(table), as.character(1:4)),
                  all(is.na(table$band_lower)), is.double(table$band_lower),
                  all(is.na(table$band_upper)), is.double(table$band_upper),
                  all(is.na(table$valid_replicates)), is.integer(table$valid_replicates),
                  identical(table$band, rep("none", 4L)), identical(table$level, rep(.9, 4L)),
                  identical(table$replicates, rep(1000L, 4L)), identical(table$seed, rep(NA_integer_, 4L)),
                  identical(graphics::par("mfrow"), layout))

# A tie split by its events pools into one estimate; logical truth is the event.
tied <- corp_plot(c(FALSE, FALSE, TRUE, FALSE, TRUE), c(.2, .2, .6, .6, .9), band = "none",
                  level = .8, replicates = 10, seed = -3)
assert_conditions(identical(tied$probability, c(.2, .6, .9)), identical(tied$estimate, c(0, .5, 1)),
                  identical(tied$count, c(2L, 2L, 1L)), identical(tied$level, rep(.8, 3L)),
                  identical(tied$replicates, rep(10L, 3L)), identical(tied$seed, rep(-3L, 3L)))
single <- corp_plot(1, .4, band = "none")
assert_conditions(identical(single$estimate, 1), identical(single$count, 1L))

# The decomposition by hand: Brier 0.122 against 0.1 for the PAV estimate
# (0, 0, 0.5, 0.5, 1) and 0.24 for the event rate 0.4.
decomposition <- withVisible(corp_decomposition_plot(c(0, 0, 1, 0, 1), c(.2, .2, .6, .6, .9)))
parts <- decomposition$value
assert_conditions(!decomposition$visible, is.data.frame(parts), identical(names(parts), components),
                  nrow(parts) == 1L,
                  isTRUE(all.equal(unlist(parts), c(brier = .122, mcb = .022, dsc = .14, unc = .24))),
                  isTRUE(all.equal(parts$brier, brier_score(c(0, 0, 1, 0, 1), c(.2, .2, .6, .6, .9)))))

# The identity holds up to rounding on continuous and discrete probabilities.
set.seed(20261020)
for (values in list(stats::runif(200), round(stats::runif(200), 1))) {
  events <- stats::rbinom(200, 1, values^2)
  parts <- corp_decomposition_plot(events, values)
  curve <- corp_plot(events, values, band = "none")
  assert_conditions(abs(parts$brier - (parts$mcb - parts$dsc + parts$unc)) < 1e-12,
                    parts$mcb > -1e-12, parts$dsc > -1e-12,
                    identical(curve$probability, sort(unique(values))),
                    sum(curve$count) == 200L, !is.unsorted(curve$estimate))
}

# Multiclass: every factor, even with two levels, evaluates each class
# against the rest on normalized rows, one panel per class in its own grid.
classes <- c("a", "b", "c")
truth <- factor(c("a", "b", "c", "a", "b", "c", "a", "a"), levels = classes)
probability <- matrix(c(.6, .2, .2, .2, .5, .3, .1, .3, .6, .4, .4, .2,
                        .3, .3, .4, .2, .2, .6, .5, .25, .25, .6, .2, .2),
                      ncol = 3, byrow = TRUE, dimnames = list(NULL, classes))
curves <- corp_plot(truth, probability, band = "none", seed = 7)
assert_conditions(identical(names(curves), c("class", columns)),
                  identical(levels(curves$class), classes), identical(unique(as.character(curves$class)), classes),
                  identical(rownames(curves), as.character(seq_len(nrow(curves)))),
                  identical(curves$seed, rep(7L, nrow(curves))),
                  identical(graphics::par("mfrow"), layout))
for (class in classes) {
  own <- curves[curves$class == class, columns]
  rownames(own) <- NULL
  assert_conditions(identical(own, corp_plot(truth == class, probability[, class], band = "none", seed = 7)))
}
parts <- corp_decomposition_plot(truth, probability)
assert_conditions(identical(names(parts), c("class", components)),
                  identical(parts$class, factor(classes, levels = classes)),
                  isTRUE(all.equal(sum(parts$brier), brier_score(truth, probability))),
                  all(abs(parts$brier - (parts$mcb - parts$dsc + parts$unc)) < 1e-12),
                  identical(graphics::par("mfrow"), layout))
for (class in classes) {
  own <- parts[parts$class == class, components]
  rownames(own) <- NULL
  assert_conditions(identical(own, corp_decomposition_plot(truth == class, probability[, class])))
}
two <- factor(c("no", "yes", "yes", "no"))
binary <- cbind(no = c(.8, .3, .4, .6), yes = c(.2, .7, .6, .4))
assert_conditions(identical(unique(as.character(corp_plot(two, binary, band = "none")$class)), c("no", "yes")),
                  identical(as.character(corp_decomposition_plot(two, binary)$class), c("no", "yes")))

# Graphical options override the defaults: a single panel leaves its
# coordinates, so `xlim` shows through `usr`.
corp_plot(c(0, 1, 1), c(.2, .5, .9), band = "none", xlim = c(.2, .8), main = "Own title")
assert_conditions(identical(graphics::par("usr")[1:2], c(.2, .8)))
corp_decomposition_plot(truth, probability, ylim = c(0, 2))
assert_conditions(identical(graphics::par("usr")[4], 2))

# A failing draw restores the grid as well.
error(corp_plot(truth, probability, band = "none", ylim = "bad"), "ylim")
assert_conditions(identical(graphics::par("mfrow"), layout))

# Inputs follow the scoring contract and its messages, without `na.rm`.
for (fun in list(function(...) corp_plot(..., band = "none"), corp_decomposition_plot)) {
  error(fun(c(0, 2), c(.1, .2)), "`truth` must be a logical or numeric vector containing only 0 and 1.")
  error(fun(c(0, NA), c(.1, .2)), "`truth` must be a logical or numeric vector containing only 0 and 1.")
  error(fun(numeric(), numeric()), "`truth` must be non-empty.")
  error(fun(c(0, 1), c(.1, 1.2)), "`probability` must be a numeric vector of finite values between 0 and 1.")
  error(fun(c(0, 1), c(.1, NA)), "`probability` must be a numeric vector of finite values between 0 and 1.")
  error(fun(c(0, 1), .1), "`probability` must have the same length as `truth`.")
  error(fun(factor(c("a", NA)), probability[1:2, ]), "`truth` must be a non-empty factor")
  error(fun(truth, unname(probability)), "`probability` column names must be non-empty, unique, and match the classes.")
  error(fun(truth, probability[-1L, ]), "`probability` must have one row per `truth` value.")
  error(fun(truth, probability * 2), "`probability` must be a numeric matrix of finite values between 0 and 1.")
}

# Every control is validated, also without a band.
for (bad in list("bootstrap", c("none", "confidence"), NA_character_, 1)) {
  error(corp_plot(c(0, 1), c(.1, .2), band = bad), "`band` must be one of: consistency, confidence, none.")
}
for (bad in list(0, 1, -.1, NA_real_, c(.8, .9), "0.9", Inf)) {
  error(corp_plot(c(0, 1), c(.1, .2), band = "none", level = bad),
        "`level` must be one number strictly between 0 and 1.")
}
for (bad in list(0, 1.5, -1, NA_real_, c(1000, 2000), "1000", Inf)) {
  error(corp_plot(c(0, 1), c(.1, .2), band = "none", replicates = bad),
        "`replicates` must be one finite whole number of at least 1.")
}
for (case in list(c(19, .9), c(1000, .999), c(1, .5))) {
  error(corp_plot(c(0, 1), c(.1, .2), band = "none", replicates = case[1L], level = case[2L]),
        "`replicates * (1 - level) / 2` must be at least 1; increase `replicates` or lower `level`.")
}
# Exactly one replicate beyond each quantile is enough.
assert_conditions(nrow(corp_plot(c(0, 1), c(.1, .2), band = "none", replicates = 20)) == 2L,
                  nrow(corp_plot(c(0, 1), c(.1, .2), band = "none", replicates = 2000, level = .999)) == 2L)
for (bad in list(1.5, NA_real_, c(1, 2), "1", .Machine$integer.max + 1, -.Machine$integer.max - 1)) {
  error(corp_plot(c(0, 1), c(.1, .2), band = "none", seed = bad),
        "`seed` must be NULL or one integer in R's seed range.")
}

# Bands: limits in [0, 1] around the estimate, one row per unique
# probability, resting on at most `replicates` replicates each.
warned <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, messages = messages)
}
bands <- c("band_lower", "band_upper", "valid_replicates")
set.seed(20261021)
values <- stats::runif(60)
events <- stats::rbinom(60, 1, values^2)
for (band in c("consistency", "confidence")) {
  for (x in list(values, round(values, 1))) {
    drawn <- withVisible(corp_plot(events, x, band = band, replicates = 199, seed = 11))
    curve <- drawn$value
    assert_conditions(!drawn$visible, identical(names(curve), columns),
                      identical(curve$probability, sort(unique(x))),
                      identical(curve$estimate, corp_plot(events, x, band = "none")$estimate),
                      is.double(curve$band_lower), is.integer(curve$valid_replicates),
                      !anyNA(curve$band_lower), !anyNA(curve$band_upper),
                      all(curve$band_lower >= 0), all(curve$band_upper <= 1),
                      all(curve$band_lower <= curve$band_upper),
                      all(curve$valid_replicates >= 1L), all(curve$valid_replicates <= 199L),
                      identical(curve$band, rep(band, nrow(curve))),
                      identical(curve$replicates, rep(199L, nrow(curve))),
                      identical(curve$seed, rep(11L, nrow(curve))))
  }
}
# The consistency band surrounds the diagonal; ties pool in every replicate,
# so a probability drawn in each replicate has all of them valid.
curve <- corp_plot(events, round(values, 1), replicates = 199, seed = 11)
inner <- curve$probability > .1 & curve$probability < .9
assert_conditions(all(curve$band_lower[inner] <= curve$probability[inner]),
                  all(curve$band_upper[inner] >= curve$probability[inner]),
                  all(curve$valid_replicates[inner] == 199L))

# The defaults: a consistency band at 0.9 from 1000 replicates.
curve <- corp_plot(events, values, seed = 3)
assert_conditions(identical(curve$band, rep("consistency", 60L)), identical(curve$level, rep(.9, 60L)),
                  identical(curve$replicates, rep(1000L, 60L)), !anyNA(curve$band_lower),
                  max(curve$valid_replicates) == 1000L)

# Inner probabilities with a confidence estimate of exactly 0 or 1 take limits
# interpolated linearly from the smallest and largest probability.
probabilities <- (1:6) / 7
curve <- corp_plot(c(0, 0, 0, 1, 1, 1), probabilities, band = "confidence", replicates = 199, seed = 5)
inner <- 2:5
assert_conditions(identical(curve$estimate, c(0, 0, 0, 1, 1, 1)),
                  isTRUE(all.equal(curve$band_lower[inner], stats::approx(probabilities[c(1, 6)],
                    curve$band_lower[c(1, 6)], xout = probabilities[inner])$y)),
                  isTRUE(all.equal(curve$band_upper[inner], stats::approx(probabilities[c(1, 6)],
                    curve$band_upper[c(1, 6)], xout = probabilities[inner])$y)))

# A single outcome omits the confidence band with a warning; the consistency
# band depends on the probabilities only.
omitted <- warned(corp_plot(c(0, 0, 0), c(.2, .5, .7), band = "confidence", replicates = 199, seed = 1))
assert_conditions(identical(omitted$messages, "Confidence band omitted: `truth` has only one outcome."),
                  all(is.na(omitted$value$band_lower)), all(is.na(omitted$value$band_upper)),
                  all(is.na(omitted$value$valid_replicates)),
                  identical(omitted$value$band, rep("confidence", 3L)))
kept <- warned(corp_plot(c(1, 1, 1), c(.2, .5, .7), replicates = 199, seed = 1))
assert_conditions(length(kept$messages) == 0L, !anyNA(kept$value$band_lower))
single <- corp_plot(1, .4, replicates = 199, seed = 1)
assert_conditions(identical(single$valid_replicates, 199L), single$band_lower <= .4, single$band_upper >= .4)
lonely <- factor(c("a", "a", "b", "b"), levels = classes)
shares <- cbind(a = c(.6, .5, .3, .2), b = c(.3, .3, .5, .6), c = c(.1, .2, .2, .2))
omitted <- warned(corp_plot(lonely, shares, band = "confidence", replicates = 199, seed = 1))
assert_conditions(identical(omitted$messages,
                            "Confidence band omitted: `truth` has only one outcome for class \"c\"."),
                  all(is.na(omitted$value$band_lower[omitted$value$class == "c"])),
                  !anyNA(omitted$value$band_lower[omitted$value$class != "c"]))

# A seed reproduces the band and leaves the caller's random state and RNG
# kind as they were; one seed covers all classes in turn, so the first class
# matches its binary call.
set.seed(99)
state <- get(".Random.seed", envir = globalenv())
kind <- RNGkind()
first <- corp_plot(events, values, replicates = 199, seed = 4)
assert_conditions(identical(get(".Random.seed", envir = globalenv()), state), identical(RNGkind(), kind),
                  identical(corp_plot(events, values, replicates = 199, seed = 4), first),
                  !identical(corp_plot(events, values, replicates = 199, seed = 5)[bands], first[bands]))
saved <- get(".Random.seed", envir = globalenv())
rm(".Random.seed", envir = globalenv())
corp_plot(events, values, replicates = 199, seed = 4)
absent <- !exists(".Random.seed", envir = globalenv(), inherits = FALSE)
assign(".Random.seed", saved, envir = globalenv())
assert_conditions(absent)
curves <- corp_plot(truth, probability, replicates = 199, seed = 4)
own <- curves[curves$class == "a", columns]
rownames(own) <- NULL
assert_conditions(identical(own, corp_plot(truth == "a", probability[, "a"], replicates = 199, seed = 4)),
                  identical(graphics::par("mfrow"), layout))

# Without a seed, the band draws from the global stream and advances it, as
# the seeded call does from set.seed().
set.seed(4)
unseeded <- corp_plot(events, values, replicates = 199)
assert_conditions(!identical(get(".Random.seed", envir = globalenv()), state),
                  identical(unseeded$seed, rep(NA_integer_, 60L)),
                  identical(unseeded[bands], first[bands]))

  testthat::expect_true(TRUE)
})

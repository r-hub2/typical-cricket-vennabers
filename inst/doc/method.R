## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", error = FALSE)

## ----library------------------------------------------------------------------
library(vennabers)

## ----ivap-by-hand-------------------------------------------------------------
set.seed(48)
k <- 40
x <- rnorm(k)
truth <- runif(k) < plogis(2 * x)
score <- plogis(0.5 * x)
new_score <- plogis(0.5 * 0.3)

# Value at the new case of the isotonic fit that includes it with `label`.
isotonic_at_new <- function(score, label, new_score, new_label) {
  scores <- c(score, new_score)
  o <- order(scores)
  fit <- stats::isoreg(scores[o], c(label, new_label)[o])
  fit$yf[o == length(scores)]
}
p0 <- isotonic_at_new(score, truth, new_score, 0)
p1 <- isotonic_at_new(score, truth, new_score, 1)

calibrator <- venn_abers(truth, score)
package <- predict(calibrator, new_score)
rbind(by_hand = c(lower = p0, upper = p1), package = unlist(package[1:2]))
stopifnot(isTRUE(all.equal(c(p0, p1), c(package$lower, package$upper))))

## ----validity-simulation------------------------------------------------------
set.seed(2)
repetitions <- 2000
selected <- outcome <- numeric(repetitions)
for (r in seq_len(repetitions)) {
  x <- rnorm(k + 1)
  y <- runif(k + 1) < plogis(2 * x)
  s <- plogis(0.5 * x)
  pair <- predict(venn_abers(y[1:k], s[1:k]), s[k + 1])
  selected[r] <- if (y[k + 1]) pair$upper else pair$lower
  outcome[r] <- y[k + 1]
}
group <- cut(selected, 0:5 / 5, include.lowest = TRUE)
mean_p <- tapply(selected, group, mean)
data.frame(cases = tabulate(group, nlevels(group)),
           mean_selected = mean_p,
           event_frequency = tapply(outcome, group, mean),
           standard_error = sqrt(mean_p * (1 - mean_p) / tabulate(group, nlevels(group))))

## ----merge--------------------------------------------------------------------
merged <- c(log = p1 / (1 - p0 + p1), brier = p1 + (p0^2 - p1^2) / 2)
package <- c(log = predict(calibrator, new_score)$probability,
             brier = predict(calibrator, new_score, loss = "brier")$probability)
rbind(by_hand = merged, package = package)
stopifnot(isTRUE(all.equal(merged, package)))

predict(calibrator, min(score) / 2)

## ----cvap-by-hand-------------------------------------------------------------
set.seed(15)
data <- data.frame(x = rnorm(150))
data$truth <- runif(150) < plogis(data$x^3)
fold <- rep(1:3, length.out = 150)
new_x <- data.frame(x = c(-1, 0.2, 1.5))

cal_score <- numeric(nrow(data))
new_score <- matrix(NA_real_, nrow(new_x), 3, dimnames = list(NULL, 1:3))
for (j in 1:3) {
  model <- glm(truth ~ x, data = data[fold != j, ], family = binomial())
  cal_score[fold == j] <- predict(model, data[fold == j, ], type = "response")
  new_score[, j] <- predict(model, new_x, type = "response")
}
cross <- cross_venn_abers(data$truth, cal_score, fold)

pairs <- lapply(1:3, function(j) {
  predict(venn_abers(data$truth[fold == j], cal_score[fold == j]),
          new_score[, j])
})
lower <- sapply(pairs, `[[`, "lower")
upper <- sapply(pairs, `[[`, "upper")
gm <- function(p) exp(rowMeans(log(p)))
by_hand <- cbind(log = gm(upper) / (gm(1 - lower) + gm(upper)),
                 brier = rowMeans(upper + (lower^2 - upper^2) / 2))
package <- cbind(log = predict(cross, new_score)$probability,
                 brier = predict(cross, new_score, loss = "brier")$probability)
by_hand
stopifnot(isTRUE(all.equal(by_hand, package)))

## ----multiclass-by-hand-------------------------------------------------------
set.seed(3)
classes <- c("a", "b", "c")
softmax <- function(z) exp(z) / rowSums(exp(z))
x <- rnorm(300)
truth <- factor(classes[apply(softmax(cbind(2 * x, -2 * x, 0)), 1,
                             function(p) sample(3, 1, prob = p))],
                levels = classes)
p <- softmax(cbind(a = x, b = -x, c = 0))
new_p <- softmax(cbind(a = c(-1, 0, 1), b = c(1, 0, -1), c = 0))

ovo <- venn_abers(truth, p)
pairs <- predict(ovo, new_p, type = "details")$pairs

rows <- truth %in% c("a", "b")
ab <- venn_abers(truth[rows] == "b", p[rows, "b"] / rowSums(p[rows, c("a", "b")]))
ab_new <- predict(ab, new_p[, "b"] / rowSums(new_p[, c("a", "b")]))
stopifnot(isTRUE(all.equal(ab_new$probability,
                           pairs$probability[pairs$non_event == "a" &
                                               pairs$event == "b"])))

r <- function(i, j) {
  forward <- pairs$non_event == j & pairs$event == i
  if (any(forward)) pairs$probability[forward] else
    1 - pairs$probability[pairs$non_event == i & pairs$event == j]
}
q <- sapply(classes, function(i) {
  1 / (Reduce(`+`, lapply(setdiff(classes, i), function(j) 1 / r(i, j))) -
         (length(classes) - 2))
})
coupled <- q / rowSums(q)
coupled
stopifnot(isTRUE(all.equal(coupled, predict(ovo, new_p),
                           check.attributes = FALSE)))

## ----ivar-by-hand-------------------------------------------------------------
set.seed(4)
k <- 25
prediction <- rnorm(k)
outcome <- 1.5 * prediction + rnorm(k)
new_prediction <- 0.4

m <- 1
sorted <- sort(outcome)
a <- sorted[m]
b <- sorted[k - m + 1]
lower_outcome <- pmin(pmax(outcome, a), sorted[k - m])
upper_outcome <- pmin(pmax(outcome, sorted[m + 1]), b)
L <- isotonic_at_new(prediction, lower_outcome, new_prediction, a)
U <- isotonic_at_new(prediction, upper_outcome, new_prediction, b)
point <- ((L - a) * L + (b - U) * U) / ((L - a) + (b - U))

by_hand <- c(prediction = point, lower = L, upper = U)
package <- unlist(predict(venn_abers_regression(outcome, prediction, m = m),
                          new_prediction))
rbind(by_hand, package)
stopifnot(isTRUE(all.equal(by_hand, package)))


## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", error = FALSE,
                      fig.width = 7, fig.height = 7)

## ----library------------------------------------------------------------------
library(vennabers)

## ----data---------------------------------------------------------------------
simulate <- function(event_probability, seed) {
  set.seed(seed)
  sizes <- c(training = 500, calibration = 300, assessment = 2000)
  x <- rnorm(sum(sizes))
  data <- data.frame(x = x, truth = runif(sum(sizes)) < event_probability(x))
  split(data, factor(rep(names(sizes), sizes), levels = names(sizes)))
}
data <- simulate(function(x) plogis(x^3), seed = 49)
sapply(data, nrow)

## ----calibrate----------------------------------------------------------------
calibrate <- function(data) {
  model <- glm(truth ~ x, data = data$training, family = binomial())
  score <- predict(model, data$calibration, type = "response")
  new_score <- predict(model, data$assessment, type = "response")
  truth <- data$calibration$truth

  platt <- glm(truth ~ qlogis(score), family = binomial())
  platt_probability <- plogis(coef(platt)[[1]] +
                                coef(platt)[[2]] * qlogis(new_score))

  o <- order(score)
  isotonic <- stats::isoreg(score[o], truth[o])
  isotonic_map <- stats::approxfun(score[o], isotonic$yf, method = "constant",
                                   rule = 2, ties = mean)

  venn_abers_probability <- predict(venn_abers(truth, score),
                                    new_score)$probability

  list(truth = data$assessment$truth,
       probability = list(Uncalibrated = new_score,
                          Platt = platt_probability,
                          Isotonic = isotonic_map(new_score),
                          `Venn-Abers` = venn_abers_probability))
}
result <- calibrate(data)

## ----scores-------------------------------------------------------------------
score_methods <- function(result) {
  t(sapply(result$probability, function(probability) {
    c(log_loss = log_loss(result$truth, probability),
      brier_score = brier_score(result$truth, probability),
      at_0_or_1 = sum(probability %in% c(0, 1)))
  }))
}
scores <- score_methods(result)
knitr::kable(scores, digits = 4)

## ----corp, fig.cap="CORP reliability diagrams on the assessment set."---------
old <- par(mfrow = c(2, 2))
for (method in names(result$probability)) {
  corp_plot(result$truth, result$probability[[method]], seed = 1, main = method)
}
par(old)

## ----extremes-----------------------------------------------------------------
isotonic <- result$probability$Isotonic
extreme <- isotonic %in% c(0, 1)
table(probability = isotonic[extreme], event = result$truth[extreme])

## ----contrast-----------------------------------------------------------------
contrast <- score_methods(calibrate(simulate(function(x) plogis(0.3 + 1.5 * x),
                                             seed = 49)))
knitr::kable(contrast, digits = 4)


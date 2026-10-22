## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", error = FALSE,
                     fig.width = 6, fig.height = 4)

## ----data---------------------------------------------------------------------
library(vennabers)
set.seed(17)
data <- data.frame(x = rnorm(900))
data$truth <- runif(nrow(data)) < plogis(0.3 + 1.5 * data$x)
training <- data[1:300, ]
calibration <- data[301:600, ]
assessment <- data[601:900, ]

model <- glm(truth ~ x, data = training, family = binomial())
calibration_p <- predict(model, calibration, type = "response")
assessment_p <- predict(model, assessment, type = "response")
calibrator <- venn_abers(calibration$truth, calibration_p)
probability <- predict(calibrator, assessment_p)$probability

## ----fit-map, fig.height=6, fig.cap="Calibration map and Venn-Abers width of an Inductive calibrator."----
map <- plot(calibrator, which = c("map", "width"))
head(map)

## ----fit-cross, fig.cap="Fold maps of a Cross calibrator."--------------------
fold <- rep(c("A", "B", "C"), length.out = nrow(calibration))
cross <- cross_venn_abers(calibration$truth, calibration_p, fold)
fold_maps <- plot(cross)
table(fold_maps$fold)

## ----fit-regression, fig.cap="Calibration map of a regression calibrator."----
set.seed(52)
x <- runif(600, 0, 10)
regression_data <- data.frame(x = x, truth = x + sin(x) + rnorm(600, sd = 0.5))
regression_model <- lm(truth ~ x, data = regression_data[1:200, ])
uncalibrated <- predict(regression_model, regression_data)
regression <- venn_abers_regression(regression_data$truth[201:400],
                                    uncalibrated[201:400])
regression_map <- plot(regression)

## ----corp, fig.cap="CORP reliability diagram with the consistency band."------
corp <- corp_plot(assessment$truth, probability, seed = 1)
head(corp)

## ----corp-confidence, fig.cap="CORP reliability diagram with the approximate confidence band."----
confidence <- corp_plot(assessment$truth, probability,
                        band = "confidence", seed = 1)

## ----decomposition, fig.cap="CORP decomposition of the Brier score."----------
parts <- corp_decomposition_plot(assessment$truth, probability)
parts
brier_score(assessment$truth, probability)

## ----regression-plot, fig.height=6, fig.cap="Observed values and residuals against calibrated regression points."----
held_out <- 401:600
point <- predict(regression, uncalibrated[held_out])
held_out_pairs <- regression_plot(regression_data$truth[held_out],
                                  point$prediction,
                                  which = c("observed", "residual"))
head(held_out_pairs)


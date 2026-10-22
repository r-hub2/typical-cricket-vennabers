## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", error = FALSE,
                     fig.width = 6, fig.height = 4)

## ----binary-data--------------------------------------------------------------
library(vennabers)
set.seed(12)
data <- data.frame(x = rnorm(240))
data$truth <- runif(nrow(data)) < plogis(data$x)
training <- data[1:120, ]
calibration <- data[121:200, ]
assessment <- data[201:240, ]

model <- glm(truth ~ x, data = training, family = binomial())
calibration_p <- predict(model, calibration, type = "response")
new_p <- predict(model, assessment["x"], type = "response")
calibrator <- venn_abers(calibration$truth, calibration_p)
result <- predict(calibrator, new_p)
head(result)

## ----binary-options-----------------------------------------------------------
head(predict(calibrator, new_p, loss = "brier"))
head(predict(calibrator, new_p, type = "class"))

## ----binary-scores------------------------------------------------------------
rbind(
  original = c(brier = brier_score(assessment$truth, new_p),
               log_loss = log_loss(assessment$truth, new_p)),
  calibrated = c(brier = brier_score(assessment$truth, result$probability),
                 log_loss = log_loss(assessment$truth, result$probability))
)

## ----binary-plot, fig.cap="Calibration on the separate assessment set (five equal-width bins)."----
calibration_plot(assessment$truth, result$probability, bins = 5)

## ----binary-cross-fit---------------------------------------------------------
development <- data[1:200, ]
set.seed(31)
fold <- sample(rep(c("A", "B", "C", "D"), length.out = nrow(development)))
models <- setNames(vector("list", 4), c("A", "B", "C", "D"))
oof_p <- numeric(nrow(development))
for (id in names(models)) {
  held <- fold == id
  models[[id]] <- glm(truth ~ x, data = development[!held, ],
                      family = binomial())
  oof_p[held] <- predict(models[[id]], development[held, ], type = "response")
}
cross <- cross_venn_abers(development$truth, oof_p, fold)

# Columns identify models; rows are the same new observations in every column.
fold_p <- vapply(models, predict, numeric(nrow(assessment)),
                 newdata = assessment["x"], type = "response")
head(predict(cross, fold_p))
cross_details <- predict(cross, fold_p, type = "details")
head(cross_details$folds)

## ----multiclass-input---------------------------------------------------------
classes <- c("cat", "dog", "bird")
truth <- factor(rep(classes, 3), levels = classes)
p <- rbind(c(.7, .2, .1), c(.2, .6, .2), c(.1, .2, .7),
           c(.5, .3, .2), c(.3, .5, .2), c(.2, .3, .5),
           c(.6, .3, .1), c(.1, .7, .2), c(.2, .2, .6))
colnames(p) <- classes
new_p <- rbind(c(.6, .3, .1), c(.1, .4, .5), c(.2, .6, .2))
colnames(new_p) <- classes
new_truth <- factor(c("cat", "bird", "dog"), levels = classes)

## ----multiclass-fit-----------------------------------------------------------
ovo <- venn_abers(truth, p, method = "ovo")
ovr <- venn_abers(truth, p, method = "ovr")
predict(ovo, new_p)
predict(ovr, new_p)
predict(ovo, new_p, type = "class")
head(predict(ovo, new_p, type = "details")$pairs)
head(predict(ovr, new_p, type = "details")$classes)
c(brier = brier_score(new_truth, predict(ovo, new_p)),
  log_loss = log_loss(new_truth, predict(ovo, new_p)))

## ----multiclass-cross---------------------------------------------------------
multi_fold <- rep(c("A", "B"), c(6, 3))
fold_new <- list(A = new_p, B = (new_p + 1 / 3) / 2)
for (method in c("ovo", "ovr")) {
  fitted <- cross_venn_abers(truth, p, multi_fold, method = method)
  print(predict(fitted, fold_new))
}

## ----regression-inductive-----------------------------------------------------
set.seed(52)
regression <- data.frame(x = rnorm(240))
regression$truth <- 2 * regression$x + rnorm(240, sd = .7)
reg_model <- lm(truth ~ x, data = regression[1:120, ])
reg_cal <- regression[121:200, ]
reg_test <- regression[201:240, ]
ivar <- venn_abers_regression(reg_cal$truth, predict(reg_model, reg_cal), m = 1)
reg_result <- predict(ivar, predict(reg_model, reg_test["x"]))
head(reg_result)
sqrt(mean((reg_test$truth - reg_result$prediction)^2))

## ----regression-cross---------------------------------------------------------
reg_dev <- regression[1:200, ]
reg_models <- setNames(vector("list", 4), names(models))
oof_numeric <- numeric(nrow(reg_dev))
for (id in names(reg_models)) {
  held <- fold == id
  reg_models[[id]] <- lm(truth ~ x, data = reg_dev[!held, ])
  oof_numeric[held] <- predict(reg_models[[id]], reg_dev[held, ])
}
cvar <- cross_venn_abers_regression(reg_dev$truth, oof_numeric, fold, m = 1)
new_numeric <- vapply(reg_models, predict, numeric(nrow(reg_test)),
                      newdata = reg_test["x"])
reg_details <- predict(cvar, new_numeric, type = "details")
head(reg_details$prediction)
head(reg_details$folds$A)


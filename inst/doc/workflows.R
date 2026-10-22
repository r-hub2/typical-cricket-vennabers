## ----setup, include=FALSE-----------------------------------------------------
# The code needs optional packages from Suggests; without them it is shown
# but not run.
required <- c("parsnip", "workflows", "hardhat", "rsample", "recipes", "rlang",
              "tibble", "withr", "generics", "tidyselect", "rpart")
unavailable <- required[!vapply(required, requireNamespace, logical(1),
                                quietly = TRUE)]
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", error = FALSE,
                      eval = length(unavailable) == 0)

## ----missing-packages, eval=length(unavailable) > 0, echo=FALSE, results="asis"----
# cat("> **Note:** The code in this vignette was not run because these",
#     "suggested packages are not installed:",
#     paste0("`", unavailable, "`", collapse = ", "), "\n")

## ----data---------------------------------------------------------------------
library(vennabers)
set.seed(12)
data <- data.frame(x = rnorm(240))
data$truth <- factor(ifelse(runif(240) < plogis(data$x), "yes", "no"),
                     levels = c("yes", "no"))
development <- data[1:200, ]
assessment <- data[201:240, ]
new_data <- assessment["x"]

## ----inductive----------------------------------------------------------------
recipe <- recipes::recipe(truth ~ x, data = development) |>
  recipes::step_impute_median(recipes::all_numeric_predictors()) |>
  recipes::step_normalize(recipes::all_numeric_predictors())
base <- workflows::workflow() |>
  workflows::add_recipe(recipe) |>
  workflows::add_model(parsnip::logistic_reg(engine = "glm"))

spec <- calibrated_workflow(base, outcome = "truth",
                            calibration_prop = .3, seed = 42)
fitted <- generics::fit(spec, development)
probability <- predict(fitted, new_data, type = "prob")
head(probability)
head(predict(fitted, new_data, type = "class"))
brier_score(assessment$truth == "yes", probability$.pred_yes)

## ----details------------------------------------------------------------------
details <- venn_abers_details(fitted, new_data)
head(details$details)
stopifnot(identical(probability, details$prediction))

## ----cross--------------------------------------------------------------------
cross_spec <- calibrated_workflow(base, "truth", method = "cross", v = 3,
                                  seed = 42)
cross_fit <- generics::fit(cross_spec, development)
head(predict(cross_fit, new_data, type = "prob"))
head(venn_abers_details(cross_fit, new_data)$details)

## ----multiclass---------------------------------------------------------------
set.seed(33)
iris_split <- rsample::initial_split(iris, prop = .8, strata = Species)
iris_dev <- rsample::analysis(iris_split)
iris_test <- rsample::assessment(iris_split)
tree <- workflows::workflow() |>
  workflows::add_formula(Species ~ .) |>
  workflows::add_model(parsnip::decision_tree(mode = "classification",
                                             engine = "rpart"))
for (classification in c("ovo", "ovr")) {
  multi_fit <- calibrated_workflow(tree, "Species",
    classification = classification, method = "cross", v = 3, seed = 33) |>
    generics::fit(iris_dev)
  print(head(predict(multi_fit, iris_test, type = "prob"), 2))
}

pair_fit <- calibrated_workflow(tree, "Species", classification = "ovo",
  model_training = "pairwise", seed = 33) |>
  generics::fit(iris_dev)
head(venn_abers_details(pair_fit, iris_test)$details)

## ----regression---------------------------------------------------------------
set.seed(52)
regression <- data.frame(x = rnorm(240))
regression$truth <- 2 * regression$x + rnorm(240, sd = .7)
reg_base <- workflows::workflow() |>
  workflows::add_formula(truth ~ x) |>
  workflows::add_model(parsnip::linear_reg(engine = "lm"))
reg_fit <- calibrated_workflow(reg_base, "truth", method = "cross", v = 3,
                               m = 1, seed = 52) |>
  generics::fit(regression[1:200, ])
head(predict(reg_fit, regression[201:240, "x", drop = FALSE]))
head(venn_abers_details(reg_fit, regression[201:240, ])$aggregate)

## ----frozen-------------------------------------------------------------------
earlier <- generics::fit(base, development[1:120, ])
frozen <- calibrate_workflow(earlier, development[121:200, ], outcome = "truth")
head(predict(frozen, new_data, type = "prob"))
head(venn_abers_details(frozen, new_data)$details)

## ----serialization------------------------------------------------------------
path <- tempfile(fileext = ".rds")
saveRDS(frozen, path)
restored <- readRDS(path)
stopifnot(identical(predict(frozen, new_data, type = "prob"),
                    predict(restored, new_data, type = "prob")))
unlink(path)

## ----postprocessing-----------------------------------------------------------
sigmoid_fit <- calibrated_workflow(base, "truth", seed = 42,
  postprocessing = list(method = "sigmoid", strategy = "holdout", prop = .3,
                        sigmoid_input = "probability")) |>
  generics::fit(development)
head(predict(sigmoid_fit, new_data, type = "prob"))

## ----early-stopping, eval=FALSE-----------------------------------------------
# spec <- parsnip::boost_tree(mode = "classification", trees = 100L, stop_iter = 10L) |>
#   parsnip::set_engine("xgboost", validation = 0.2)
# boosted <- workflows::workflow() |>
#   workflows::add_formula(truth ~ .) |>
#   workflows::add_model(spec)
# fitted <- calibrated_workflow(boosted, "truth", seed = 42) |>
#   generics::fit(development)


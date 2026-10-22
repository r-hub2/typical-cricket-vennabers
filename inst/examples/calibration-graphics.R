# After installation: source(system.file("examples", "calibration-graphics.R",
#                                       package = "vennabers"))
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

# Fit diagnoses draw the stored state of a calibrator and need no new truth
# values. Assessment graphics need held-out pairs of truth and prediction.
set.seed(20261021)
observations <- data.frame(x = rnorm(900L))
observations$truth <- rbinom(900L, 1L, stats::plogis(0.3 + 1.5 * observations$x))
training <- observations[1:300, ]
calibration <- observations[301:600, ]
assessment <- observations[601:900, ]

model <- stats::glm(truth ~ x, data = training, family = stats::binomial())
calibration_probability <- predict(model, calibration, type = "response")
assessment_probability <- predict(model, assessment, type = "response")
calibrator <- venn_abers(calibration$truth, calibration_probability)
probability <- predict(calibrator, assessment_probability)$probability

# A temporary PDF also works without an interactive display. Close and remove
# it even on error. For interactive plots, call the functions directly.
local({
  plot_file <- tempfile(fileext = ".pdf")
  on.exit(unlink(plot_file), add = TRUE)
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)

  # Fit diagnosis: the calibration map and the Venn-Abers width upper - lower
  # at every stored knot. The width is neither a confidence level nor a
  # nominal coverage.
  map <- plot(calibrator, which = c("map", "width"))
  stored <- predict(calibrator, map$score)
  assert_conditions(
    identical(names(map), c("score", "lower", "upper", "probability", "width")),
    identical(map$width, map$upper - map$lower),
    identical(map$probability, stored$probability),
    identical(range(map$score), c(0, 1))
  )
  print(head(map))

  # Cross fits overlay the map of every fold calibrator; the table keeps one
  # block per fold, not the aggregated Cross probability.
  fold <- rep(c("A", "B", "C"), length.out = nrow(calibration))
  cross <- cross_venn_abers(calibration$truth, calibration_probability, fold)
  fold_maps <- plot(cross, loss = "brier")
  stopifnot(identical(unique(fold_maps$fold), c("A", "B", "C")))

  # Assessment graphic: the CORP reliability diagram of the held-out pairs.
  # The default consistency band shows where the estimate of perfectly
  # calibrated probabilities would fall; a seed makes it reproducible without
  # changing the caller's random stream.
  state <- .Random.seed
  corp <- corp_plot(assessment$truth, probability, seed = 1)
  assert_conditions(
    identical(.Random.seed, state),
    identical(corp, corp_plot(assessment$truth, probability, seed = 1)),
    identical(corp$probability, sort(unique(probability))),
    sum(corp$count) == nrow(assessment),
    all(diff(corp$estimate) >= 0),
    all(corp$band_lower >= 0 & corp$band_upper <= 1)
  )
  print(head(corp))

  # The decomposition splits the Brier score of the same pairs into
  # miscalibration, discrimination and uncertainty.
  parts <- corp_decomposition_plot(assessment$truth, probability)
  assert_conditions(
    abs(parts$brier - brier_score(assessment$truth, probability)) < 1e-12,
    abs(parts$brier - (parts$mcb - parts$dsc + parts$unc)) < 1e-12
  )
  print(parts)

  # Regression: the map over the uncalibrated prediction with reference lines
  # at the moderated truths a and b, then held-out pairs against the
  # calibrated points.
  x <- runif(600L, 0, 10)
  regression_data <- data.frame(x = x, truth = x + sin(x) + rnorm(600L, sd = 0.5))
  regression_model <- stats::lm(truth ~ x, data = regression_data[1:200, ])
  uncalibrated <- predict(regression_model, regression_data)
  regression <- venn_abers_regression(regression_data$truth[201:400],
                                      uncalibrated[201:400])
  regression_map <- plot(regression, which = c("map", "width"))
  stopifnot(identical(regression_map$width,
                      regression_map$upper - regression_map$lower))
  point <- predict(regression, uncalibrated[401:600])$prediction
  held_out_truth <- regression_data$truth[401:600]
  held_out_pairs <- regression_plot(held_out_truth, point,
                                    which = c("observed", "residual"))
  stopifnot(identical(held_out_pairs$residual, held_out_truth - point))
})

# CORP bands are pointwise resampling bands for the assessment pairs given.
# They are not Venn-Abers lower and upper probabilities, and neither graphic
# makes a statement about Venn-Abers validity.

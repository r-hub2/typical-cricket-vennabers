\name{regression_plot}
\alias{regression_plot}
\title{Assess Regression Predictions on Held-Out Pairs}
\usage{
regression_plot(truth, prediction, which = "observed", ...)
}
\arguments{
  \item{truth}{A non-empty numeric vector of finite observed values.}
  \item{prediction}{A numeric vector of finite predicted values, the same
  length as \code{truth}, such as the \code{prediction} column of
  \code{\link{predict.venn_abers_regression}}.}
  \item{which}{One or more distinct values of \code{"observed"} (default)
  and \code{"residual"}: the diagrams to draw.}
  \item{...}{Base R graphical options passed to \code{plot()}, taking
  precedence over the defaults in every panel.}
}
\value{
Invisibly, a data frame with one row per pair in input order and the
columns \code{prediction}, \code{truth}, and \code{residual}
(\code{truth - prediction}), the same whatever \code{which} holds.
}
\description{
Assesses regression predictions on evaluation pairs of truth and
prediction, held out from fitting the calibrator, such as an assessment set
or a test set. Unlike \code{\link{plot.venn_abers_regression}}, which shows
the stored state of a fit, this graphic needs new truth values.

\code{which = "observed"} draws the observed against the predicted values
with the dashed diagonal of equality, titled \code{"Observed vs
predicted"}. \code{which = "residual"} draws the residuals, observed minus
predicted, against the predicted values with a dashed horizontal line at
zero, titled \code{"Residuals vs predicted"}.
}
\details{
The pairs are drawn as points, without a smoother. Choosing both diagrams
draws them in this order, one above the other over the same predictions, in
their own grid, and restores the previous graphical parameters on exit, also
after errors.

Invalid pairs are errors; there is no \code{na.rm}. The graphic describes
the evaluation pairs given and makes no statement about the Venn-Abers
regression bounds \code{lower} and \code{upper}, which it does not use.
}
\seealso{
\code{\link{venn_abers_regression}} for regression calibrators,
\code{\link{plot.venn_abers_regression}} for their fit diagnoses, and
\code{\link{corp_plot}} and \code{\link{calibration_plot}} for held-out
classification probabilities.
\code{vignette("diagnostics", package = "vennabers")} compares both kinds of
graphics.
}
\examples{
set.seed(1)
score <- runif(300, 0, 10)
truth <- score + sin(score) + rnorm(300, sd = .5)
fitting <- 1:150
calibrator <- venn_abers_regression(truth[fitting], score[fitting])
assessment <- predict(calibrator, score[-fitting])$prediction

regression_plot(truth[-fitting], assessment)
regression_plot(truth[-fitting], assessment, which = c("observed", "residual"))
}

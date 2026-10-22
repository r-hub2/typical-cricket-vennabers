\name{calibration_plot}
\alias{calibration_plot}
\title{Plot Classification Probability Calibration}
\usage{
calibration_plot(truth, probability, bins = 10, ...)
}
\arguments{
  \item{truth}{A non-empty logical or numeric vector containing only
  \code{FALSE}/\code{0} and \code{TRUE}/\code{1}. \code{TRUE}/\code{1}
  is the event class. Alternatively, factor truth under the scoring contract
  in \code{\link{brier_score}}; assessment may omit observations of some levels.}
  \item{probability}{A numeric vector of regularized event probabilities,
  the same length as \code{truth}, with finite values in \code{[0, 1]}.
  For factor truth, the final class matrix under the input contract of
  \code{\link{venn_abers}}, never the internal pair probabilities.}
  \item{bins}{One finite whole number of equal-width bins.}
  \item{...}{Base R graphical options passed to \code{plot()}.}
}
\value{
Invisibly, a data frame with columns \code{mean_probability},
\code{event_rate}, and \code{count} for each non-empty bin.
Multiclass adds a leading factor column \code{class} with the truth levels,
grouping rows by class order and then increasing bin order. Each class is
evaluated against the rest using its final probability. Empty bins are omitted.
}
\description{
Draws a dependency-free reliability diagram of mean regularized probability
against observed event rate.
Multiclass draws one panel per class in a rectangular layout and restores
the previous graphical parameters on exit. Graphical options apply to every
panel; by default its title identifies the class. Bins are equal-width,
left-closed and right-open, except that the last bin includes one.
}
\seealso{
\code{\link{corp_plot}} for CORP reliability diagrams without bins and with
uncertainty bands, \code{\link{regression_plot}} for held-out regression
predictions, and \code{\link{plot.venn_abers}} for fit diagnoses of the
stored calibration map. \code{vignette("diagnostics", package = "vennabers")}
compares both kinds of graphics.
}
\examples{
truth <- c(FALSE, TRUE, FALSE, TRUE)
uncalibrated_probability <- c(0.1, 0.4, 0.35, 0.8)
calibrator <- venn_abers(truth, uncalibrated_probability)
prediction <- predict(calibrator, uncalibrated_probability)

calibration_plot(truth, prediction$probability)
}

\name{brier_score}
\alias{brier_score}
\alias{log_loss}
\title{Score Classification Probabilities}
\usage{
brier_score(truth, probability)

log_loss(truth, probability)
}
\arguments{
  \item{truth}{A non-empty logical or numeric vector containing only
  \code{FALSE}/\code{0} and \code{TRUE}/\code{1}. \code{TRUE}/\code{1}
  is the event class. Alternatively a non-empty factor with at least two
  non-empty levels and no missing values. Assessment subsets may omit
  observations of some levels; retain the full fitted class schema.}
  \item{probability}{A numeric vector of regularized event probabilities,
  the same length as \code{truth}, with finite values in \code{[0, 1]}.
  For factor truth, a numeric class matrix following the class-name,
  value, and row-sum contract of \code{\link{venn_abers}}.}
}
\value{
One unnamed double. \code{brier_score()} returns the binary mean squared
error. \code{log_loss()} returns the mean negative natural log-likelihood;
an incorrect probability of zero or one produces \code{Inf}.
For multiclass, Brier is the mean of the row sums of squared errors against
one-hot truth, without division by two or class count. A two-column Brier score
is twice the binary vector score for the same probabilities. Multiclass log
loss is the mean negative natural logarithm of the observed class probability;
zero observed probability gives \code{Inf}. Neither score clips probabilities.
}
\description{
Computes the Brier Score or Log Loss from final binary or multiclass
probabilities, regardless of the loss chosen during calibration prediction.
}
\examples{
truth <- c(FALSE, TRUE, FALSE, TRUE)
uncalibrated_probability <- c(0.1, 0.4, 0.35, 0.8)
calibrator <- venn_abers(truth, uncalibrated_probability)
prediction <- predict(calibrator, uncalibrated_probability)

brier_score(truth, prediction$probability)
log_loss(truth, prediction$probability)
}

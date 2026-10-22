\name{venn_abers_regression}
\alias{venn_abers_regression}
\alias{predict.venn_abers_regression}
\title{Moderated Inductive Venn-Abers Regression}
\usage{
venn_abers_regression(truth, prediction, ..., m = 1, epsilon = NULL)

\method{predict}{venn_abers_regression}(object, newdata, ...)
}
\arguments{
  \item{truth}{A nonempty numeric vector of finite real-valued outcomes.}
  \item{prediction}{A numeric vector of finite uncalibrated regression predictions,
  the same length as \code{truth}, from an independently fixed model.}
  \item{m}{One finite numeric whole number satisfying
  \eqn{1 \le m \le \lfloor(k-1)/2\rfloor}, for calibration size \eqn{k}.
  The default is 1, requiring at least three observations.}
  \item{epsilon}{Alternatively, one finite numeric scalar strictly between 0 and 1.
  Computes \eqn{m=\lfloor\epsilon(k+1)/2\rfloor}, then checks the same domain.
  Supplying both \code{m} and \code{epsilon} explicitly is an error; omitting
  \code{m} permits \code{epsilon} alone. Explicit \code{NULL} is invalid.}
  \item{object}{A fitted \code{venn_abers_regression} object.}
  \item{newdata}{A numeric vector of finite new predictions from the same fixed
  model. Zero or one observation is allowed.}
  \item{...}{Must be empty. Unknown options are errors.}
}
\value{
Fit returns an object of class \code{venn_abers_regression}. Its documented
metadata are \code{observations}, \code{requested} (an empty list for the
implicit default, otherwise a named list containing \code{m} or \code{epsilon}),
\code{m} (effective integer-valued moderation), \code{epsilon} (effective level
\eqn{2m/(k+1)}), and \code{a}, \code{b} (outer bounds in outcome units).
Other fields are internal. The fit can be saved with \code{saveRDS()} and
restored with \code{readRDS()} in a fresh R process with this package loaded.

Predict returns a data frame with double columns \code{prediction},
\code{lower}, \code{upper}, in input order, including for zero or one row.
All results are finite and the point lies between the bounds up to numerical
roundoff. Prediction only looks up fitted endpoints and merges them; it trains nothing.
\code{\link{plot.venn_abers_regression}} draws the stored calibration map and
returns it invisibly as a data frame.
}
\description{
Calibrates the numeric predictions of an already fitted regression model with
an Inductive Venn-Abers regressor (IVAR), returning a calibrated point
prediction between a lower and an upper bound. It needs the observed outcomes
and the model's predictions for a separate calibration set; the caller trains
and retains the model and chooses all data roles. Choose it when one held-out
calibration set is affordable; use \code{\link{cross_venn_abers_regression}}
to calibrate with retained fold models, or \code{\link{calibrated_workflow}} to
let the package fit the model and split the data.
}
\details{
The calibrator applies asymmetric moderation and isotonic regression to the
supplied real-valued predictions, then combines the resulting bounds into a
point estimate.

All vectors align by position, regardless of names. Negative values and values
outside \eqn{[0,1]} are valid; numeric 0/1 truth remains regression. Factors,
logical values, matrices, lists, classed vectors, missing/nonfinite values,
length mismatches and invalid options are rejected. Inputs are never rounded,
clipped or silently removed. Tied predictions and new predictions beyond the
calibration prediction range are supported.

For sorted outcomes \eqn{y_{(1)},\ldots,y_{(k)}}, set \eqn{a=y_{(m)}} and
\eqn{b=y_{(k-m+1)}}. The lower calibration replaces outcomes below \eqn{a}
with \eqn{a} and those above \eqn{y_{(k-m)}} with that value. The upper
calibration replaces outcomes below \eqn{y_{(m+1)}} with that value and
those above \eqn{b} with \eqn{b}. Every observation is retained. Lower and
upper isotonic endpoints correspond to inserting the new score with outcome
\eqn{a} and \eqn{b}, respectively; tied scores are pooled.

For endpoints \eqn{L,U}, use weights \eqn{w_L=L-a}, \eqn{w_U=b-U} and point
\eqn{(w_LL+w_UU)/(w_L+w_U)}. Only when both weights are exactly zero is the
midpoint used, including for nonconstant outcomes. Small positive weights use
the weighted merger. Constant outcomes and coincident endpoints return their
common finite value. Arithmetic is scaled to avoid avoidable overflow.

The theoretical setting assumes IID observations and a base model determined
independently of calibration observations, including preprocessing and tuning.
Assessment outcomes must be excluded from both model training and calibration.
The bounds concern a regularized conditional expectation; they are not
confidence or prediction intervals for individual future outcomes. No
individual coverage, universal improvement in RMSE/MAE, or general selector
guarantee for this particular point merger is asserted.

This manual R entry point is an R extension, not a public manual Python wrapper
API. Offline fixtures use valid regression components from Python venn-abers
1.5.4, commit \code{0e7b4bc575ce230134576ba1c76688d46bf99ac7}. Reference lower-bound
defects are checked with reflected upper components and independent isotonic
calculations; undefined reference point results use the explicit contract above.
See \file{tests/testthat/extdata/regression-reference.md} in the source package
for provenance.
Package-managed training, automatic splits, new regression score/plot functions,
exact minimax, observation weights and rounding options are not provided.
}
\seealso{\code{vignette("method", package = "vennabers")} for the theory and
its limits, \code{vignette("manual-calibration", package = "vennabers")}
for a worked introduction, and \link{vennabers_examples} for annotated
installed scripts.

\code{\link{cross_venn_abers_regression}}}
\examples{
# The installed example fixes disjoint training/calibration/assessment roles,
# fits lm() only on training data, and reports Base R RMSE and MAE.
source(system.file("examples", "inductive-regression.R", package = "vennabers"))

fit <- venn_abers_regression(c(-10, 0, 2, 6, 20), rep(0, 5))
fit$requested
fit$m
fit$epsilon
predict(fit, c(-1, 0, 1))
}

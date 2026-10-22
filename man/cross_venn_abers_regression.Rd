\name{cross_venn_abers_regression}
\alias{cross_venn_abers_regression}
\alias{predict.cross_venn_abers_regression}
\title{Manual Cross Venn-Abers Regression with Matched Fold Models}
\usage{
cross_venn_abers_regression(truth, prediction, fold, ..., m = 1,
                          epsilon = NULL)

\method{predict}{cross_venn_abers_regression}(object, newdata, ...,
                                           type = "response")
}
\arguments{
  \item{truth}{A nonempty numeric vector of finite real-valued outcomes.}
  \item{prediction}{A numeric vector of finite out-of-fold regression predictions,
  the same length as \code{truth}. Each observation must have been excluded from
  the entire learned procedure supplying its prediction.}
  \item{fold}{One fold ID per observation: nonempty character labels or finite
  positive whole-number numeric IDs. At least two represented folds are required;
  unequal sizes are allowed. Factors, logical IDs, missing values and invalid
  lengths are errors. Labels match literally, without trimming or case folding;
  numeric IDs use Base R \code{as.character()} names.}
  \item{m}{One finite positive whole number, default 1. Every fold of size
  \eqn{k_f} must satisfy \eqn{1 \le m \le \lfloor(k_f-1)/2\rfloor}.
  Even the default requires at least three observations per fold.}
  \item{epsilon}{Alternatively, one finite numeric scalar strictly between 0 and 1.
  Each fold uses \eqn{m_f=\lfloor\epsilon(k_f+1)/2\rfloor}, checked against the
  same domain. Explicitly supplying both \code{m} and \code{epsilon} is an error;
  omitting \code{m} permits \code{epsilon} alone. Explicit \code{NULL} is invalid.}
  \item{object}{A fitted \code{cross_venn_abers_regression} object.}
  \item{newdata}{A finite real-valued numeric matrix with exactly one uniquely
  named column per stored fold. Each column comes from that fold's retained
  model; every row refers to the same new observation across columns.
  Columns are matched by name and may be reordered. Missing, additional,
  duplicated, empty, missing-value or absent column names are errors.}
  \item{type}{Exactly \code{"response"} (default) or \code{"details"}.}
  \item{...}{Must be empty. Unknown options are errors.}
}
\value{
Fit returns an object of class \code{cross_venn_abers_regression}. Documented
fields are \code{folds} (character IDs in first-observed order),
\code{observations} (total calibration size), and \code{calibrators}, a named
list of complete \code{venn_abers_regression} fits, one per fold. Each records
\code{observations}, \code{requested} (empty for the implicit default, otherwise
the explicit \code{m} or \code{epsilon}), effective \code{m}, effective
\code{epsilon} equal to \eqn{2m_f/(k_f+1)}, and outer bounds \code{a}, \code{b}.
Unequal folds can have different effective moderation levels. Any invalid
fold fails the entire fit with fold context; nothing is capped or omitted.
Fits can be saved with \code{saveRDS()} and restored with \code{readRDS()} in
a fresh R process with the installed package loaded.

Predict returns a data frame with double columns \code{prediction},
\code{lower}, \code{upper}, in supplied row order. With \code{type = "details"},
it returns a list containing \code{prediction} (that standard data frame) and
\code{folds} (a named list of fold data frames with the same three columns,
in stored fold order). Zero and one input row retain these schemas; empty
details retain every fold ID. Outputs are finite and the point lies between
the bounds up to checked numerical roundoff. Row names do not align observations.
\code{\link{plot.cross_venn_abers_regression}} draws each fold's stored
calibration map and returns the maps invisibly as a data frame.
}
\description{
Calibrates the numeric predictions of regression models with a Cross
Venn-Abers regressor (CVAR): one moderated Inductive Venn-Abers regression
calibrator per fold, with complete fold results averaged with equal weights.
It needs the observed outcomes, the out-of-fold predictions of one retained
model per fold and the fold labels; the caller trains and retains the models
and chooses all data roles. Choose it over \code{\link{venn_abers_regression}}
when a single held-out calibration set would cost too much training data, or
use \code{\link{calibrated_workflow}} with \code{method = "cross"} to let the
package fit the fold models.
}
\details{
Each fold follows the asymmetric moderation, isotonic endpoints and Formula 6
weighted point merger documented in \code{\link{venn_abers_regression}}.
The midpoint is used only when both merger weights vanish, including with
nonconstant truth; small positive weights retain the weighted merger. Constant
truth and coincident endpoints return their common finite value.

For each new observation, first compute the complete point, lower and upper
IVAR results for every fold. Then average each of those three quantities
separately, with equal fold weights even for unequal fold sizes. No pooling,
fold-size weighting, or new point merger on averaged bounds is performed.
The three aggregates can be reconstructed by averaging the corresponding
columns from \code{details$folds}. Means use scaling and compensated summation
to avoid intermediate overflow and preserve cancellation residuals, including
on platforms without extended floating-point precision.

Truth and prediction vectors align by position, ignoring names. Negative
values and values outside \eqn{[0,1]} are valid; numeric 0/1 truth remains
regression. Factors, logical values, lists, classed vectors, matrices at fit,
nonfinite values and mismatched lengths are rejected. No values are rounded,
clipped or silently dropped. Ties and new scores outside calibration ranges
are supported. Predict trains nothing and accepts no new moderation options.

The theoretical setting assumes IID observations and a base procedure fixed
independently of its calibration observations. Every observation in a fold
must have been excluded from its corresponding model's entire learned flow,
including preprocessing and tuning. That same retained model supplies the
fold's new prediction column. Assessment outcomes stay outside training and
calibration. The caller is responsible for model provenance and semantic row
identity; matching column names cannot prove these. A single model retrained
on all development observations cannot replace the matching fold models.

The aggregated regression bounds are descriptive summaries of fold bounds,
without any additional coverage or selector guarantee. They are not confidence
intervals or prediction intervals for individual future outcomes. Neither the
Cross point nor these bounds inherit a general selector/autocalibration claim,
and improvement in RMSE or MAE is not promised. Epsilon is a moderation level,
not a prediction-interval error rate.

This manual Cross entry point is an R extension. Offline reference fixtures
combine valid IVAR components from Python venn-abers 1.5.4 and independently
computed arithmetic means; no public manual Python CVAR wrapper is claimed.
See \file{tests/testthat/extdata/cross-regression-reference.md} in the source package for provenance,
retained reference defects and fixed scale-aware tolerances.
}
\seealso{\code{vignette("method", package = "vennabers")} for the theory and
its limits, \code{vignette("manual-calibration", package = "vennabers")}
for a worked introduction, and \link{vennabers_examples} for annotated
installed scripts.

\code{\link{venn_abers_regression}}}
\examples{
# Fits retained lm() fold models with the held-out observations excluded,
# reconstructs the aggregates, and scores on separate assessment data.
source(system.file("examples", "cross-regression.R", package = "vennabers"))

fit <- cross_venn_abers_regression(
  c(-10, 0, 2, 6, 20, 10, 10, 10), rep(0, 8),
  rep(c("A", "B"), c(5, 3))
)
predict(fit, cbind(B = c(1, 0), A = c(-1, 0)), type = "details")
}

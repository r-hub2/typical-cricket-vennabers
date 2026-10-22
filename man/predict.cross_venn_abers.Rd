\name{predict.cross_venn_abers}
\alias{predict.cross_venn_abers}
\title{Predict with Matched Fold-Model Probabilities}
\usage{
\method{predict}{cross_venn_abers}(object, newdata, ..., loss = "log", type = "response")
}
\arguments{
  \item{object}{A fitted \code{\link{cross_venn_abers}} object.}
  \item{newdata}{A numeric/integer matrix of finite event probabilities in
  \code{[0, 1]}, one row per new observation and one column per corresponding
  fold model. Column names must be complete, non-empty, non-missing, unique,
  and match the fitted fold identities exactly. Any column order is accepted.
  A correctly named zero-row matrix is accepted. For multiclass, a list
  named uniquely and completely by fitted fold IDs, each element a numeric
  class matrix with the same number of new rows and all fitted class columns.
  Classes and folds are aligned by name. Zero-row matrices are accepted.}
  \item{...}{Must be empty; unsupported arguments are errors.}
  \item{loss}{Exactly \code{"log"} (default) or \code{"brier"}.}
  \item{type}{Exactly \code{"response"} (default), \code{"class"}, or
  \code{"details"}. Supply both options by their full argument names.}
}
\value{
For binary calibration, \code{type = "response"} returns a data frame with exactly one double column, \code{probability}, containing
one finite regularized event probability in \code{[0, 1]} per input row,
in input row order. A zero-row input returns the same empty schema.

For \code{type = "class"}, an unnamed logical vector: \code{TRUE} when the
combined event probability is at least 0.5, otherwise \code{FALSE}. This uses
the combined probability, not a vote among folds. Zero rows return \code{logical(0)}.

For \code{type = "details"}, a list with \code{prediction} (the standard
data frame) and \code{folds} (a data frame with columns \code{row}, \code{fold},
\code{lower}, \code{upper}, \code{probability}). The integer \code{row} is the
one-based position in this call's new data, independent of input row names.
The character \code{fold} is the stored fold ID. The three double columns
contain that fold's IVAP endpoints and its probability for the chosen loss.
Rows are grouped by stored fold order, then new-row order. There is exactly
one record per new row and fold: zero input rows give an empty fold table,
and one input row gives one record per fold. These are new-observation
results, not out-of-fold calibration probabilities.

For multiclass, the standard output is a double class-probability matrix in
stored factor-level order, with row names from the first stored fold's matrix.
Class output is a factor with those levels, choosing the highest final
probability and the first level on exact ties (also for two classes).
Details contain \code{prediction} plus \code{pairs} for OvO or \code{classes}
for OvR. The table columns are integer \code{row}, character \code{non_event}
(OvO only), character \code{event}, character \code{fold}, and double
\code{lower}, \code{upper}, \code{probability}. Events are the later pair
class for OvO or the literal class against the rest for OvR. Records follow
pair/class order, then fold first-appearance order within that subproblem's
calibration rows, then new-row order. Each row/subproblem/fold appears once.
The points use the chosen loss; endpoints do not change with loss. Zero or
one new row retains matrix, factor, and detail-table types and mappings.
}
\description{
Applies a fitted \code{\link{cross_venn_abers}} calibrator to new observations
and returns calibrated probabilities, predicted classes or per-fold details.
It needs the new probabilities from every retained fold model, as one named
column per fold (binary) or one named class matrix per fold (multiclass), with
the same observation in each row. Use \code{\link{predict.venn_abers}} for an
Inductive calibrator, which needs only the probabilities of a single model.
}
\details{
New-observation columns are matched to their fold calibrators by identity, and
the resulting IVAP results are combined using the chosen equal-weight
aggregation.

Character fold identities match literally. Numeric fitting IDs match their
Base R character representations. Missing, extra, duplicate, empty, or absent
column names fail the entire call, as do invalid probability values and
unsupported binary types or shapes (including vectors, lists, data frames, logical
matrices and character matrices). A column for an unobserved fold is an error.
Values are never dropped, clamped or rounded. For multiclass, invalid list
names, unequal row counts, or an invalid class matrix in any fold fail the
whole call. The class value, name, and row-sum rules match fit time; accepted
row sums within tolerance are normalized. Row names are not alignment keys.

Use each retained fold model that supplied that fold's calibration
probabilities; predictions from one all-data refit cannot substitute for the
matched columns. The caller ensures the same observation occupies a given
row in every column. Row names do not establish that alignment, and valid
column names do not prove model provenance or exclusion of calibration folds
from training. See \code{\link{cross_venn_abers}} for pipeline responsibilities.

For log loss, let \eqn{p_0^k, p_1^k} be the IVAP endpoints from fold
\eqn{k}. If \eqn{A} is the geometric mean of \eqn{1-p_0^k} and \eqn{B} is the
geometric mean of \eqn{p_1^k}, the result is \eqn{B/(A+B)}. Evaluation uses
\code{plogis(mean(log(upper)) - mean(log1p(-lower)))} in log space, with equal
fold weight and no endpoint clipping. It does not pool calibration rows,
weight folds by size, or average log-loss fold point probabilities.

For Brier, each fold probability is \code{upper + (lower^2 - upper^2) / 2}.
The combined probability is the arithmetic mean of these fold probabilities,
with equal weights even for unequal fold sizes. Both aggregates can be
reconstructed from \code{type = "details"}. Loss changes neither the fold
endpoints nor stored learning state and requires no new fit.

No aggregate endpoints are returned. The geometric pair
\eqn{(1-A, B)} can cross; it is not an ordered interval or confidence interval.
Neither that pair nor the combined point inherits the IVAP perfect-calibration
guarantee. There is no universal score-improvement guarantee. The output
feeds the existing scores and Base R calibration plot directly. For
multiclass, apply this binary merger within each pair/class first, then OvO
coupling or OvR normalization as described in \code{\link{cross_venn_abers}}.
Reconstruct the final matrix from these per-fold diagnostics; there is no
average of already coupled/normalized multiclass fold predictions.
}
\seealso{\code{vignette("method", package = "vennabers")} explains the Cross
merger and the limits of its guarantees;
\code{vignette("manual-calibration", package = "vennabers")} shows matched
fold-model columns in a worked example.

\code{\link{cross_venn_abers}}, \code{\link{venn_abers_one_hot}},
\code{\link{brier_score}}, \code{\link{log_loss}}}
\examples{
calibrator <- cross_venn_abers(c(0, 1, 0, 1), c(.1, .8, .2, .9), c(1, 1, 2, 2))
newdata <- cbind("1" = c(.2, .8), "2" = c(.3, .7))
predict(calibrator, newdata)
stopifnot(identical(predict(calibrator, newdata), predict(calibrator, newdata[, 2:1])))
predict(calibrator, newdata[FALSE, ])
details <- predict(calibrator, newdata, loss = "brier", type = "details")
details$folds
reconstructed <- vapply(seq_len(nrow(newdata)), function(i) {
  mean(details$folds$probability[details$folds$row == i])
}, 0)
stopifnot(max(abs(reconstructed - details$prediction$probability)) < 1e-10)
predict(calibrator, newdata, loss = "brier", type = "class")
}

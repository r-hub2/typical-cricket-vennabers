\name{predict.venn_abers}
\alias{predict.venn_abers}
\title{Predict with a Venn-Abers Calibrator}
\usage{
\method{predict}{venn_abers}(object, newdata, ..., loss = "log", type = "response")
}
\arguments{
  \item{object}{A fitted \code{venn_abers} object.}
  \item{newdata}{A numeric vector of uncalibrated event probabilities with
  finite values in \code{[0, 1]}. \code{numeric(0)} is accepted for binary
  objects. For multiclass, a numeric matrix under the fit-time class,
  name, range, and row-sum contract; zero rows are accepted.}
  \item{...}{Must be empty; unsupported arguments are errors.}
  \item{loss}{Exactly \code{"log"} (default) or \code{"brier"}.}
  \item{type}{Exactly \code{"response"} (default) or \code{"class"}.
  Multiclass also accepts \code{"details"}.
  Supply both options by their full argument names.}
}
\value{
For \code{type = "response"}, a data frame with one row per value in \code{newdata} and three double
columns: \code{lower} and \code{upper} are the Venn-Abers interval bounds,
and \code{probability} is the regularized event probability
using the chosen loss. For \code{type = "class"}, an unnamed logical vector,
\code{TRUE} when the final event probability is at least 0.5 and \code{FALSE}
otherwise. Zero observations return the same empty data-frame schema or
\code{logical(0)}; one observation retains its data-frame/vector shape.

For multiclass, \code{"response"} returns a numeric class-probability matrix
in fitted factor-level order, retaining new-row order and matrix row names.
\code{"class"} returns an unnamed factor with the stored levels; the first
stored class wins exact ties. This also applies with exactly two classes,
whereas binary vector ties select the event. Classes use the final probabilities
of the fitted method. \code{\link{venn_abers_one_hot}} represents binary or
multiclass class decisions as a one-hot matrix. For OvO, \code{"details"} returns a list
with \code{prediction} (the response matrix) and \code{pairs}, a data frame
with \code{row}, \code{non_event}, \code{event}, \code{lower}, \code{upper},
and \code{probability}. Row is the integer one-based new-row position;
the two class names are character columns. The remaining columns are doubles.
Pairs follow factor-level combination order, then new-row order.
For OvR, the list contains \code{prediction} and \code{classes}, a data frame
with integer \code{row}, character \code{event}, and double \code{lower},
\code{upper}, \code{probability}. The event is the named class against all
other classes; no artificial rest label is introduced. Records follow stored
class order, then new-row order. Zero or one
new rows retain these matrix, factor, list, and table schemas.
}
\description{
Applies a fitted \code{\link{venn_abers}} calibrator to new uncalibrated
probabilities from the same model and returns calibrated probabilities or
predicted classes, preserving the order of \code{newdata}. For binary
calibration it returns the Venn-Abers interval of two probabilities together
with the regularized event probability that merges them according to
\code{loss}; multiclass calibration returns a class-probability matrix. Use
\code{\link{predict.cross_venn_abers}} for a Cross calibrator, which needs one
probability column per retained fold model.
}
\details{
Log-loss regularization is \code{upper / (1 - lower + upper)}. Brier
regularization is \code{upper + (lower^2 - upper^2) / 2}. Changing loss does
not change the interval endpoints or mutate the fitted calibrator. New
probabilities are never rounded, even if calibration used \code{precision}.
Input names are not alignment keys. The endpoints are Venn-Abers bounds,
not confidence intervals. The probability column can be passed to the
existing scores and calibration plot with either loss choice.

For OvO, class and pair normalization precedes prediction, including
the zero-pair replacement described in \code{\link{venn_abers}}. New pair
probabilities remain unrounded. After binary regularization, for class j let
q[j,i] be its probability in its pair with class i. Reciprocal coupling is
\code{w[j] = 1 / (sum(1 / q[j,i]) - (K - 2))}, summed over the other K-1
classes. The response is \code{w / sum(w)}. For each diagnostic record,
\code{probability} is the event's q and \code{1 - probability} is the
non-event's q.

For OvR, input class normalization precedes prediction. Each class probability
enters its binary calibrator unrounded. After log or Brier regularization,
divide the class points by their sum for each new row to obtain the response.
The \code{classes} table supplies those pre-normalization points, so details
reconstruct the final matrix for either method and loss.
The endpoints describe only the binary subproblems, not multiclass confidence
intervals. Binary Brier regularization does not imply optimality for the final
multiclass Brier score or guarantee improvement in assessment scores.
}
\seealso{\code{vignette("method", package = "vennabers")} explains the
interval and the merger into one probability;
\code{vignette("manual-calibration", package = "vennabers")} shows prediction
in a worked example.

\code{\link{venn_abers}}, \code{\link{venn_abers_one_hot}},
\code{\link{brier_score}}, \code{\link{log_loss}}}
\examples{
calibrator <- venn_abers(c(0, 1, 0), c(0.1, 0.2, 0.3))
predict(calibrator, c(0.05, 0.2, 0.4))
predict(calibrator, c(0.05, 0.2, 0.4), loss = "brier")
predict(calibrator, c(0.05, 0.2, 0.4), loss = "brier", type = "class")
}

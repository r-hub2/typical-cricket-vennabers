\name{venn_abers}
\alias{venn_abers}
\title{Fit an Inductive Venn-Abers Calibrator}
\usage{
venn_abers(truth, probability, precision = NULL, method = "ovo")
}
\arguments{
  \item{truth}{A non-empty logical or numeric vector containing only
  \code{FALSE}/\code{0} and \code{TRUE}/\code{1}. \code{TRUE}/\code{1}
  is the event class. For multiclass, a non-empty factor without missing
  values, with at least two non-empty levels. Levels fix class order;
  every level must be observed in the calibration set.}
  \item{probability}{A numeric vector of uncalibrated event probabilities,
  the same length as \code{truth}, with finite values in \code{[0, 1]}.
  For factor truth, a numeric matrix with one row per truth and complete,
  unique, non-empty column names matching the factor levels. Columns are
  reordered by name. Entries must be finite and in \code{[0, 1]}; row sums
  must be within \code{1e-8} of one and are normalized within that tolerance.
  No clipping, name repair, or dropped rows.}
  \item{precision}{\code{NULL} (default) disables rounding. Otherwise, one
  finite numeric whole number, including zero or negative values, passed to
  native R \code{round()} for calibration probabilities only. Logical,
  missing, nonfinite, fractional, and multiple values are errors.}
  \item{method}{Exactly \code{"ovo"} (default) or \code{"ovr"}. Explicitly
  supplying this argument requires factor truth and a class matrix.}
}
\value{
An S3 object of class \code{"venn_abers"}: a list holding the fitted
calibration state, from which new probabilities are calibrated without
refitting. For binary truth it records the sorted unique calibration
probabilities with their lower and upper Venn-Abers values and the
calibration and event counts; for factor truth, the class levels, the method
and one binary calibrator per class pair (OvO) or class (OvR). These list
elements are internal and not public API. Use the methods instead:
\itemize{
\item \code{\link{predict.venn_abers}} returns the lower and upper values with
the regularized event probability (binary), class probabilities (multiclass),
class decisions or binary subproblem details, with the loss chosen at
prediction time.
\item \code{\link{plot.venn_abers}} draws the stored calibration map and
returns it invisibly as a data frame.
\item \code{\link{print.venn_abers}} summarizes the calibration size, events and
probability knots, or the method, calibration size, classes and number of
subproblems.
}
No predictive model or new data are stored.
}
\description{
Turns the uncalibrated probabilities of an already fitted binary or multiclass
classifier into calibrated probabilities with an Inductive Venn-Abers predictor
(IVAP). It needs the true classes and the model's probabilities for a separate
calibration set; it creates neither predictive models nor data splits. Choose
it when one held-out calibration set is affordable; use
\code{\link{cross_venn_abers}} when data are scarce and one model per fold can
be retained, or \code{\link{calibrated_workflow}} to let the package fit the
model and split the data.
}
\details{
The calibration set must be separate from the fixed predictive model's
training, including learned preprocessing and tuning. Single-class and
singleton binary calibration sets are valid. Rows align by position, not names.

One-vs-One (\code{method = "ovo"}) uses each pair in factor-level order and only that pair's
calibration observations. The later class is the binary event. Its probability
is divided by the pair sum before optional rounding. An exactly zero pair sum
uses event probability 0.5, including during prediction. This fit replacement
is also rounded: at \code{precision = 0}, native R rounds 0.5 to zero.
Only the event probability is rounded, not its complement independently.
There are no separate trained pair models. All global classes must occur,
so an Inductive OvO pair always has both observed classes. The underlying
binary calibrator still accepts single-class subsets.

One-vs-Rest (\code{method = "ovr"}) uses all calibration observations for
each class. That class is the binary event; its class probability is the
calibration probability, optionally rounded at fit time. The binary points
are normalized across classes only after the chosen regularization.

Original probabilities are validated before rounding; an out-of-range input
cannot become valid through rounding. Exact ties after rounding are combined.
New probabilities supplied to \code{predict()} remain unrounded, and changing
precision requires a new fit. Native R rounding can differ from NumPy at
floating-point boundaries; NumPy rounding is not emulated.

Choose log-loss (default) or Brier regularization at prediction time without
refitting. See \code{\link{predict.venn_abers}} for probabilities, classes,
pair coupling or class normalization, and details. Manual Brier for both
multiclass methods is a deliberate R extension over the examined Python
1.5.4 manual path, which supports log regularization and precision.
Exactly two factor levels also select the
multiclass matrix contract, with its distinct tie rule and Brier score scale.
}
\seealso{\code{vignette("method", package = "vennabers")} for the theory and
its guarantees, \code{vignette("manual-calibration", package = "vennabers")}
for a worked introduction, and \link{vennabers_examples} for annotated
installed scripts.

\code{\link{predict.venn_abers}}, \code{\link{cross_venn_abers}},
\code{\link{venn_abers_regression}}}
\examples{
calibrator <- venn_abers(
  truth = c(FALSE, TRUE, FALSE, TRUE),
  probability = c(0.1, 0.4, 0.35, 0.8)
)
calibrator

predict(calibrator, c(0.2, 0.5, 0.9))
rounded <- venn_abers(c(0, 1, 0, 1), c(.501, .502, .509, .510), precision = 2)
predict(rounded, c(.501, .505, .51), loss = "brier")

# Fixed-model probabilities on a separate calibration set.
classes <- c("cat", "dog", "bird")
truth <- factor(rep(classes, 2), levels = classes)
p <- rbind(c(.7, .2, .1), c(.2, .6, .2), c(.1, .2, .7),
           c(.5, .3, .2), c(.3, .5, .2), c(.2, .3, .5))
colnames(p) <- classes
fit <- venn_abers(truth, p, precision = 2)
new_p <- rbind(c(.6, .3, .1), c(.1, .4, .5))
colnames(new_p) <- classes
assessment <- factor(c("cat", "bird"), levels = classes)
details <- predict(fit, new_p, loss = "brier", type = "details")
details$pairs
predict(fit, new_p, loss = "brier", type = "class")
brier_score(assessment, details$prediction)
log_loss(assessment, details$prediction)
calibration_plot(assessment, details$prediction)

ovr <- venn_abers(truth, p, method = "ovr", precision = 2)
predict(ovr, new_p)
ovr_details <- predict(ovr, new_p, loss = "brier", type = "details")
ovr_details$classes
predict(ovr, new_p, loss = "brier", type = "class")
brier_score(assessment, ovr_details$prediction)
log_loss(assessment, ovr_details$prediction)
calibration_plot(assessment, ovr_details$prediction)

# Installed standalone examples for binary options (precision, loss, event
# selection, Cross) and both multiclass methods:
system.file("examples", c("binary-options.R", "inductive-ovo.R", "inductive-ovr.R"),
  package = "vennabers")
}

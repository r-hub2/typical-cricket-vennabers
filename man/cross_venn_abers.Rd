\name{cross_venn_abers}
\alias{cross_venn_abers}
\title{Fit a Binary or Multiclass Cross Venn-Abers Calibrator}
\usage{
cross_venn_abers(truth, probability, fold, precision = NULL, method = "ovo")
}
\arguments{
  \item{truth}{A non-empty logical or numeric/integer vector containing only
  \code{FALSE}/\code{0} and \code{TRUE}/\code{1}. The event class is always
  \code{TRUE}/\code{1}. Alternatively, a factor with at least two levels for
  multiclass calibration; every level must be globally observed.}
  \item{probability}{A numeric/integer vector of out-of-fold event
  probabilities, finite and in \code{[0, 1]}, aligned with \code{truth}. For
  factor truth, a numeric matrix with one row per observation and complete,
  unique class column names matching its levels; see \code{\link{venn_abers}}.}
  \item{fold}{An aligned vector of non-empty, non-missing character labels
  or finite positive whole-number numeric/integer IDs. Character matching is
  literal, without trimming or case folding. Numeric IDs match their Base R
  \code{as.character} representations. Repeated IDs indicate membership in
  the same fold. Factors and logical IDs are not supported.}
  \item{precision}{\code{NULL} (default) disables rounding. Otherwise, one
  finite numeric whole number, including zero or negative values, used for
  native R rounding of calibration probabilities in every fold. See
  \code{\link{venn_abers}} for validation and rounding semantics.}
  \item{method}{For factor truth, exactly \code{"ovo"} (default) or
  \code{"ovr"}. Explicit use with binary vector truth is an error.}
}
\value{
An S3 object of class \code{"cross_venn_abers"}: a list retaining fold
identities and calibration state for repeated predictions from the matching
fold models. For binary truth it records the fold IDs in first-observed order,
the total calibration size and one Inductive \code{venn_abers} calibrator per
fold; for factor truth, the class levels, the method, the fold IDs and one
Cross calibrator per class pair (OvO) or class (OvR). These list elements are
internal and not public API. Use the methods instead:
\itemize{
\item \code{\link{predict.cross_venn_abers}} combines the fold calibrators on
new probabilities from each fold model into regularized probabilities, class
decisions or per-fold details, with the loss chosen at prediction time.
\item \code{\link{plot.cross_venn_abers}} draws each fold's stored calibration
map and returns the maps invisibly as a data frame.
\item \code{\link{print.cross_venn_abers}} summarizes the number of folds and
calibration observations, and the classes for multiclass.
}
Predictive models are retained separately by the caller.
}
\description{
Calibrates binary or multiclass classifier probabilities with a Cross
Venn-Abers predictor (CVAP): every fold calibrates the model that was trained
without it, and prediction averages the fold results. It needs the true
classes, the out-of-fold probabilities of one retained model per fold and the
fold labels; it fits no models itself. Choose it over
\code{\link{venn_abers}} when a single held-out calibration set would cost too
much training data, or use \code{\link{calibrated_workflow}} with
\code{method = "cross"} to let the package fit the fold models.
}
\details{
The function fits Inductive Venn-Abers calibrators per calibration fold and, for
multiclass, per class pair or class against the rest. Prediction
combines their results with equal-weight log-loss (default) or Brier aggregation
to produce a regularized event probability.

For binary calibration, all three inputs must be vectors of equal, nonzero length with no missing
values. Names attached to vector elements are not alignment keys: positions
must refer to the same observation. Invalid types, shapes, nonbinary truth,
nonfinite values, or out-of-range probabilities cause an error for the entire
call. No rows are dropped and probabilities are not clamped or coerced from
other representations. Original probabilities are checked before optional
rounding; exact ties after rounding are combined. New prediction values are
never rounded. Changing precision requires refitting, while loss is chosen
at prediction time. R and NumPy rounding can differ at floating-point boundaries.

Supply one non-overlapping partition with at least two represented non-empty
folds. Unequal sizes, singleton folds, and single-class folds are supported;
each fold receives equal aggregation weight regardless of its size. There is
no separate empty-fold declaration. Use \code{venn_abers()} for one calibration
set.

Each fold model must exclude its calibration fold from the entire learned
training pipeline, including preprocessing and tuning. That same model must
produce both its fold's out-of-fold probabilities and its corresponding named
column of new-observation probabilities. Retain those models for subsequent
prediction: a single model refitted on all data is not a replacement. Keep
evaluation truth outside model fitting and calibration.

Callers ensure each underlying calibration observation appears once and that
new-observation rows align across model columns. The supplied arrays cannot
prove model provenance, training exclusion, absence of leakage, unique
underlying observations, or semantic row alignment.

Cross output defaults to a regularized probability. Classes and individual
fold results are available via \code{\link{predict.cross_venn_abers}}. There is
no aggregate interval. Geometric aggregate endpoints can cross.
Neither those endpoints nor the combined point inherit the IVAP
perfect-calibration guarantee. They are not confidence intervals, and score
improvement is not guaranteed on every data set.
}
\section{Unsupported inputs and features}{
Package-managed predictive-model fitting, fold generation, preprocessing,
tuning or resampling; repeated or overlapping partitions; regression calibration; separate
pairwise model training or pair-specific fold plans; observation weights; alternate event-class selection;
automatic coercion or missing-value handling; fold-size weighting or pooling;
public aggregate endpoints; formula or data-frame cross
interfaces; tidymodels/yardstick integrations; and new metrics, smoothing or
plotting features are not supported. Only Base R and standard R distribution
packages are needed. Installation and acceptance require neither Python nor
networking; Python reference-fixture regeneration is optional and separate.
}
\section{Multiclass with a shared fold plan}{
Factor levels fix the global class order. Supply one out-of-fold class matrix
and a shared fold vector. All classes must be observed globally; individual
folds may omit classes. Each OvO pair must have at least one observation in
every fold. Empty pairs fail with their class names and fold ID; singleton
and single-class binary subcalibrations are allowed. OvR uses every row in
each fold. No pair-specific training plans are introduced.

Class matrices must have finite entries in \code{[0, 1]} and row sums within
\code{1e-8} of one; accepted rows are normalized. In OvO, pair normalization
precedes calibration rounding, with zero pair sums replaced by an event
probability of 0.5 (also subject to fit rounding). The later class in factor
order is the event; in OvR, the selected class is the event against the rest.
New prediction probabilities remain unrounded.

Predict with a list of class matrices named by fold ID, from the corresponding
retained models, for the same new observations in the same order. Binary
Cross aggregation happens separately for each subproblem before reciprocal
OvO coupling or OvR normalization. Log loss uses the stable equal-weight
endpoint merger; Brier uses the equal-weight mean of binary fold points.
No pooling, size weighting, or averaging of finished multiclass fold outputs
is performed. Scores and plots use the final class matrix.

This manual shared-fold interface is an R extension. Python 1.5.4 has no
corresponding public manual Cross-multiclass whole-procedure API; numerical
checks use binary reference components and independent combination. Binary
diagnostics are not multiclass or aggregate Cross confidence intervals.
}
\seealso{\code{vignette("method", package = "vennabers")} for the theory and
its guarantees, \code{vignette("manual-calibration", package = "vennabers")}
for a worked introduction, and \link{vennabers_examples} for annotated
installed scripts.

\code{\link{predict.cross_venn_abers}},
\code{\link{print.cross_venn_abers}}, \code{\link{brier_score}},
\code{\link{log_loss}}, \code{\link{calibration_plot}}}
\examples{
# Precomputed probabilities from three retained fold models. Each model
# excluded its corresponding calibration fold from its learned pipeline.
calibrator <- cross_venn_abers(
  truth = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE),
  probability = c(.1, .8, .2, .9, .3, .7),
  fold = c("A", "A", "B", "B", "C", "C")
)
print(calibrator)

# Each row is the same new observation across all corresponding fold models.
newdata <- cbind(A = c(.2, .8), B = c(.3, .7), C = c(.1, .9))
result <- predict(calibrator, newdata)
result
test_truth <- c(FALSE, TRUE)  # Untouched by fitting and calibration.
brier_score(test_truth, result$probability)
log_loss(test_truth, result$probability)
calibration_plot(test_truth, result$probability)

# After installation, find the standalone example with:
# system.file("examples", "cross-multiclass.R", package = "vennabers")
library(vennabers)

# Precomputed out-of-fold probabilities from two retained multiclass models.
# Model A excluded rows 1:2; model B excluded rows 3:6, including preprocessing
# and tuning. Retain both models yourself. No model is trained by this package.
classes <- c("cat", "dog", "bird")
truth <- factor(c("cat", "dog", "cat", "bird", "dog", "bird"), levels = classes)
fold <- c("A", "A", "B", "B", "B", "B")
p <- rbind(c(.7, .2, .1), c(.2, .6, .2), c(.5, .3, .2),
           c(.1, .2, .7), c(.3, .5, .2), c(.2, .3, .5))
colnames(p) <- classes
# The same three untouched assessment observations, in the same row order,
# predicted by each corresponding retained model. List order is immaterial.
new_p <- list(B = rbind(c(.6, .3, .1), c(.1, .4, .5), c(.1, .8, .1)),
              A = rbind(c(.5, .4, .1), c(.2, .2, .6), c(.2, .7, .1)))
new_p <- lapply(new_p, function(x) { colnames(x) <- classes; x })
assessment <- factor(c("cat", "bird", "dog"), levels = classes)

local({
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({ grDevices::dev.off(); unlink(plot_file) })
  for (method in c("ovo", "ovr")) {
    fit <- cross_venn_abers(truth, p, fold, method = method, precision = 2)
    print(fit)
    for (loss in c("log", "brier")) {
      details <- predict(fit, new_p, loss = loss, type = "details")
      probability <- predict(fit, new_p, loss = loss)
      predicted_class <- predict(fit, new_p, loss = loss, type = "class")
      # Reconstruct from the fold endpoints/points, then combine subproblems.
      parts <- details[[2L]]
      problems <- unique(parts[c(if (method == "ovo") "non_event", "event")])
      merged <- vapply(seq_len(nrow(problems)), function(i) {
        selected <- parts$event == problems$event[i]
        if (method == "ovo") selected <- selected & parts$non_event == problems$non_event[i]
        vapply(seq_len(nrow(probability)), function(row) {
          f <- parts[selected & parts$row == row, ]
          if (loss == "log") {
            plogis(mean(log(f$upper)) - mean(log1p(-f$lower)))
          } else mean(f$probability)
        }, 0)
      }, numeric(nrow(probability)))
      recovered <- if (method == "ovr") merged else vapply(classes, function(class) {
        q <- cbind(1 - merged[, problems$non_event == class, drop = FALSE],
                   merged[, problems$event == class, drop = FALSE])
        1 / (1 + rowSums((1 - q) / q))
      }, numeric(nrow(probability)))
      recovered <- recovered / rowSums(recovered)
      stopifnot(max(abs(recovered - probability)) <= 1e-10,
                identical(probability, details$prediction),
                identical(predicted_class,
                          factor(classes[max.col(probability, ties.method = "first")],
                                 levels = classes)))
      print(details)
      print(predicted_class)
      print(c(brier = brier_score(assessment, probability),
              log = log_loss(assessment, probability)))
      groups <- calibration_plot(assessment, probability, bins = 1)
      stopifnot(max(abs(groups$mean_probability - colMeans(probability))) <= 1e-10,
                all(groups$count == length(assessment)), all(groups$event_rate == 1 / 3))
      print(groups)
    }
  }
})
# One shared fold plan serves all pairs/classes. Fold A lacks bird but every
# OvO pair remains nonempty. Empty pairs fail with the pair and fold identified.
# Endpoints describe binary subproblems, with no aggregate confidence interval.
# Manual Cross multiclass is an R extension; Python has no matching public
# manual whole-procedure API. Numerical tests reconstruct binary components.
}

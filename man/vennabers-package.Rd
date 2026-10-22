\name{vennabers-package}
\alias{vennabers}
\alias{vennabers-package}
\docType{package}
\title{Venn-Abers Classification and Regression Calibration}
\description{
Calibrate binary and multiclass classification probabilities or real-valued
regression predictions using Inductive or Cross Venn-Abers methods.
Supply predictions from an existing model, or use optional tidymodels
workflows to train preprocessing, models and calibrators together.

Start with \code{vignette("method", package = "vennabers")}: it explains what
Venn-Abers predictors compute, which guarantees the literature gives and which
it does not, and checks the formulas against the package.
}
\section{Calibrating supplied predictions}{
The manual interfaces require no additional R packages. Use
\code{\link{venn_abers}} for Inductive classification calibration and
\code{\link{cross_venn_abers}} for matched out-of-fold probabilities.
Both support binary inputs and multiclass one-vs-one or one-vs-rest calibration.
For numeric predictions, use \code{\link{venn_abers_regression}} or
\code{\link{cross_venn_abers_regression}}.
\code{\link{venn_abers_one_hot}} represents class decisions as a 0/1 matrix
with one column per class.

Continue with \code{vignette("manual-calibration", package = "vennabers")}.
}
\section{Training and calibrating workflows}{
Use \code{\link{calibrated_workflow}} with an unfitted tidymodels workflow,
or \code{\link{calibrate_workflow}} to calibrate an already fitted workflow
on new data. These interfaces require the optional packages in DESCRIPTION.
\code{\link{model_adapter}} turns your own fit and predict functions into a
workflow model without registering a parsnip engine.
Follow \code{vignette("workflows", package = "vennabers")} for a complete
example, then \code{vignette("tuning", package = "vennabers")} for Grid/Bayes
selection and independent assessment.

\link{workflow_training} documents splits and prediction models;
\link{workflow_postprocessing} covers learned postprocessors.
Use \code{\link{venn_abers_details}} for prediction diagnostics and
\code{\link{venn_abers_tuning}} for search configuration and results.
}
\section{Example scripts}{
\link{vennabers_examples} indexes the installed standalone examples and their
annotated articles. Each collection has a contents list for navigating scripts:
\code{vignette("examples-manual", package = "vennabers")},
\code{vignette("examples-workflows", package = "vennabers")},
\code{vignette("examples-postprocessing", package = "vennabers")} and
\code{vignette("examples-tuning", package = "vennabers")}.
Use \code{system.file("examples", package = "vennabers")} to locate the scripts
for inspection or execution with \code{source()}.
}
\section{Assessment and interpretation}{
Score final classification probabilities with \code{\link{brier_score}} and
\code{\link{log_loss}}, and assess them on held-out pairs with
\code{\link{corp_plot}}, \code{\link{corp_decomposition_plot}} or
\code{\link{calibration_plot}}; assess regression predictions with
\code{\link{regression_plot}}. \code{\link[=plot.venn_abers]{plot()}} draws
the fit diagnosis of a fitted calibrator from its stored state, without new
truth values. \code{vignette("diagnostics", package = "vennabers")} compares
both kinds of graphics.
Keep assessment data separate from model training, selection and calibration.
Manual Cross calibration requires predictions from the corresponding retained
fold models.

Classification bounds are Venn-Abers probability bounds, not confidence
intervals. Regression bounds concern regularized conditional expectations,
not individual outcome coverage; Cross aggregate bounds are descriptive.
The installed \code{system.file("PARITY.md", package = "vennabers")} separates
numerical parity with the Python reference from R extensions and limitations.
}
\seealso{
\code{browseVignettes("vennabers")}, \code{help(package = "vennabers")}.

Cite the package and the method papers with \code{citation("vennabers")}.

Project: \url{https://github.com/RMKruse/vennabers}

Bug reports: \url{https://github.com/RMKruse/vennabers/issues}
}
\examples{
calibrator <- venn_abers(
  truth = c(FALSE, TRUE, FALSE, TRUE),
  probability = c(0.1, 0.4, 0.35, 0.8)
)
prediction <- predict(calibrator, c(0.2, 0.5, 0.9))
prediction
brier_score(c(FALSE, TRUE, TRUE), prediction$probability)
}
\keyword{package}

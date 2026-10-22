\name{venn_abers_details}
\alias{venn_abers_details}
\title{Inspect Calibrated Workflow Predictions}
\description{
Return native predictions together with calibration bounds, fold identities
and source-model diagnostics from a fitted trainable or frozen workflow.
Prediction uses the retained models and calibrators without learning again.
}
\usage{
venn_abers_details(object, new_data, ..., loss = NULL)
}
\arguments{
\item{object}{A fitted workflow returned by \code{\link{calibrated_workflow}}
and \code{generics::fit()}, or by \code{\link{calibrate_workflow}}.}
\item{new_data}{A data frame with required raw predictors. The same blueprint,
schema and inner recipe/engine checks apply as for \code{predict()}.}
\item{...}{Must be empty. Unknown, irrelevant and partially matched names are errors.}
\item{loss}{For classification, omit to use the stored default, or supply
\code{"log"} or \code{"brier"} to override it without refitting or changing
stored state. Explicit NULL is invalid. Regression rejects explicit loss.}
}
\value{
\code{venn_abers_details()} returns a list with \code{prediction} (the same native
probability table) and \code{details}: integer \code{row}, character \code{fold}
(\code{"Inductive1"} for Inductive; the supplied fold IDs for Cross), character \code{event_class}, and double
\code{lower}, \code{upper}, \code{probability}. Cross details contain one row per
new observation and fold; their probability is the fold point, while
\code{prediction} contains the final aggregate.
Trainable workflow details also contain character \code{prediction_model} and
\code{model_source}. The source is the fold ID in matched mode, or \code{"Full1"}
for every fold in full mode. These are model metadata, separate from original
class identities and the \code{fold} identifying the held-out calibrator;
\code{venn_abers_decode()} preserves them without decoding them as classes.
Trainable OvO/OvR details include \code{model_training}.
OvO/OvR details additionally contain character \code{method} and either
\code{class_1}/\code{class_2} (OvO) or \code{class} (OvR). Each row identifies
one new observation, fold and binary subproblem. \code{event_class} is the later class
in an OvO pair or the considered OvR class. These binary point probabilities
must first be merged across folds for Cross (Log/Brier), then coupled (OvO)
or normalized (OvR) to reconstruct \code{prediction}.
Regression returns the native \code{.pred} table as \code{prediction}, and
\code{details} with integer \code{row}, character \code{fold}, and double
\code{prediction}, \code{lower}, \code{upper}. Cross additionally returns
\code{aggregate} with \code{row}, \code{lower}, \code{upper}: descriptive
separate equal-weight means of fold bounds. Trainable regression aggregates also
contain character \code{prediction_model} and a \code{model_source} list column:
each row lists all contributing source model IDs (only \code{"Full1"} in full mode).
Zero/one rows retain these schemas. Frozen regression outputs retain their existing schema.
}
\section{Postprocessing diagnostics}{
For learned postprocessing, the additional \code{postprocessing} element
records local data memberships, fitted mappings and fit counts. Pairwise OvO
uses unique \code{pair_N/VA_fold} keys; ensemble prediction details retain
inner-fold names and sizes in list-columns. See \link{workflow_postprocessing}
for the complete record layout and method-specific diagnostics.
}
\section{Interpreting bounds}{
Classification bounds are Venn-Abers probability intervals, not confidence
intervals. Multiclass details describe binary subproblems, not multiclass
intervals. Regression bounds concern regularized conditional expectations,
not individual outcome coverage; Cross aggregate bounds are descriptive.
}
\seealso{
\code{\link{calibrated_workflow}}, \code{\link{calibrate_workflow}},
\code{\link{venn_abers_decode}}, \code{\link{venn_abers_tuning}}, and
\link{workflow_training} for the fitted engine's stored state.
For worked diagnostics, see \code{vignette("workflows", package = "vennabers")}.
}
\examples{
packages <- c("parsnip", "workflows", "recipes", "rsample", "hardhat",
              "generics", "rlang", "tibble", "tidyselect", "withr")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  development <- data.frame(
    x = sin(1:80),
    truth = factor(rep(c("yes", "no", "no", "yes"), 20),
                   levels = c("yes", "no"))
  )
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  specification <- calibrated_workflow(base, outcome = "truth",
                                       method = "cross", v = 3, seed = 42)
  fitted <- generics::fit(specification, development)
  new_data <- data.frame(x = c(-0.5, 0.5))

  # The final probabilities, then one record per new row and Cross fold.
  result <- venn_abers_details(fitted, new_data)
  print(result$prediction)
  result$details
}
}

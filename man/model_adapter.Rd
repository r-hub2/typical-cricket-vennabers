\name{model_adapter}
\alias{model_adapter}
\title{Own Fit and Predict Functions as a Workflow Model}
\description{Turn a pair of your own functions into an unfitted parsnip model
specification for \code{workflows} and \code{\link{calibrated_workflow}}. You
do not register a parsnip engine. The workflow keeps preprocessing, data
roles, calibration splits and calibration, and calls your functions on each
local training and prediction step.}
\usage{
model_adapter(mode, fit, predict)
}
\arguments{
\item{mode}{\code{"classification"} or \code{"regression"}.}
\item{fit}{A \code{function(x, y)} that trains a fresh model. \code{x} is a base
data frame of the predictors after the workflow's preprocessing, for the local
training rows only. \code{y} is the matching outcome: a factor with all stored
levels for classification, or a numeric vector for regression. Return any
object that holds everything \code{predict} needs.}
\item{predict}{A \code{function(object, new_data)}. \code{object} is the value
returned by \code{fit}. \code{new_data} is a base data frame with the same
preprocessing applied. For classification, return a numeric matrix with one
column per level of \code{y}, named by the original class names. Binary
classification also needs both columns. Rows must keep the input order. Values
must be finite, in [0, 1] and sum to one within 1e-8 in every row. For
regression, return a plain numeric vector with one finite value per row.}
}
\value{An unfitted \code{parsnip} model specification of class
\code{vennabers_adapter} with engine \code{"callbacks"}. Add it to a workflow
with \code{workflows::add_model()}. After fitting,
\code{workflows::extract_fit_engine()} returns a \code{vennabers_adapter_fit}
whose \code{state} element is the object your \code{fit} returned.}
\section{Lifecycle}{
Every local fit calls \code{fit} once, on the rows of that fit, after the
workflow learned its preprocessing on the same rows: each Inductive or Cross
source model, the \code{Full1} refit of \code{prediction_model = "full"}, and
each analysis set of an outer resample. Calibration rows never reach
\code{fit}. A fitted adapter workflow is an ordinary fitted workflow and can be
adopted frozen with \code{\link{calibrate_workflow}}.

Supported tasks are binary classification, shared multiclass OvO/OvR and
regression, with Inductive or Cross calibration, matched or full prediction
models and outer resampling. For shared multiclass, the callbacks must learn
and predict all local classes. Encodings and transformations belong in the
workflow; the adapter does not guess model-specific preprocessing.

The existing compositions work unchanged:
\itemize{
\item \code{model_training = "pairwise"}: a binary callback is fitted on each
class pair's rows, with \code{y} holding that pair's two original levels.
This does not make it a shared multiclass model.
\item \code{postprocessing}: \code{fit} learns on B only, for a holdout and for
every inner ensemble member; P and C never reach it. Sigmoid on
probabilities, Tailor logistic and multinomial on class probabilities, and
Tailor linear or numeric isotonic on regression predictions keep their task,
method and convergence rules.
\item \code{tuning}: Grid and Bayes searches tune calibration and declared
recipe parameters, each local search fitting fresh callbacks.
}

The callbacks receive no case weights, additional data or tuning parameters.
Bind fixed model options inside the functions. An adapter has no \code{tune()}
parameters: a search grid column for it is rejected, and so is
\code{update()} with any argument. \code{postprocessing$sigmoid_input = "glm_link"}
requires the native GLM engine, also when a callback calls \code{stats::glm()}.
}
\section{Predictions and errors}{
Class probability columns are aligned to the stored class order by name, never
by position, and returned as \code{.pred_<class>}. Class predictions of a plain
adapter workflow take the first maximum in level order. Calibrated workflows
apply their own class rules. Regression returns \code{.pred}. Empty
\code{new_data} returns an empty result without calling \code{predict}.

Malformed output is an error, never repaired: a non-matrix or unnamed result,
missing, duplicate or unknown class columns, a wrong row count, or invalid
probabilities. For regression, a matrix, a wrong length or non-finite values
are errors. Errors from your functions start with
\code{model_adapter() `fit` callback:} or
\code{model_adapter() `predict` callback:}. Calibrated workflows add the fold
or refit, for example \code{Inductive1}, \code{Fold2} or \code{Full1}, and the
class pair (\code{pair_2 Fold1}) or inner postprocessing fold
(\code{Inductive1 / Fold2 sigmoid}).
}
\section{Reproducibility and saving}{
The local \code{seed} of \code{\link{calibrated_workflow}} also applies to
callbacks that use R's random number generator. Callbacks must seed other
native generators themselves.

A saved fit stores the functions themselves, not their names. To use it in a
fresh R process, define the functions in a factory that binds their helpers
and fixed settings locally, and call package functions as \code{pkg::fun}.
Values looked up later in the global environment, or state shared and changed
between fits, are not supported. Models whose state needs their own save and
restore procedure, such as external pointers, are not covered.
}
\examples{
packages <- c("parsnip", "workflows", "recipes", "rsample", "withr")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  glm_adapter <- model_adapter("classification",
    fit = function(x, y) list(classes = levels(y),
      model = stats::glm(y ~ ., data = cbind(x, y = y), family = stats::binomial())),
    predict = function(object, new_data) {
      p <- stats::predict(object$model, new_data, type = "response")
      matrix(c(1 - p, p), ncol = 2L, dimnames = list(NULL, object$classes))
    })
  development <- data.frame(
    x = sin(1:80) + rep(c(1, 0, 0, 1), 20),
    truth = factor(rep(c("yes", "no", "no", "yes"), 20), levels = c("yes", "no"))
  )
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ x) |>
    workflows::add_model(glm_adapter)
  fitted <- generics::fit(calibrated_workflow(base, outcome = "truth", seed = 42), development)
  predict(fitted, data.frame(x = c(-0.5, 0.5)), type = "prob")
}
}
\seealso{\file{model-adapter.R} in \link{vennabers_examples} for GLM, multinomial
and linear regression adapters with saved fits, pairwise OvO, Sigmoid
postprocessing and a grid search.

\code{\link{calibrated_workflow}}, \code{\link{calibrate_workflow}},
\code{\link{workflow_training}}.
}

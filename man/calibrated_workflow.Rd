\name{calibrated_workflow}
\alias{calibrated_workflow}
\alias{predict._vennabers_fit}
\alias{predict.vennabers_workflow}
\title{Train a Calibrated Classification or Regression Workflow}
\description{
Wraps an unfitted tidymodels workflow for classification or regression so that
\code{generics::fit()} trains preprocessing, models and a Venn-Abers calibrator
on separate local splits of the same data. It needs the workflow and the
outcome column; the fitted workflow is then used with \code{predict()} on raw
predictors. Choose it when the package should fit the model and keep the
calibration rows out of training; use \code{\link{calibrate_workflow}} for a
workflow that is already fitted, or \code{\link{venn_abers}} and
\code{\link{venn_abers_regression}} for predictions computed elsewhere.
}
\usage{
calibrated_workflow(workflow, outcome, ..., predictors = NULL,
  method = "inductive", classification = "binary", event_level = "first",
  calibration_prop = 0.2, training_size = NULL, calibration_size = NULL,
  v = 5L, strata = NULL, splits = NULL, seed = NULL,
  precision = NULL, loss = "log", m = NULL, epsilon = NULL, schema = NULL,
  prediction_model = "matched", model_training = "shared", tuning = NULL,
  postprocessing = NULL)
\method{predict}{`_vennabers_fit`}(object, new_data, type = NULL,
  opts = list(), ..., loss = NULL)
}
\arguments{
\item{workflow}{An unfitted workflows object with a classification or regression model and
recipe, formula, or variable selection. Prepped or partly trained recipes,
unresolved tune placeholders without \code{tuning}, and postprocessors without an
explicit matching \code{postprocessing} request are rejected.}
\item{outcome}{One original outcome column name, or \code{schema$encoded_outcome}
when a schema is supplied. Fit requires a nonmissing factor
with exactly two levels for binary or at least two for OvO/OvR. All stored levels
must be represented in the development data. Regression requires finite numeric
truth; numeric 0/1 outcomes remain regression.}
\item{...}{Must be empty. Unknown, irrelevant and partially matched names are errors.}
\item{predictors}{Optional character vector of raw predictor names, excluding
outcome. NULL selects all other supplied columns. Their stored order is canonical;
unused selected columns still have to be supplied for prediction.}
\item{method}{Exactly \code{"inductive"} (default) or \code{"cross"}.}
\item{classification}{Exactly \code{"binary"} (default), \code{"ovo"}, or
\code{"ovr"}. OvO/OvR explicitly use the
factor/matrix core, including with two classes. Mode comes from the base model.
Regression rejects explicit classification, event_level, precision and loss options,
including their defaults; implicit classification defaults do not apply.}
\item{event_level}{Exactly \code{"first"} or \code{"second"}, identifying the
original factor level whose probability is calibrated. Only for binary;
explicitly supplying this option for OvO/OvR is an error.}
\item{calibration_prop}{A number strictly between zero and one.
\code{rsample::initial_split()} uses training proportion \code{1-calibration_prop}.
Explicitly supplying this option requires Inductive.}
\item{training_size, calibration_size}{Independent Inductive training and
calibration sizes, like sklearn's \code{train_size}/\code{test_size}. Both NULL
(default) keep \code{calibration_prop} and its draw unchanged. A value strictly
between zero and one is a proportion of the local rows; a whole number of at
least one is a row count. Proportions round training down and calibration up;
row counts are exact. A missing size is the complement. Rows left over when both
are supplied form the unassigned set U: the matched source model, its B/P
postprocessing and the calibrator never use them, while
\code{prediction_model = "full"} trains its additional workflow on all local
rows, including calibration rows and U. Classification draws are stratified by
outcome with exact totals (largest-remainder class allocation, ties in level
order) unless \code{strata = FALSE}; regression is unstratified. Sizes are
checked against the local rows and never shrunk: more rows than available, or an
empty training or calibration set, is an error naming both sizes. They cannot be
tuned and cannot be combined with \code{method = "cross"}, explicit
\code{calibration_prop} or \code{splits}; the default \code{calibration_prop}
is then ignored. With \code{model_training = "pairwise"} the sizes resolve on
each class pair's local rows, with U per pair; after calibration selection they
resolve on the final-fit F. Calibration tuning and joint model tuning fit
Venn-Abers candidates that resolve the same sizes on their resample rows, so
they accept only proportions; two-stage model tuning also accepts row counts.
See \file{inductive-sizes.R} in \code{\link{vennabers_examples}}.}
\item{v}{Cross fold count, a whole number of at least two, default five.
Automatic \code{rsample::vfold_cv()} uses one repetition and never reduces this
count; it must not exceed the local row count. Explicitly supplying \code{v}
requires Cross.}
\item{strata}{NULL stratifies classification by outcome with \code{pool=0}; regression
is unstratified. FALSE disables
stratification. Other values are rejected. No retry or automatic split adjustment.}
\item{splits}{NULL generates a local split. Otherwise a mixed or ordered holdout
\code{rsplit} from \code{rsample::initial_split()} or \code{rsample::make_splits()},
for Inductive; a simple \code{rset} from \code{rsample::vfold_cv()} or
\code{rsample::manual_rset()} for Cross; or a serializable function receiving only
selected local raw predictors followed by outcome and returning such a plan. Explicit plans cannot be combined with
explicit \code{calibration_prop}, \code{v} or \code{strata}, even their defaults.
Use \code{\link{local_calibration_splits}} for prepared assignments.
Pairwise accepts a function of the current local pair data, or for direct fits a
list named by exactly \code{pair_1}, \code{pair_2}, etc. A single global
multiclass plan is not a pair plan. Automatic holdouts/partitions are generated after filtering.
Explicit Cross plans may have different fold counts per pair.}
\item{seed}{NULL uses the caller's RNG stream. A finite whole number in R's integer
seed range sets a local seed for splitting and base fitting, restoring the caller's
state and RNG kind afterwards, including on error. The scope includes the plan
callback and model learning. Seed the outer plan separately. Reproducibility also
depends on RNG kind, engine, package versions, thread and parallel configuration.}
\item{precision}{NULL or a finite numeric whole number, including negative numbers.
Only calibration probabilities are rounded, following \code{venn_abers()}.}
\item{loss}{Constructor: stored prediction default, exactly \code{"log"} or
\code{"brier"}. Prediction/details: omit to use that default; an explicit choice
overrides it without refitting or changing stored state. Explicit NULL is invalid.
Regression rejects explicit loss at construction and prediction/details.}
\item{m, epsilon}{Regression moderation: both NULL means effective \code{m=1}.
Only one non-NULL option may be supplied and only that option reaches the core.
\code{m} is a positive numeric whole number; \code{epsilon} is strictly between
zero and one. Before any inner training, each fold of calibration size k must
satisfy \code{1 <= m <= floor((k-1)/2)}, with
\code{m = floor(epsilon*(k+1)/2)} for an epsilon request. Invalid folds fail the
entire fit with fold/size context; there is no clipping or fold-count adjustment.
Classification rejects explicit m/epsilon, including NULL.}
\item{schema}{NULL retains the direct naming contract. A schema from
\code{venn_abers_schema()} fixes the encoded input names, mode, class order and binary event.
Omitted classification/event options inherit the schema; contradictory explicit
options fail before learning. Build the base formula/recipe and resamples on
\code{venn_abers_encode()} data. All schema predictors must be selected.}
\item{prediction_model}{Exactly \code{"matched"} (default) or \code{"full"}.
Available for binary, shared OvO/OvR, pairwise OvO classification and regression.
This is a fixed fit option; predict cannot change it.
Matched uses each calibration source model for new observations. Full applies
the unchanged held-out calibrators to predictions from one additional workflow
trained on the entire permitted local fit domain, including calibration rows.
For pairwise training there is one such full workflow per local class pair.
With postprocessing, full means an additional complete B/P holdout or ensemble
on that domain, with its own inner splits; each engine still learns only on B.
Matching-model calibration guarantees do not automatically transfer to full
models; see \link{workflow_training}.}
\item{model_training}{Exactly \code{"shared"} (default) or \code{"pairwise"}.
Explicit values are rejected for binary and regression workflows. A shared model
produces one complete class matrix; OvO splits calibration into pairs without
training pair models. Pairwise requires \code{classification="ovo"},
and a workflow
that can train on a two-level factor. Each pair has its own models, recipes and
holdout or Cross partition. The later global level is its event; no extra shared/OvR model is fitted.}
\item{tuning}{NULL keeps fixed settings. Otherwise, a named list with
\code{strategy = "joint"} or \code{"two_stage"}, and \code{model} and/or
\code{calibration} Grid or Bayes searches matching the \code{tune()} placeholders.
See \code{\link{venn_abers_tuning}} for the configuration, selection splits,
pairwise search rules and diagnostics; start with
\code{vignette("tuning", package = "vennabers")} for examples.}
\item{postprocessing}{NULL leaves the unprocessed training path unchanged;
it does not remove a supplied postprocessor. Otherwise, a named list selecting
\code{method} (\code{"sigmoid"}, \code{"tailor_logistic"},
\code{"tailor_multinomial"}, \code{"tailor_linear"} or \code{"isotonic_numeric"})
and \code{strategy} (\code{"holdout"} or \code{"ensemble"}). Holdouts need
\code{prop} or local \code{splits}; ensembles need \code{v} or local \code{splits}.
Sigmoid also requires an explicit \code{sigmoid_input}.
See \link{workflow_postprocessing} for compatible methods, data roles and the
complete configuration.}
\item{object}{A fitted calibrated workflow. The predict method is normally
dispatched by \code{predict()} on that workflow.}
\item{new_data}{A data frame with required raw predictors. Normal workflow blueprint
and \code{hardhat::forge()} checks apply equally to prediction and details. Novel
factor levels reach the inner workflow, whose recipe/engine handles them or errors.}
\item{type}{Classification: \code{"prob"} or \code{"class"}, default class. Regression:
\code{"numeric"}, also the default, returns double \code{.pred} in outcome units.
Neither \code{"conf_int"} nor \code{"pred_int"} is registered.}
\item{opts}{The standard parsnip prediction options list; normally empty.}
}
\value{The constructor returns an unfitted ordinary workflow. After
\code{generics::fit()}, native probability predictions are tibbles of double
\code{.pred_<original class>} columns in original level order. Class predictions
contain \code{.pred_class}, a factor retaining all levels. Binary favors its event
at exactly 0.5; OvO/OvR select the first maximum in stored level order from the
final probabilities, even with two classes. Empty predictions keep the schema without calling the inner engine.
\code{\link{venn_abers_one_hot}} represents these class decisions as a one-hot matrix.

Regression predictions contain a double \code{.pred} column in outcome units.
Use \code{\link{venn_abers_details}} for bounds, fold identities and source models,
and \code{\link{venn_abers_tuning}} for selection diagnostics.
The fitted engine, \code{workflows::extract_fit_engine(fitted)}, is a
\code{vennabers_fit}: printing it shows source mappings and base-fit counts, and
\code{\link{plot.vennabers_fit}} draws its stored calibration map.
}
\section{Learning path}{
Start with the binary example below, then follow
\code{vignette("workflows", package = "vennabers")} through learned preprocessing,
Cross calibration, multiclass, regression, frozen workflows and postprocessing.
\code{vignette("tuning", package = "vennabers")} continues with parameter
selection, outer resampling and final assessment.

For individual rules, use these reference topics:
\describe{
\item{Training and evaluation}{\link{workflow_training}: dependencies, split
validation, stored state, pairwise models, full-model transfer and resampling.}
\item{Prediction diagnostics}{\code{\link{venn_abers_details}}: returned
columns, bounds and fold/source identities.}
\item{Tuning}{\code{\link{venn_abers_tuning}}: search configuration, selection
data, local searches and diagnostics.}
\item{Postprocessing}{\link{workflow_postprocessing}: supported methods,
B/P/C data roles, holdouts, ensembles and full compositions.}
}
For a model that is already fitted, use \code{\link{calibrate_workflow}}.
}
\section{Data roles and limitations}{
The split isolates new base learning from
calibration observations. An unfitted workflow does not prove that earlier tuning,
feature choices, or external transformations excluded these observations. Keep
assessment observations separate from all training, selection and calibration.
Calibration is not an independent evaluation. Cross model parameters must be
fixed using independent data or selected through the configured local
\code{tuning} searches inside the allowed training portions.
The statistical setting assumes IID observations and a predictive procedure fixed
independently of its calibration data. The bounds are Venn-Abers probability
intervals, not confidence intervals. Regression bounds concern regularized
conditional expectations, not individual outcome coverage; Cross aggregate bounds
are descriptive. Improved scores on every data set are not
promised. The engine cannot validate an outer partition it never receives.
Validate outer analysis/assessment disjointness before calling fit_resamples;
assessment sets from different outer resamples may overlap.
}
\examples{
packages <- c("parsnip", "workflows", "recipes", "rsample", "hardhat",
              "generics", "rlang", "tibble", "tidyselect", "withr")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  # Development rows are available for training and calibration.
  # Reserve separate assessment data before choosing this procedure.
  development <- data.frame(
    x = sin(1:80),
    truth = factor(rep(c("yes", "no", "no", "yes"), 20),
                   levels = c("yes", "no"))
  )
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))

  specification <- calibrated_workflow(base, outcome = "truth", seed = 42)
  fitted <- generics::fit(specification, development)
  new_data <- data.frame(x = c(-0.5, 0.5))
  predict(fitted, new_data, type = "prob")
  venn_abers_details(fitted, new_data)$details
}
}
\seealso{For a step-by-step introduction, start with
\code{vignette("workflows", package = "vennabers")};
\code{vignette("method", package = "vennabers")} explains the theory and its
guarantees. \link{vennabers_examples} lists annotated installed scripts.

\code{\link{venn_abers}}, \code{\link{cross_venn_abers}},
\code{\link{venn_abers_schema}}, \code{\link{calibrate_workflow}},
\code{\link{model_adapter}}.
}

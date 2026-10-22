\name{calibrate_workflow}
\alias{calibrate_workflow}
\alias{predict.vennabers_frozen}
\alias{print.vennabers_frozen}
\title{Calibrate a Frozen Workflow}
\description{
Adds an Inductive Venn-Abers calibrator to an already fitted classification or
regression workflow and fits only that calibrator, on the workflow's final
predictions; the complete workflow, including learned preprocessing and any
learned postprocessor, is retained unchanged. It needs the fitted workflow and
new calibration data that played no part in its training, tuning or selection.
Choose it when the workflow is final or expensive to refit; use
\code{\link{calibrated_workflow}} when the package should train the workflow
and create the calibration split itself.
}
\usage{
calibrate_workflow(fitted_workflow, data_calibration, outcome, ...,
  classification = "binary", event_level = "first", precision = NULL,
  loss = "log", m = NULL, epsilon = NULL, schema = NULL)
\method{predict}{vennabers_frozen}(object, new_data, type = NULL, ...,
  loss = NULL)
\method{print}{vennabers_frozen}(x, ...)
}
\arguments{
\item{fitted_workflow}{A fully fitted classification or regression workflow. Single
parsnip fits, extracted engines and unfitted workflows are rejected.}
\item{data_calibration}{A nonempty data frame in the original learned input
format, including the outcome. Do not pre-bake the data. Missing predictor
values are allowed only if the stored prediction path handles them.}
\item{outcome}{One column name matching the trained outcome. Calibration truth
must be a nonmissing factor with exactly the trained class universe. Factor
levels may be reordered; the stored workflow fixes the output class order.
Binary calibration may observe only one class. OvO/OvR must observe every
stored class globally. Regression truth must be a finite numeric vector,
including when all values are 0/1; the workflow fixes the mode.}
\item{classification}{Exactly \code{"binary"} (default), \code{"ovo"} or \code{"ovr"}.
Choose OvO/OvR explicitly or through the training schema, including for a two-class matrix problem.}
\item{event_level}{Exactly \code{"first"} or \code{"second"}, selecting the
event by the trained class order. Only accepted explicitly for binary calibration.}
\item{precision}{NULL or a finite numeric whole number. Only calibration
probabilities are rounded, following \code{\link{venn_abers}}. New probabilities
are never rounded. Changing this option or the event requires a new calibrator,
but no new base workflow fit.}
\item{loss}{Constructor: exactly \code{"log"} or \code{"brier"}, saved as the
prediction default. Predict and \code{venn_abers_details()}: omit to retain that
default, or override without refitting or changing stored state. Explicit NULL
is invalid at prediction. Only classification accepts loss or precision.
Regression rejects explicit classification, event_level, precision and loss
options, including explicit defaults or NULL.}
\item{m, epsilon}{Regression moderation. With both NULL, use m = 1. Otherwise
supply exactly one non-NULL value: a whole number with
1 <= m <= floor((k-1)/2), or 0 < epsilon < 1 giving
m = floor(epsilon*(k+1)/2) within that same range, where k is the calibration
size. At least three observations are needed. Requested values and effective
m and epsilon = 2*m/(k+1) are retained. Classification rejects these options.}
\item{schema}{NULL for an originally named workflow, or the fixed
\code{\link{venn_abers_schema}} used to encode the original training data.
With a schema, explicitly encode calibration and prediction data using
\code{venn_abers_encode(data, schema)} and use \code{schema$encoded_outcome}.
Omitted classification and binary event options inherit the schema; explicit
options must agree. The workflow mode, input columns, class order and ordered
factor status must match. A schema does not rename an existing workflow.}
\item{object, x}{A fitted predictor returned by \code{calibrate_workflow()}.}
\item{new_data}{A data frame in the learned input format. The outcome is not
needed and, if present, is removed before the original workflow predicts.}
\item{type}{NULL selects \code{"prob"} for classification and \code{"numeric"} for
regression. Classification also accepts \code{"class"}. Regression returns
only numeric points; \code{"conf_int"} and \code{"pred_int"} are unsupported.}
\item{...}{Must be empty. Constructor argument names must be exact. Unknown,
training, split, Cross and full-model options are rejected.}
}
\value{An already fitted \code{vennabers_frozen} predictor, also inheriting
\code{vennabers_fit}. It is not a trainable workflow or a specification for
\code{fit_resamples()}. The object retains complete workflows in \code{models},
the \code{calibrator}, \code{schema}, \code{class_levels}, \code{event_class}, input column
metadata, calibration size, requested/effective \code{parameters}, and relevant
R, package and engine \code{versions}. No additional complete copy of new raw
calibration data is stored; data already retained by the original workflow remain.
Prediction returns a tibble with named \code{.pred_<class>} probabilities in
trained class order, or a factor \code{.pred_class}, ordered exactly when the
trained outcome is ordered. The event wins at probability
0.5 in binary calibration. OvO/OvR ties choose the first stored level.
Regression returns one double \code{.pred} column.
\code{venn_abers_details()} returns the same prediction table and binary
lower/upper bounds and points with row, Inductive fold and event identity;
OvO/OvR add class-pair or class-against-rest identities. Regression details
contain row, fold, prediction, lower and upper. IVAR bounds concern regularized
conditional expectations, not confidence or prediction intervals.
Zero and one row preserve types and columns; zero rows do not call the base engine.
Print reports the adoption without base fits, the one retained workflow, the
number of calibration rows and the Inductive fold-to-model mapping.
\code{\link{plot.vennabers_fit}} draws the stored calibration map.}
\details{
Calibration uses the complete public workflow \code{predict(type = "prob")}
or \code{predict(type = "numeric")} path, including its existing postprocessor
exactly once. Missing, extra,
duplicate or unknown class columns, nonfinite/out-of-range probabilities,
invalid row sums and changed row counts cause errors. Multiclass matrices
are normalized only within the existing 1e-8 row-sum tolerance, by the manual
core. OvO zero pair sums, coupling, OvR normalization and calibration-only
rounding use those same cores. Regression requires finite final numeric values. Failures retain original
class/input names and prediction-stage context. No learned state is renamed,
refitted or bypassed through an extracted engine or a raw fallback.

The user is responsible for keeping the new calibration data outside all earlier
model, preprocessing, postprocessing and parameter-selection data. Assessment
data must also be separate from calibration and selection. This history is not
automatically verified. Frozen adoption is an R extension, neither Cross
calibration nor full-model transfer. \code{\link{calibrated_workflow}} is the
separate path that trains new base workflows on its own split.

Native \code{fit_resamples()} and \code{last_fit()} evaluate the base workflow
itself. Their fitted workflows, from \code{tune::extract_workflow()} or saved
extracts, can be adopted with calibration rows outside both their analysis and
assessment data; the adoption then learns only the calibrator. The returned
predictor is not itself a \code{tune} input.

To choose between the Log and Brier rules after adoption, compare
\code{predict(fitted, validation, type = "prob", loss = ...)} with yardstick
metrics on separate, previously unused validation rows, then pass the chosen
\code{loss} explicitly when predicting. Neither the workflow nor the calibrator
is refitted, and no further helper is needed. Validation rows remain selection
data; final assessment uses other rows. Automatic search with new training
belongs to \code{calibrated_workflow()} tuning.

\code{saveRDS()} and \code{readRDS()} preserve prediction and diagnostics without
learning for the tested GLM, multinomial and linear workflows, including recipes
and Tailor, in fresh R processes and parallel workers. Required
packages must be available in the new R process. Engines with external resources
retain their own storage contracts; no universal portability is promised.
The integration dependencies are optional; manual cores work without them.
The installed \code{examples/frozen-binary.R} demonstrates earlier model/recipe
and postprocessor learning, separate calibration, the loss-only choice on
validation rows and independent assessment.
\code{examples/frozen-multiclass-regression.R} demonstrates learned multinomial
and linear postprocessing with OvO/OvR and IVAR. The original workflow must
support the requested batch size: for example, the tested probably 1.2.0
multinomial smooth=FALSE path fails on one row; smooth=TRUE supports it.
Original prediction failures propagate without modifying the workflow.

For encoded training, retain the original schema together with the fitted
workflow. No schema is learned from later calibration data. Inputs must carry
the encoded names; valid reordered input columns and factor levels are matched
by identity. Predictions and details carry the schema binding and remain
encoded. \code{venn_abers_decode()} returns the existing roles with original
names and classes, including diagnostics; it invents neither outer identities
nor Cross aggregates for these Inductive methods. See the installed
\code{examples/frozen-schema.R} for all four paths with learned postprocessing.

Bindings retained on the workflow, its original recipe template, or incoming
data are checked, along with learned metadata. Workflows and table operations
can discard attributes: matching detached codes alone cannot prove which schema
was used during training. The caller must supply the actually associated schema.
An originally named workflow must continue receiving its original inputs with
\code{schema = NULL}. Creating a schema afterward does not authorize recipe
rewriting, refitting or encoded input substitution. An old workflow whose
complete prediction path fails on reserved names reports the stage, names and
underlying cause. The limited binary GLM raw-response probe without a
postprocessor is not a product fallback; no engine or postprocessor is bypassed.
}
\seealso{\code{vignette("workflows", package = "vennabers")} for frozen
workflows in context, \code{vignette("method", package = "vennabers")} for the
theory and its guarantees, and \link{vennabers_examples} for annotated
installed scripts.

\code{\link{calibrated_workflow}}, \code{\link{venn_abers_details}},
\code{\link{venn_abers_schema}}, \code{\link{venn_abers}}, \code{\link{venn_abers_regression}}}
\examples{
if (all(vapply(c("parsnip", "workflows", "hardhat", "rsample", "recipes",
                "rlang", "tibble", "withr", "generics", "tidyselect"),
               requireNamespace, logical(1), quietly = TRUE))) {
  d <- data.frame(x = sin(1:80), truth = factor(rep(c("yes", "no", "no", "yes"), 20)))
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  earlier <- generics::fit(base, d[1:40, ])
  calibrated <- calibrate_workflow(earlier, d[41:60, ], "truth")
  predict(calibrated, d[61:80, "x", drop = FALSE])
}
}

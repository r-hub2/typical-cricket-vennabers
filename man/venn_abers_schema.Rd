\name{venn_abers_schema}
\alias{venn_abers_schema}
\alias{venn_abers_encode}
\alias{venn_abers_decode}
\title{Reversible Names for Calibrated Workflows}
\description{Explicitly encode input columns and outcome classes before
constructing a base workflow, and decode native results into separate roles.
These helpers learn no preprocessing, statistics, or calibration.}
\usage{
venn_abers_schema(data, outcome, classification = "binary", event_level = "first")
venn_abers_encode(data, schema)
venn_abers_decode(x, schema)
}
\arguments{
\item{data}{A data frame containing exactly the selected input columns. At schema
construction the outcome must be a factor with exactly two declared levels for
\code{"binary"}, at least two for \code{"ovo"}/\code{"ovr"}, or a finite numeric
vector for regression (including numeric 0/1).
Encoding accepts new data without the outcome, but requires every predictor.}
\item{outcome}{One exact original outcome name, including reserved or
nonsyntactic names.}
\item{classification}{Exactly \code{"binary"}, \code{"ovo"}, or \code{"ovr"}
for factor outcomes. Must be omitted for numeric regression outcomes.}
\item{event_level}{Exactly \code{"first"} or \code{"second"}, selecting the event
from the original factor levels in their stored order. Only binary classification
accepts an explicit event option; omit it for OvO/OvR and regression.}
\item{schema}{The fixed object returned by \code{venn_abers_schema()}. Pass the
same schema to encoding, \code{calibrated_workflow()} or
\code{\link{calibrate_workflow}}, and decoding.}
\item{x}{A native probability, class or numeric prediction data frame, an unsummarized
\code{tune::collect_predictions()} table, or a \code{venn_abers_details()} result.}
}
\value{
\code{venn_abers_schema()} returns a \code{vennabers_schema} metadata list.
\code{version} is the format version; \code{outcome} and \code{encoded_outcome}
are the original and technical target names. \code{columns} maps original column
names to codes, in \code{original_columns} order. \code{classes} maps original
class names to codes, in \code{class_levels} order. \code{ordered} records whether
the original target is ordered. \code{mode}, \code{classification},
\code{event_level}, and original \code{event_class} define the task. Multiclass schemas have NULL event fields; regression
has NULL class levels, classification and event fields, and an empty class map.
No observation values or additional calibration data are retained in the schema.

\code{venn_abers_encode()} returns a data frame in schema column order with
encoded column names and, when present, encoded target levels. Predictor values,
types, attributes, factor levels, row order, and row names are preserved.
The target is aligned to the saved class order even if incoming levels are
reordered or known levels are absent. Missing factor target values stay missing; numeric targets must be finite;
the fitting boundary rejects missing truth.

\code{venn_abers_decode()} returns a list with \code{outcome} (the original
outcome name as metadata), and only the result roles available in \code{x}:
\code{probability} (a numeric matrix with original class columns in schema order),
\code{class} (an original-level factor), \code{truth} (an original-level factor
or numeric regression truth), \code{prediction} (numeric regression points),
\code{identity} (available \code{.row}, \code{id}, \code{id2}, further numbered
\code{id} columns, and \code{.config}), \code{details} (existing binary, multiclass or regression diagnostics), and
\code{aggregate} (existing Cross regression mean bounds). Multiclass details
decode \code{class_1}/\code{class_2} for OvO, \code{class} for OvR, and
\code{event_class}; fold, row and subproblem identities remain unchanged. It preserves numbers, dimensions, row order,
and identities, including zero and one row. No class is inferred from probabilities,
and no truth or details are inferred from a prediction table. Regression never
adds class or probability roles. Multiclass ties favor the first maximum in the
stored class order, also with two classes. The decoder preserves decisions
already made by the predictor; it does not recompute them.
}
\details{
Select the intended columns before creating the schema. All selected names and
all declared outcome levels get deterministic, collision-free technical codes.
Mapping uses exact names, never name repair, Unicode normalization, case folding,
or frequency sorting. Even original names that resemble codes remain distinct.
Empty, missing or duplicate names, unknown columns or classes, missing predictors,
and damaged or incomplete schema mappings are errors. The schema is fixed across
folds, prediction batches and processes; do not edit it after construction.

Build formulas, recipes, variable selections, and outer resamples explicitly on
the encoded data. Use \code{schema$encoded_outcome} in the workflow constructor.
Existing recipes are not rewritten. To retain prescribed outer memberships,
rebuild their splits on the encoded data with the same row positions and IDs.
Inner plan functions receive local encoded data; the exact local plan/data
comparison and all training, calibration and assessment separation rules remain.

For a frozen workflow originally trained on encoded data, pass that original
schema to \code{calibrate_workflow()}, with encoded calibration data and
\code{schema$encoded_outcome}. Encode new prediction data with the same schema.
The stored workflow is never renamed or refitted. Without encoded training,
keep original input names and \code{schema = NULL}; a new schema cannot repair
an old workflow's failed complete prediction path. See
\code{examples/frozen-schema.R} for the four frozen Inductive methods.

Native predictions remain encoded. Separate decoded roles avoid collisions
between an original class \code{class}, outcome \code{.pred_class}, and external
\code{id}/\code{id2}. The decoder does not flatten or repair these names.
Probability columns must cover every schema class exactly once; named columns
are reordered to schema order. Unknown encoded output names or class codes fail.

Encoded data, native predictions, and details carry a \code{vennabers_schema}
attribute. Existing foreign bindings are rejected. Ordinary table operations or
collectors may drop attributes: detached tables therefore require the actually
associated schema supplied by the caller. Matching codes alone cannot establish
provenance. Decode unsummarized predictions to preserve repetitions; the decoder
cannot recover rows already summarized or dropped. Collect notes and unsummarized
metrics separately, audit completeness, and explicitly extract fitted states
when inner models or diagnostics are needed. Match \code{control_resamples()}'s
\code{event_level} to binary schemas; no control is silently changed.

Schema helpers use base R and do not load optional integration packages.
The schema and fitted pure R workflow state survive RDS restoration in a fresh
process with the required engine namespaces loaded. Outer workers require
\code{pkgs = c("vennabers", "stats")} for GLM/LM examples (add \code{"nnet"}
for the installed multinomial example),
\code{parallel_over = "resamples"}, and consistent seeds and thread settings.
}
\examples{
source(system.file("examples", "encoded-binary.R", package = "vennabers"))
source(system.file("examples", "encoded-multiclass-regression.R", package = "vennabers"))
}
\seealso{\code{\link{calibrated_workflow}}, \code{\link{local_calibration_splits}}}

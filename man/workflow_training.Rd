\name{workflow_training}
\alias{workflow_training}
\title{Calibrated Workflow Training Reference}
\description{
Split validation, stored models, outer resampling, pairwise training and
full-model transfer for \code{\link{calibrated_workflow}}. For a worked
introduction, start with \code{vignette("workflows", package = "vennabers")}.
These are reference rules for the trainable workflow, not a separate constructor.
}
\section{Supported scope and dependencies}{
All eight combinations of binary, OvO, OvR and regression with Inductive or
Cross training are supported, including complete outer resampling.
The integration packages in DESCRIPTION are optional Suggests and must be
installed explicitly; manual calibration needs none of them. Use the declared
minimum versions or newer supported versions. Tuning additionally needs
\code{tune}, \code{dials} and \code{yardstick}; Tailor postprocessing needs
\code{tailor} and \code{probably}. Each path checks its packages and minimum
versions at construction and again at fit, and stops with the missing package.
Basis engines also need their own packages: the examples use \code{stats} for
GLM/LM, \code{nnet} for a shared multinomial model and \code{rpart} for
decision trees. Fresh parallel workers must load \code{vennabers} and those
engine namespaces, for example through \code{control_resamples(pkgs = ...)}.
Reversible name encoding for binary, OvO/OvR and regression is available through \code{venn_abers_schema()}.
Frozen classification/regression workflow adoption is available through \code{\link{calibrate_workflow}}.
Own fit and predict functions enter as a workflow model through \code{\link{model_adapter}}.
Ordered reference holdouts and binary/shared multiclass full-model transfer are supported.
Pairwise OvO training supports matched/full Inductive and Cross calibration.
Binary, shared OvO/OvR, pairwise OvO and regression paths support joint or
two-stage Grid and Bayes selection, and learned postprocessing. Shared multiclass
training does not establish complete Python training parity.
\code{system.file("PARITY.md", package = "vennabers")} separates numerical
parity, training/data parity, R extensions and reference deviations;
the installed \code{NEWS.md} states the 1.x compatibility promise.
}
\section{Preprocessing and calibration splits}{
All learned preprocessing belongs in the base workflow; the outer Variables
preprocessor only carries raw data. Training and calibration positions must form a
complete nonempty disjoint partition. Binary allows singleton and one-class calibration sets.
Inductive OvO/OvR require every stored class in the holdout. Cross requires every
class globally in the out-of-fold calibration; individual folds may omit classes.
Nonempty one-class or singleton binary subcalibrations remain valid. An empty
OvO pair in any fold fails before any base training, naming the pair and fold.
Splits are never repaired. Trainability of the base engine is a separate requirement, and engine
errors identify the training or probability-prediction stage.
}
\section{Engine options and model-internal validation}{
Fixed, data-independent model options belong in the base workflow's model
specification: main arguments and \code{parsnip::set_engine()} arguments.
They travel with the base workflow, so every local fit uses them on its own
rows. \code{calibrated_workflow()} and \code{fit()} accept no extra fit
arguments; there is no interface for separately supplied validation data.
An engine's internal validation, such as native XGBoost early stopping, splits
its share from the prepared rows of each local fit, after the base recipe was
learned on all of them. It is not an independent assessment of the procedure.
\tabular{ll}{
  \strong{Local fit} \tab \strong{Rows for learning and internal validation} \cr
  Inductive/Cross source \tab its analysis rows, or their B with postprocessing; never C \cr
  Class pair \tab the analysis rows of that pair with its two original levels \cr
  Full1 \tab all local fit rows (F), or its own B with postprocessing \cr
  Search candidate \tab its local training rows; held-out rows stay out \cr
  Frozen workflow \tab none: no base training, only calibration \cr
}
P trains only the postprocessor and outer assessment rows never reach a fit.
Invalid engine settings fail the whole fit with its fold, pair or B/P context;
there is no fallback. Validation of engine-specific values stays with the
engine. No support is promised for observation weights, offsets or an external
\code{eval_set}. Not run, because xgboost is not a vennabers dependency:
\preformatted{spec <- parsnip::boost_tree(mode = "classification",
  trees = 100L, stop_iter = 10L) |>
  parsnip::set_engine("xgboost", validation = 0.2)
base <- workflows::workflow() |>
  workflows::add_formula(truth ~ .) |>
  workflows::add_model(spec)
fitted <- calibrated_workflow(base, outcome = "truth", seed = 42) |>
  generics::fit(data = training_data)
}
}
\section{Shared multiclass probabilities}{
With \code{model_training = "shared"}, multiclass uses one trained base workflow per split complement. Its probability columns must
be exactly \code{.pred_} plus each original class name, unique and complete;
column order is matched by name. Values must be finite in [0,1], with row sums
within 1e-8 of one. The existing matrix core normalizes within this tolerance.
OvO normalizes each pair and uses a neutral 0.5 score for zero-sum pairs before
optional native fit rounding (so 0.5 becomes 0 at precision 0). OvR uses all
calibration observations for each class against the rest. New probabilities
remain unrounded. Coupling, normalization and Log/Brier regularization are the
same as in the manual core; no separate pair models are trained.
}
\section{Names and schemas}{
Input and class names must be unique, nonmissing and nonempty. Class \code{class}
and names starting with \code{.pred_} are reserved. Outcome names \code{.pred},
\code{.pred_*}, \code{.row}, \code{.config}, and \code{id} or \code{id} followed
by a positive integer are reserved. Without a schema, names are never repaired; distinct names such
as \code{a b} and \code{a.b} require a compatible base engine. Raw predictor
\code{id} remains valid. Details do not join raw observations to predictions.
}
\section{Stored state and prediction models}{
\code{workflows::extract_fit_engine()} returns a \code{vennabers_fit} with:
\code{schema} (NULL or the explicitly supplied metadata),
\code{parameters} (\code{requested}, logical \code{explicit} flags and
\code{effective}), \code{outcome}, \code{predictors}, \code{columns} (predictors
then outcome), \code{class_levels}, \code{event_class}, \code{loss},
\code{folds} (by fold ID, with integer \code{analysis}/\code{assessment} positions and
\code{training_size}/\code{calibration_size}; with requested independent sizes,
\code{Inductive1} also stores integer \code{unassigned} positions, and
\code{effective} holds the resolved \code{training_size}, \code{calibration_size}
and \code{unassigned_size} row counts; fits without the field have no
unassigned rows), \code{models} (complete fitted
prediction sources named by the same IDs, or only \code{Full1} in full mode;
learned postprocessing uses package-owned compositions),
\code{calibrator}, \code{calibration} (for binary: a list named by fold ID, each
with original local integer \code{row}, logical event \code{truth} and unrounded
double held-out \code{probability}; for multiclass: one list per fold with
integer \code{row}, factor \code{truth} retaining all levels, and the complete
unrounded named class matrix \code{probability}; for regression: a data frame per
fold with integer \code{row}, numeric \code{truth} and held-out double
\code{prediction}), \code{versions}, and \code{rng_kind}.
Full mode fits the source recipes/models only on each analysis portion, then
fits a fresh complete recipe/model on all local fit rows. No source recipe
statistics are reused. Only after success are temporary source workflows
discarded. Held-out values and the separate calibrators never change; there is
no pooled OOF calibration or in-sample replacement. Training failure invalidates
the whole fit; there is no fallback. A one-row full-model prediction validates
its output contract during fit and is discarded without calibration.
Print identifies the mode and fold-to-model
mapping, and the training, calibration and unassigned row counts when
independent sizes were requested. Saved RDS objects retain the complete prediction procedure.
The cost is one/two base fits for matched/full Inductive and K/K+1 for Cross.
Full prediction calls its one model once per nonempty batch, reusing that output
across calibrators within the call. There is no persistent prediction cache.
The existing Log/Brier mergers and prediction-time loss override are unchanged.
For regression, \code{parameters$requested$mode} is \code{"regression"};
\code{effective$m} and \code{effective$epsilon} are numeric vectors named by
fold ID. Each \code{folds[[id]]$moderation} stores \code{requested}, effective
\code{m} and effective \code{epsilon=2*m/(k+1)}, alongside the fold sizes.
Implicit classification options are absent from effective settings; class levels,
event class and loss are NULL. The calibrator is the unchanged IVAR or CVAR core.
IVAR uses weighted points, including the midpoint when both weights vanish.
CVAR averages complete fold points, lower bounds and upper bounds separately
with equal weights, regardless of fold size; it never merges averaged bounds.
Matched mode retains each source workflow/recipe. Full mode retains one fresh
full workflow and applies its numeric prediction to every held-out calibrator.
No rounding or probability constraint is applied to regression predictions.
For OvO/OvR, the global \code{event_class} and effective \code{event_level} are
NULL; events belong to subproblems in details. The calibrator stores the original
classes and chosen multiclass method.
For Cross, \code{calibrator} is a \code{cross_venn_abers} core with
\code{folds}. Binary has a correspondingly named \code{calibrators} list;
multiclass has one binary Cross core per stored pair or class, each retaining
its fold IDs and corresponding ordered Inductive calibrators.
Inspect learned recipes with the public workflow extractors; predict with full
workflows. RDS restoration for pure R engines works in a fresh process after
loading vennabers, workflows, and required engine packages. External-pointer engines
may require their own serialization support.
}
\section{Explicit split validation}{
For an explicit plan, selected raw columns (in canonical predictor order, then
outcome), values, column types, factor levels and row order must match local fit
data. Container class and row names are ignored; equal observations still occupy
distinct positions. Training and calibration indices must be finite, whole,
unique, in range, disjoint, nonempty and cover all local rows. The single internal
fold ID is always \code{Inductive1}; rsample's creation/outer IDs are metadata,
not configurable inner fold IDs. Cross requires unique nonempty character IDs,
at least two folds, exactly one assessment assignment for every local row, and
the exact complementary analysis set in each fold. All folds are checked before
any learning starts. Repeated inner folds and bootstrap duplicates are rejected.
Ordinary and classification-stratified partitions may be mixed or non-shuffled.
For an ordered Inductive holdout use \code{make_splits()} with the first
\code{n - ceiling(n * q)} local rows for training and the remaining
\code{ceiling(n * q)} for calibration, requiring both parts to be nonempty.
This example convention rounds the calibration size upwards: 11 rows and
\code{q = .25} give training 1--8 and calibration 9--11. Explicit plans are
used as supplied, without sorting, mixing or stratification. Order is local to
each outer analysis set. Non-shuffled Cross folds retain exact training
complements, including later rows, and are not forward validation.
Group, time and other specialized split classes
are outside this partition contract. An unbound global positions plan cannot be reused
on different outer analysis data.
See installed \code{examples/reference-splits.R} for all eight training paths,
known ordinary/stratified memberships, data roles and assessment scores.
Equal R/Python seeds alone do not establish equal splits or training; numerical
comparisons require identical model outputs, assignments and options.
}
\section{Cross aggregation}{
Matched Cross retains exactly one full base workflow per fold and all corresponding
binary subproblem calibrators, with no
additional full-data model. New data traverse every matching learned recipe/model
path. A failed fold fails the entire fit. Log loss uses the core's stable logarithmic
aggregation of endpoints; Brier uses the equal-weight mean of fold Brier points.
Neither merger weights by fold size or pools calibration data. Precision rounds
only calibration probabilities, never new model probabilities. Fold details
reconstruct the final probability using these mergers. For multiclass, merge each
binary subproblem across folds first, then couple OvO or normalize OvR. Averaging
already coupled/normalized multiclass fold matrices is a different procedure.
This shared R fold plan does not establish Python pair-model training parity.
}
\section{Outer resampling and assessment}{
Ordinary \code{tune::fit_resamples()} fits the complete predictor independently
on every outer analysis set. The stored loss controls native predictions even
when tune calls parsnip directly; there is no loss argument to fit_resamples and
no implicit control modification. For binary, match \code{event_level} in the constructor and
\code{tune::control_resamples()}. Yardstick's binary Brier is the binary vector
score. For OvO/OvR, yardstick Brier is exactly half the manual matrix score,
including the two-class matrix contract. Its log-loss clipping affects scoring
only, never stored predictions. Regression uses normal yardstick RMSE/MAE
(\code{.estimator = "standard"}) in outcome units, scored on the complete
procedure's held-out predictions.

Use \code{save_pred = TRUE} and all public collectors, including
\code{collect_notes()}, to check complete unsummarized resample, repeat,
configuration and assessment-row identities and finite metrics. Inner failures,
metric errors, NA removal or missing rows mean incomplete evaluation. All fits
failing may make a collector itself error; a mean of surviving fits is insufficient.
The installed \code{examples/calibrated-resampling.R} shows binary Inductive/Cross
and multiclass Inductive/Cross OvO/OvR, plus both regression methods with an
executable completeness check. Five outer by five inner Cross folds require exactly
25 base fits; no final fit runs automatically. Keeping extracted Cross states retains
all 25 workflows and their recipes/calibrators, so memory scales with the number
of retained base fits. A separate final development fit costs another five fits.

Use \code{control_resamples(extract = workflows::extract_fit_engine)} for all
fitted engine states, \code{extract = identity} for full fitted workflows, or a
function selecting compact diagnostics. \code{save_workflow = TRUE} alone does
not preserve each fitted composite. Inner training is serial; use
\code{parallel_over = "resamples"} for outer parallelism. Fresh workers must load
vennabers and required base-engine namespaces through \code{pkgs}. Reproducibility
is checked for the same versions, seeds, threads, strategy and deterministic
engine; changing strategies or external engine RNGs is not covered.
}
\section{Optional package registration}{
Integration packages are optional Suggests and version-checked on use. Neither
core loading nor manual calibration requires or loads parsnip. Registration uses
public parsnip APIs, is idempotent in either load order, and rejects foreign model
registrations under the name \code{vennabers}. Packages are never auto-installed.
}
\section{Pairwise model training}{
Global factor levels fix the class order before decomposition. Each stable
\code{pair_id} names its two classes through
\code{pair_models[[pair_id]]$class_levels}; the second is always the event.
Each complete binary workflow
learns exclusively from its pair's training rows, including all recipe estimates.
Only its held-out event probabilities fit the existing binary calibrator.
For Cross, the existing binary core first aggregates within each pair using stable
log-loss aggregation or equal-weight Brier points. Only then are pair points
combined with the same reciprocal OvO merger as shared training; coupled
multiclass fold predictions are never averaged.
The global schema remains unchanged; only the unfitted recipe's outcome type is
adapted to the two-level local domain, retaining steps, roles and blueprint options.

\code{pair_models[[pair_id]]} retains \code{pair_id}, \code{class_levels},
\code{event_class}, \code{rows} (positions in the current global fit domain),
\code{folds} (positions in the local pair, plus sizes and, with requested
independent sizes, the pair's \code{unassigned} positions; the resolved counts
in \code{effective} are then named by pair ID), \code{models},
\code{calibrator}, and \code{calibration}. Folds, models and calibration records
are named by local fold ID: \code{Inductive1} or the supplied Cross IDs.
The Cross calibrator retains its ordered, named binary fold calibrators.
Prediction re-derives the pair order from the global class levels and rejects
components whose \code{pair_id}, \code{class_levels} or fold names disagree,
as well as models whose fitted two-class outcome mold is not the pair's.
These replace the shared model/fold/calibrator fields. Calibration records retain
local row positions, binary truths and unrounded probabilities. Details add
\code{pair_id}, \code{model_training="pairwise"},
\code{prediction_model}, and \code{model_source} equal to the pair ID in matched mode.
The pair ID together with \code{fold} identifies the source workflow
and calibrator in matched mode. Full mode instead names each new model source
\code{pair_N/Full1}, separately from the unchanged pair/fold calibrator. Decoding changes class identities while preserving pair/fold/source
IDs. Zero-row outputs retain all columns without calling any base model.

All pair plans and data bindings are validated before the first base fit.
For Inductive, both classes must occur in each pair's training and calibration
portions. Cross permits nonempty one-sided calibration folds when the training
complement supports a trainable engine returning correct two-class probabilities.
Within each pair, assessment folds partition every local row exactly once and
training is the exact complement. Repeated/overlapping or incomplete plans,
foreign data and missing/duplicate IDs fail before training. Automatic Cross uses
the requested \code{v} without reduction; explicit valid partitions may have
different fold counts across pairs. Equal fold IDs in different pairs have no
shared meaning. Ordinary/stratified and non-shuffled reference partitions follow
the same local IID split rules as shared training.
Observations may recur across pairs. Explicit plans retain their local position
order; an ordered pair holdout does not imply a global time cutoff. With a fixed
\code{seed}, each pair resolver/automatic plan starts from that same seed
in a restored local scope, so plans do not depend on preceding fits or call counters.
Base fits then consume the enclosing fit's seeded stream in stable pair order.
With NULL, splitting and fitting consume the caller's stream. Prepared assignments
use \code{local_calibration_splits()} over every local pair, which rejects
conflicting plans for identical data. Resampling callbacks see only the current
outer analysis pair; development-wide plans cannot be replayed on different data.

For K classes, P=K*(K-1)/2 pairs cost P base fits for Inductive and the sum of
local fold counts for Cross (three classes with five folds cost 15 fits).
Each nonempty prediction call uses each retained model exactly once. Recipes, calibrators and models remain stored; prediction never learns.
With \code{prediction_model="full"}, all source workflows and calibrators are
built first, then one extra full workflow learns its recipe/model on each complete
local two-class domain. This order preserves source RNG draws and held-out values.
After every full workflow succeeds, \code{models} retains only \code{Full1} per
pair; all calibrators, held-out records and identities remain. No in-sample
recalibration or OOF pooling occurs. Full Inductive costs 2P fits; full Cross costs
the sum of local fold counts plus P (18 for three classes and five folds).
Each nonempty prediction uses each pair full model once and shares its named event
probability across that pair's calibrators before binary aggregation and global
coupling. Outer assessment/test rows never enter these fits.
Any failed component fails the whole fit. See installed
\code{examples/pairwise-inductive.R} for separate binary references, explicit
pair plans, independent coupling, decoding and complete outer evaluation.
The installed \code{examples/pairwise-cross.R} adds explicit unequal Cross
partitions, reconstruction, native/decoded diagnostics and outer fit costs.
Both examples compare matched/full predictions, costs and held-out scores;
matching-model guarantees do not automatically transfer and better quality is not promised.
}
\section{Full-model transfer}{
For \code{prediction_model="full"}, the split isolates
calibration observations from source learning only. The extra full model includes
them. Existing matching-model calibration or coverage guarantees do not
automatically transfer to this differently trained model; improved quality is
not guaranteed. Binary, shared multiclass and regression grids give Full1 its own local model/recipe search;
the training method and prediction-model mode remain fixed per run.
\code{fit_resamples()} rebuilds the complete procedure inside each outer analysis
domain. \code{last_fit()} trains solely on its analysis domain, never analysis
plus test. \code{calibrate_workflow()} remains a frozen-model path and does not
accept this training option. Shared multiclass full-model transfer is an R
extension, not Python pair-model training parity. It uses the same full class
matrix for every subproblem and fold, with binary aggregation before multiclass
combination. See the installed examples \code{full-model-workflow.R} and
\code{full-multiclass-workflow.R} under \code{system.file("examples", package="vennabers")};
the latter shows native/decoded diagnostics, repeated outer assessment,
\code{last_fit()} and the unchanged yardstick Brier factor of two.
\code{full-regression-workflow.R} demonstrates IVAR/CVAR transfer, RMSE/MAE,
decoded numeric predictions and reconstruction of descriptive bounds. The Python
1.5.4 full-model regression option belongs to the exported lower
\code{VennAbersCV} class, not its regression facade.
}
\section{Printed costs}{
The fitted engine prints source mappings, the number of base
fits (including temporary source fits discarded in full mode), models used per
nonempty prediction, and stored calibration rows containing held-out truth and
model predictions. Pairwise calibration rows count membership in each pair.
}
\seealso{
\code{\link{calibrated_workflow}} for arguments and data-role limitations,
\code{\link{venn_abers_details}} for prediction diagnostics,
\code{\link{local_calibration_splits}} for prepared plans,
\link{workflow_postprocessing} for B/P compositions, and
\code{\link{venn_abers_tuning}} for local parameter searches.
}

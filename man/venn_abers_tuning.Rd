\name{venn_abers_tuning}
\alias{venn_abers_tuning}
\title{Inspect Classification or Regression Grid or Bayes Selection}
\usage{
venn_abers_tuning(object)
}
\arguments{
\item{object}{A fitted binary, shared multiclass, pairwise OvO, or regression \code{calibrated_workflow()} with grid or Bayes tuning.}
}
\description{
Returns the stored configuration and selection diagnostics without refitting.
This page also documents the \code{tuning} configuration passed to
\code{\link{calibrated_workflow}}. For a worked progression from Grid to
Bayes and outer assessment, see \code{vignette("tuning", package = "vennabers")}.
}
\value{
A list with \code{config}, \code{calibration}, \code{models}, \code{seed},
class/schema and method metadata, final fit roles and model source mappings.
Each search records candidate IDs and grid values, eligibility and exclusion reasons, per-resample
main and report scores, native \code{tune_grid()} or \code{tune_bayes()} results (including prediction
rows and notes), the stable candidate-order tie break, selected values, local plan positions,
gross base-fit budget estimate, attempted evaluation count and observed base-fit
attempts, including failures. Bayes diagnostics additionally retain each public one-step
continuation, its stop reason, elapsed time, and the global no-improvement count. The estimate is unavailable when a provider
function determines a nested plan size at fit time; its count factors remain
visible. In Bayes candidate tables, parameter IDs that collide with
\code{status}, \code{estimate} or \code{reason} receive unique
numeric suffixes separated by an underscore (for example, \code{status_1}).
Selected values and native results keep the original parameter IDs. The ID
\code{id} itself is unsupported because native tune uses it for resample
identities.
Calibration selection also records the
disjoint selection (S) and final fit (F) row positions. In full-model mode,
\code{models$Full1} is its separate local search. With postprocessing, each
\code{models} entry is a list of local B searches named by inner fold, with
\code{fold_id}, \code{inner_fold}, B positions and \code{search_seed} (one draw per B from the fit stream, so each search has its own RNG stream); budgets report \code{inner_pairs}.
For regression calibration searches, saved native predictions also carry the
requested \code{m} or \code{epsilon} value (\code{va_m} or \code{va_epsilon}
when the outcome itself is named \code{m} or \code{epsilon}). Actual candidate
values, including locally generated grids, are stored separately in \code{grid}.
Effective per-fold moderation is
recorded in \code{fit_roles[[fold]]$moderation}; epsilon may imply different
integer m values on different fold sizes. These are moderation parameters,
not confidence or prediction interval levels.
}
\section{Tuning configuration}{
NULL keeps fixed settings. For binary, shared OvO/OvR or pairwise OvO classification and regression, a named list
with \code{strategy = "joint"} or \code{"two_stage"} and \code{model} and/or
\code{calibration} searches matching the present \code{tune()} placeholders.
Each search has \code{method = "grid"} or \code{"bayes"}. Grid uses a concrete grid data frame with columns
named by placeholder IDs or a positive integer grid size with dials
\code{param_info}. Bayes uses positive whole \code{initial} and \code{iter} budgets plus
\code{param_info}; \code{initial} must exceed the number of tuned parameters. Its
\code{control_bayes()} must save predictions, use \code{uncertain = Inf}, and supplies
the seed, \code{no_improve}, and optional time limit. Local model searches and calibration searches require
\code{resamples = function(data) ...} (or a bound rsample rset); pairwise
\code{joint} model grids use the global plan described below. Optional \code{primary_metric} defaults to
\code{yardstick::mn_log_loss} for classification or \code{yardstick::rmse} for regression; \code{report_metrics} are independent yardstick
metrics. Explicit \code{control} must be the matching \code{tune::control_grid()} or \code{tune::control_bayes()}
with \code{save_pred=TRUE}. Bayes calls the public \code{tune::tune_bayes()} one
iteration at a time, filtering incomplete candidates before the next call. Failed,
dropped and repeated proposals spend one step of the \code{iter} budget and count
toward \code{no_improve}; the time limit is checked between complete steps. Too few
complete initial candidates stop the fit with their reasons. Inner parallelism is disabled, and a custom
\code{extract} callback is unavailable because selection uses extraction for
diagnostics.
Model/recipe and calibration IDs must be unique. The ID \code{id} is unsupported
because native tune uses it for resample identities. Calibration search requires
exactly one of \code{selection_prop} (the proportion reserved for selection S)
or \code{selection_split} (rsplit or function, with analysis=final fit F and
assessment=selection S). S is never included in F. Automatic calibration
selection requires \code{splits=NULL} or a local plan provider function.
In pairwise OvO \code{two_stage}, model searches and resamples run on each local class pair. Calibration
selection scores the coupled complete Multiclass output on S; its final fit
uses only F. Pair plans must be local provider functions when calibration is tuned.
For pairwise \code{joint} Grid model tuning, candidate choices form a Cartesian grid
over class pairs (and separate full-model choices when requested). Bayes model
searches instead advance independently on each local pair and source/full domain;
they require local \code{model$resamples}. Global calibration candidates use
complete coupled predictions. Model-only joint Bayes evaluates the complete
locally selected procedure on \code{joint_resamples}; model-only joint Grid
uses that plan to select its Cartesian candidates. Omit \code{model$resamples}
only for the Cartesian Grid route. Global plans are never used as pair plans.
A joint grid above 4096
combinations is rejected before fitting; reduce the grids. With
\code{postprocessing}, pairwise \code{joint} Grid follows the Bayes route: each
local pair B searches itself and requires \code{model$resamples}; there is no
Cartesian grid.
Loss candidates are \code{"log"}/\code{"brier"}; precision candidates are
\code{"none"}, \code{"0"}, \code{"2"}, and \code{"-2"}. Regression searches exactly one of
\code{m = tune::tune("moderation")} or \code{epsilon = tune::tune("moderation")};
classification options are rejected. Every inner calibration fold must satisfy
\code{1 <= m <= floor((k-1)/2)}, with \code{m = floor(epsilon*(k+1)/2)}
when epsilon is supplied. Invalid local candidates are excluded without clipping.
Joint model selection scores complete IVAR/CVAR predictions; two-stage selection
scores the uncalibrated numeric predictions first. Fixed model/recipe or
moderation options stay fixed. Each Cross source and Full1 searches only its
own training domain; Full1 never averages fold choices. The return value above describes selections and diagnostics.
}
\section{Tuned postprocessing}{
With \code{postprocessing}, every engine/recipe search runs inside its own B:
each VA source Q, each inner ensemble pair and each B of a full F composition,
per class pair for pairwise OvO. The search resamples, and for \code{joint} the
VA and B/P splits of every candidate composition, are formed within that B.
Only then are the chosen settings fitted on B and the postprocessor on its P;
P, C, outer assessment A and a reserved selection share S never enter a search,
so mutating P leaves its search, recipe and engine unchanged. Inner pairs,
VA folds and F may choose different settings; none are averaged or copied.

\code{joint} scores complete calibrated compositions with fixed candidate
settings (no nested search). The first \code{two_stage} step scores the raw
engine/recipe without postprocessor or VA; calibration options are then chosen
on S as before, and S never returns to final B/P fits. A candidate is eligible
only with valid finite main scores on every expected resample, inner pair and row;
a failed final composition fit aborts without a fallback candidate.

A static \code{postprocessing$splits} plan cannot bind to candidate data, so
joint model selection and calibration selection require a plan function or
prop/v. \code{venn_abers_tuning(object)$models[[fold]]} (and, for pairwise,
\code{$models[[pair]][[fold]]}) is a list of local searches named by inner fold
(\code{Holdout1} or the ensemble fold IDs); each records \code{fold_id},
\code{inner_fold}, its B positions and the usual search diagnostics. Each
postprocessing record also stores its \code{selected} values. Budgets multiply
by the inner pair count J when it is known in advance; \code{cost$observed_base_fits}
counts every search and final B fit. See installed \code{examples/postprocessing-grid.R}.

Bayes searches follow the same B-local rule. Only candidates complete over
every resample, inner pair and row with a valid main score enter the next GP; a
failed inner model/postprocessor pair, a missing row or an invalid main score
excludes the whole candidate. Excluded candidates, including proposals that
\code{tune::tune_bayes()} drops after failing on every resample, keep their
global candidate IDs and reasons; \code{steps}, \code{notes},
\code{stop_reason}, \code{attempts}, \code{bayes_seed} and \code{search_seed}
record each local search. There is no switch to Grid or another method. See
installed \code{examples/postprocessing-bayes.R}.
}
\section{Related training rules}{
See \link{workflow_training} for local split validation and full-model
transfer, and \link{workflow_postprocessing} for B/P data roles, methods and
compositions. These choices remain fixed while the configured parameters are
searched.
}
\seealso{\link{vennabers_examples} for annotated installed scripts.

\code{\link{calibrated_workflow}}. Installed Grid examples
\code{binary-grid.R}, \code{shared-multiclass-grid.R}, \code{pairwise-grid.R},
\code{regression-grid.R} and \code{postprocessing-grid.R}, and Bayes examples
\code{shared-multiclass-bayes.R}, \code{pairwise-bayes.R}, \code{regression-bayes.R}
and \code{postprocessing-bayes.R}, are under
\code{system.file("examples", package = "vennabers")}.}
\examples{
packages <- c("parsnip", "workflows", "recipes", "rsample", "hardhat",
              "generics", "rlang", "tibble", "tidyselect", "withr",
              "tune", "dials", "yardstick")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  development <- data.frame(
    x = sin(1:80),
    truth = factor(rep(c("yes", "no", "no", "yes"), 20),
                   levels = c("yes", "no"))
  )
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))

  # Choose the regularization loss on a reserved selection share S; the
  # final calibrator is fitted on the remaining rows F.
  calibration <- list(method = "grid",
                      grid = data.frame(va_loss = c("log", "brier")),
                      resamples = function(d) rsample::vfold_cv(d, v = 2))
  tuning <- list(strategy = "two_stage", selection_prop = 0.3,
                 calibration = calibration)
  specification <- calibrated_workflow(base, outcome = "truth", seed = 42,
                                       loss = tune::tune("va_loss"),
                                       tuning = tuning)
  fitted <- generics::fit(specification, development)

  choices <- venn_abers_tuning(fitted)
  print(choices$calibration$candidates)
  print(choices$calibration$selected)
  # Selection rows never enter the final calibration fit.
  length(intersect(choices$calibration$selection_rows,
                   choices$calibration$final_rows))
}
}

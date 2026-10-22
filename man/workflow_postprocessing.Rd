\name{workflow_postprocessing}
\alias{workflow_postprocessing}
\title{Postprocessing Before Workflow Calibration}
\description{
Configuration and reference rules for the \code{postprocessing} argument of
\code{\link{calibrated_workflow}}. Start with the worked Sigmoid example in
\code{vignette("workflows", package = "vennabers")} before using this reference.
These options learn a postprocessor before the Venn-Abers calibrator.
}
\section{Data roles}{
Each VA source training portion Q is split into B for recipe and engine
learning and P for postprocessor learning. Held-out C fits the VA calibrator.
Outer assessment A is excluded from all learning. When calibration is tuned,
selection rows S are excluded from the final fit domain F. With
\code{prediction_model = "full"}, F gets its own complete B/P composition.
An inner pair is a model/postprocessor pair; an OvO class pair is a separate
two-class training domain.
}
\section{Configuration}{
NULL preserves the existing unprocessed training path and
does not remove a supplied postprocessor. Otherwise a named list with explicit
\code{method = "tailor_logistic"}, \code{"tailor_multinomial"}, \code{"tailor_linear"}, \code{"sigmoid"}, or \code{"isotonic_numeric"},
\code{strategy = "holdout"}, and either
\code{prop} (the P fraction, strictly between zero and one) or \code{splits}
(a local rsplit or serializable function of local Q data). Analysis is B;
assessment is P. A supplied plan cannot be combined with prop or v. Unknown,
duplicate, irrelevant and conflicting fields are errors. Currently binary, shared OvO/OvR and pairwise OvO
Inductive/Cross VA and regression IVAR/CVAR with matched or full models are implemented,
with fixed engine/recipe settings or local selection (see \code{\link{venn_abers_tuning}});
Binary and pairwise OvO \code{tailor_logistic}/\code{sigmoid}, shared OvO/OvR
\code{tailor_logistic} (two classes)/\code{tailor_multinomial} (three or more), and
regression \code{tailor_linear}/\code{isotonic_numeric} support
\code{strategy = "ensemble"}, with an explicit integer \code{v >= 2} or
\code{splits} (a local unrepeated rset or serializable function of Q).
Ensembles reject prop; explicit splits reject v.
Method, Sigmoid input, strategy, prop/v and prepared plans stay fixed per run.
}
\section{Supported compositions and costs}{
The postprocessing scope is a finite matrix of 80 fixed compositions:
method x holdout/inner ensemble x Inductive/Cross x matched/full over binary
(Sigmoid, Tailor logistic), pairwise OvO (Sigmoid, Tailor logistic), shared
OvO/OvR with three or more classes (Tailor multinomial), shared two-class
OvO/OvR (Tailor logistic) and regression (reference isotonic, Tailor linear).
Every composition also accepts joint or two-stage Grid/Bayes selection of
engine/recipe settings inside each B. B trains recipe and engine, P the
postprocessor, C the unchanged VA calibrator; F is the permitted final fit
data of a full composition. Without tuning, J inner pairs cost J (Inductive
matched), 2J (Inductive full), KJ (Cross matched with K folds) or (K+1)J (Cross full)
engine and postprocessor fits per domain, retaining J, J, KJ or J pairs;
class pairs add up. Out of scope: out-of-fold postprocessing followed by a
complete engine refit, foreign engines beyond \code{\link{model_adapter}}
callbacks, Tailor step sequences, further postprocessing methods, and other new
public constructors. Frozen workflows
(\code{\link{calibrate_workflow}}) keep their existing contract.
}
\section{Binary Sigmoid holdouts}{
For binary matched or full Inductive/Cross workflows with fixed parameters, set
\code{postprocessing = list(method = "sigmoid", strategy = "holdout", prop = .3,
  sigmoid_input = "probability")} or explicitly select \code{sigmoid_input = "glm_link"}.
The base workflow must have no postprocessor. Recipe and engine learn only on B;
Sigmoid learns only on P; the unchanged VA core learns on held-out C. Both B and
P must contain both factor levels, and P needs at least two rows. Constant P
inputs are valid; coefficients are not identifiable in that case.

\code{probability} uses the final named probability of the second stored factor
level from the complete public base workflow. \code{glm_link} requires a
\code{logistic_reg(engine = "glm")} binomial model with logit link and calls
\code{predict(workflow, new_data, type = "raw", opts = list(type = "link"))},
including recipe processing. It never reconstructs links from probabilities.
Both VA event choices are supported by named probability selection after Sigmoid.

The mapping minimizes unregularized logistic loss with targets
\eqn{(N_1+1)/(N_1+2)} for the second class and \eqn{1/(N_0+2)} for the first.
It follows scikit-learn's Platt calibration, including its start values,
L-BFGS-B tolerances (relative loss reduction of 64 machine epsilons, gradient
1e-6) and SciPy's ten correction pairs. These tolerances stop short of the exact
optimum, so the optimizer path determines the result. Stored \code{coefficients} are
\code{a, b} on the stored input \code{scale}: the second-class output is
\code{plogis(-(a * (x / scale) + b))}. The scale is the maximum absolute P
input when at least 30, otherwise one. Constant inputs use the exact zero-slope
solution (the mean target). This is an R choice, not Python parity:
scikit-learn's iterative estimate stops within its gradient tolerance and, for
very small P, can differ from it by more than 1e-7. A line-search stop is kept
only when its mean gradient per P row is at most 1e-6 in each coefficient, as at
a machine-precision optimum.
Nonfinite or incomplete inputs, every other unsuccessful optimization, and incomplete
probability outputs fail with method, role and local fold context.

The fitted state retains the complete base workflow, mapping and package
versions. Details expose B/P/C memberships, \code{input}, \code{sigmoid_class},
scale, coefficients, convergence, warnings and fit counts. Empty batches retain their schema without engine calls.
No Tailor, probably or Python installation is required for Sigmoid. Python is
used only by the separate checked-in reference generator; regular tests read
its fixed CSV outputs. See \code{system.file("examples", "sigmoid-holdout.R",
package = "vennabers")} for both inputs and native outer evaluation.
}
\section{Tailor classification holdouts}{
Each VA source training portion Q is partitioned into disjoint B and P. The
recipe and engine learn only on B through the public workflow fit. A
package-owned composition validates its P predictions and fits an unchanged
canonical Tailor through its public fit, using
\code{tailor::adjust_probability_calibration(method = "logistic", smooth = FALSE)}
for two classes, or \code{method = "multinomial", smooth = FALSE} for three or
more classes. The requested method must match the class count; no automatic
exchange is performed. Shared two-class OvO/OvR retains multiclass class and
tie rules, not the binary event contract.
One exactly matching untrained Tailor adjustment may already be configured;
it is adopted once. Foreign, trained or multiple adjustments are rejected.
The stored composition predicts held-out VA data C and new predictors through
the base workflow, input validation and the public Tailor prediction path.
The VA core, event alignment, loss and rounding rules are unchanged.

Both B and P must contain all K model classes. P needs at least K+1 rows and
a full-rank design: the logistic intercept and one probability for two classes,
or the multinomial intercept and K-1 probabilities for K classes. Rank and
convergence are checked on the actual learned model after the necessary B
predictions. Failed fits, absent learned state and method switches to none
are errors. Tailor multinomial uses the \code{nnet::multinom()} defaults of
\pkg{probably} (at most 100 iterations), and the fixed Tailor arguments offer
no iteration setting: weakly informative class probabilities can stop there
unconverged even when the B engine converges. Setting \code{maxit} on the B
engine does not change the P postprocessor's limit. The composition fails
with fold, inner-pair and B/P context, without retry, method switch or partial
ensemble. This is an accepted limit of the 1.x contract, not a pending
option. This R composition limit does not change
numerical VA parity on identical valid inputs. Warnings from successful fits are retained with source/fold context.
The existing VA rules for missing C classes remain unchanged (including valid
one-class Cross subproblems); B/P requirements do not apply to C. Use \code{venn_abers_schema()} for nonsyntactic outcome
or class names that the underlying Tailor/probably formula interface cannot use.
}
\section{Tailor linear holdouts}{
Linear regression postprocessing requires \code{method = "tailor_linear"} explicitly and uses
\code{tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE)}
through the same package-owned composition. Recipe/engine learning uses B;
the intercept and slope are learned only on P. Numeric 0/1 targets remain
regression. Automatic B/P splits are unstratified. P needs at least three rows,
finite predictions and truth, full actual design rank two, finite coefficients
and convergence; constant predictions, absent fits and switches to none fail.
The complete composition's numeric output enters the unchanged IVAR/CVAR core
with the requested valid \code{m} or \code{epsilon}. CVAR averages fold points
and both bounds separately. No confidence or prediction intervals are added.
See installed \code{examples/tailor-regression-holdout.R} for both paths,
encoded names, complete native evaluation and data-role checks.
}
\section{Numeric isotonic holdouts}{
\code{method = "isotonic_numeric"} accepts a regression base workflow without
a postprocessor. The recipe and engine learn on B, and unweighted increasing
least-squares isotonic regression learns on finite P predictions and truths.
At least two P rows are needed; tied predictions retain their frequencies.
Fitted support values are linearly interpolated, extended constantly beyond
the support, and are not restricted to [0,1]. A single support predicts the
mean P truth everywhere. No full-rank or Ridge condition applies. The complete
saved B workflow and isotonic fit enter the same IVAR/CVAR cores. See installed
\code{examples/isotonic-regression-holdout.R} for both native evaluation paths.
}
\section{Local plans and seeds}{
All source B/P plans are validated before any source is trained. Static plans must
match local selected columns in order, values, types, levels and row order.
Automatic classification splits use the local outcome with \code{pool = 0}, without retries
or size adjustments; ordered valid explicit holdouts are allowed. The enclosing
seed scope covers VA splitting, B/P splitting and fitting in stored fold order;
NULL consumes the caller's stream and an explicit seed restores it.
}
\section{Inner ensembles}{
For inner ensembles, P folds partition Q exactly once and each B is its
own P complement. B/P disjointness is local to a pair: a P row can belong to
another pair's B. Gaps, overlap, bootstrap duplicates, repeated folds and
incorrect local data bindings fail before learning. Each pair uses the same
complete holdout fitting rules. Aligned final probability vectors are averaged
with equal pair weights before VA, without inner VA calibrators or size weights.
Shared OvO/OvR averages complete class vectors in stored class order before
binary VA aggregation and subsequent OvO coupling or OvR normalization.
Regression averages finished numeric predictions from the retained complete
model/postprocessor pairs before IVAR/CVAR. It does not average coefficients,
isotonic supports or separately VA-calibrated inner results. Each pair retains
its method's P-size and rank rules, including valid constant isotonic inputs.
CVAR's outer fold bounds remain descriptive VA bounds, not confidence or
prediction intervals. See installed \code{examples/regression-ensemble.R}
for both numeric methods, native evaluation, schema and held-out data checks.
All plans are resolved first, followed by fits in stored VA/inner-fold order in
one RNG stream; the seed is never restarted for each inner B/P pair. Local plan functions
must be pure functions of their supplied data, not call counters or global row
positions. Different Q domains may have different inner fold counts.
}
\section{Pairwise OvO compositions}{
For pairwise OvO, restrict the data and outcome levels to each class pair before
resolving any VA or B/P plan. Each pair supports \code{tailor_logistic} or
\code{sigmoid} with either explicit input type; no foreign class reaches its
recipe, model, postprocessor or plan function. Pair-local complete predictions
are averaged before the binary VA core, and only the resulting pair points enter
global OvO coupling. Different pairs and VA folds can have different local sizes
and inner fold counts. Explicit B/P plans must bind to their own Q, using a pure
local function or \code{local_calibration_splits()} for multiple prepared domains.
Global position plans are rejected. Existing isolated VA pair-plan seeds remain
unchanged; all B/P plans are then resolved in class-pair/VA-fold order, followed
by fits in class-pair/VA-fold/inner-fold order in the enclosing RNG stream.

Pairwise diagnostics are stored at
\code{extract_fit_engine(fitted)$pair_models[[pair_id]]$postprocessing}.
\code{venn_abers_details()$postprocessing} exposes these records under unique
\code{pair_N/VA_fold} keys. Each record includes \code{pair_id} and
\code{pair_rows} (positions in the enclosing fit data); Q/B/P/C positions are
relative to that class pair. Inner-fold names identify model/postprocessor pairs,
not class pairs. The decoder preserves the two-class outcome domain and translates
its names. Costs sum the actual inner counts across every class pair and VA fold;
all required paths remain stored. Missing members and failed fits/predictions
abort the complete result with pair/fold context. See installed
\code{examples/pairwise-postprocessing.R} for native resampling and last_fit.
}
\section{Shared diagnostics and serialization}{
For shared models, \code{extract_fit_engine(fitted)$postprocessing} and the \code{postprocessing}
element of \code{venn_abers_details()} store actual Q/B/P/C row memberships,
Q-local analysis/assessment positions, sizes, class/event and schema bindings,
RNG settings, method-specific fit state, warnings and observed fit counts. Prediction
details add the inner fold, method, input and B/P sizes. Package versions are
stored in the fitted engine's \code{versions}. Each nonempty prediction calls
each retained base workflow and Tailor once; empty batches call none. Predict needs
no truth or plan resolution and never learns. RDS retains the full trained state.
Ensemble records contain a named \code{inner_folds} list of complete pair
diagnostics and summed observed fit counts. Their prediction details use
list-columns for \code{inner_fold}, \code{B_size} and \code{P_size}, retaining
one VA result row and the ordered identities/sizes of every contributing pair.
An Inductive ensemble fits J engines and J postprocessors; a K-fold Cross
ensemble with J pairs per fold fits K*J of each. Missing or failed pairs abort
the whole composition with fold context. See installed
\code{examples/binary-ensemble.R} and \code{examples/multiclass-ensemble.R}
for complete native evaluation.
Tailor and probably are optional and required only for this path.
See installed \code{examples/tailor-holdout.R} for native fit, complete repeated
outer resampling, last_fit, diagnostics and checks of held-out data roles.
\code{examples/tailor-multiclass-holdout.R} adds shared two-, three- and four-class
OvO/OvR, encoded names and the yardstick Brier factor of one half. Single-row
multinomial predictions use a duplicated prediction row at the composition
boundary, then discard that extra result. The base workflow is predicted once;
the canonical Tailor and its learned nnet model remain unchanged. Learned-state
compatibility checks are confined to the package's Tailor boundary because
Tailor has no public learned-model extractor.
}
\section{Full postprocessing compositions}{
All fixed-parameter holdout and inner-ensemble methods also accept
\code{prediction_model = "full"}. First every original Q composition predicts
its held-out C and builds its VA calibrator. Only then does each permitted F
domain resolve its own B/P plan and fit an additional complete composition.
This ordering preserves the source RNG stream, C values and calibrators between
matched and full fits with the same seed. F excludes outer assessment A and
any reserved selection S; for pairwise OvO it contains only the local class pair.
C can influence this extra composition, never its own original Q composition.
Full refers to the whole composition: no inner engine is trained on all F rows.

After successful fitting and prediction validation, only the F composition
remains in \code{models$Full1}. Source workflows are discarded, while their
held C values, VA calibrators, plans and postprocessing diagnostics remain.
The additional \code{postprocessing$Full1} record uses \code{Q} for the F row
membership, has empty \code{C}, and records its own B/P memberships and fit counts.
Pairwise details expose it as \code{pair_N/Full1}. Prediction-row inner-fold names
and sizes describe the retained F composition; source diagnostic records remain
available without claiming new predictions from discarded workflows.
If a caller's VA fold already uses \code{Full1}, the final ID gets the next free
numeric suffix (\code{Full1.1}, etc.); \code{sources} records the actual ID.

For constant inner count J (J=1 for holdout), full Inductive performs 2J engine
and postprocessor fits; full Cross with K VA folds performs (K+1)J. Only J pairs
are retained per domain. Different local counts are summed, including F's own J.
Three class pairs with K=5 and J=3 perform 54 engine and 54 postprocessor fits,
retaining nine pairs. Each nonempty batch predicts those retained pairs once,
then reuses the composed F prediction for every original VA calibrator. Empty
batches call no engine. A failed final plan, fit or prediction aborts the result;
there is no fallback, retention switch or prediction-time model-mode switch.

RDS retains complete F workflows and all diagnostics. Prediction needs neither
truth, plan functions nor additional learning. See installed
\code{examples/full-postprocessing.R} for all 80 fixed matched/full combinations,
including shared two-class OvO/OvR, native outer evaluation and last_fit.
The full-model transfer limitations in \link{workflow_training} apply to these
complete compositions as well.
}
\seealso{\link{vennabers_examples} for annotated installed scripts.


\code{\link{calibrated_workflow}} for the constructor,
\code{\link{venn_abers_details}} for returned diagnostics,
\code{\link{venn_abers_tuning}} for searches inside each B, and
\link{workflow_training} for VA splits and full-model transfer.
}

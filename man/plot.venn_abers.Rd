\name{plot.venn_abers}
\alias{plot.venn_abers}
\alias{plot.cross_venn_abers}
\alias{plot.venn_abers_regression}
\alias{plot.cross_venn_abers_regression}
\alias{plot.vennabers_fit}
\title{Plot the Calibration Map of a Fitted Venn-Abers Calibrator}
\usage{
\method{plot}{venn_abers}(x, which = "map", loss = "log", ...)

\method{plot}{cross_venn_abers}(x, which = "map", loss = "log", ...)

\method{plot}{venn_abers_regression}(x, which = "map", ...)

\method{plot}{cross_venn_abers_regression}(x, which = "map", ...)

\method{plot}{vennabers_fit}(x, which = "map", loss = x$loss, ...)
}
\arguments{
  \item{x}{A fitted binary or multiclass \code{venn_abers} or
  \code{cross_venn_abers} object, or a fitted \code{venn_abers_regression} or
  \code{cross_venn_abers_regression} object. For \code{plot.vennabers_fit},
  a frozen fit from \code{\link{calibrate_workflow}} or the engine fit of a
  trained \code{\link{calibrated_workflow}}, extracted with
  \code{workflows::extract_fit_engine()}; see the Workflows section.}
  \item{which}{A character vector of distinct diagrams from \code{"map"}
  (default), the calibration map, and \code{"width"}, the Venn-Abers width
  \code{upper - lower}. Both draw the map before the width, whatever their
  order.}
  \item{loss}{Exactly \code{"log"} (default) or \code{"brier"}: the
  regularization of the drawn and returned probability, as in
  \code{\link{predict.venn_abers}}. Classification only: regression methods
  take no \code{loss}, and passing one is an error. For a workflow fit,
  the default is the \code{loss} stored at fitting.}
  \item{...}{Base R graphical options passed to \code{plot()}, taking
  precedence over the defaults in every panel. For regression, \code{xlim}
  also sets the range of the returned table.}
}
\value{
Invisibly, a data frame with the same columns whatever \code{which} holds:
\code{score}, \code{lower}, \code{upper}, \code{probability}, and
\code{width}. Rows hold each stored probability knot plus the range ends 0
and 1 in increasing order. Values are those of \code{predict()} at
\code{score} with the chosen \code{loss}, and \code{width} is
\code{upper - lower}.

For \code{cross_venn_abers}, a leading character column \code{fold} holds the
stored fold IDs. Rows group by fold in stored order, each with its own knots
and the values of its own fold calibrator, not the aggregated Cross
probability.

For multiclass fits, the tables of all subproblems are stacked in stored
order, with identity columns named as in \code{\link{venn_abers_details}}
after \code{fold} and before the values: \code{class} for one-vs-rest,
\code{class_1} and \code{class_2} for one-vs-one. Within a Cross
subproblem, folds follow the stored fold IDs of the fit.

For regression, the columns are \code{score}, \code{lower}, \code{upper},
\code{prediction}, and \code{width}, after \code{fold} for
\code{cross_venn_abers_regression}, with the values of \code{predict()} at
\code{score}. The range ends replace 0 and 1; see the Regression section.

For a workflow fit, the table is that of its stored Venn-Abers calibrator.
A pairwise fit stacks the tables of its pair components in stored order,
led by the identity columns \code{pair_id}, then \code{fold} for Cross
fits, then \code{class_1} and \code{class_2}, as in
\code{\link{venn_abers_details}}.

Between neighboring rows, the map holds the \code{lower} value of the left
row and the \code{upper} value of the right row, so the width between knots
exceeds the width at them. The drawn curves follow these values exactly.
}
\description{
Draws a fit diagnosis from the stored state of a fitted calibrator, without
new truth values. For classification, the calibration map shows the lower
and upper Venn-Abers values as steps and the regularized probability as a
dashed line over the uncalibrated probability on \code{[0, 1]}. The width diagram shows
\code{upper - lower} over the same scores; it is neither a confidence level
nor a nominal coverage.

Cross fits overlay the map of every fold calibrator in one panel,
distinguished by color and legend. Comparing fold maps at common scores does
not describe the spread of an individual Cross prediction.

Choosing both diagrams draws them in their own grid of panels and restores
the previous graphical parameters on exit, also after errors. A single
diagram uses the current layout. Fits without a stored calibration map
(\code{knots}, \code{lower}, \code{upper}) are errors; refit them with
vennabers 1.0.0.
}
\section{Multiclass}{
A multiclass fit draws one panel per subproblem and diagram, titled
\code{"Calibration map: <class>"} for one-vs-rest and
\code{"Calibration map: <class_1> vs <class_2>"} for one-vs-one (likewise
\code{"VA width"}). Each panel shows the map of the subproblem's own binary
calibrator over its score: the probability of the class for one-vs-rest, and
the probability of \code{class_2} renormalized within the pair for
one-vs-one. The final coupling into class probabilities is not drawn. Cross
fits overlay the fold maps of each subproblem in its panel, each fold in the
same color across panels.

The panels fill their own grid row by row, nearly square as in
\code{\link{calibration_plot}}. With both diagrams, each row of maps is
followed by the widths of the same subproblems. A single panel, as for one
one-vs-one pair, uses the current layout. Options in \code{...} apply to
every panel.
}
\section{Regression}{
A regression fit draws the lower and upper Venn-Abers values as steps and the
point prediction as a dashed line over the uncalibrated prediction, with
gray dashed reference lines at the moderated truth values \code{a} and
\code{b} that bound every Venn-Abers value. The width diagram shows
\code{upper - lower} in units of the outcome.

The drawn range spans the stored knots, widened on each side by 5 \% of
their span; a single knot widens by 5 \% of its absolute value, or 0.05 at
zero. A Cross fit spans the knots of all folds, overlays its fold maps with
a reference line at each fold's \code{a} and \code{b}, and returns every
fold map over this shared range. A finite \code{xlim} in \code{...}
replaces the range: the table then holds its ends and the knots within.
}
\section{Workflows}{
\code{plot()} on a frozen fit from \code{\link{calibrate_workflow}}, and on
\code{workflows::extract_fit_engine(fit)} for a trained
\code{\link{calibrated_workflow}}, draws the fit diagnosis of the stored
Venn-Abers calibrator: binary, multiclass, or regression, inductive or
Cross, as described above for that calibrator, with the same table. There is
no method for the workflow itself. Postprocessing models (sigmoid, tailor,
isotonic, or their ensembles) are source models before the calibrator: the
map shows the calibrator over their output, not the postprocessing step.

A pairwise fit (\code{model_training = "pairwise"}) draws one panel per pair
component and diagram, titled \code{"Calibration map: <class_1> vs <class_2>"},
in a grid as for multiclass fits. Each panel shows the binary calibrator of
its pair over the probability of \code{class_2} from that pair's own model.

A fit without a stored Venn-Abers calibrator is an error.
}
\seealso{
Assessment graphics need held-out truth values instead of the stored state:
\code{\link{corp_plot}} and \code{\link{calibration_plot}} for
classification probabilities, \code{\link{regression_plot}} for regression
predictions. \code{vignette("diagnostics", package = "vennabers")} compares
both kinds of graphics.
}
\examples{
truth <- c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE)
uncalibrated_probability <- c(0.1, 0.2, 0.3, 0.4, 0.4, 0.6, 0.7, 0.9)
calibrator <- venn_abers(truth, uncalibrated_probability)

map <- plot(calibrator, which = c("map", "width"))
head(map)

cross <- cross_venn_abers(truth, uncalibrated_probability, fold = rep(1:2, 4))
plot(cross, loss = "brier")

classes <- c("a", "b", "c")
class_truth <- factor(rep(classes, 4), levels = classes)
class_probability <- prop.table(matrix(c(
  6, 3, 1, 2, 5, 3, 1, 2, 7, 5, 1, 4, 3, 4, 3, 2, 2, 6,
  4, 4, 2, 1, 6, 3, 3, 3, 4, 5, 3, 2, 2, 7, 1, 2, 2, 6
), 12, 3, byrow = TRUE, dimnames = list(NULL, classes)), 1)
ovr <- venn_abers(class_truth, class_probability, method = "ovr")
head(plot(ovr, which = c("map", "width")))

outcome <- c(-4, 1, 0, 3, 2, 6, 5, 9, 7, 12)
uncalibrated_prediction <- c(-2, 0, 0, 1, 2, 3, 3, 5, 6, 8)
regression <- venn_abers_regression(outcome, uncalibrated_prediction)
head(plot(regression, which = c("map", "width")))
plot(regression, xlim = c(0, 5))
if (all(vapply(c("parsnip", "workflows", "hardhat", "rsample", "recipes",
                "rlang", "tibble", "withr", "generics", "tidyselect"),
               requireNamespace, logical(1), quietly = TRUE))) {
  d <- data.frame(x = sin(1:80), truth = factor(rep(c("yes", "no", "no", "yes"), 20)))
  base <- workflows::workflow() |>
    workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  frozen <- calibrate_workflow(generics::fit(base, d[1:40, ]), d[41:80, ], "truth")
  print(head(plot(frozen)))

  trained <- generics::fit(calibrated_workflow(base, "truth", seed = 42), d)
  plot(workflows::extract_fit_engine(trained), which = "width")
}
}

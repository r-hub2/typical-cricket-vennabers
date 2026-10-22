\name{venn_abers_one_hot}
\alias{venn_abers_one_hot}
\title{One-Hot Matrices of Class Decisions}
\description{Represent existing class decisions as a one-hot matrix with one
column per class. The helper converts the decision a predictor has already made;
it never chooses a class from probabilities.}
\usage{
venn_abers_one_hot(x)
}
\arguments{
\item{x}{A class decision: the factor returned by multiclass
\code{predict(..., type = "class")} or by \code{venn_abers_decode()$class}; the
logical vector returned by binary \code{predict(..., type = "class")}; or a data
frame with a factor \code{.pred_class} column, as returned by workflow, Grid,
parsnip and frozen \code{predict(..., type = "class")}. Other columns of a data
frame are ignored.}
}
\value{
A base \code{matrix} of type double with values 0 and 1, one row per element of
\code{x} in the same order and no row names. Column names are the factor levels in
their stored order, or \code{"FALSE"} and \code{"TRUE"} for a logical vector, so
the binary event is the \code{TRUE} column. Every class keeps its column, also when
no row predicts it. Zero rows give a 0 x K matrix and one row a 1 x K matrix; the
result is never reduced to a vector. A missing decision gives a row of \code{NA}.
}
\details{
The ones follow the existing tie rules: a binary probability of exactly 0.5 is
the event, and a multiclass tie keeps the first stored level, also with exactly
two classes. The Python reference facade breaks a binary tie towards the first
sorted class instead; \code{system.file("PARITY.md", package = "vennabers")}
records this difference.

Ordered factors only set the column order. Schema-encoded \code{.pred_class}
columns keep their encoded classes as column names; pass
\code{venn_abers_decode(prediction, schema)$class} to obtain original class names.
Probability matrices or tables, numeric or character vectors, data frames without
a factor \code{.pred_class}, and factors without levels are errors.
}
\examples{
source(system.file("examples", "one-hot-classes.R", package = "vennabers"))
}
\seealso{\code{\link{predict.venn_abers}}, \code{\link{calibrated_workflow}},
\code{\link{venn_abers_decode}}}

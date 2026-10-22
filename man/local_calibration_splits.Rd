\name{local_calibration_splits}
\alias{local_calibration_splits}
\title{Resolve Prepared Local Calibration Splits}
\description{Validate all prepared assignments and return a serializable function
mapping selected local raw data to one mixed or ordered Inductive holdout or Cross partition.}
\usage{local_calibration_splits(splits, columns)}
\arguments{
\item{splits}{A nonempty list of data-bound holdouts from
\code{rsample::initial_split()} or \code{rsample::make_splits()}, or simple Cross
\code{rset} partitions from \code{rsample::vfold_cv()} or
\code{rsample::manual_rset()}. Each split's
stored data is its local fitting domain, typically one outer analysis set.}
\item{columns}{Unique selected raw predictor names in carrier order followed by
the outcome name. These columns, including outcome, define data identity. They
must exactly match the canonical columns received by the resolver; an omitted
outcome or selected predictor is an error.}
}
\value{A function of local raw \code{data}, suitable for the \code{splits} argument
of \code{calibrated_workflow()}. It returns a data-bound \code{rsplit} or Cross \code{rset}.}
\details{Every prepared split is validated before the function is returned.
Identity includes column values, types, factor levels and row order. Container
class and row names are ignored. Extra unselected columns are ignored; duplicate
observations remain distinct positions. Inductive uses the fixed single fold ID
\code{Inductive1}; Cross retains the ordered unique fold IDs. Repeated, overlapping,
incomplete and noncomplementary inner Cross partitions are rejected before fitting.
Identical data require the same method, ordered fold IDs and ordered
analysis/assessment indices. Equivalent rset containers normalize to the same plan.
Duplicate identical plans are allowed; conflicting entries are errors identifying
their positions. A lookup without a matching snapshot is an error. No call counter,
first conflicting match, or outer assessment access is used. Construction compares
pairs of prepared entries; lookup scans the saved snapshots. Large collections may
need a more efficient user resolver.

The resolver supports \code{saveRDS()} and \code{readRDS()} and arbitrary lookup
order. An arbitrary user callback is trusted R code and must obey the same data
and plan contract; this helper does not prove the absence of hidden state in other
callbacks. Keep model selection and all other learning outside assessment data.
An ordered holdout can be constructed with \code{rsample::make_splits()}: reserve
\code{ceiling(n * q)} trailing calibration rows and the preceding
\code{n - ceiling(n * q)} training rows, requiring both parts to be nonempty.
At \code{n = 11, q = .25}, these are training positions 1--8 and calibration
positions 9--11, with no sorting, mixing or stratification. The order is local
to each outer analysis set. Cross accepts mixed or non-shuffled ordinary or
classification-stratified complete partitions. Exact training complements are
required even for non-shuffled folds; these are not forward validation.
See the executable installed \code{examples/calibrated-resampling.R} and
\code{examples/reference-splits.R}.}
\examples{
if (all(vapply(c("parsnip", "workflows", "hardhat", "rsample", "recipes",
                "rlang", "tibble", "withr", "generics", "tidyselect"),
               requireNamespace, logical(1), quietly = TRUE))) {
  data <- data.frame(x = 1:20, truth = factor(rep(c("yes", "no"), 10)))
  split <- withr::with_seed(42, rsample::initial_split(data))
  resolve <- local_calibration_splits(list(split, split), c("x", "truth"))
  stopifnot(identical(as.integer(resolve(data)), as.integer(split)))
}
}
\seealso{\code{\link{calibrated_workflow}}}

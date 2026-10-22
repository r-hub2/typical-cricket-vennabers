\name{print.cross_venn_abers}
\alias{print.cross_venn_abers}
\title{Print a Cross Venn-Abers Calibrator}
\usage{
\method{print}{cross_venn_abers}(x, ...)
}
\arguments{
  \item{x}{A fitted \code{\link{cross_venn_abers}} object.}
  \item{...}{Unused.}
}
\value{
\code{x}, unchanged and invisibly.
}
\description{
Prints three lines: \code{<cross_venn_abers>}, \code{Folds: K}, and
\code{Calibration observations: N}, with the actual counts substituted.
For multiclass, also identifies the method and class names in stored order.
It exposes neither calibration arrays nor predictive models.
}
\examples{
calibrator <- cross_venn_abers(c(0, 1, 0, 1, 0, 1), c(.1, .8, .2, .9, .3, .7),
                               c("A", "A", "B", "B", "C", "C"))
print(calibrator)
}

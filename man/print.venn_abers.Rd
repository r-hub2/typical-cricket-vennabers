\name{print.venn_abers}
\alias{print.venn_abers}
\title{Print a Venn-Abers Calibrator}
\usage{
\method{print}{venn_abers}(x, ...)
}
\arguments{
  \item{x}{A fitted \code{venn_abers} object.}
  \item{...}{Unused.}
}
\value{
\code{x}, invisibly.
}
\description{
Prints the calibrator class and counts of calibration observations, event
observations, and unique probability knots.
For multiclass objects, prints the method, calibration count, class names in
stored order, and number of class pairs or classes against rest instead.
}
\examples{
venn_abers(c(0, 1, 0), c(0.1, 0.2, 0.3))
}

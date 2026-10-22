\name{vennabers-engine}
\alias{fit_vennabers}
\alias{fit_vennabers_proxy}
\alias{predict_vennabers_prob}
\alias{predict_vennabers_class}
\alias{predict_vennabers_numeric}
\title{Parsnip Engine Entry Points}
\description{Exported entry points for the optional parsnip registration. Users
should construct and predict with \code{\link{calibrated_workflow}} or \code{\link{model_adapter}} instead.}
\usage{
fit_vennabers(x, y, base_workflow, outcome, parameters, ...)
fit_vennabers_proxy(x, y, outcome, parameters, ...)
predict_vennabers_prob(object, new_data)
predict_vennabers_class(object, new_data)
predict_vennabers_numeric(object, new_data)
}
\arguments{
\item{x, y}{Raw predictor data frame and factor classification or finite numeric regression truth from parsnip.}
\item{base_workflow}{An unfitted configured base workflow.}
\item{outcome}{Original outcome name.}
\item{parameters}{Constructor-validated requested options and explicit flags, internal proxy candidate plans, or the validated callbacks of \code{\link{model_adapter}} (then \code{base_workflow} and \code{outcome} are unused).}
\item{...}{Must be empty for \code{fit_vennabers}; an internal grid candidate ID or Bayes candidate values for \code{fit_vennabers_proxy}.}
\item{object}{The fitted \code{vennabers_fit}, grid proxy or \code{vennabers_adapter_fit} engine state.}
\item{new_data}{Forged raw predictors.}
}
\value{The public fit entry point returns a \code{vennabers_fit}, or a \code{vennabers_adapter_fit} for adapter callbacks; the proxy returns a \code{vennabers_grid_fit}. Prediction entry points
return a named probability tibble, factor of original classes, or double regression points.}
\keyword{internal}

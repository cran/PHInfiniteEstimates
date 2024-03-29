% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/newllk.R
\name{newllk}
\alias{newllk}
\title{Proportional hazards partial likelihood, using Breslow method for ties, excluding some observations.}
\usage{
newllk(
  beta,
  fit,
  exclude = NULL,
  minus = FALSE,
  keeponly = NULL,
  justd0 = FALSE,
  cc1 = 0
)
}
\arguments{
\item{beta}{parameter vector.}

\item{fit}{Output from a Cox PH regression, with x=TRUE and y=TRUE}

\item{exclude}{data set with stratum and patient number to exclude.}

\item{minus}{logical flag to change sign of partial likelyhood}

\item{keeponly}{variables to retain.  Keep all if this is null or NA.}

\item{justd0}{logical variable, indicating whether to calculate only the function value and skip derivatives.}

\item{cc1}{Continuity correction for first component of the score.}
}
\value{
a list with components
\itemize{
  \item d0 partial likelihood
  \item d1 first derivative vector
  \item d2 second derivative matrix
}
}
\description{
This function implements the approximate conditional inferential approach of \insertCite{kz19;textual}{PHInfiniteEstimates} to proportional hazards regression.
}
\references{
\insertRef{kz19}{PHInfiniteEstimates}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertstoml.R
\name{convertstoml}
\alias{convertstoml}
\title{Convert a proportional hazards regression to a multinomial regression.}
\usage{
convertstoml(survobj, covmat)
}
\arguments{
\item{survobj}{A survival object, with potentially right censoring.}

\item{covmat}{a matrix of covariates.}
}
\value{
a data set on which to apply conditional multinomial regression, corresponding to the proportional hazards regression analysis.
In order to run the line commented out below, you would need this:
# @importFrom mlogit mlogit.data
}
\description{
Convert a proportional hazards regression to a multinomial regression.
}
\details{
Implements version of \insertCite{kz19}{PHInfiniteEstimates}.  
The proportional hazards regression is converted to a multinomial regression logistic regression, and methods of \insertCite{kolassa16}{PHInfiniteEstimates} may be applied.
This function is intended to produce intermediate results to be passed to \code{convertmtol}, and then to \code{reduceLR} of \insertCite{kolassa97}{PHInfiniteEstimates}.  See examples in the \code{convertmtol} documentation.
}
\references{
\insertRef{kolassa97}{PHInfiniteEstimates}

\insertRef{kolassa16}{PHInfiniteEstimates}

\insertRef{kz19}{PHInfiniteEstimates}
}

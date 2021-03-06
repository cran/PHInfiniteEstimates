% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aalenjohansen.R
\name{aalenjohansen}
\alias{aalenjohansen}
\title{Calculate the Aalen-Johansen (1978) estimate in the Competing risk context.  See Aalen, Odd O., and Søren Johansen. "An Empirical Transition Matrix for Non-Homogeneous Markov Chains Based on Censored Observations." Scandinavian Journal of Statistics 5, no. 3 (1978): 141-50. Accessed January 15, 2021. http://www.jstor.org/stable/4615704.}
\usage{
aalenjohansen(times, causes)
}
\arguments{
\item{times}{Event times.}

\item{causes}{Causes, with 0 coded as censored, 1 as cause of interest, other for competing.}
}
\value{
a list with components
\itemize{
  \item times Unique times
  \item surv Aalen-Johansen estimator for cause 1.
}
}
\description{
Calculate the Aalen-Johansen (1978) estimate in the Competing risk context.  See Aalen, Odd O., and Søren Johansen. "An Empirical Transition Matrix for Non-Homogeneous Markov Chains Based on Censored Observations." Scandinavian Journal of Statistics 5, no. 3 (1978): 141-50. Accessed January 15, 2021. http://www.jstor.org/stable/4615704.
}

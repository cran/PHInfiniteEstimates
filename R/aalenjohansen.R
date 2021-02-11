#' Calculate the Aalen-Johansen (1978) estimate in the Competing risk context.  See Aalen, Odd O., and Søren Johansen. "An Empirical Transition Matrix for Non-Homogeneous Markov Chains Based on Censored Observations." Scandinavian Journal of Statistics 5, no. 3 (1978): 141-50. Accessed January 15, 2021. http://www.jstor.org/stable/4615704.
#'
#' @param times Event times.
#' @param causes Causes, with 0 coded as censored, 1 as cause of interest, other for competing.
#' @return a list with components
#' \itemize{
#'   \item times Unique times
#'   \item surv Aalen-Johansen estimator for cause 1.
#' }
#' @export
aalenjohansen<-function(times,causes){
  utimes<-sort(unique(times))
  fix<-kmall<-rep(NA,length(utimes))
  for(jj in seq(length(utimes))){
     kmall[jj]<-1-sum(causes[times==utimes[jj]]>0)/sum(times>=utimes[jj])
     fix[jj]<-sum(causes[times==utimes[jj]]==1)/sum(times>=utimes[jj])
  }
  kmall<-c(1,cumprod(kmall))[-length(utimes)]
  fhat<-c(0,cumsum(kmall*fix))
  utimes<-c(0,utimes)
  return(list(time=utimes,surv=1-fhat))
}

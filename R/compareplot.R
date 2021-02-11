#' Plot resuts of simcode
#'
#' @param simresults the result of simcode
#' @importFrom graphics axis legend plot points segments
#' @importFrom stats as.formula coefficients optimize rweibull uniroot
#' @return nothing.
#' @importFrom graphics plot lines legend abline
#' @export
compareplot<-function(simresults){
   nsim<-dim(simresults$out)[1]
   plot(c(0,1),c(0,1),type="n",xlab="True level",ylab="Nominal level")
   abline(0,1)
   lines((1:nsim)/nsim,sort(simresults$out[,"Waldpv"]),lty=2)
   lines((1:nsim)/nsim,sort(simresults$out[,"LRTpv"]),lty=3)
   legend(.18,.23,lty=1:3,legend=c("Ideal",
                           "Wald Test from Penalized Regression",
                           "Signed Root of the Likelihood Ratio Test"))
   plot(simresults$out[,"Waldpv"],
      simresults$out[,"LRTpv"]-simresults$out[,"Waldpv"],
      main="Comparison of Wald and SRLRT p-values",
      xlab="Wald p value",ylab="LRT p value - Wald p value")
}

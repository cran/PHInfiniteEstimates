#' Check how censoring impacts sampling properties of KM fit and log rank test.
#'
#' @param nsamp Number of MC samples
#' @param nobs Number of observations
#' @return biases of fits.
#' @importFrom stats rexp
#' @importFrom graphics hist
#' @importFrom survival survdiff
#' @export
checkcensor<-function(nsamp=1000,nobs=1000){
   out<-array(NA,c(nsamp,3))
   for(jj in seq(dim(out)[1])){
      dataset<-data.frame(t=rexp(nobs),g=rep(c(1,2),nobs/2),c=2*rexp(nobs))
      dataset$c<-dataset$c*dataset$g
      dataset$x<-pmin(dataset$t,dataset$c)
      dataset$delta<-(dataset$t<dataset$c)+0
      out[jj,1]<-1-pchisq(survdiff(Surv(x,delta)~g,data=dataset)$chisq,1)
      plot(a<-survfit(Surv(t,delta)~g,data=dataset),col=c(1,2),xlim=c(0,6))
      lines(seq(60)/10,exp(-seq(60)/10),col=3)
      out[jj,2:3]<-summary(a,times=3)$surv
   }
   hist(out[,1])
   return(apply(out[,-1],2,mean)-exp(-4))
}

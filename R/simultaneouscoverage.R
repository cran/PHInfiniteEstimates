#' Calculate simultaneous coverage of pointwise confidence intervals.
#'
#' Simulate exponential event times with expecation 1.  Simulate censoring times with expectation 2.  Calculate confidence intervals and check simultaneous coverage.
#'
#' @param nsamp Number of Monte Carlo samples.
#' @param nobs Number of observations per sample.
#' @return Simultaneous coverage proportion.
#' @export
#' @importFrom stats update
#' @examples
#' simultaneouscoverage(1000,20)
simultaneouscoverage<-function(nsamp,nobs){
   topok<-botok<-rep(NA,nsamp)
   for(i in 1:nsamp){
      t1<--log(runif(nobs))
      t2<--2* log(runif(nobs))
      tt<-apply(cbind(t1,t2),1,max)
      status<-t1<t2
      uuuu<-survfit(Surv(tt,status)~1)
      tail<-exp(-uuuu$time)
      topok[i]<-all((uuuu$upper[-nobs]>tail[-1])&(uuuu$upper[-nobs]>tail[-nobs]))
      botok[i]<-all((uuuu$lower[-nobs]>tail[-1])&(uuuu$lower[-nobs]>tail[-nobs]))
   }
   return(mean(topok&botok))
}

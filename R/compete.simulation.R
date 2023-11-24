#' Simulate from a competing risk model with correlated log normal errors, and plot various estimates.
#'
#' @param ncr Number of competing risks.
#' @param sig correlation among competing risks.
#' @param ns number of observations.
#' @importFrom stats median
#' @importFrom stats rnorm
#' @importFrom survival survfit
#' @export
compete.simulation<-function(ncr=4,sig=.8,ns=1000){
  common<-rnorm(ns)
  sep<-array(rnorm(ns*ncr),c(ns,ncr))
  obs<-exp(cbind(rnorm(ns),sig*outer(common,rep(1,ncr))+sqrt(1-sig^2)*sep))
  aggregate<-apply(obs,1,min)
  whichmin<-apply(obs,1,order)[1,]-1
  c1<-obs[,2]
  net<-survfit(Surv(aggregate,whichmin==1)~1)
  cruderesp<-aggregate  
  cruderesp[whichmin!=1]<-Inf 
  plot(c(0,max(aggregate)),c(0,1),type="n", xlab="Survival",ylab="Time",
     sub=paste("Lognormal Frailty Model, ",ncr," hazards, sigma=",sig,sep=""),
     main="Behavior of Measures for Competing Risk")
  lines(sort(aggregate),rev(seq(length(aggregate))/length(aggregate)),type="s")
  lines(sort(cruderesp),rev(seq(length(cruderesp))/length(aggregate)),
      lty=2,type="s")
  lines(net,conf.int=FALSE,lty=3)
  aj<-aalenjohansen(aggregate,whichmin)
  lines(aj$time,aj$surv,lty=4,col=2)
  legend(
#    median(range(aggregate)),1,
     0,.4,
     lty=1:4,col=c(1,1,1,2), 
     legend=c("Aggregate","Crude","Net","Aalen-Johansen"))
}

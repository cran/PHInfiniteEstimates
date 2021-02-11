#' Test size of asymptotic Cox tests.
#'
#' @param nsamp Number of MC samples
#' @param nobs Number of observations
#' @param ncov Number of covariates
#' @param randomcov Indicator of whether to draw random covariates.
#' @return level of two-sided test of nominal size 0.05.
#' @importFrom stats rexp
#' @export
testcox<-function(nsamp=1000,nobs=50,ncov=5,randomcov=TRUE){
   out<-rep(NA,nsamp)
   temp<-vector("list",length=ncov+2)
   names(temp)<-c("tt","dd",paste("x",1:ncov,sep=""))
   for(j in seq(length(out))){
      if((10*j/nsamp)==floor(10*j/nsamp)) cat("Progress")
      yy<-rexp(nobs)
      cc<-rexp(nobs)
      temp$tt <- pmin(yy,cc)
      temp$dd <- (yy>cc)+0
      if(randomcov){
         for(k in seq(ncov)) temp[[paste("x",k,sep="")]]<-(runif(nobs)>.5)+0
      }else{
         nn<-nobs
         for(k in seq(ncov)){
             temp[[paste("x",k,sep="")]]<-(2*((0:(nobs-1))%%nn))>nn
             nn<-floor(nn/2)
         }
      }
      out[j]<-summary(coxph(Surv(tt,dd)~x1+x2+x3+x4+x5,data=temp))$coef[1,5]}
   plot(seq(length(out))/length(out),sort(out),type="l",
      xlab="True p value",ylab="Nominal p value",
      main="Operating Characteristics of Proportional Hazards Regression",
      sub=paste( c(nobs,5,length(out),.5),
         c("Observatons","Covarites", "MC Samples","Censoring") ,collapse=","))
   abline(0,1,lty=2)
   return(mean(out>0.05))
}


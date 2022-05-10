#' Simulate operating characteristics of repaired Cox regression and competitors.
#'
#'
#' This function is intended to verify the operating characteristics of the approximate conditional inferential approach of \insertCite{kz19;textual}{PHInfiniteEstimates} to proportional hazards regression.  An exponential regression model, corresponding to the proportional hazards regression model, is fit to the data, and new data sets are simulated from this model.  P-values are calculated for these new data sets, and their empirical distribution is compared to the theoretical uniform distribution.
#' @param nobs number of observations in simulated data set.
#' @param k number of covariates in simulated data set.  Each covariate is dochotomous.
#' @param B odds of 1 vs. 0 in dichotomous variables.
#' @param c censoring proportion.
#' @param nsamp number of samples.
#' @param beta regression parameters, all zeros if null, and all the same value if a scalar.
#' @param add partial simulation results to be added to, or NULL if de novo.
#' @param half does nothing; provided for compatabilitity with simcode.
#' @param verbose Triggers verbose messages.
#' @param smoothfirst Triggers normal rather than dichotomous interest covariate.
#' @return a list with components
#' \itemize{
#'   \item out matrix with columns corresponding to p-values.
#' }
#' @importFrom stats runif
#' @export
heinzeschemper<-function(nobs=50,k=5,B=1,c=0,nsamp=1000,beta=NULL,add=NULL,half=NULL,verbose=FALSE,smoothfirst=FALSE){
   if (is.null(add)) {
      set.seed(202043125)
      start <- 0
   } else {
      outout <- rbind(add$out, array(NA, c(nsamp, dim(add$out)[2])))
      start <- dim(add$out)[1]
      set.seed(add$seed)
   }
   if(is.null(beta)) beta<-rep(0,k)
   if(length(beta)==1) beta<-rep(beta,k)
   gg<-as.formula(paste("Surv(times,delta)~",paste("x",seq(k),sep="",collapse="+")))
   hh<-as.formula(paste("Surv(times,delta)~",paste("x",(2:k),sep="",collapse="+")))
   d1 <- Sys.time()
   cenp<-rep(NA,nsamp)
   for(kk in seq(nsamp)){
      if (verbose) {
         d2 <- Sys.time()
         message("kk=",kk," of ",nsamp,".  Completion time ",(d2 - d1) * (nsamp - kk)/kk + d2)
      }
      randdat<-if(smoothfirst) cbind(rnorm(nobs),as.data.frame(array(runif(nobs*(k-1))>(B/(1+B)),c(nobs,k-1)))+0) else as.data.frame(array(runif(nobs*k)>(B/(1+B)),c(nobs,k)))+0
      names(randdat)<-paste("x",seq(k),sep="")
      randdat$x<-as.matrix(randdat)
      randdat$times<--log(runif(nobs))/exp(randdat$x%*%beta)
      randdat$delta<-runif(nobs) > c
      cenp[kk]<-mean(randdat$delta)
      randdat$y<-Surv(randdat$t,randdat$delta)
#     cat("About to run fixcoxph\n")
      repairedfit<-fixcoxph(randdat,randdat$x,"x1")
      penalizedout<-coxphf(gg,randdat,maxit=400,maxstep=0.05)
      penalizedoutsmaller<-coxphf(hh,randdat,maxit=400,maxstep=0.05)
      myout<-summarizefits(repairedfit,penalizedout,penalizedoutsmaller,"x1")
      if((start+kk)==1){
         outout<-array(NA,c(nsamp,length(myout)))
         dimnames(outout)<-list(NULL,names(myout))
      }
      outout[start+kk,]<-myout
   }
   return(list(out=outout,seed=.Random.seed,settings=list(nobs=nobs,k=k,B=B,c=c,nsamp=nsamp,beta=beta,half=half,verbose=verbose),cenp=cenp))
}

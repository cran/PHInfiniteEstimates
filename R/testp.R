#' Simulate Weibull survival data from a model, perform reduction to remove infinite estimates, and calculate p values.
#'
#' Operating characteristics for the approximate conditional inferential approach to proportional hazards.
#'
#' This function is intended to verify the operating characteristics of the approximate conditional inferential approach of \insertCite{kz19;textual}{PHInfiniteEstimates} to proportional hazards regression.  A Weibull regression model, corresponding to the proportional hazards regression model, is fit to the data, and new data sets are simulated from this model.  P-values are calculated for these new data sets, and their empirical distribution is compared to the theoretical uniform distribution.
#' @param dataset the data set to use
#' @param myformula the formula for the Cox regression
#' @param iv name of the variable of interest, as a character string
#' @param ctime fixed censoring time
#' @param nsamp number of samples.
#' @param add preliminary results, if any.
#' @param nobs number of observations in target models, if different from that of dataset.
#' @param half logical flag triggering a less extreme simulation by dividing the Weibull regression parameters in half.
#' @param verbose logical flag triggering intermediate messaging.
#' @return a list with components
#' \itemize{
#'   \item out matrix with columns corresponding to p-values.
#'   \item seed random seed
#'   \item bad unused.
#'   \item srreg parametric lifetime regression
#' }
#' @importFrom survival Surv coxph survreg
#' @importFrom stats pnorm
#' @export
#' @references 
#' \insertRef{kz19}{PHInfiniteEstimates}
#' @examples
#' data(breast)
#'\donttest{
#' breattestp<-testp(breast,Surv(TIME,CENS)~ T+ N+ G+ CD,"T",72,nsamp=100,verbose=TRUE)
#'}
testp<-function(dataset,myformula,iv,ctime,nsamp=10000,add=NULL,nobs=NA,half=FALSE,verbose=FALSE){
   bad<-NULL
   gg<-as.formula(paste("Surv(newt,nc)~",deparse(myformula[[3]])))
   sr1<-survreg(formula=myformula,data=dataset,x=TRUE,y=TRUE)
   if(is.na(nobs)) nobs<-dim(sr1$x)[1]
   use<-seq(nobs)
   csr1<-coefficients(sr1)
   ivl<-(names(csr1)==iv)[-1]
   if(half) csr1<-csr1/2
   csr1[iv]<-0
   randdat<-cbind(as.data.frame(list(newt=rep(NA,nobs))),sr1$x[use,])
   if(is.null(add)){
      outout<-array(NA,c(nsamp,6))
      dimnames(outout)<-list(NULL,c("Waldp","LRTp","Est","SE","SRLRT","dim"))
      start<-0
      set.seed(194837485)
   }else{
      outout<-rbind(add$out,array(NA,c(nsamp,6)))
      start<-dim(add$out)[1]
      set.seed(add$seed)
   }
   d1<-Sys.time()
   for(kk in seq(nsamp)){
      if(verbose){
         d2<-Sys.time()
         if(verbose) message("kk=",kk," of ",nsamp,".  Completion time ",(d2-d1)*(nsamp-kk)/kk+d2)
         
      }
      if(nobs<dim(sr1$x)[1]) use<-sample(dim(sr1$x)[1],nobs)
      for(j in seq(nobs)){
         elpred<-exp(sr1$x[use,]%*%csr1)
         randdat$newt[j]<-rweibull(1,shape=1/sr1$scale,scale=elpred[j])
      }
      randdat$nc<-randdat$newt<ctime
      randdat$newt[randdat$newt>=ctime]<-ctime
      randdat$y<-Surv(randdat$newt,randdat$nc)
      xxx <- sr1$x[use,-1]# Remove intercept
      randdat$x<-xxx
      repairedfit <- fixcoxph(randdat,xxx,iv,verbose=verbose)
      if(is.atomic(repairedfit)){
#        cat("Branch a\n")
         if(repairedfit==-1){
            outout[start+kk,1]<-1
         }else{
            if(verbose) message("Error 2")
         }
      }else{
#        cat("Branch b\n")
         if(repairedfit$iter>0){
            print(repairedfit$keepme)
            pvs2<-summary(repairedfit)$coefficients[iv,5]
            km<-repairedfit$keepme
            km[iv]<-TRUE
            v1<-try(det(repairedfit$var[km,km,drop=FALSE]))
            if(inherits(v1,"try-error")){
               message("Error from v1")
               browser()
            }
            km[iv]<-FALSE
            v2<-try(if(any(km)) det(repairedfit$svar[km,km,drop=FALSE]) else 1)
            if(inherits(class(v2),"try-error")){
               message("Error from v2")
               browser()
            }
            ssd<-sqrt(abs(v1/v2))
#           if(is.na(ssd)){
#              message("NA in ssd")
#              browser()
#           }
            est<-repairedfit$coefficients[iv]
            pvs<-2*pnorm(-abs(repairedfit$coefficients/sqrt(diag(repairedfit$var))))
            outout[start+kk,1]<-pvs2
            outout[start+kk,5]<-sqrt(2*(repairedfit$loglik[2]-repairedfit$dropone))*sign(est)
            pvs3<-1-pchisq(2*(repairedfit$loglik[2]-repairedfit$dropone),1)
            outout[start+kk,2]<-pvs3
            outout[start+kk,3]<-est
            outout[start+kk,4]<-ssd
            outout[start+kk,6]<-sum(!is.na(repairedfit$coefficients))
#           browser()
#           cat("End branch b1\n")
         }else{
            if(verbose) message("Error 2")
         }
      }
#     cat("Ping\n")
      if(verbose) message(outout[start+kk])
      if(is.na(outout[start+kk])){
         dump("randdat")
#        browser()
      }
   }
   return(list(out=outout,seed=.Random.seed,bad=bad,srreg=sr1))
}

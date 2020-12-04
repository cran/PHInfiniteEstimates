#' Proportional hazards partial likelihood, using Breslow method for ties, excluding some observations.
#'
#' This function implements the approximate conditional inferential approach of \insertCite{kz19;textual}{PHInfiniteEstimates} to proportional hazards regression.
#' @param beta parameter vector.
#' @param fit Output from a Cox PH regression, with x=TRUE and y=TRUE
#' @param exclude data set with stratum and patient number to exclude.
#' @param minus logical flag to change sign of partial likelyhood
#' @param keeponly variables to retain.  Keep all if this is null or NA.
#' @param justd0 logical variable, indicating whether to calculate only the function value and skip derivatives.
#' @return a list with components
#' \itemize{
#'   \item d0 partial likelihood
#'   \item d1 first derivative vector
#'   \item d2 second derivative matrix
#' }
#' @export
#' @importFrom stats update
#' @references
#' \insertRef{kz19}{PHInfiniteEstimates}
newllk<-function(beta,fit,exclude=NULL,minus=FALSE,keeponly=NULL,justd0=FALSE){
#  cat("On entry to newllk beta",beta,"keeponly",keeponly,"\n")
#  print("fit");print(fit)
   d0<-0
   if(is.null(keeponly)) keeponly<-rep(TRUE,length(beta))
   if(is.na(keeponly[1])) keeponly<-rep(TRUE,length(beta))
   if(!justd0){
      d1<-rep(0,sum(keeponly))
      d2<-array(0,rep(sum(keeponly),2))
   }else{
      d1<-rep(NA,sum(keeponly))
      d2<-array(NA,rep(sum(keeponly),2))
   }
   uuu<-sort(unique(fit$y[fit$y[,2]==1,1]))
   if(!is.array(fit$x)) fit$x<-matrix(fit$x,ncol=1)
   for(mm in seq(length(uuu))){
      thistime<-uuu[mm]
      keep<-rep(TRUE,dim(fit$x)[1])
      keep[as.numeric(exclude$obs[exclude$str==mm])]<-FALSE
      eventset<-(fit$y[,1]==thistime)&(fit$y[,2]==1)&keep
      riskset<-(fit$y[,1]>=thistime)&keep
#     cat("mm",mm,"riskset",riskset,"\n")
#     print(fit$x[riskset,])
      if(sum(riskset)>1){
         xxx<-fit$x[eventset,keeponly,drop=F]
         xxxr<-fit$x[riskset,keeponly,drop=F]
#        cat("dim(xxxr)",dim(xxxr),"\n")
#        cat("length(beta)",length(beta),"\n")
         ddd<-sum(eventset)
         numers<-as.vector(xxxr%*%beta[keeponly])
         denom<-sum(exp(numers))
         pp<-exp(numers)/denom
         d0<-d0+sum(xxx%*%beta[keeponly])- ddd*log(denom)
         if(!justd0){
            mmm<-t(xxxr)%*%pp
            d1<-d1+apply(xxx,2,sum)-ddd*mmm
#           message("dim(xxxr)",dim(xxxr))
            dd2<-try(ddd*(t(apply(xxxr,2,"*",pp))%*%xxxr-mmm%*%t(mmm)))
            if(inherits(dd2,"try-error")){
               message("Non-conformable array")
            }else{
               d2<-d2-ddd*(t(apply(xxxr,2,"*",pp))%*%xxxr-mmm%*%t(mmm))
            }
         }
      }
   }
   out<-if(minus){-d0}else{list(d0=d0,d1=d1,d2=d2)}
   if(minus){
      out<- -d0
#     attr(out,'gradient')<- -d1
#     attr(out,'hessian')<- -d2
   }else{
      out<-list(d0=d0,d1=d1,d2=d2)
   }
#  cat("On exit from newllk out")
#  print(out)
   return(out)
}

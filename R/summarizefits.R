#' Summarize proportional hazards fits
#'
#' @param repairedfit coxph fit
#' @param penalizedout coxphf fit
#' @param penalizedoutsmaller smaller coxphf fit
#' @param iv name of the variable of interest, as a character string
#' @param verbose logical flag triggering intermediate messaging.
#' @return a vector with components
#' \itemize{
#'     \item Wald p-value from the Cox regression fit.
#'     \item partial likelihood ratio p-value from Cox regression fit.
#'     \item parameter estimate from the Cox regression fit.
#'     \item standard error from the Cox regression fit.
#'     \item Conditional Skovgaard standard error from the Cox regression fit.
#'     \item Signed root of the partial likelihood ratio statistic from Cox regression fit.
#'     \item partial likelihood ratio statistic p-value from coxphf
#'     \item Wald p-value from coxphf
#'     \item parameter estimate from coxphf
#'     \item standard error from coxphf
#'     \item number of parameters in reduced fit.
#' }
#' @importFrom survival Surv coxph survreg
#' @importFrom stats pnorm
#' @importFrom coxphf coxphf
#' @export
#' @references 
#' \insertRef{kz19}{PHInfiniteEstimates}
#'
summarizefits<-function(repairedfit,penalizedout,penalizedoutsmaller,iv,verbose=TRUE){
#  print(repairedfit$var)
   nn<- c("Waldpv", "LRTpv", "Est","SE", "Skov SE", "SRLRT", "Waldpvf","LRTpvf","Estf", "SEf","npar","nspar")
   mindiff<-rep(0,2)
   names(mindiff)<-c("Unpenalized","Penalized")
   out<-rep(NA,length(nn))
   names(out)<- nn
   if (is.atomic(repairedfit)) {
      if (repairedfit == -1) {
         out["Waldpv"]<-1
         if (verbose) message("Error 1")
      }else{
         if (verbose) message("Error 2")
         browser()
      }
   } else {
      if(!is.na(repairedfit$var[iv,iv])){
         if(repairedfit$var[iv,iv]>=0){
            myse<-sqrt(repairedfit$var[iv,iv])
         }else{
            message("There should not be a negative variance")
            myse<-NA
         }
      }else{
         myse<-NA
      }
#     message("repairedfit$coefficients",paste(repairedfit$coefficients,collapse=","))
      if(length(repairedfit$coefficients)==0){#This branch gets triggered if repairedfit gets set to a scalar.
         out["Waldpv"]<-1
      }else{
         if(is.na(repairedfit$coefficients[iv])){
            out["Waldpv"]<-1
         }else{
            if(abs(repairedfit$coefficients[iv])==Inf) {
               out["Waldpv"]<-0
#              message("Infinite interest parameter estimate")
            }else{
#              message("repairedfit$wald.test",repairedfit$wald.test)
               out["Waldpv"] <- summary(repairedfit)$coefficients[iv, 5]
            }
         }
     }
      km <- repairedfit$keepme
      km[iv] <- TRUE
      v1 <- try(det(repairedfit$var[km, km, drop = FALSE]),silent=TRUE)
         if(inherits(v1, "try-error")) {
            message("Error from v1")
            browser()
         }
      km[iv] <- FALSE
      v2 <- try(if(any(km)) det(repairedfit$svar[km,km,drop=FALSE]) else 1,silent=TRUE)
      if(inherits(class(v2), "try-error")) {
         message("Error from v2")
         browser()
      }
      out["Skov SE"] <- sqrt(abs(v1/v2))
      out["Est"] <- repairedfit$coefficients[iv]
      out["LRTpv"] <- 1-pchisq(2 * (repairedfit$loglik[2] - repairedfit$dropone), 1)
      dlrt<-repairedfit$loglik[2] - repairedfit$dropone
      mindiff[1]<-min(mindiff[1],dlrt)
      dlrt=max(dlrt,0)
      out["SRLRT"]<-sqrt(2 * dlrt) * sign(out["Est"])
      out["Waldpvf"]<-penalizedout$p[iv]
      temp<-diff(penalizedout$loglik)-diff(penalizedoutsmaller$loglik)
      mindiff[2]<-min(mindiff[2],temp)
      temp<-max(temp,0)
      out["LRTpvf"] <- 1-pchisq(2*max(temp,0),1)
      out["Estf"]<-penalizedout$coefficients[iv]
      out["SEf"]<-sqrt(penalizedout$var[iv,iv])
      out["SE"]<-sqrt(repairedfit$var[iv,iv])
      out["npar"]<-sum(repairedfit$keepme)
      out["nspar"]<-repairedfit$nssmall
      if(any(is.na(out[c("Waldpv", "LRTpv", "Est","SRLRT","LRTpvf","Estf", "SEf","npar","nspar")]))){
         message("NA among output")
#        print("out"); print(out)
#        browser()
      }
   }
   return(out)
}

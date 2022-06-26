#' Newton Raphson Fitter for partial likelihood
#'
#' This function implements the approximate conditional inferential approach of \insertCite{kz19;textual}{PHInfiniteEstimates} to proportional hazards regression.
#' @param fit Output from a Cox PH regression, with x=TRUE and y=TRUE
#' @param exclude data set with stratum and patient number to exclude.
#' @param start Starting value
#' @param touse columns of the design matrix to use.
#' @param usecc Logical variable indicating whether to use a continuity correction.
#' @examples
#' bfit<-coxph(Surv(TIME,CENS)~T+N+CD,data=breast,x=TRUE)
#' bestbeta(bfit)
#' bestbeta(bfit,usecc=TRUE)
#' @export
#' @importFrom stats update
#' @importFrom stats nlm
#' @return Fitted survival analysis regression parameter of class coxph
#' @references
#' \insertRef{kz19}{PHInfiniteEstimates}
bestbeta<-function (fit, exclude = NULL, start = NULL,touse=NA,usecc=FALSE){
#   cat("Just entered bestbeta touse",touse,"\n")
#   browser()
    if(is.null(start)) start <- rep(0, dim(fit$x)[2])
#   message("In bestbeta before raw call to newllk touse",touse)
    llkout <- newllk(start, fit, exclude = exclude,keeponly=touse,justd0=FALSE)
#   browser()
    cc1<-if(usecc) sign(llkout$d1[1])/2 else 0.0
#   cat("In bestbeta after raw call to newllk\n")
    if(is.na(touse[1])){
       touse <- rep(TRUE, dim(fit$x)[2])
       names(touse)<-dimnames(fit$x)[[2]]
    }
#   browser()
    newvar<-array(NA,rep(length(touse),2))
    dimnames(newvar)<-list(names(touse),names(touse))
#   print(names(touse))
#   print.default(newvar)
    newcoef<-rep(NA,length(touse))
    if(any(touse)){
#      cat("In bestbeta about to call nlmo touse",touse,"\n")
       aaa<-fit
       aaa$x<-fit$x[,touse,drop=FALSE]
       nlmo<-nlm(newllk,rep(0,sum(touse)),hessian=TRUE,fit=aaa,exclude=exclude,minus=TRUE,cc1=cc1)
#      browser()
       newcoef[touse]<-nlmo$estimate
       names(newcoef)<-names(touse)
#      message("In bestbeta Hessian diagonal",paste(diag(nlmo$hessian),collapse=" "))
#      message("In bestbeta inv Hessian diagonal",paste(diag(solve(nlmo$hessian)),collapse=" "))
#      if(any(diag(nlmo$hessian)<0)){
#         message("In bestbeta Hessian should be negative definite")
#         browser()
#      }
       temp<-try(solve(nlmo$hessian),silent=TRUE)
       if(inherits(temp,"try-error")){
          message("Inverting hessian failed")
       }else{
          newvar[touse,touse]<-temp
#         cat("temp"); print(temp); cat("newvar"); print(newvar)
       }
       score<-NA
#      print("touse");print(touse)
#      print("newcoef");print(newcoef)
#      print("dim(fit$x)");print(dim(fit$x))
       linear.predictors<-try(as.vector(fit$x[,touse,drop=FALSE]%*%newcoef[touse])-as.vector(apply(fit$x[,touse,drop=FALSE],2,mean)%*%newcoef[touse]),silent=TRUE)
       if(inherits(linear.predictors,"try-error")){
          message("Error in linear predictor"); browser()
       }
       vv<-try(solve(newvar[touse,touse,drop=FALSE],newcoef[touse]),silent=TRUE)
       if(inherits(vv,"try-error")){
          message("Inverting matrix for wald test failed")
          wald.test<-NA
       }else{
          wald.test<-newcoef[touse]%*%vv
#         cat("temp"); print(temp); cat("newvar"); print(newvar)
       }
    }else{
       var<-array(NA,c(length(newcoef),length(newcoef)))
       nlmo<-list(minimum=0,iterations=0)
       score=0
       wald.test=0
       linear.predictors<-rep(0,dim(fit$x)[1])
    }
    out<-list(coefficients=newcoef,var=newvar,loglik=c(llkout$d0,-nlmo$minimum),score=score, wald.test=wald.test, linear.predictors=linear.predictors,iterations=nlmo$iterations)
#   print(out)
    class(out)<-"coxph"
    return(out)
}

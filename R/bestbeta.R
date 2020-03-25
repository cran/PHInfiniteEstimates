#' Newton Raphson Fitter for partial likelihood
#'
#' This function implements the approximate conditional inferential approach of \insertCite{kz19;textual}{PHInfiniteEstimates} to proportional hazards regression.
#' @param fit Output from a Cox PH regression, with x=TRUE and y=TRUE
#' @param exclude data set with stratum and patient number to exclude.
#' @param start Starting value
#' @param touse columns of the design matrix to use.
#' @export
#' @importFrom stats update
#' @importFrom stats nlm
#' @return Fitted survival analysis regression parameter of class coxph
#' @references
#' \insertRef{kz19}{PHInfiniteEstimates}
bestbeta<-function (fit, exclude = NULL, start = NULL,touse=NA){
#   cat("Just entered bestbeta touse",touse,"\n")
#   browser()
    start <- rep(0, dim(fit$x)[2])
#   cat("In bestbeta before raw call to newllk touse",touse,".\n")
    llkout <- newllk(start, fit, exclude = exclude,keeponly=touse,justd0=TRUE)
#   cat("In bestbeta after raw call to newllk\n")
    if(is.na(touse[1])){
       touse <- rep(TRUE, dim(fit$x)[2])
    }
#   browser()
    newvar<-array(NA,rep(length(touse),2))
    if(any(touse)){
#      cat("In bestbeta about to call nlmo touse",touse,"\n")
       aaa<-fit
       aaa$x<-fit$x[,touse,drop=FALSE]
       nlmo<-nlm(newllk,rep(0,sum(touse)),hessian=TRUE,fit=aaa,exclude=exclude,minus=TRUE)
#      browser()
       newcoef<-rep(NA,length(touse))
       newcoef[touse]<-nlmo$estimate
       names(newcoef)<-names(touse)
       temp<-try(solve(nlmo$hessian))
       if(inherits(temp,"try-error")){
          message("Inverting hessian failed")
       }else{
          newvar[touse,touse]<-temp
       }
       score<-NA
#      print("touse");print(touse)
#      print("newcoef");print(newcoef)
#      print("dim(fit$x)");print(dim(fit$x))
       linear.predictors<-try(as.vector(fit$x[,touse,drop=FALSE]%*%newcoef[touse])-as.vector(apply(fit$x[,touse,drop=FALSE],2,mean)%*%newcoef[touse]),silent=TRUE)
       if(inherits(linear.predictors,"try-error")){
          message("Error in linear predictor"); browser()
       }
       wald.test<-try(newcoef[touse]%*%solve(newvar[touse,touse,drop=FALSE],newcoef[touse]),silent=TRUE)
    }else{
       newcoef<-NULL
       var<-array(0,c(0,0))
       loglik<-list(d0=0)
       nlmo<-list(minimum=0,iterations=0)
       score=0
       wald.test=0
       linear.predictors<-rep(0,dim(fit$x)[1])
    }
    out<-list(coefficients=newcoef,var=newvar,loglik=c(llkout$d0,-nlmo$minimum),score=score, wald.test=wald.test, linear.predictors=linear.predictors,iterations=nlmo$iterations)
    class(out)<-"coxph"
    return(out)
}

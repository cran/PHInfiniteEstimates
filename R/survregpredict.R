#' Fit survival probabilties from a survreg object.
#' @param fit a survreg object.  This should not contain strata().  It also must use the log transformation.
#' @param newdata a new data set with covariates from the fit.
#' @param time a time value (on the original, and not log, scale).
#' @importFrom stats delete.response model.frame model.matrix predict terms
#' @examples
#' #Fit the survival probability for an individual with extent 1 and
#' #differentiation 2 at 700 days from a Weibull regression using the 
#' #colon cancer data set distributed as part of the survival package.
#' fit<-survreg(Surv(time,status)~factor(extent)+differ,data=colon)
#' survregpredict(fit,data.frame(extent=1,differ=2),700)
#' @export
survregpredict<-function(fit,newdata,time){
# Give the survival fit and standard error using survival pack-
# age tools acting on the fit.  
# I patterned this from predict.lm .  delete.response doesn't
# seem to do anyting for these survival models, but I'll try
# anyway.
   Terms<-delete.response(terms(fit))
# Build the model frame using the list of terms, and  chanage 
# the frame into a design matrix, respecting, for example, 
# levels of factors.  
   z<-model.matrix(Terms, model.frame(Terms,newdata,xlev=fit$xlevels),
      contrasts.arg = fit$contrasts)
# z is now a matrix with a column for the intercept, and 
# one row.  Change this to a vector without the intercept.
   z<-z[1,-1,drop=TRUE]
# The distribution is coded into the fitted structure; this in 
# turn codes the rescaled distribution.  It also codes the 
# transformation, but since we're always using the log, we ig-
# nore this.  Next function will return a vector with four 
# entries; the first is CDF and the second is a density.
   dist<-survival::survreg.distributions[[
      survival::survreg.distributions[[fit$dist]]$dist]]
   out<-predict(fit,type="lp",newdata=newdata,se.fit=TRUE)
   ord<-(log(time)-out$fit)/fit$scale
   do<-dist$density(ord)[1:2]
   out<-c(1-do[1],NA)
   v<-do[2]*c(c(1,z)/fit$scale,ord)
   out[2]<-sqrt(v%*%fit$var%*%v)
   names(out)<-c("fit","se")
   return(out)
}

convertmtol<-function(xmat,str,yvec,subjects){
#' Convert a polytomous regression to a conditional logistic regression.
#'
#' @param xmat regression matrix
#' @param str stratum label
#' @param yvec vector of responses
#' @param subjects vector of subject labels passed directly to the output.
#' @return a data set on which to apply conditional logistic regression, corresponding to the multinomial regression model.
#' @export
#' @details
#' Implements version of \insertCite{kolassa16}{PHInfiniteEstimates}.  
#' The multinomial regression is converted to a conditional logistic regression, and methods of \insertCite{kolassa97}{PHInfiniteEstimates} may be applied.
#' This function differs from \code{convertbaselineltolr} of this package in that the former treats the richer data structure of package \code{mlogit}, and this function treats a less complicated structure.
#' Data in the example is the breast cancer data set \code{breast} of package \code{coxphf}.
#' @references 
#' \insertRef{kolassa97}{PHInfiniteEstimates}
#'
#' \insertRef{kolassa16}{PHInfiniteEstimates}
#' @examples
#' #Uses data set breast from package coxphf.
#' data(breast)
#' out<-convertstoml(Surv(breast$TIME,breast$CENS),breast[,c("T","N","G","CD")])
#' out1<-convertmtol(out[,c("T","N","G","CD")],out[,"chid"],out[,"choice"],
#'    out[,"patients"])
#' glmout<-glm.fit(out1$xmat,out1$y,family=binomial())
#' #In many practice examples, the following line shows which observations to retain
#' #in the logistic regression example.
#' moderate<-(fitted(glmout)<1-1.0e-8)&(fitted(glmout)>1.0e-8)
#' # Proportional hazards fit illustrating infinite estimates.
#' coxph(Surv(TIME,CENS)~ T+ N+ G+ CD,data=breast)
#' # Wrong analysis naively removing covariate with infinite estimate
#' coxph(Surv(TIME,CENS)~ T+ N+ CD,data=breast)
#' summary(glm((CENS>22)~T+N+G+CD,family=binomial,data=breast))
#'\donttest{
#' out2<-reduceLR(out1$xmat,yvec=out1$y,keep="CD")
#' bestcoxout<-coxph(Surv(TIME,CENS)~ T+ N+ G+ CD,data=breast,
#'    subset=as.numeric(unique(out1$subjects[out2$moderate])))
#'}
   stru<-unique(str)
   z<-array(0,c(length(yvec),length(stru)))
   nvec<-rep(NA,length(yvec))
   dimnames(z)<-list(NULL,paste("str",stru,sep=""))
   for(j in seq(length(stru))) {
      z[str==stru[j],j]<-1
      nvec[str==stru[j]]<-sum(yvec[str==stru[j]])
   }
   return(list(xmat=cbind(z,as.matrix(xmat)),n=nvec,y=yvec,str=str,subjects=subjects,
      use=dim(z)[2]+seq(dim(xmat)[2])))
}

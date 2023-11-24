#' Perform inference on conditional sample space.
#'
#' This function performs classical frequentist statistical inference to a discrete multivariate canonical exponential family.  
#' It produces the maximum likelihood estimator, one- and two-sided p-values for the test that model parameters are zero, and providing confidence intervals for the parameters.  The discrete probability model is given by a set of possible values of the random vectors, and null weights for these vectors.  
#' Such a discrete probability model arises in logistic regression, and this function is envisioned to be applied to the results of a network algorithm for conditional logistic regression.  
#' Examples apply this to data from \insertCite{mehtapatel;textual}{PHInfiniteEstimates}, 
#' citing \insertCite{goorinetal87;textual}{PHInfiniteEstimates}.
#' @param netout List of the sort provided by network.
#' @param alpha Test level, or 1- confidence level.
#' @param rng Range of possible parameter values.
#' @param alternative String indicating two- or one-sided alternative, and, if one-sided, direction.
#' @return List of outputs, including
#' \itemize{
#'   \item ospv Observed one-sided p values
#'   \item tspv Observed two-sided p value.
#'   \item ci confidence interval.
#'   \item estimate Maximum conditional likelihood estimator.
#'   \item null.value Value of parameter under null hypothesis.
#'   \item data.name Name of data set
#'   \item method Method used to generate test.
#'   \item statistic sufficient statistic value for inference variable.
#'   \item p.value p.value
#'   \item conf.int confidence interval.
#'   \item alternative String indicating two- or one-sided alternative, and, if one-sided, direction.
#' }
#' and including standard stats:::orint.htest components, and of class htest.
#' @examples
#' #Columns in table are:
#' # Lymphocytic Infiltration (1=low, 0=high)
#' # Sex (1=male, 0=female)
#' # Any Ostioid Pathology (1=yes, 0=no)
#' # Number in LI-Sex-AOP group
#' # Number in LI-Sex-AOP group with disease free interval greater than 3 y
#' goorin<-data.frame(LI=c(0,0,0,0,1,1,1,1),Sex=c(0,0,1,1,0,0,1,1),
#'    AOP=c(0,1,0,1,0,1,0,1),N=c(3,2,4,1,5,5,9,17),Y=c(3,2,4,1,5,3,5,6))
#'\donttest{
#' netout<-network(goorin[,1:3],goorin[,4],conditionon=1:3,resp=goorin[,5])
#' inference(netout)
#'}
#' @references
#' \insertRef{mehtapatel}{PHInfiniteEstimates}
#'
#' \insertRef{goorinetal87}{PHInfiniteEstimates}
#' @export
inference<-function(netout,alpha=0.05,rng=c(-5,5),
   alternative = c("two.sided", "less", "greater")){
   alternative <- match.arg(alternative)
   nv<-dim(netout$possible)[2]
   ci<-ospv<-array(NA,c(2,nv))
   estimate<-rep(NA,nv)
   for(jj in seq(nv)){
      ospv[,jj]<-c(sum(netout$count[netout$possible[,jj]<=netout$obsd[jj]]),
         sum(netout$count[netout$possible[,jj]>=netout$obsd[jj]]))/sum(netout$count)
      ll<-function(psi,netout,jj){
         rrr<-sum(netout$count[netout$possible[,jj]==netout$obsd[jj]]*
            exp(psi*netout$obsd[jj]))/sum(netout$count*exp(netout$possible[,jj]*psi))
         return(rrr)
      }
#     browser()
      estimate[jj]<-optimize(ll,rng,netout,jj,maximum=TRUE)$maximum
      ff<-function(psi,netout,alpha){
         probs<-netout$count*exp(psi*netout$possible[,jj])
         pd<-sum(probs[netout$possible[,jj]>=netout$obsd[jj]])/sum(probs)-(alpha/2)
         return(pd)
      }
      gg<-function(psi,netout,alpha){
         probs<-netout$count*exp(psi*netout$possible[,jj])
         pd<-sum(probs[netout$possible[,jj]<=netout$obsd[jj]])/sum(probs)-(alpha/2)
         return(pd)
      }
      ci[,jj]<-c(uniroot(ff,rng,netout,alpha/2)$root, 
         uniroot(gg,rng,netout,alpha/2)$root)
   }
#  browser()
   tspv<-2*apply(ospv,2,min)
# Borrowed and modified from t.test.default
   pval <- if (alternative == "less") ospv[1] else if (alternative == "greater") ospv[2] else tspv 
   names(netout$obsd)<-paste("Logistic Regression Sufficient Statistic for",netout$comparison)
   nulllor<-0
   names(nulllor)<-"Log odds ratio"
   cat("conditionon",netout$conditionon,"\n")
   infout<-list(ospv=ospv,tspv=tspv,ci=ci,estimate=estimate,null.value=0,
      data.name=netout$data.name,
      method=paste("Exact Logistic Regression Conditioned on\n",
         paste(netout$conditioned,collapse=" ")),
      statistic=netout$obsd,p.value=pval,conf.int=ci, alternative=alternative)
   if(nv==1){
      class(infout)<-"htest"
   }
   return(infout)
}

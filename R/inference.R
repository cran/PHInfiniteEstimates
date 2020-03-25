#' Perform inference on conditional sample space.
#'
#' This function performs classical frequentist statistical inference to a discrete multivariate canonical exponential family.  
#' It produces the maximum likelihood estimator, one- and two-sided p-values for the test that model parameters are zero, and providing confidence intervals for the parameters.  The discrete probability model is given by a set of possible values of the random vectors, and null weights for these vectors.  
#' Such a discrete probability model arises in logistic regression, and this function is envisioned to be applied to the results of a network algorithm for conditional logistic regression.  
#' Examples apply this to data from \insertCite{mehtapatel;textual}{PHInfiniteEstimates}, 
#' citing \insertCite{goorinetal87;textual}{PHInfiniteEstimates}.
#' @param out List of the sort provided by network.
#' \itemize{
#'   \item possible   matrix with vectors of possible unconditioned values of the sufficient statistic.
#'   \item count   count of entries in the conditional distribution.
#'   \item obsd   Observed value of unconditioned sufficient statistics.
#' }
#' @param alpha Test level, or 1- confidence level.
#' @param rng Range of possible parameter values.
#' @return List with components:
#' \itemize{
#'   \item ospv Observed one-sided p values
#'   \item tspv Observed two-sided p value.
#'   \item ci confidence interval.
#'   \item mle Maximum likelihood estimator.
#' }
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
#' out<-network(goorin[,1:3],goorin[,4],conditionon=1:3,resp=goorin[,5])
#' inference(out)
#'}
#' @references
#' \insertRef{mehtapatel}{PHInfiniteEstimates}
#'
#' \insertRef{goorinetal87}{PHInfiniteEstimates}
#' @export
inference<-function(out,alpha=0.05,rng=c(-5,5)){
   nv<-dim(out$possible)[2]
   ci<-ospv<-array(NA,c(2,nv))
   mle<-rep(NA,nv)
   for(jj in seq(nv)){
      ospv[,jj]<-c(sum(out$count[out$possible[,jj]<=out$obsd[jj]]),
         sum(out$count[out$possible[,jj]>=out$obsd[jj]]))/sum(out$count)
      ll<-function(psi,out,jj){
         rrr<-sum(out$count[out$possible[,jj]==out$obsd[jj]]*
            exp(psi*out$obsd[jj]))/sum(out$count*exp(out$possible[,jj]*psi))
         return(rrr)
      }
#     browser()
      mle[jj]<-optimize(ll,rng,out,jj,maximum=TRUE)$maximum
      ff<-function(psi,out,alpha){
         probs<-out$count*exp(psi*out$possible[,jj])
         pd<-sum(probs[out$possible[,jj]>=out$obsd[jj]])/sum(probs)-(alpha/2)
         return(pd)
      }
      gg<-function(psi,out,alpha){
         probs<-out$count*exp(psi*out$possible[,jj])
         pd<-sum(probs[out$possible[,jj]<=out$obsd[jj]])/sum(probs)-(alpha/2)
         return(pd)
      }
      ci[,jj]<-c(uniroot(ff,rng,out,alpha/2)$root, 
         uniroot(gg,rng,out,alpha/2)$root)
   }
#  browser()
   tspv<-2*apply(ospv,2,min)
   return(list(ospv=ospv,tspv=tspv,ci=ci,mle=mle))
}

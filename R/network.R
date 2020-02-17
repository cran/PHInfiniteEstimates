#' This function enumerates conditional sample spaces associated with logistic regression,
#'
#' This function uses a network algorithm to enumerate conditional sample spaces associated with logistic regression,
#' using a minimal version of the algorithm of \insertCite{mehtapatel;textual}{PHInfiniteEstimates}.
#'
#' Examples apply this to data from \insertCite{mehtapatel;textual}{PHInfiniteEstimates}, 
#' citing \insertCite{goorinetal87;textual}{PHInfiniteEstimates}.
#'
#' @param dm matrix of covariates
#' @param n Vector of number of trials.  If null, make them all ones.
#' @param resp vector of successes.  Used only to calculate the sufficient statistics, unless
#' sufficient statistics are entered directly.  Either resp or sst must be provided.
#' @param conditionon indices of covariate matrix indicating sufficient statistics to be conditioned on.
#' @param sst sufficient statistic vector, if input directly.  Otherwise, recomputed from resp.
#' @param addint logical, true if a column of 1s must be added to the covariate matrix.
#' @param verbose logical; if true, print intermediate results.
#' @return For a successful run, a list with components:
#' \itemize{
#'   \item possible   matrix with vectors of possible unconditioned values of the sufficient statistic.
#'   \item count   count of entries in the conditional distribution.
#'   \item obsd   Observed value of unconditioned sufficient statistics.
#' }
#' For an unsuccessful run (because of input inconsistencies) NA
#' @export
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
network<-function(dm,n=NULL,resp=NULL,conditionon=NULL,sst=NULL,addint=TRUE,verbose=FALSE){
   out<-0#out will be used to track error conditons.  NA indicates errors; anything else indicates ok so far.
   dm<-as.matrix(dm)
   if(addint){
      dm<-cbind(1,dm)
      if(verbose) message("Added intercept to design matrix.\n")
   }
   if(!is.null(conditionon)&is.null(sst)){
      if(!is.null(resp)){
         sst<-(t(dm)%*%resp)
         if(verbose) message("Created SST",sst,"\n")
      }else{
         out<-NA
      }
   }
   if(!is.na(out)){
      dd<-dim(dm)[2]
      mm<-dim(dm)[1]
      if(is.null(n)) n<-rep(1,mm)
      possible<-array(0,c(1,dd))
      count<-c(1)
      for(i in seq(mm)){
         if(verbose) message("i",i," of ",mm,"\n")
         oldfound<-dim(possible)[1]
         np<-possible[rep(seq(dim(possible)[1]),n[i]+1),]
         nc<-c(count,rep(NA,length(count)*n[i]))
         for(y in 1:n[i]){
            for(j in seq(oldfound)){
               np[oldfound*y+j,]<- possible[j,]+as.matrix(dm[i,]*y)
               nc[oldfound*y+j]<-count[j]*choose(n[i],y)
            }
          }
          if(verbose) message("sum(nc>0)",sum(nc>0),"\n")
          if(!is.null(conditionon)){
             for(j in seq(length(conditionon))){
                bnds<- if(i<mm) c(0,sum(abs(dm[(i+1):mm,conditionon[j]])*n[(i+1):mm])) else c(0,0)
                if(verbose) message("bnds",bnds,"\n")
                for(k in seq(length(nc))){
                   if((np[k,conditionon[j]]+bnds[1])>sst[conditionon[j]]) nc[k]<-0
                   if((np[k,conditionon[j]]+bnds[2])<sst[conditionon[j]]) nc[k]<-0
                }
             }
          }
          if(verbose) message("sum(nc>0)",sum(nc>0),"\n")
          for(m in seq(length(nc)-1)){
             if(nc[m]>0) for(q in (m+1):length(nc)) if(nc[q]>0) {
                if(all(np[m,]==np[q,])){
                   nc[m]<-nc[q]+nc[m]
                   nc[q]<-0
                }
             }
          }
          possible<-np[nc>0,,drop=FALSE]
          count<-nc[nc>0]
          if(verbose) message("Length count",length(count),"\n")
      }
      possible<-possible[,-conditionon,drop=FALSE]
#     oo<-order(possible)
#     possible<-possible[oo,]
#     count<-count[oo]
      out<-list(possible=possible,count=count,obsd=sst[-conditionon])
   }
   return(out)
}

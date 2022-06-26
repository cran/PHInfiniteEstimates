#' Convert a proportional hazards regression to a multinomial regression.
#'
#' @param survobj A survival object, with potentially right censoring.
#' @param covmat a matrix of covariates.
#' @return a data set on which to apply conditional multinomial regression, corresponding to the proportional hazards regression analysis.
#' In order to run the line commented out below, you would need this:
#' # @importFrom mlogit mlogit.data
#' @export
#' @details
#' Implements version of \insertCite{kz19}{PHInfiniteEstimates}.  
#' The proportional hazards regression is converted to a multinomial regression logistic regression, and methods of \insertCite{kolassa16}{PHInfiniteEstimates} may be applied.
#' This function is intended to produce intermediate results to be passed to \code{convertmtol}, and then to \code{reduceLR} of \insertCite{kolassa97}{PHInfiniteEstimates}.  See examples in the \code{convertmtol} documentation.
#' @references 
#' \insertRef{kolassa97}{PHInfiniteEstimates}
#'
#' \insertRef{kolassa16}{PHInfiniteEstimates}
#'
#' \insertRef{kz19}{PHInfiniteEstimates}
convertstoml<-function(survobj,covmat){
   class(survobj)<-NULL
   out<-NULL
   count<-0
   id<-seq(dim(survobj)[1])
   for(tt in unique(sort(survobj[,1]))){
      count<-count+1
      if(any((survobj[,1]==tt)&(survobj[,2]==1))){
         atrisk<-(survobj[,1]>=tt)
         if(sum(atrisk)>1){
            new1<-data.frame(list(chid=as.character(rep(count,sum(atrisk))),
              patients=as.character(id[atrisk]),choice=rep(FALSE,sum(atrisk))),
              stringsAsFactors=FALSE)
            new1<-as.data.frame(new1)
            new1[(survobj[atrisk,2]==1)&(survobj[atrisk,1]==tt),3]<-TRUE
            out<-rbind(out,cbind(new1,covmat[atrisk,,drop=FALSE]))
         }
      }
   }
#  out<-mlogit.data(as.data.frame(out),choice="choice",id.var="chid",alt="alt")
   return(out)
}

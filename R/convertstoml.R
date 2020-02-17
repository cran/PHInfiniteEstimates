#' Convert a proportional hazards regression to a multinomial regression.
#'
#' @param survobj A survival object, with potentially right censoring.
#' @param covmat a matrix of covariates.
#' @return a data set on which to apply conditional multinomial regression, corresponding to the proportional hazards regression analysis.
#' @importFrom mlogit mlogit.data
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
   oo<-order(survobj[,1],1-survobj[,2])
   id<-seq(dim(survobj)[1])
   for(ii in id){
      if(survobj[oo[ii],2]==1){
         atrisk<-(survobj[,1]>=survobj[oo[ii],1])
         if(sum(atrisk)>1){
            new1<-data.frame(list(chid=as.character(rep(ii,sum(atrisk))),
              patients=as.character(id[atrisk]),choice=rep(FALSE,sum(atrisk))),
              stringsAsFactors=FALSE)
            new1<-as.data.frame(new1)
            new1[id[atrisk]==oo[ii],3]<-TRUE
            out<-rbind(out,cbind(new1,covmat[atrisk,,drop=FALSE]))
         }
      }
   }
#  out<-mlogit.data(as.data.frame(out),choice="choice",id.var="chid",alt="alt")
   return(out)
}

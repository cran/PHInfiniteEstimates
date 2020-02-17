reduceLR<-function(Z,nvec=NULL,yvec=NULL,keep,sst=NULL){
#' Reduce a logistic regression with monotone likelihood to a conditional regression with
#' double descending likelihood.
#'
#' @param Z regression matrix
#' @param nvec vector of sample sizes
#' @param yvec vector of responses
#' @param keep vector of variable names to block from consideration for removal.
#' @param sst vector of sufficient statistics
#' @return a list with components
#' \itemize{
#'   \item keepme indicators of which variables are retained in the reduced data set
#'   \item moderate indicatiors of which observations are retained in the reduced data set
#'   \item extreme indicators of which observations are removed in the reduced data set
#'   \item toosmall indicator of whether resulting data set is too small to fit the proportional hazards regression
#' }
#' @details
#' This function implements version of \insertCite{kolassa97;textual}{PHInfiniteEstimates}.  
#' It is intended for use with extensions to multinomial regression as in \insertCite{kolassa97;textual}{PHInfiniteEstimates} and to survival analysis as in \insertCite{kz19;textual}{PHInfiniteEstimates}.
#' The method involves linear optimization that is potentially repeated.  Initial calculations were done using a proprietary coding of the simplex, in a way that allowed for later iterations to be restarted from earlier iterations; this computational advantage is not employed here, in favor of computational tools in the public domain and included in the R package lpSolve.  
#' Furthermore, \insertCite{kolassa97;textual}{PHInfiniteEstimates} removed regressors that became linearly dependent using orthogonalization, but on further reflection this computation is unnecessary.
#' Data in the examples are from \insertCite{mehtapatel;textual}{PHInfiniteEstimates}, 
#' citing \insertCite{goorinetal87;textual}{PHInfiniteEstimates}.
#' @importFrom lpSolve lp
#' @export
#' @references 
#' \insertRef{mehtapatel}{PHInfiniteEstimates}
#'
#' \insertRef{goorinetal87}{PHInfiniteEstimates}
#'
#' \insertRef{kolassa97}{PHInfiniteEstimates}
#'
#' \insertRef{kolassa16}{PHInfiniteEstimates}
#'
#' \insertRef{kz19}{PHInfiniteEstimates}
#' @examples
#' #Cancer Data
#' Z<-cbind(rep(1,8),c(rep(0,4),rep(1,4)),rep(c(0,0,1,1),2),rep(c(0,1),4))
#' dimnames(Z)<-list(NULL,c("1","LI","SEX","AOP"))
#' nvec<-c(3,2,4,1,5,5,9,17); yvec<-c(3,2,4,1,5,3,5,6)
#' reduceLR(Z,nvec,yvec,c("SEX","AOP"))
#' #CD4, CD8 data
#' Z<-cbind(1,c(0,0,1,1,0,0,1,0),c(0,0,0,0,1,1,0,1),c(0,0,0,0,0,1,1,0),c(0,1,0,1,0,0,0,1))
#' dimnames(Z)<-list(NULL,c("1","CD41","CD42","CD81","CD82"))
#' nvec<-c(7,1,7,2,2,13,12,3); yvec<-c(4,1,2,2,0,0,4,1)
#' reduceLR(Z,nvec,yvec,"CD41")
#  message("Entering reduceLR")
   if(is.null(nvec)) nvec<-rep(1,dim(Z)[1])
   explore<-dimnames(Z)[[2]][-match(keep,dimnames(Z)[[2]])]
   moderate<-seq(length(nvec))
   extreme<-rep(0,length(nvec))
   oldkeepme<-rep(TRUE,dim(Z)[2])
#  message("Initial number of parameters",length(oldkeepme))
   names(oldkeepme)<-dimnames(Z)[[2]]
   done<-FALSE
   if(is.null(sst)) sst<-t(Z)%*%yvec
   if(any(is.na(sst))) message("NA in SST mark 1")
   origsst<-sst
   tempz<-Z
   toosmall<-FALSE
   if(class(try(solve(t(tempz)%*%tempz),silent=T))[1]!="try-error"){
      done<-FALSE
      toosmall<-FALSE
   }else{
      done<-TRUE
      toosmall<-NA
      oldkeepme<-NULL
      moderate<-NULL
      extreme<-NULL
   }
   while(!done){
      im<-solve(t(tempz)%*%tempz)%*%t(tempz)
      if(any(is.na(sst))){
#        browser()
         toosmall<-TRUE
         done<-TRUE
         message("NA appears in sst")
      }
      constraints<-t(sst)%*%im
      im2<-diag(rep(1,length(moderate)))-(tempz%*%im)
      constraints<-rbind(c(constraints-nvec[moderate],-constraints),
         cbind(im2,-im2),1)
# Set up a linear program with twice as many variables as there are model
# parameters.  Constraints correspond to logistic observations, and a final
# constraint bounding the sum of the indicators.
      xx<-lp("max",rep(1,2*length(moderate)),constraints,
         c(rep("=",length(moderate)+1),"<="), c(rep(0,length(moderate)+1),1))
#     browser()
      fo<-xx$solution[seq(length(moderate))]==0
      fz<-xx$solution[-seq(length(moderate))]==0
      extreme[moderate[!fo]]<-1
      extreme[moderate[!fz]]<--1
      done<-all(c(fz,fo))
      if(!done){
         atone<-moderate[!fo]
         moderate<-moderate[fo&fz]
#        browser()
#        message("Length oldkeepme",length(oldkeepme)," Lenth of sst",length(sst))
         keepme<-rep(FALSE,length(oldkeepme))
         names(keepme)<-names(oldkeepme)
         keepme[keep]<-TRUE
         if(kappa(t(Z[moderate,keep])%*%Z[moderate,keep])>1.0e10){
            toosmall<-TRUE
            done<-TRUE
#           message("Exitting because of poor condition")
         }
#        cat("explore",explore)
#        for(j in seq(length(explore))){
         for(j in explore){
            if(keepme[j]==FALSE){
               keepme[j]<-TRUE
               if(kappa(t(Z[moderate,keepme])%*%Z[moderate,keepme])>1.0e10){
#                 cat("keepme",keepme)
#                 cat("names(keepme)",names(keepme))
                  keepme[j]<-FALSE
#                 cat("Dropping variable",j,"\n")
               }
            }
         }
#        message("Length keepme",length(keepme))
         if(!toosmall){
            newsst<-if(length(atone)>0){origsst[keepme]-as.vector(t(Z[atone,keepme,drop=FALSE])%*%nvec[atone])}else{origsst[keepme]}
            if(any(is.na(newsst))){
                message("NA in SST mark 2; keepme",keepme)
#               browser()
            }
            sst<-newsst
            oldkeepme<-keepme
            tempz<-Z[moderate,oldkeepme]
         }else{
            done<-TRUE
         }
      }
   }
   return(list(keepme=oldkeepme,moderate=moderate,extreme=extreme,toosmall=toosmall))
}

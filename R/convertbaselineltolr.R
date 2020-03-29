#' Convert a baseline logit model data set, formatted in the long form as described in the documentation for mlogit.data from mlogit package, to a conditional logistic regression.
#'
#' @param dataset in formatted as in the output from mlogit.data of the mlogit packages
#' @param choice name of variable in dataset representing choice, a logical variable indicating whether this choice is actually chosen.
#' @param covs vector of names of covariates
#' @param strs name of variable in data set indicating independent subject
#' @param alt name of variable in data set indicating potential choice.
#' @return a data set on which to apply conditional logistic regression, corresponding to the baseline logit model.
#' @export
#' @details
#' This function implements version of \insertCite{kolassa16}{PHInfiniteEstimates}.
#' The multinomial regression is converted to a conditional logistic regression, and methods of \insertCite{kolassa97}{PHInfiniteEstimates} may be applied.  
#' This function differs from \code{convertmtol} of this package in that \code{convertmtol} treats a less-rich data structure, and this function treats the richer data structure that is an output of \code{mlogit.data} from package \code{mlogit}.
#' Data in the example is from \insertCite{bes;textual}{PHInfiniteEstimates}.
#' @references 
#' \insertRef{bes}{PHInfiniteEstimates}
#'
#' \insertRef{kolassa97}{PHInfiniteEstimates}
#'
#' \insertRef{kolassa16}{PHInfiniteEstimates}
#' @examples
#' data(voter.ml)
#' covs<-c("Labor","Liberal.Democrat","education")
#' #Fit the multinomial regression model, for comparison purposes.
#' ## Lines beginning ## give mlogit syntax that has been made obsolete.
#' #Add the index attribute to the data set, giving the index of choice made and the index of the 
#' #alternative, and a boolean variable giving choice.
#' ##attributes(voter.ml)$index<-voter.ml[,c("chid","alt")]
#' ##attributes(voter.ml)$choice<-"voter"
#' ##mlogit(voter~1|Labor+Liberal.Democrat+education,data=voter.ml)
#' mlogit(voter~1|Labor+Liberal.Democrat+education,data=voter.ml, 
#'    chid.var = "chid", alt.var = "alt")
#' #Convert to a data set allowing treatment as the equivalent conditional logistic regression.  
#' #This result will be processed using reduceLR of this package to give an equivalent conditional
#' # regression model avoiding infinite estimates.
#' out<-convertbaselineltolr(voter.ml,"voter",c("Labor","Liberal.Democrat","education"))
#' #Fit the associated unconditional logistic regression for comparison purposes.
#' glm(out[,"y"]~out[,1:75],family=binomial)
convertbaselineltolr<-function(dataset,choice,covs,strs="chid",alt="alt"){
   nobs<-dim(dataset)[1]
   mxstr<-0
   ustrs<- unique(dataset[[strs]])
   for(str in ustrs){
      mxstr<-max(mxstr,sum(!is.na(match(dataset[[strs]],str))))
   }
   altnames<-unique(dataset[[alt]])
   out<-array(0,c(nobs,length(ustrs)+(length(covs)+1)*(mxstr-1)+1))
   dimnames(out)<-list(NULL,c(paste("s",ustrs,sep=""),
      as.vector(t(outer(altnames[-1],c("Intercept",covs),paste,sep=":"))),"y"))
   for(str in ustrs){
      out[str==dataset[[strs]], paste("s",str,sep="") ]<-1
   }
   for(jj in seq(nobs)){
      for(ii in 2:length(altnames)){
         if(dataset[[alt]][jj]==altnames[ii]){
            out[jj,paste(altnames[ii],c("Intercept",covs),sep=":")]<-c(1,as.matrix(dataset[jj,covs]))
         }
      }
   }
   out[,"y"]<-0+dataset[,choice]
   return(out)
}

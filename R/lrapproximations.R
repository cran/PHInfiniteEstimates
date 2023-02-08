#' Assess the accuracy of the log rank statistic approximation to the true value, in the case without censoring.  Provides plots of statistics, and empirical test level.
#'
#' @param nobs number of observations in each group.  This currently supports only equal group size data sets.
#' @param ratio Ratio of group means; use 1 for null.
#' @param nsamp Monte Carlo sample size.
#' @return a vector of empirical test sizes.
#' @export
#' @importFrom graphics pairs
#' @importFrom nph logrank.test
#' @examples
#' lrapproximations(nsamp=100)
lrapproximations<-function(nobs=10,ratio=1,nsamp=1000){
   out<-array(NA,c(nsamp,4))
   myranksa<--(1+log(1-(1:(2*nobs))/(2*nobs+1)))
   myranks<-cumsum(1/(2*nobs+1-seq(2*nobs)))-1
#  message("myranks",myranks)
#Kolassa (2020) Eq 3.11
   m1<-nobs*mean(myranks)
   v1<-nobs*nobs*(mean(myranks^2)-mean(myranks)^2)/(2*nobs-1)
   mymarks<-c(rep(0,nobs),rep(1,nobs))
   for(jj in seq(dim(out)[1])){
      x<-rexp(nobs)
      y<-rexp(nobs)/ratio
      oo<-order(c(x,y))
      out[jj,1]<-(sum(myranks*mymarks[oo])-m1)/sqrt(v1)
      out[jj,2]<-(sum(myranksa*mymarks[oo])-m1)/sqrt(v1)
# Put the x and y into a single vector, and tack on an indicator of which group.
      mydata<-rbind(cbind(x,0),cbind(y,1))
# Order the data, from largest to smallest.
      mydata<-mydata[order(mydata[,1],decreasing=TRUE),]
# Calculate the vector Y of individuals at risk
      Y<-cumsum(rep(1,2*nobs))
# Calculate the vector Y1 of idividuals at risk in group 1.
      Y1<-cumsum(mydata[,2])
      v2<-sum((1*(Y-1)*(Y-Y1)*Y1/(Y^2*(Y-1)))[Y>1])
#     cat("v2",v2,"\n")
      out[jj,3]<-sum(mydata[,2]-Y1*1/Y)/sqrt(v2)
# From package nph
      out[jj,4]<-logrank.test(c(x,y),rep(1,2*nobs),mymarks)$test$z           
   }
   dimnames(out)<-list(NULL,c("True Log Rank Statistic", 
         "Approximation as Sums of Inverses",
         "Survival Analysis Approximation","Result from nph"))
   pairs(out[,],main="Two-Sample Rank Statistics",
      labels=dimnames(out)[[2]],
      sub=paste(c("Observations per group","Relative Risk","Number of Samples"),
         c(nobs,ratio,nsamp),sep="=",collapse=","))
   return(apply(out>=1.96,2,mean))
}

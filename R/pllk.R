#' Partial likelihood for proportional hazards
#'
#' @param beta parameter vector
#' @param xmat regression matrix
#' @param ind censoring indicator, 1 for event and any other value otherwise.
#' @param cc Continuity correction for sum of x vectors with multiple occurrences in risk set.  For binary covariates is half.  Default a vector of zeros.
#' @return a list with components
#' \itemize{
#'   \item d0 partial likelihood
#'   \item d1 first derivative vector
#'   \item d2 second derivative matrix
#' }
#' @export
#' @examples
#' #Uses data set breast from package coxphf.
#' data(breast)
#'  xmat<-as.matrix(breast)[order(breast$TIME),c("T","N")]
#'  ind<-breast$CENS[order(breast$TIME)]
#'  short<-coxph(Surv(TIME,CENS)~ T+ N,data=breast)
#'  pllk(as.vector(coef(short)),xmat,ind)
pllk<-function(beta,xmat,ind,cc=NULL){
    if(is.null(cc)) cc<-rep(0,length(beta))
    i<-0; d0<-0; d1<-cc;d2<-array(0,rep(dim(xmat)[2],2))
    ip<-xmat%*%beta; n<-length(ip)
    for(i in seq(n-1)){
       if(ind[i]==1) {
          we<-exp(ip[i:n])
          psum<-sum(we)
          d0<-d0+ip[i]-log(psum)
          meanvv<-as.vector(we%*%xmat[i:n,])/psum
          d1<-d1+xmat[i,]-meanvv
          wx<-apply(xmat[i:n,],2,"*",we)
          d2<-d2-t(xmat[i:n,])%*%wx/psum+outer(meanvv,meanvv)
       }
    }
    return(list(d0=d0,d1=d1,d2=d2))
}

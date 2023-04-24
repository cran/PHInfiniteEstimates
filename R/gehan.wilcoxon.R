#' Perform Gehan's application to the Wilcoxon test for multiple samples, testing for equivalance of survival curve.  See Klein and Moeschberger (1997) Survival Analysis (7.3.3) and pp. 193-194.
#'
#' @param myformula Proportional hazards formula appropriate for survfit
#' @param data the data set
#' @param gehan logical flag triggering the Wilcoxon test (gehan=TRUE), with weights equal to total at risk, or the log rank test (gehan=FALSE) with weights all 1.
#' @param plot logical flag triggering plotting.
#' @param alpha Nominal test level for plotting on graph
#' @param subset Apply to a subset of the data.
#' @return An htest-like object with the chi-square version of the test.
#' @examples
#' data(breast)#From package coxphf
#' gehan.wilcoxon.test(Surv(TIME,CENS)~G,data=breast)
#' @importFrom stats qnorm
#' @export
gehan.wilcoxon.test<-function(myformula,data,gehan=TRUE,plot=FALSE,alpha=0.05,subset=NULL){
   a<-survfit(myformula,data=data)
   cn<-cumsum(a$strata)
   ttt<-unique(sort(a$time))
   grps<-rep(NA,sum(a$strata))
   outmat<-array(NA,c(length(ttt),length(a$n),3))
   for(ii in seq(length(grps))){
      grps[ii]<-sum(ii>cn)+1
      myrow<-(ttt==a$time[ii])
      outmat[ttt==a$time[ii],grps[ii],]<-c(a$n.risk[ii],a$n.event[ii],a$n.censor[ii])
   }
   outmat[1,,1]<-a$n
   for(jj in seq(length(a$n))){
      if(is.na(outmat[1,jj,2])) outmat[1,jj,2]<-0
      if(is.na(outmat[1,jj,3])) outmat[1,jj,3]<-0
      for(ii in (2:(length(ttt)))) {
         if(is.na(outmat[ii,jj,2])) outmat[ii,jj,2]<-0
         if(is.na(outmat[ii,jj,3])) outmat[ii,jj,3]<-0
         if(is.na(outmat[ii,jj,1])) outmat[ii,jj,1]<-outmat[ii-1,jj,1]-sum(outmat[ii-1,jj,2:3])
      }
   }
   outmat<-outmat[apply(outmat[,,1],1,sum)>1,,]
   Ytot<-apply(outmat[,,1],1,sum)
   Dtot<-apply(outmat[,,2],1,sum)
   stat<-rep(NA,length(a$strata))
   v<-array(NA,rep(length(a$strata),2))
   w<-if(gehan) Ytot else rep(1,length(Ytot))
   for(jj in seq(length(a$n))){
      stat[jj]<-sum((outmat[,jj,2]-Dtot*outmat[,jj,1]/Ytot)*w)
      v[jj,jj]<-sum((outmat[,jj,1]*(Ytot-outmat[,jj,1])*(Ytot-Dtot)*Dtot/(Ytot^2*(Ytot-1)))*w^2)
      if(jj<length(a$n)) for(kk in (jj+1):length(a$n)){
         v[jj,kk]<--sum((outmat[,jj,1]*outmat[,kk,1]*(Ytot-Dtot)*Dtot/ (Ytot^2*(Ytot-1)))*w^2)
         v[kk,jj]<-v[jj,kk]
      }
   }
   
   yvec<-cumsum(c((outmat[,2,2]-Dtot*outmat[,2,1]/Ytot)*w,0))/sqrt(v[2,2])
   cv<--qnorm(alpha/2)
   if(plot){
      if(any(abs(yvec)>cv)){
         cutt<-min(ttt[abs(yvec)>cv])
         sgnme<-sign(yvec[ttt==cutt])
         flipt<-ttt[ttt>=cutt]
         flipy<-2*cv*sgnme-yvec[ttt>=cutt]
      }else{
          flipt<-NULL; flipy<-NULL
      }
      plot(range(c(0,ttt)),range(c(cv,-cv,flipy,yvec)),type="n", 
         xlab="Time",ylab="Statistic Value", 
         main=paste(if(gehan) "Gehan Wilcoxon" else "Log Rank", "and Reyni tests"),
         sub=paste(c("One-Sided Log Rank Level","One-Sided Reyni Level"),
             c(alpha/2,alpha),sep="=",collapse=","))
      lines(c(0,ttt),c(0,yvec),type="s")
      abline(h=cv,lty=2)
      abline(h=-cv,lty=2)
      lines(flipt,flipy,col=2)
   }
   sq<-(stat[-1,drop=FALSE]%*%solve(v[-1,-1,drop=FALSE],stat[-1,drop=FALSE]))[1,1]
   testout<-list(method=if(gehan) "Gehan-Wilcoxon" else "log rank",alternative="two-sided",
      statistic=sq,p.value=1-pchisq(sq,length(stat)-1),renyi=max(abs(yvec)),
      both=any((yvec>cv)&any(yvec<(-cv))))
   class(testout)<-"htest"
   return(testout)
}

#' Examine the potential role of treatment in treatment in 
#' a model already including sex.
#' Straight lines that are not 45 degrees
#' indicate the appropriateness of new variable as a linear
#' effect.
#' @param formulastring A formula for a coxph fit.
#' @param time The name of the time variable
#' @param stratifier The name of the stratifier variable
#' @param status The name of the status variable
#' @param mydata The data frame.
arjasplot<-function(formulastring,time,stratifier,status,mydata){
   fit <- coxph(as.formula(formulastring),mydata)
   bbb<-survfit(fit,newdata=mydata)$cumhaz
   new<-as.integer(as.factor(mydata[,stratifier]))
   mmm<-max(new)
   utime<-sort(unique(mydata[,time]))
   totdel<-totgmat<-array(0,c(mmm,length(utime)))
   for(jj in seq(dim(mydata)[1])){
      for(kk in seq(length(utime))){
         ll<-min(kk,sum(utime<=mydata[jj,time]))
         totgmat[new[jj],kk]<-totgmat[new[jj],kk]+bbb[ll,jj]
         totdel[new[jj],kk]<-totdel[new[jj],kk]+
            mydata[jj,status]*(mydata[jj,time]<=utime[kk])
      }
   }
   plot(range(totdel),range(totgmat),main="Arjas Plot",
      xlab="Number of failures in stratum",
         ylab="Estimated Cumulative Hazard",type="n")
   for(ll in seq(mmm)) lines(totdel[ll,],totgmat[ll,],lty=ll)
   legend(0,max(totgmat),legend=paste(stratifier,
      names(table(mydata[,stratifier]))),lty=seq(mmm))
#  return(invisible(out))
}

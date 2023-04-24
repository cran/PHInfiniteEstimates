#' Plot hazards for two strata for each time.  At times with an event in one but not the other group, the fitted hazard remains constant, and so the plot is a step function.  If hazards are proportional between strata, then the plot should be close to a straight line.
#' @param fit A coxph fit with a stratification term.
#' @export
#' @importFrom stats lm
#' @importFrom survival basehaz
andersenplot<-function(fit){
  out<-basehaz(fit)
  out<-out[order(out$time),]
  nstrat<-length(levels(out[,3]))
  out1<-rbind(0,cbind(out$time,array(NA,c(dim(out)[1],nstrat))))
  for(j in seq(dim(out)[1])){
     out1[j+1,as.numeric(out[j,3])+1]<-out[j,1]
     for(i in seq(nstrat)) if(is.na(out1[j+1,i+1]))
       out1[j+1,i+1]<-out1[j,i+1]
  }
  plot(out1[,2],out1[,3],xlab="Group 1 Hazard",
     ylab="Group 2 Hazard",main="Andersen Plot")
  abline(lm(out1[,3]~-1+out1[,2]))
}

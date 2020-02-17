#' Produce a graphical assessment of Monte Carlo experiment on fidelity of proportional hazards regression to the uniform ideal.
#'
#' This function draws a quantile plot for Monte Carlo assessments of fit to the corrected proportional hazards fit.
#' @param regnsimulation A matrix with six columns and as many rows as there MC samples.
#' @export
##' @importFrom stats update
#' @importFrom graphics lines
#' @importFrom stats dnorm pchisq
#' @return A list with components of consisting of simulated Wald p-values, likelihood ratio p-values, and corrected likelihood ratio p-values.
checkresults<-function(regnsimulation){
   hw<-regnsimulation[,"SRLRT"]
   hz<-regnsimulation[,"Est"]/regnsimulation[,"SE"]
   ospvb<-pnorm(hw)+(1/hw-1/hz)*dnorm(hw)
   flip<-ospvb<0
   flip[is.na(flip)]<-FALSE
   ospvb[flip]<-0
   flip<-ospvb>1
   flip[is.na(flip)]<-FALSE
   ospvb[flip]<-1
   tspv<-ospvb
   flip<-tspv>.5
   flip[is.na(flip)]<-FALSE
   tspv[flip]<-1-tspv[flip]
   tspv<-2*tspv
   uu<-sort(tspv)
   ss<-sort(regnsimulation[,1])
   tt<-sort(regnsimulation[,2])
   use<-seq(length(ss)/10)
   plot(use/length(ss),ss[use],type="l")
   use<-seq(length(tt)/10)
   lines(use/length(tt),tt[use],lty=2)
   use<-seq(length(uu)/10)
   lines(seq(use)/length(uu),uu[use],lty=3)
   legend(0,max(ss[use]), lty=1:3,
      legend=c("Wald Test","Likelihood Ratio Test","Saddlepoint Corrected"))
   return(list(waldpv=ss,lrpval=tt,lrcorrectedpval=uu))
}

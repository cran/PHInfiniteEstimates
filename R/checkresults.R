#' Produce a graphical assessment of Monte Carlo experiment on fidelity of proportional hazards regression to the uniform ideal.
#'
#' This function draws a quantile plot for Monte Carlo assessments of fit to the corrected proportional hazards fit.
#' @param regnsimulation A structure with a component out, matrix with columns representing definitions of p-values and as many rows as there MC samples.
#' @param frac Proportion for bottom of distribution to be assessed.
#' @export
#' @importFrom graphics lines.default plot.default
#' @importFrom stats dnorm pchisq
#' @return A list with components of consisting of simulated Wald p-values, likelihood ratio p-values, and corrected likelihood ratio p-values.
checkresults<-function(regnsimulation,frac=0.1){
# Block to return both one and two-sided Skovgaard p values, constrained to [0,1].
   out<-regnsimulation$out
   settings<-regnsimulation$settings
   firth<-FALSE; saddlepoint<-FALSE
   if(!firth){
      hw<-out[,"SRLRT"]
      hz<-out[,"Est"]/out[,"Skov SE"]
      if(saddlepoint){
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
      }
   }
## Skovgaard p-value complete.
   yyy<-xxx<-vector(mode="list",length=4)
   linenames<-c("Waldpv","LRTpv","Waldpvf", "LRTpvf")
   for(jj in 1:4){
      xxx[[jj]]<-sort(out[,linenames[jj]])
      use<-seq(length(xxx[[jj]]))<=(frac*length(xxx[[jj]]))
      yyy[[jj]]<-(seq(length(xxx[[jj]]))/length(xxx[[jj]]))[use]
      xxx[[jj]]<-xxx[[jj]][use]
#     message("Length xxx",length(xxx[[jj]]))
#     message("Length yyy",length(yyy[[jj]]))
   }
   my<-max(unlist(xxx))
   plot.default(c(0,frac),c(0,my),type="n",main="P-value comparisons",
      xlab="True quantile p value", ylab="Nominal p value",
      sub=paste(c("k","B","c","nsamp","beta[1]"),c(settings$k,settings$B,settings$c,settings$nsamp,settings$beta[1]),
         collapse=",",sep="="))
   abline(0,1)
   for(jj in 1:4) lines.default(yyy[[jj]],xxx[[jj]],lty=jj+1)
   legend(0,my, lty=1:5,
      legend=c("Target", "Wald Test","Likelihood Ratio Test",
         "Firth Wald Test","Firth Likelihood Ratio Test"))
}

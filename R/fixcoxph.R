#' Remove observations from a proportional hazards regression, and return the fit of the reduced model.
#'
#' This function implements the approximate conditional inferential approach of \insertCite{kz19;textual}{PHInfiniteEstimates} to proportional hazards regression.
#' @param randdat A list with at least the component y, representing the Surv() object.  I expect that this will be output from an initial non-convergent regression.
#' @param xxx a design matrix for the gregression.  I expect that this will be the $x component of the output from an initial non-convergent regression, run with x=TRUE .
#' @param iv name of the variable of interest, as a character string
#' @param verbose logical flag governing printing.
#' @return Fitted survival analysis regression parameter of class coxph, fitted form data set with observations forcing infinite estimation removed.
#' @export
#' @importFrom stats update
#' @references 
#' \insertRef{kz19}{PHInfiniteEstimates}
#' @examples
#' data(breast) # From library coxphf
#' bcfit<-coxph(Surv(TIME,CENS)~ T+ N+ G+ CD,data=breast,x=TRUE)
#'\donttest{
#' fixcoxph(bcfit,bcfit$x,"T",Surv(TIME,CENS)~ T+ N+ G+ CD)
#'}
#' data(HgFish) # From library NADA
#' HgFish$PropWetland<-HgFish$PctWetland/100
#' OneSpecies<-HgFish[HgFish$Species=="LargemouthBass",]
#' OneSpecies$SedMeHgH<-(OneSpecies$SedMeHg>median(OneSpecies$SedMeHg))*1
#' OneSpecies$SedAVSH<-(OneSpecies$SedAVS>median(OneSpecies$SedAVS))*1
#' OneSpecies$SedLOIH<-(OneSpecies$SedLOI>median(OneSpecies$SedLOI))*1
#' largemouthbassfit<-coxph(Surv(1/Hg,!HgCen)~SedAVSH+SedLOIH+SedMeHgH+
#'   SedMeHgH*SedAVSH+SedMeHgH*SedLOIH+PropWetland ,data=OneSpecies,x=TRUE)
#' \donttest{
#' fixcoxph(largemouthbassfit,largemouthbassfit$x,"PropWetland")
#' }
fixcoxph<-function(randdat,xxx,iv,verbose=FALSE){
   bad<-NULL
# Test first to make sure there is at least one event.  Otherwise return an error.
   if(sum(randdat$y[,2])>0){
      out<-convertstoml(randdat$y,xxx)
      if(is.null(out[,"chid"])){
         message("Why is chid component of out null?")
         browser()
      }
      out1<-convertmtol(out[,dimnames(xxx)[[2]]],out[,"chid"],out[,"choice"],out[,"patients"])
# Model here is glm(out1$y~out1$xmat-1,family=binomial)
      out2<-reduceLR(out1$xmat,yvec=out1$y,keep=iv)
      nv<-dim(randdat$x)[2]
      nstr<-length(out2$keepme)-nv
      keepme<-out2$keepme[nstr+seq(nv)]
      out2$toosmall<-FALSE
   }else{
      out2<-list(toosmall=TRUE)
   }
   if(!out2$toosmall){
      exclude<-as.list(as.data.frame(list(str=out$chid,obs=out$patients),stringsAsFactors=FALSE)[out2$extreme!=0,])
      larger<-bestbeta(randdat,exclude=exclude,touse=keepme)
      keepme[iv]<-FALSE
      smaller<-bestbeta(randdat,exclude=exclude,touse=keepme)
      larger$svar<-smaller$var
      larger$dropone<-smaller$loglik[2]
      larger$keepme<-keepme
#     browser()
      if(any(abs(larger$linear.predictor)>8)){
         message("Big linear predictor")
      }
   }else{
#     browser()
      nd<-dim(randdat$x)[2]
      no<-dim(randdat$x)[1]
      larger<-list(coefficients=rep(NA,nd),var=array(NA,c(nd,nd)),loglik=c(0,0),
         score=NA,iter=0,linear.predictors=rep(0,no),residuals=rep(0,no),
         dropone=smaller$loglik[2],keepme=keepme)
      class(larger)<-"coxph"
   }
#  if(verbose) message("Exiting fixcoxph")
   return(larger)
}

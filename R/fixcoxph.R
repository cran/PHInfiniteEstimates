#' Remove observations from a proportional hazards regression, and return the fit of the reduced model.
#'
#' This function implements the approximate conditional inferential approach of \insertCite{kz19;textual}{PHInfiniteEstimates} to proportional hazards regression.
#' @param randdat A list with at least the component y, representing the Surv() object.  I expect that this will be output from an initial non-convergent regression.
#' @param xxx a design matrix for the regression.  I expect that this will be the $x component of the output from an initial non-convergent regression, run with x=TRUE .
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
#'testdat2 <- data.frame(Time=c(4,3,1,1,2,2,3),
#'   Cen=c(1,1,1,0,0,0,0), Primary=c(0,2,1,1,1,0,0), Sex=c(0,0,0,0,1,1,1))
#'(bcfit<-coxph(Surv(Time,Cen)~Primary + Sex, testdat2, x=TRUE, ties="breslow"))
#'fixcoxph(bcfit,bcfit$x,"Primary")
fixcoxph<-function(randdat,xxx,iv,verbose=FALSE){
#  message("Mark a")
   bad<-NULL
   nnn<-dim(xxx)[1]
# Test first to make sure there are at least two events.  Otherwise return an error.
   if(sum(randdat$y[,2])>1){
      out<-convertstoml(randdat$y,xxx)
      if(is.null(out[,"chid"])){
         message("Why is chid component of out null?")
         browser()
      }
      out1<-convertmtol(out[,dimnames(xxx)[[2]]],out[,"chid"],out[,"choice"],
         out[,"patients"])
# Model here is glm(out1$y~out1$xmat-1,family=binomial)
      out2<-reduceLR(out1$xmat,yvec=out1$y,keep=iv)
      smallerxmat<-out1$xmat[,-match(iv,dimnames(out1$xmat)[[2]])]
      out3<-reduceLR(smallerxmat,yvec=out1$y,keep=NULL)
      e3<-length(out3$extreme)
      if(e3>0) e3<-if(any(is.na(out3$extreme))){-1}else{sum(abs(out3$extreme))}
      e2<-length(out2$extreme)
      if(e2>0) e2<-if(any(is.na(out2$extreme))){-1}else{sum(abs(out2$extreme))}
#     message("Before checking for infinite estimate")
      ivinf<-(e3!=e2)*sign(newllk(rep(0,dim(randdat$x)[2]),randdat)$d1[iv,1])
      if(is.na(ivinf)){
         message("ivinf na")
         browser()
      }
#     message("After checking for infinite estimate ivinf",ivinf)
#     message("sum(abs(out3$extreme))",sum(abs(out3$extreme)))
#     message("ivinf",ivinf)
      nv<-dim(randdat$x)[2]
      nstr<-length(out2$keepme)-nv
      keepme<-out2$keepme[nstr+seq(nv)]
      names(keepme)<-dimnames(randdat$x)[[2]]
      out2$toosmall<-FALSE
   }else{
      out2<-list(toosmall=TRUE)
   }
   if(!out2$toosmall){
      exclude<-as.list(as.data.frame(list(str=out$chid,obs=out$patients),stringsAsFactors=FALSE)[out2$extreme!=0,])
      if(all(!keepme)) message("Why is larger model empty?")
#     message("Before bestbeta")
      larger<-bestbeta(randdat,exclude=exclude,touse=keepme)
#     message("larger");print.default(larger)
      if(ivinf!=0) {
         larger$coefficients[iv] <- ivinf*Inf
         larger$var[iv,iv] <- 0
         larger$wald.test <-Inf
      }
      smallkeep<-keepme
      smallkeep[iv]<-FALSE
#     message("Checking contents of smaller model")
      if(all(!smallkeep)){
#        message("Nothing in smaller model")
         larger$dropone<-larger$loglik[1]
         larger$svar<-NULL
         larger$nssmall<-0
      }else{
         smaller<-bestbeta(randdat,exclude=exclude,touse=smallkeep)
         larger$dropone<-smaller$loglik[2]
         larger$svar<-smaller$var
         larger$nssmall<-sum(smallkeep)
      }
#     message("Change in log likelihood is ",larger$loglik[1]-larger$dropone)
      larger$keepme<-keepme
#     browser()
#     if(any(abs(larger$linear.predictor)>8)){
#        message("Big linear predictor")
#     }
#     message("Mark b, diag(larger$var)",paste(diag(larger$var),collapse=" "))
   }else{#Too small branch
      message("In too small branch")
#     browser()
      nd<-dim(randdat$x)[2]
      no<-dim(randdat$x)[1]
      keepme<-rep(FALSE,nd)
      names(keepme)<-dimnames(randdat$x)[[2]]
      message("End branch before bestbeta, nd and no",nd,no)
      vari<-array(NA,c(nd,nd))
      dimnames(vari)<-list(dimnames(randdat$x)[[2]],dimnames(randdat$x)[[2]])
      larger<-list(coefficients=rep(NA,nd),var=vari,loglik=c(0,0),
         score=NA,iter=0,linear.predictors=rep(0,no),residuals=rep(0,no),
         dropone=0,keepme=keepme,nssmall=0)
      class(larger)<-"coxph"
   }
#  if(verbose) message("Exiting fixcoxph")
   larger$n<-nnn
   return(larger)
}

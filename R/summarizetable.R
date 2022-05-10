#' Summarize the results of simulations investigating operating conditions for the data reduction method to avoid monotone likelihood.  Files are of form "hsxxx", for xxx numerals.
#' @importFrom stats quantile pchisq
#' @export
summarizetable<-function(){
   fl<-system("ls hsout*",intern=TRUE)
   biglist<-vector(mode="list",length=length(fl))
   names(biglist)<-fl
   for(fn in fl){source(fn); biglist[[fn]]<-eval(parse(text=fn))}
   settings<-array(NA,c(length(fl),5)) 
   dimnames(settings)<-list(fl,c("bb","nobs","k","c","B"))
   for(fn in fl){
      for(sn in dimnames(settings)[[2]][-1]){
         settings[fn,sn]<-biglist[[fn]]$settings[[sn]]
      }
      settings[fn,"bb"]<-exp(biglist[[fn]]$settings$beta[1])
   }
   tt<-table(as.data.frame(settings))
   temp<-dimnames(tt)
   temp$compare<-c("Proposed","HS")
   temp$assess<-c("Estimate Median Bias","Wald CV","LR CV","Wald Test Level",
     "LR Test Level")
   out<-array(NA,c(dim(tt),length(temp$compare),length(temp$assess)))
   out2<-array(NA,c(dim(tt),2))
   dimnames(out)<-temp
   dimnames(out2)<-temp[seq(length(temp)-1)]
   dimnames(out2)[[length(dim(tt))+1]]<-c("Infinite","BadFirth")
   index<-array(NA,c(1,length(dim(out))))
   index1<-array(NA,c(1,length(dim(out2))))
   dimnames(index)<-list(NULL,names(dimnames(out)))
   dimnames(index1)<-list(NULL,names(dimnames(out2)))
   nnn<-dimnames(out)[[length(index)]]
   mmm<-dimnames(out)[[length(index)-1]]
   for(fn in fl){
      for(kk in seq(length(index)-2)) index[1,kk]<-match(settings[fn,kk],
         dimnames(out)[[kk]])
#  Bad
      index1<-index[,seq(length(index)-1),drop=FALSE]
      index1[1,length(index)-1]<-1
#     browser()
      out2[index1]<-mean(biglist[[fn]]$out[,"npar"]<4)
      index1[1,length(index)-1]<-2
      out2[index1]<-mean(is.na(biglist[[fn]]$out[,"npar"]))
      print("Final"); print(out2[index1])
#Median bias
      target<-log(as.numeric(settings[fn,"bb"]))
#  HS
      index[1,length(index)-1]<-match("HS",mmm)
      index[1,length(index)]<-match("Estimate Median Bias",nnn)
      out[index]<-median(biglist[[fn]]$out[,"Estf"])-target
#  Proposed
      index[1,length(index)-1]<-match("Proposed",mmm)
      index[1,length(index)]<-match("Estimate Median Bias",nnn)
      out[index]<-median(biglist[[fn]]$out[,"Est"])-target
#Wald test level/power
#  HS
      index[1,length(index)-1]<-match("HS",mmm)
      index[1,length(index)]<-match("Wald Test Level",nnn)
      indexb<-index; indexb[1,"bb"]<-1; 
      pp<-if(settings[fn,"bb"]==1){0.05}else{out[indexb]}
      out[index]<-mean(biglist[[fn]]$out[,"Waldpvf"]<=pp  ,na.rm=TRUE)
#  Proposed
      index[1,length(index)-1]<-match("Proposed",mmm)
      index[1,length(index)]<-match("Wald Test Level",nnn)
      indexb<-index; indexb[1,"bb"]<-1; 
      pp<-if(settings[fn,"bb"]==1){0.05}else{out[indexb]}
      out[index]<-mean(biglist[[fn]]$out[,"Waldpv"]<=pp  ,na.rm=TRUE)
#LR test level/power
#  HS
      index[1,length(index)-1]<-match("HS",mmm)
      index[1,length(index)]<-match("LR Test Level",nnn)
      indexb<-index; indexb[1,"bb"]<-1; 
      pp<-if(settings[fn,"bb"]==1){0.05}else{out[indexb]}
      out[index]<-mean(biglist[[fn]]$out[,"LRTpvf"]<=pp  ,na.rm=TRUE)
#  Proposed
      index[1,length(index)-1]<-match("Proposed",mmm)
      index[1,length(index)]<-match("LR Test Level",nnn)
      indexb<-index; indexb[1,"bb"]<-1; 
      pp<-if(settings[fn,"bb"]==1){0.05}else{out[indexb]}
      out[index]<-mean(biglist[[fn]]$out[,"LRTpv"]<=pp  ,na.rm=TRUE)
#Wald test CV
#  HS
      index[1,length(index)-1]<-match("HS",mmm)
      index[1,length(index)]<-match("Wald CV",nnn)
      out[index]<-quantile(biglist[[fn]]$out[,"Waldpvf"],p=0.05,na.rm=TRUE)
#  Proposed
      index[1,length(index)-1]<-match("Proposed",mmm)
      index[1,length(index)]<-match("Wald CV",nnn)
      out[index]<-quantile(biglist[[fn]]$out[,"Waldpv"],p=0.05,na.rm=TRUE)
#LR test CV
#  HS
      index[1,length(index)-1]<-match("HS",mmm)
      index[1,length(index)]<-match("LR CV",nnn)
      out[index]<-quantile(biglist[[fn]]$out[,"LRTpvf"],p=0.05,na.rm=TRUE)
#  Proposed
      index[1,length(index)-1]<-match("Proposed",mmm)
      index[1,length(index)]<-match("LR CV",nnn)
      out[index]<-quantile(biglist[[fn]]$out[,"LRTpv"],p=0.05,na.rm=TRUE)
   }
   outa<-out[ ,1, , , , , -grep("CV",dimnames(out)[[7]]),drop=FALSE]
   temp<-dimnames(outa)
   for(jj in seq(length(temp))) temp[[jj]]<-paste(names(temp)[jj],temp[[jj]],sep="=")
   dimnames(outa)<-temp
   outb<-outa
   for(s0 in dimnames(outa)[["assess"]]){
   for(s1 in dimnames(outa)[["bb"]]) for(s2 in dimnames(outa)[["c"]]) 
      for(s3 in dimnames(outa)[["B"]]){
         if(s0=="Estimate Median Bias"){
            testme<-sign(diff(outa[s1,1,1,s2,s3,,s0]))
         }else{
            if(s1=="1"){
               testme<-sign(diff(abs(outa[s1,1,1,s2,s3,,s0]-0.05)))
            }else{
               testme<-sign(diff(outa[s1,1,1,s2,s3,,s0]))
            }
         }
         outb[s1,1,1,s2,s3,,s0]<-c(testme,-testme)
      }
   }
   outb[is.na(outb)]<-0
   outc<-array(paste("xx",format(outa,digits=2),"yy",outb,"zz",sep=""),
      dim(outa),dimnames=dimnames(outa))
   bad<-NA
   return(list(outc=outc,out2=out2))
#  library(R.utils) ; library(xtable)
#  out<-summarizetable()
#  print(xtable(wrap(out$outc[,1,1,,,,2:3],map=list(c(4,5,1),c(2,3)))),file="Survival/Tex/tabout.tex")
#  print(xtable(wrap(out$outc[,1,1,,,,1],map=list(c(4,1),c(2,3)))),file="Survival/Tex/tabouta.tex")
#  system("sed -i -e 's/xx/\\\\cfmacro{/g' -e 's/yy/}{/g' -e 's/zz/}/g' -e 's/compare=//' -e 's/.assess=/ /' -e 's/.bb=/ R=/g' -e 's/Test//' Survival/Tex/tabout.tex; sed -e 's/Level R=2/Power R=2/' -e 's/Level R=4/Power R=4/' -i Survival/Tex/tabout.tex")
#  system("sed -i -e 's/xx/\\\\cfmacro{/g' -e 's/yy/}{/g' -e 's/zz/}/g' -e 's/compare=//' -e 's/.assess=/ /' -e 's/.bb=/ R=/g' -e 's/Test//' Survival/Tex/tabouta.tex; sed -e 's/Level R=2/Power R=2/' -e 's/Level R=4/Power R=4/' -i Survival/Tex/tabouta.tex")
}

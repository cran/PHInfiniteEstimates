#' Draw diagram for toy PH example.
#'
#' @importFrom graphics axis legend plot points segments
#' @importFrom stats as.formula coefficients optimize rweibull uniroot
#' @return nothing.
#' @export
drawdiagram<-function(){
   plot(c(0,1),c(1,5),type="n",yaxt="n",xlab="Time",xaxt="n",ylab="Group")
   ttt<-c(.1,.15,.4,.7,.95)
   segments(c(0,0,0,0,0),1:5,ttt                ,1:5)
   axis(2,at=1:5,tick=TRUE,labels=c("A","B","A","B","A"))
   points(ttt,1:5,pch=5-c(1,1,4,4,1))
   legend(.5,2,pch=c(4,1),legend=c("Event","Censored"))
}

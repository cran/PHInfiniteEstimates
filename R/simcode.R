#' Simulate Weibull survival data from a model, perform reduction to remove infinite estimates, and calculate p values.
#'
#' Operating characteristics for the approximate conditional inferential approach to proportional hazards.
#'
#' This function is intended to verify the operating characteristics of the approximate conditional inferential approach of \insertCite{kz19;textual}{PHInfiniteEstimates} to proportional hazards regression.  A Weibull regression model, corresponding to the proportional hazards regression model, is fit to the data, and new data sets are simulated from this model.  P-values are calculated for these new data sets, and their empirical distribution is compared to the theoretical uniform distribution.
#' @param dataset the data set to use
#' @param myformula the formula for the Cox regression
#' @param iv name of the variable of interest, as a character string
#' @param ctime fixed censoring time
#' @param nsamp number of samples.
#' @param add preliminary results, if any.
#' @param nobs number of observations in target models, if different from that of dataset.
#' @param half logical flag triggering a less extreme simulation by dividing the Weibull regression parameters in half.
#' @param verbose logical flag triggering intermediate messaging.
#' @return a list with components
#' \itemize{
#'   \item out matrix with columns corresponding to p-values.
#'   \item seed random seed
#'   \item bad unused.
#'   \item srreg parametric lifetime regression
#' }
#' @importFrom survival Surv coxph survreg
#' @importFrom stats pnorm
#' @importFrom coxphf coxphf
#' @export
#' @references 
#' \insertRef{kz19}{PHInfiniteEstimates}
#' @examples
#' data(breast)
#'\donttest{
#' breasttestp<-simcode(breast,Surv(TIME,CENS)~ T+ N+ G+ CD,"T",72,nsamp=100,verbose=TRUE)
#'}
simcode<-function (dataset, myformula, iv, ctime, nsamp = 10000, add = NULL, 
    nobs = NA, half = FALSE, verbose = FALSE) {   
 bad <- NULL
    gg <- as.formula(paste("Surv(newt,nc)~", deparse(myformula[[3]])))
    sr1 <- survreg(formula = myformula, data = dataset, x = TRUE, y = TRUE)
    if (is.na(nobs)) nobs <- dim(sr1$x)[1]
    use <- seq(nobs)
    csr1 <- coefficients(sr1)
    ivl <- (names(csr1) == iv)[-1]
    if (half) csr1 <- csr1/2
#Replace the para of interest as 0
    csr1[iv] <- 0
    randdat <- cbind(as.data.frame(list(newt = rep(NA, nobs))), sr1$x[use, ])
    start<-0
    if (is.null(add)) {
#       set.seed(185928396)
#       set.seed(194837485)
        set.seed(202043125)
    }   else {
        outout <- rbind(add$out, array(NA, c(nsamp, dim(add$out)[2])))
        start <- dim(add$out)[1]
        set.seed(add$seed)
    }
    d1 <- Sys.time()
    for (kk in seq(nsamp)) {
        if (verbose) {
            d2 <- Sys.time()
            if (verbose) 
                message("kk=", kk, " of ", nsamp, 
                  ".  Completion time ", (d2 - d1) * (nsamp - 
                    kk)/kk + d2)
        }
        if (nobs < dim(sr1$x)[1]) 
            use <- sample(dim(sr1$x)[1], nobs)
        for (j in seq(nobs)) {
            elpred <- exp(sr1$x[use, ] %*% csr1)
            randdat$newt[j] <- rweibull(1, shape = 1/sr1$scale, 
                scale = elpred[j])
        }
        randdat$nc <- randdat$newt < ctime
        randdat$newt[randdat$newt >= ctime] <- ctime
        randdat$y <- Surv(randdat$newt, randdat$nc)
        xxx <- sr1$x[use, -1]
        randdat$x <- xxx
        repairedfit <- fixcoxph(randdat, xxx, iv, verbose = verbose)
        penalizedout <- coxphf(gg, randdat, maxit = 200, maxstep = 0.1)
        hh<-update(as.formula(gg),paste("~. -",iv))
        penalizedoutsmaller <- coxphf(hh, randdat, maxit = 200, maxstep = 0.1)
        myout<-summarizefits(repairedfit,penalizedout,penalizedoutsmaller,iv)
        if((start+kk)==1){
           outout <- array(NA, c(nsamp, length(myout)))
           dimnames(outout)<-list(NULL,names(myout))
        }
        outout[start + kk,]<-myout
#       if (verbose) message("outout",outout[start + kk,])
#       if (is.na(outout[start + kk])) {
#           dump("randdat")
#       }
    }
    return(list(out = outout, seed = .Random.seed, bad = bad, srreg = sr1))
}

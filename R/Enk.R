###############################################################################
## package 'secr'
## Enk.R
## 2022-11-18 cf pdot()
## 2023-04-28 improve nkpointcpp (allows multi), block unsupported detector types
## 2023-04-28 simulation option: average nrepl simulations
## 2023-05-13 nk()
###############################################################################

nk <- function(capthist) {
    if (ms(capthist)) {
        lapply(capthist, nk)
    }
    else {
        apply(apply(abs(capthist),c(1,3),sum)>0, 2, sum)
    }
}

Enk <- function (D, mask, traps, detectfn = 0, 
    detectpar = list(g0 = 0.2, sigma = 25, z = 1),
    noccasions = NULL, binomN = NULL, userdist = NULL, 
    ncores = NULL, nrepl = NULL) {
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)
    if ((detectfn > 9) & (detectfn<14) & is.null(detectpar$cutval))
        stop ("requires 'cutval' for detectfn 10:13")
    if (ms(traps) || ms(mask))
        stop ("requires single-session traps and mask")
    
    truncate <- ifelse(is.null(detectpar$truncate), 1e+10, detectpar$truncate)
    
    detectpars <- unlist(detectpar[parnames(detectfn)])
    if ((detectfn>9) & (detectfn<14))  detectpars <- c(detectpars, detectpar$cutval)
    if (length(detectpars)<3) detectpars <- c(detectpars,0)
    miscparm <- numeric(4);   ## dummy
    
    if (!is.null(usage(traps))) {
        usge <- usage(traps)
        if (is.null(noccasions)) {
            noccasions <- ncol(usage(traps))
        }
        else {
            if (noccasions < ncol(usage(traps))) {
                warning ("specified noccasions less than ncol of usage matrix")
            }
            if (noccasions > ncol(usage(traps)))
                stop ("specified noccasions exceeds ncol of usage matrix")
        }
    }
    else {
        if (is.null(noccasions))
            stop("must specify noccasions when traps does not have usage attribute")
        usge <- matrix(1, ndetector(traps), noccasions)
    }
    dettype <- detectorcode(traps, noccasions = noccasions)
    binomN <- getbinomN (binomN, detector(traps))
    markocc <- markocc(traps)
    
    if (is.null(markocc)) markocc <- rep(1,noccasions)
    if (!inherits(mask, 'mask')) {
        stop("mask input should be a mask object")
    }
    if (!is.null(nrepl)) {
        # simulation option
        onesimnk <- function (r) {
            pop <- sim.popn(D, core = mask, model2D = 'IHP')
            ch <- sim.capthist(
                traps      = traps, 
                popn       = pop, 
                detectfn   = detectfn, 
                detectpar  = detectpar, 
                noccasions = noccasions, 
                nsessions  = 1, 
                binomN     = binomN)
            nk(ch)  # individuals per detector
        }
        out <- sapply(1:nrepl, onesimnk)
        apply(out,1,mean)
    }
    else if (!all(detector(traps) %in% c('multi','proximity','count'))) {
        stop("Enk formula available only for multi, proximity and count detectors; try simulation")
    }
    else {
        D <- rep(D * getcellsize(mask), length.out = nrow(mask))  # per cell; includes linear
        distmat2 <- getuserdist (traps, mask, userdist, sessnum = NA, NULL, NULL, miscparm, detectfn == 20)
        ncores <- setNumThreads(ncores)
        grain <- if (ncores==1) 0 else 1
        nkpointcpp(
            as.double  (D),
            as.matrix  (distmat2),
            as.integer (dettype),
            as.matrix  (usge),
            as.integer (markocc),
            as.integer (detectfn),
            as.double  (detectpars),
            as.double  (miscparm),
            as.double  (truncate^2),
            as.integer (expandbinomN(binomN, dettype)),
            as.integer (grain),
            as.integer (ncores)
        )
    }
}
############################################################################################


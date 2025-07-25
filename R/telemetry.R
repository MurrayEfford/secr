###############################################################################
## package 'secr'
## telemetry.R
## fitting detection function to telemetry data
## 2012-10-19,20,21,22 2013-06-07
## 2012-06-15 read.telemetry
## 2013-11-15 telemloglik renamed telemetry.LC
## 2013-11-15 telemetryloglik renamed telemetry.LT
## 2013-11-16 telemetry.LCmask for smuggling telemetry centres into C code

## 2017-01-07 REDUNDANT?
###############################################################################

distances <- function (X, Y) {
    ## X and Y are 2-column matrices of coordinates
    onerow <- function (xy) {
        d <- function(xy2) {
            sqrt(sum((xy2 - xy)^2))
        }
        apply(Y, 1, d)
    }
    t(apply(X, 1, onerow))
}
###############################################################################

telemetry.LT <- function(CH, detectfn, realparval, PIA,
                         nmix = 1, knownclass = 1, uppersigma = 20) {
    ## likelihood component L_T for locations of telemetered animals
    ## no spatial covariates in this version
    ## does not use habitat mask, or truncate at boundary
    normalize <- function (parm) {
        ## Why normalize?
        par3 <- length(secr_parnames(detectfn)) == 3
        sigma <- parm['sigma']
        z <- if (par3) parm[secr_parnames(detectfn)[3]] else 1
        if (detectfn %in% c(14,16))
            lambda0 <- 10000 / (sigma^2 * 2 * pi)
        else {
            rdfn <- function (r, pars) r * dfn(r, pars, 0)
            lambda0 <- 10000/integrate (rdfn, 0, sigma * uppersigma, pars=c(1, sigma, z))$value
        }
        out <- c(lambda0=lambda0, sigma=sigma)
        if (par3) out <- c(out,z=z)
        out
    }
    Lx <- function(x) {
        ## return Pr(obs) given animal belongs to class x
        pind <- as.numeric(PIA[cbind(df$i, df$s, df$k, rep(x,n))])
        param <- Nrealparval[pind,, drop = FALSE]
        sapply(1:n, function(i) dfn(df$d[i], param[i,], 0))
    }
    px <- function(x) {
        ## return naive Pr(animal in class x) (i.e. ignoring knownclass?)
        indexmatrix <- cbind(df$i, df$s, df$k, rep(x,n))
        pind <- as.numeric(PIA[indexmatrix])

        realparval[pind,'pmix', drop = FALSE]
    }

    detectfn <- secr_valid.detectfn(detectfn, 14:19)
    dfn <- secr_getdfn (detectfn)
    Nrealparval <- t(apply(realparval, 1, normalize))  ## adjust lambda0, complete pars

    ## df with one row per detection
    ## i is integer row number within CH
    df <- data.frame(xy(CH), 
        i = animalID(CH, names = FALSE, sortorder = 'snk'), 
        s = occasion(CH, sortorder = 'snk'),
        k = trap(CH, names = FALSE, sortorder = 'snk')
    )
    df <- df[order(df$i, df$s, df$k),]              ## sort by animal, occasion, detector
    ni <- table(df$i)                               ## frequency of ith animal
    df$cx <- rep(tapply(df$x, df$i, mean),ni)       ## centre x
    df$cy <- rep(tapply(df$y, df$i, mean),ni)       ## centre y
    df$d <- sqrt((df$x-df$cx)^2 + (df$y-df$cy)^2)   ## distance from centre
    n <- nrow(df)                                   ## total number of detections
    df$knownclass <-  if (length(knownclass) == 1)  ## repeat knownclass for each detection
        rep(knownclass, n) else rep(knownclass, ni)

    L <- sapply(1:nmix, Lx)                         ## for each class, Pr(obs)
                                                 ## n x nmix matrix
    if (nmix>1) {
        pmix <- sapply(1:nmix, px)                  ## should be n x nmix matrix
        pmix[df$knownclass > 1,] <- 0
        pmix[cbind(1:n,df$knownclass-1)] <- 1       ## relies on zero index not changing
        L <- apply(L * pmix, 1, sum)
    }
    ## likelihood component for hcov binomial proportion
    if (any (knownclass != 1)) {
        pmix <- realparval[,'pmix']
        obsmix <- tabulate(knownclass)[2:(length(pmix)+1)]
        Lknown <- sum(obsmix * log(pmix))
    }
    else Lknown <- 0

    if (any(L<=0))
        list(value = -1e10, resultcode = 9)
    else
        list(value = sum(log(L)) + Lknown, resultcode = 0)
}

###############################################################################

telemetry.LC <- function(CH, detectfn, detectpar, mask, bvn = TRUE) {
    ## OBSOLETE

    ## likelihood component L_C for detection histories of telemetered animals only
    ## (including possible zero CH)
    ## experimental - doesn't allow covariates or work with traps or count detectors
    if (! all(detector(traps(CH)) %in% c('proximity')))
        stop ("requires proximity CH")
    n <- dim(CH)[1]
    J <- dim(CH)[2]
    traps <- traps(CH)
    K <- secr_ndetector(traps)
    detectfn <- secr_valid.detectfn(detectfn, 14:19)
    detectpar <- detectpar[c('lambda0','sigma','z')] ## ensure order correct
    g <- secr_getdfn (detectfn)
    xylist <- telemetryxy(CH)
    centrexy <- t(sapply(xylist, apply, 2, mean))

    ## animal x mask
    dfn <- function (xy) {
        ## xy is matrix of telemetry coordinates for one animal
        centres <- matrix(apply(xy, 2, mean), ncol = 2)
        if (bvn) {
            vcv <- var(xy)/nrow(xy)
            detS <- det(vcv)^0.5  ## sqrt(generalised variance)
            tempmask <- sweep (mask, MARGIN = 2, STATS = centres, FUN = '-')
            close <- apply(tempmask,1, function(x) sum(x^2)) < (detS*30)
            if (sum(close)<1) {
                close <- nearesttrap(matrix(0, ncol = 2, nrow = 1), tempmask)
            }
            tempmask <- tempmask[close,, drop = FALSE]
            # if (!requireNamespace(mvtnorm, quietly = TRUE)) stop()
            # dbvn <- apply(tempmask, 1, mvtnorm::dmvnorm, mean=c(0,0), sigma=vcv, log=FALSE)
            # mymvn is faster and produces identical result
            invS <- solve(vcv)
            mymvn <- function(XY) exp(-(XY %*% invS %*% XY)/2) / 2/pi/det(vcv)
            dbvn <- apply(tempmask,1,mymvn)
            dbvn <- dbvn / sum(dbvn)
            list(dbvn=dbvn, mask = mask[close,, drop=FALSE])
        }
        else {
            list(dbvn = 1.0, mask = centres)
        }
    }
    pmask <-  lapply(xylist, dfn)
    loglik <- 0
    for (id in names(xylist)) {
        m <- nrow(pmask[[id]]$mask)
        ## dtrap <- distances (traps, pmask[[id]]$mask)
        dtrap <- edist (traps, pmask[[id]]$mask)
        if (m==1) dtrap <- t(dtrap)  ## kludge to deal with vector
        gkm <- g(dtrap, unlist(detectpar), 0)   # 3rd arg is dummy for cutval
        prw <- matrix(1, nrow = K, ncol = m)
         for (j in 1:J) {
             wij <- CH[id,j,]  ## detection sites on occasion j, assume 0/1
             prw <- prw *
                 (sweep(gkm, STATS=wij, MARGIN=1, FUN = '*') +     # detected
                 sweep(1-gkm, STATS=1-wij, MARGIN=1, FUN = '*'))   # not detected
         }
        prwm <- apply(prw,2,prod)            ## product over detectors
        prwm <- prwm * pmask[[id]]$dbvn      ## m-length vectors
        loglik <- loglik + log(sum(prwm))    ## sum prob over mask_i points
    }
    loglik
}
###############################################################################
telemetry.LCmask <- function(CH, mask, bvn = TRUE) {
## the pdf the sampling distribution is normalized with the sum of on-mask points
## and implicitly assumed to be on mask
    dfn <- function (xy) {
        f <- numeric(mm)
        if (nrow(xy) == 0)
            f[] <- 1/mm
        else {
            ## xy is matrix of telemetry coordinates for one animal
            centres <- matrix(apply(xy, 2, mean), ncol = 2)
            if (bvn) {
                vcv <- var(xy)/nrow(xy)
                detS <- det(vcv)^0.5  ## sqrt(generalised variance)
                tempmask <- sweep (mask, MARGIN = 2, STATS = centres, FUN = '-')
                close <- apply(tempmask,1, function(x) sum(x^2)) < (detS*30)
                if (sum(close)<1) {
                    close <- nearesttrap(matrix(0, ncol = 2, nrow = 1), tempmask)
                }
                tempmask <- tempmask[close,, drop = FALSE]
                # if (!requireNamespace(mvtnorm, quietly = TRUE)) stop()
                # dbvn <- apply(tempmask, 1, mvtnorm::dmvnorm, mean=c(0,0), sigma=vcv, log=FALSE)
                # mymvn is faster and produces identical result
                invS <- solve(vcv)
                mymvn <- function(XY) exp(-(XY %*% invS %*% XY)/2) / 2/pi/det(vcv)
                dbvn <- apply(tempmask,1,mymvn)
                dbvn <- dbvn / sum(dbvn)
                f[close] <- dbvn
            }
            else {
                cell <- nearesttrap(centres, mask)
                f[cell] <- 1
            }
        }
        f
    }
    mm <- nrow(mask)
    xylist <- telemetryxy(CH)
    ## expand for untelemetered animals in CH
    nullxy <- rep(list(matrix(nrow=0, ncol=2)), nrow(CH))
    names(nullxy) <- rownames(CH)
    xylist <- c(xylist, nullxy[!(rownames(CH) %in% names(xylist))])
    ## compute pdf of centre for each animal
    sapply(xylist, dfn)[,rownames(CH)]  ## take care to maintain animal order
}
###############################################################################

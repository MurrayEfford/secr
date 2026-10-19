#----------------------------------------------------------------------
secr_gethazard <- function (m, binomN, cc, hk, PIA, usge) {
    nmix <- dim(PIA)[5]
    if (any(binomN == -2)) {   ## multi-catch trap
        nc <- dim(PIA)[2]
        s <- dim(PIA)[3]
        k <- dim(PIA)[4]
        haztemp <- gethcpp(
            as.integer(nc),
            as.integer(cc),
            as.integer(nmix),
            as.integer(k),
            as.integer(s),
            as.integer(m),
            as.integer(PIA),
            as.matrix(usge),
            as.double(hk))
        haztemp$h <- array(haztemp$h, dim = c(nmix, m, max(haztemp$hindex)+1))
        haztemp
    }
    else {
        list(h = array(-1, dim=c(nmix,1,1)), hindex = matrix(-1))
    }
}
#--------------------------------------------------------------------------------
## mixture proportions by animal        
## assume dim(PIA)[1] == 1
secr_getpmix <- function(knownclass, PIA, realparval)
{
    nc <- dim(PIA)[2]
    # not needed nc <- length(knownclass)   ## 2020-11-04
    k <- dim(PIA)[4]
    nmix <- dim(PIA)[5]
    pmixn <- matrix(1, nrow = nmix, ncol = nc)
    pmix <- numeric(nmix)
    if (nmix>1) {
        # index of first non-missing occasion s and detector k
        fsk <- sapply(1:nc, function(i) secr_firstsk(PIA[1,i,,,1, drop = FALSE]))
        ## fsk indexes the (s,k) array in column order: s varies fastest
        S <- dim(PIA)[3]
        sc <- as.vector((fsk-1) %% S + 1)
        kc <- as.vector((fsk-1) %/% S + 1)
        for (x in 1:nmix) {
            c <- PIA[cbind(1,1:nc,sc,kc,x)]
            pmixx <- realparval[c, 'pmix']    ## NOT CONSISTENT WITH pmix numeric(nmix)
            ## knownclass=2 maps to x=1 
            pmixn[x,] <- ifelse (knownclass > 1,
                                 ifelse (knownclass == (x+1), 1, 0),
                                 pmixx)
        }
        ## need pmix for each group... not ready yet
        attr(pmixn, 'pmix') <-  realparval[PIA[cbind(1,1,sc[1],kc[1],1:nmix)],'pmix']
    }
    pmixn
}
#--------------------------------------------------------------------------------
secr_getpID <- function(PIA, realparval, MRdata)
{
    ss <- dim(PIA)[3]
    nmix <- dim(PIA)[5]
    pID <- matrix(1, nrow = ss, ncol = nmix)
    if ('pID' %in% colnames(realparval)) {
        nc <- dim(PIA)[2]
        if (!is.null(MRdata$Tm) || !is.null(MRdata$Tn)) {
            for (s in 1:ss) {
                if (MRdata$markocc[s]<1)
                    for (x in 1:nmix) {
                        k <- match(TRUE, PIA[1,1,s,,x]>0)
                        c <- PIA[1,1,s,k,x]
                        pID[s,x] <- realparval[c, 'pID'] 
                    }
            }
        }
    }
    pID
}

#--------------------------------------------------------------------------------
secr_getmiscparm <- function(miscparm, detectfn, beta, parindx, cutval) {
    ## miscparm is used to package beta parameters that are not modelled
    ## and hence do not have a beta index specified by parindx.
    ## This includes the signal threshold and the mean and sd of noise.
    ## miscparm is passed to the C++ likelihood code, and 
    ## also as mask attribute to userdistfn
    nmiscparm <- length(miscparm)
    miscparm <- numeric(max(4, nmiscparm)) 
    if (detectfn %in% c(12,13))            ## experimental signal-noise
        miscparm[1:3] <- c(cutval, beta[max(unlist(parindx))+(1:2)])   ## fudge: last 2
    else if (detectfn %in% c(10,11))        ## Dawson & Efford 2009 models
        miscparm[1] <- cutval
    else if (nmiscparm > 0)
        miscparm[1:nmiscparm] <- beta[max(unlist(parindx)) + (1:nmiscparm)]
    miscparm
}
#--------------------------------------------------------------------------------
# returns array [nmask, ngroup, nsession] of real parameter values
# TRUE if all session masks have the same cells and covariates (or there is only one mask)
# called once per fit by secr.fit for userDfn (result stored as attribute of the function)
secr_samemask <- function (mask) {
    if (!ms(mask) || length(mask) == 1) return(TRUE)
    m1 <- mask[[1]]
    c1 <- covariates(m1)
    all(vapply(mask[-1], function (m) {
        nrow(m) == nrow(m1) && identical(m$x, m1$x) && identical(m$y, m1$y) &&
            identical(covariates(m), c1)
    }, logical(1)))
}
#--------------------------------------------------------------------------------
# Call a user-defined density function (details$userDfn) and check the result.
# The function has arguments (Dbeta, mask, ngroup, nsession). If it declares
# argument 'sessioncov' (data frame, one row per session in the call, with automatic
# covariates 'session' (factor) and 'Session' (0:(nsessions-1))) it is also supplied.
# Returns an array [nmask, ngroup, nsession].
# A result for a single session (legacy functions) is repeated for all sessions.
secr_calluserDfn <- function (f, Dbeta, mask, ngroup, nsession, sessioncov,
                              sessionlevels, which = 1:length(sessionlevels)) {
    fm <- names(formals(f))
    args <- list(Dbeta, mask, ngroup, nsession)
    if ('sessioncov' %in% fm) {
        sc <- if (is.null(sessioncov)) data.frame(row.names = seq_along(sessionlevels)) else sessioncov
        ## automatic covariates as in makeNewData: factor 'session', integer 'Session'
        if (!('session' %in% names(sc))) sc$session <- factor(sessionlevels, levels = sessionlevels)
        if (!('Session' %in% names(sc))) sc$Session <- seq_along(sessionlevels) - 1L
        args$sessioncov <- sc[which, , drop = FALSE]
    }
    D <- do.call(f, args)
    n <- nrow(mask) * ngroup
    if (length(D) == n) {
        D <- array(D, dim = c(nrow(mask), ngroup, nsession))   # same for all sessions
    }
    else if (length(D) != n * nsession) {
        stop ("details$userDfn returned an array of length ", length(D),
              "; expected dim = c(nmask, ngroup, nsession) = c(", nrow(mask), ", ", ngroup,
              ", ", nsession, ")")
    }
    array(D, dim = c(nrow(mask), ngroup, nsession))
}
#--------------------------------------------------------------------------------
secr_getD <- function (designD, beta, mask, parindx, link, fixed,
                  grouplevels, sessionlevels, parameter = 'D', sessioncov = NULL) {
    ## for .localstuff$spatialparametersD 'D', 'noneuc', 'sigmaxy', 'lambda0xy', 'a0'
    whichbeta <- beta[parindx[[parameter]]]
    if (length(whichbeta) == 0) return(NULL)
    if (!is.function(designD)) {
        if ((is.null(designD) || nrow(designD)==0) && (is.null(fixed[[parameter]]))) return(NULL)
    }
    if (ms(mask))
        nmask <- max(sapply(mask, nrow))
    else
        nmask <- nrow(mask)
    ngroup <- max(1, length(grouplevels))
    nsession <- length(sessionlevels)
    D <- array(dim = c(nmask, ngroup, nsession))
    dimnames(D) <- list(1:nrow(D), grouplevels, sessionlevels)
    if (!is.null(fixed[[parameter]])) {
        D[,,] <- fixed[[parameter]]
    }
    else {
        if (is.function(designD)) {
            ## user-defined density function; see secr_calluserDfn
            if (is.null(sessioncov)) sessioncov <- attr(designD, 'sessioncov')
            same <- attr(designD, 'samemask')
            if (!is.logical(same)) same <- secr_samemask(mask)
            if (same) {
                ## one call covering all sessions
                if (ms(mask)) mask <- mask[[1]]
                m <- nrow(mask)
                D[1:m,,] <- secr_calluserDfn(designD, whichbeta, mask, ngroup, nsession,
                                             sessioncov, sessionlevels)
            }
            else {
                ## session-specific masks differ (fixed 2021-12-10, 2022-05-24)
                for (session in 1:nsession) {
                    m <- nrow(mask[[session]])
                    D[1:m,,session] <- secr_calluserDfn(designD, whichbeta, mask[[session]],
                                                        ngroup, 1, sessioncov, sessionlevels, session)
                }
            }
        }
        else {
            Dfn <- attr(designD, 'Dfn')
            if (is.function(Dfn)) {
                    D[,,] <- Dfn(designD, whichbeta)
            }
            else {
                D[,,] <- designD %*% whichbeta  # linear predictor
            }
            D[,,] <- untransform (D, link[[parameter]])
        }
        # silently truncate D at zero
        # allow non-positive noneuc
        if (parameter == 'D')
            D[D<0] <- 0
    }
    D
}
#--------------------------------------------------------------------------------
secr_getmaskpar <- function(OK, D, m, sessnum, unmash, nmash) {
    if (!OK) {
        NULL   ## not in model
    }
    else {
        sessnum <- min(dim(D)[3],sessnum)
        density <- matrix(D[1:m,,sessnum], nrow = m)
        if (!all(is.finite(density))) {
            cat ('densities :', head(density), '\n')
            ## 2022-03-20
            warning ("bad densities in 'secrloglikfn' ",
                "(try different optimisation method, link, or model?)")
            # stop ("bad densities in 'secrloglikfn' ",
            #     "(try different optimisation method, link, or model?)")
        }
        ## optional scaling by session-specific number of clusters
        if (unmash) {
            if (!is.null(nmash))
                density <- density * length(nmash)
        }
        density
    }
}
#--------------------------------------------------------------------------------
getchat <- function (cc0, nc, n.distrib, group, usge, pmixn, pID,
                     cellsize, gkhk, pi.density, sumD,  PIA0, binomN, MRdata, miscparm, 
                     nsim, grain,ncores, debug) {
    mm <- nrow(pi.density)
    ngroup <- length(levels(group))   ## uNUSED
    ## note: should pass pi.mask as pi.density for known distribution all-sighting
    seed <- sample.int(.Machine$integer.max, 1)
    temp <- sightingchatcpp (
        as.integer(mm), 
        as.integer(nc), 
        as.integer(cc0), 
        as.integer(grain), 
        as.integer(ncores),
        as.integer(nsim),        ## number of replicate simulations for chat 
        as.integer(MRdata$sightmodel),  ## 5 allsighting known n0, 6 allsighting unknown n0, 7 telemetry marking
        as.double(sumD),
        as.double(cellsize),
        as.integer(n.distrib),
        as.integer(MRdata$anytelemetry), 
        as.integer(binomN),      ## detector -2 multi, -1 proximity 0 Poisson count 1 Binomial from usage, 2...etc. 
        as.integer(MRdata$markocc), 
        as.matrix(pID), 
        as.integer(group),      ## UNUSED? group number for 0<=n<*nc   [full likelihood only] 
        as.double(gkhk$gk), 
        as.double(gkhk$hk), 
        as.matrix(pi.density),        ## relative density - sums to 1.0
        as.integer(PIA0), 
        as.matrix(usge),         ## nk x s usage matrix 
        as.numeric(pmixn[,1]),
        as.numeric(seed),
        as.logical(debug>1)
    ) 
    if (temp$resultcode==0) {
        sumchat <- temp$chat
    }
    else {
        warning ("chat calculation failed, resultcode = ", temp$resultcode)
        sumchat <- (rep(0,3))
    }
    names(sumchat) <- c('Tu','Tm','Tn')
    if (!is.null(temp$chatmat)) {
        attr(sumchat, 'chatmat') <- temp$chatmat
    }
    return (sumchat)
}

secr_makegk <- function(dettype, detectfn, trps, mask, details, sessnum, NElist, D, 
                   miscparm, realparval, grain, ncores) {
    ## precompute gk, hk for point detectors
    if (all(dettype %in% c(0,1,2,5,8,13))) {
        distmat2 <- secr_getuserdist(trps, mask, details$userdist, sessnum, NElist, D,
                                     miscparm = miscparm, detectfn = detectfn)
        gkhk <- makegkPointcpp (
            as.integer(detectfn), 
            as.integer(grain),
            as.integer(ncores),
            as.matrix(realparval), 
            as.matrix(distmat2), 
            miscparm)
        if (any(dettype==8)) {   ## capped adjustment Not checked 2019-09-08
            gkhk <- cappedgkhkcpp (
                as.integer(nrow(realparval)),
                as.integer(nrow(trps)),
                as.double(attr(mask, "area")),
                as.double(D),
                as.double(gkhk$gk), as.double(gkhk$hk))  
        }
    }
    ## precompute gk, hk for polygon and transect detectors
    else if (all(dettype %in% c(3,4,6,7))) {
        ## k-1 because we have zero-terminated these vectors
        k <- secr_getk(trps)
        K <- if (length(k)>1) length(k)-1 else k
        cumk <- cumsum(c(0,k))[1:length(k)]
        dimension <- (dettype[1] %in% c(3,6)) + 1   ## 1 = 1D, 2 = 2D
        convexpolygon <- is.null(details$convexpolygon) || details$convexpolygon
        gkhk <- secrfunc::makegkPolygoncpp (
            as.integer(detectfn), 
            as.integer(dimension), 
            as.logical(convexpolygon), 
            as.integer(grain), 
            as.integer(ncores), 
            as.matrix(realparval), 
            as.integer(cumk),
            as.matrix(trps), 
            as.matrix(mask))
    }
    gkhk
}

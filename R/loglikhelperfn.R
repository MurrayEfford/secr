#----------------------------------------------------------------------
gethazard <- function (m, binomN, cc, hk, PIA, usge) {
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
getpmix <- function(knownclass, PIA, realparval)
{
    nc <- dim(PIA)[2]
    # not needed nc <- length(knownclass)   ## 2020-11-04
    k <- dim(PIA)[4]
    nmix <- dim(PIA)[5]
    pmixn <- matrix(1, nrow = nmix, ncol = nc)
    pmix <- numeric(nmix)
    if (nmix>1) {
        # index of first non-missing occasion s and detector k
        fsk <- sapply(1:nc, function(i) firstsk(PIA[1,i,,,1, drop = FALSE]))
        kc <- as.vector((fsk-1) %/% k + 1)
        sc <- as.vector((fsk-1) %/% k + 1)
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
getpID <- function(PIA, realparval, MRdata)
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
getmiscparm <- function(miscparm, detectfn, beta, parindx, cutval) {
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
getuserdist <- function (traps, mask, userdist, sessnum, noneuc, density, miscparm, HPX) {
    ## Apply user-provided distance function or basic distance function getdistmat2()
    if (is.null(userdist)) {
        getdistmat2(traps, mask, NULL, HPX)
    }
    # 
    # if (any(detector(traps) %in% .localstuff$polydetectors)) {
    #     matrix(0, nrow = nrow(traps), ncol = m)
    # }
    # else if (is.null(userdist)) {
    #     edist2cpp(as.matrix(traps), as.matrix(mask))
    # }
    else {
        userdistnames <- getuserdistnames(userdist)
        m <- nrow(mask)
        if (is.null(covariates(mask)))
            covariates(mask) <- data.frame(row.names = 1:m)
        if (('noneuc' %in% userdistnames) && !is.null(noneuc))
            covariates(mask)$noneuc <- noneuc  ## NE[1:m,,min(dim(NE)[3],sessnum)]
        if (('D' %in% userdistnames) && !is.null(density))
            covariates(mask)$D <- density
        ## pass miscellaneous unmodelled parameter(s)
        attr(mask, 'miscparm') <- miscparm
        distmat2 <- valid.userdist (userdist,
                                    detector(traps),
                                    xy1 = traps,
                                    xy2 = mask,
                                    mask = mask,
                                    sessnum = sessnum)^2
        baddist <- (!is.finite(distmat2)) | (distmat2<0) | is.na(distmat2)
        if (any(baddist)) {
            warning ("replacing infinite, negative and NA userdist values with 1e10")
            distmat2[baddist] <- 1e10
        }
        distmat2
    }
}
#--------------------------------------------------------------------------------
getD <- function (designD, beta, mask, parindx, link, fixed,
                  grouplevels, sessionlevels, parameter = 'D') {
    ## apply to either 'D' or 'noneuc'
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
            ## fixed 2021-12-10, 2022-05-24
            if (ms(mask)) {
                for (session in 1:nsession) {
                    m <- nrow(mask[[session]])
                    D[1:m,,session] <- designD(beta[parindx[[parameter]]], mask[[session]], ngroup, 1)
                }
            } 
            else {
                m <- nrow(mask)
                D[1:m,,1] <- designD(beta[parindx[[parameter]]], mask, ngroup, 1)
            }
        }
        else {
            ## 2021-12-09 f(x)
            beta <- beta[parindx[[parameter]]]
            f <- attr(designD, 'f')
            fcovname <- attr(designD, 'fcovname')
            if (is.function(f) && fcovname %in% dimnames(designD)[[2]]) {
                D[,,] <- f(designD[,fcovname, drop = FALSE], beta, dim(D))
            }
            else {
                D[,,] <- designD %*% beta  # linear predictor
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
getmaskpar <- function(OK, D, m, sessnum, unmash, nmash) {
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
                     nsim, grain,ncores) {
    kk <- nrow(usge)
    ss <- ncol(usge)
    mm <- nrow(pi.density)
    ngroup <- length(levels(group))   ## uNUSED
    ## note: should pass pi.mask as pi.density for known distribution all-sighting
    temp <- sightingchatcpp (
        as.integer(mm), 
        as.integer(nc), 
        as.integer(cc0), 
        as.integer(grain), 
        as.integer(ncores),
        as.integer(nsim),        ## number of replicate simulations for chat 
        as.integer(MRdata$sightmodel),  ## 5 allsighting known n0, 6 allsighting unknown n0
        as.double(sumD),
        as.double(cellsize),
        as.integer(n.distrib),
        as.integer(binomN),      ## detector -2 multi, -1 proximity 0 Poisson count 1 Binomial from usage, 2...etc. 
        as.integer(MRdata$markocc), 
        as.matrix(pID), 
        as.integer(group),      ## UNUSED? group number for 0<=n<*nc   [full likelihood only] 
        as.double(gkhk$gk), 
        as.double(gkhk$hk), 
        as.matrix(pi.density),        ## relative density - sums to 1.0
        as.integer(PIA0), 
        as.matrix(usge),         ## nk x s usage matrix 
        as.numeric(pmixn[,1])
    ) 
    if (temp$resultcode==0) {
        sumchat <- temp$chat
    }
    else {
        warning ("chat calculation failed, resultcode = ", temp$resultcode)
        sumchat <- (rep(0,3))
    }
    names(sumchat) <- c('Tu','Tm','Tn')
    return (sumchat)
}

makegk <- function(dettype, detectfn, trps, mask, details, sessnum, noneuc, D, miscparm, realparval, grain, ncores) {
    ## precompute gk, hk for point detectors
    if (all(dettype %in% c(0,1,2,5,8,13))) {
        distmat2 <- getuserdist(trps, mask, details$userdist, sessnum, noneuc, D, miscparm, detectfn == 20)
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
        k <- getk(trps)
        K <- if (length(k)>1) length(k)-1 else k
        cumk <- cumsum(c(0,k))[1:length(k)]
        dimension <- (dettype[1] %in% c(3,6)) + 1   ## 1 = 1D, 2 = 2D
        convexpolygon <- is.null(details$convexpolygon) || details$convexpolygon
        gkhk <- makegkPolygoncpp (
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

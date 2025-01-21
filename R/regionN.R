############################################################################################
## package 'secr'
## regionN.R
## population size in an arbitrary region
## 2018-08-22 all calls of sumDpdot use only first element, ignoring individual variation
## 2020-05-13 allow individual variation in sumpdot object$design0$individual
## 2020-11-04 reliance on dim(PIA0)[2] for number of individuals failed in sumDpdot
##            because PIA0 not trimmed to session-specific n; fixed 2020-11-04
## 2020-11-04 session names used in output for multi-session data
## 2022-01-22 fixed bug in sumDpdot that ignored details$ignoreusage
## 2022-10-08 fixed bug in region.N when region is mask other than original
## 2023-02-06 sumDpdot tweaked for robustness to varying detectors/session
## 2024-12-15 sumDpdot ignored varying relative density
## 2025-01-09 relativeD tweaks

################################################################################

region.N.secrlist <- function (
        object, region = NULL, spacing = NULL, session = NULL,
        group = NULL, se.N = TRUE, alpha = 0.05, loginterval = TRUE,
        keep.region = FALSE, nlowerbound = TRUE, RN.method = 'poisson',
        pooled.RN = FALSE, ncores = NULL, ...) {
    lapply(object, region.N, region, spacing, session, group, se.N, alpha, 
           loginterval,keep.region, nlowerbound, RN.method, pooled.RN, ncores)
}
################################################################################

region.N.secr <- function (object, region = NULL, spacing = NULL, session = NULL,
    group = NULL, se.N = TRUE, alpha = 0.05, loginterval = TRUE,
    keep.region = FALSE, nlowerbound = TRUE, RN.method = 'poisson',
    pooled.RN = FALSE, ncores = NULL, ...) {

    ## Notes
    ## se.N = FALSE returns scalar N
    ###########################################################
    ## for gradient of E.N wrt density betas
    betaEN <- function (betaD, object, regionmask, group, session) {
        ## regionmask is a mask (no need for spacing)
        ## assume single session
        ## indx identifies beta parameters for density D
        object$fit$par[indx] <- betaD
        if (object$CL) {
            object <- completeDbeta (object, sessnum)
            object$CL <- FALSE  # for this call
        }
        region.N(object, regionmask, spacing = NULL, session = session,
            group = group, se.N = FALSE, keep.region = FALSE,
            pooled.RN = FALSE)
    }
    ###########################################################
    ## for gradient of R.N wrt all betas
    betaRN <- function (beta, object, regionmask) {
        ## regionmask is a mask (no need for spacing)
        ## assume single session
        ## n, cellsize, sessnum global
        object$fit$par <- beta
        if (object$CL) {
            D <- covariates(derivedDsurface(object, regionmask, sessnum))$D.0
        }
        else {
            D <- predictD(object, regionmask, group, session, parameter = 'D')
        }
        noneuc <- predictD(object, regionmask, group, session, parameter = 'noneuc')
        n + sumDpdot (object, sessnum, regionmask, D, noneuc, cellsize,
                  constant = FALSE, oneminus = TRUE, pooled = pooled.RN, ncores = ncores)[1]
    }
    ###########################################################
    if (is.null(region)) {
        region <- object$mask  # using original mask as region
    }
    if (!all(session %in% session(object$capthist)))
        stop ("session incompatible with object ")

    if (is.null(group))
        group <- 1

    if (is.null(session))
        session <- session(object$capthist)

    if (!ms(object))
        pooled.RN <- FALSE

    if (pooled.RN)
        ## use only model for first session
        session <- session[1]

    ####################################################################
    ## if N requested for multiple sessions,
    ## call region.N recursively for each session
    nsess <- length(session)
    details <- object$details
    if (nsess > 1) {
        ## predict for each session
        out <- vector('list')
        for (sess in 1:nsess) {
            if (ms(region)) {
                tempregion <- region[[sess]]
            }
            else
                tempregion <- region
            detail <- details
            if (ms(detail$userdist))
                detail$userdist <- detail$userdist[[sess]]
            object$details <- detail
            out[[sess]] <- region.N (object, region = tempregion,
                spacing = spacing, session = session[sess], group = group,
                se.N = se.N, alpha = alpha, loginterval = loginterval,
                keep.region = keep.region, nlowerbound = nlowerbound,
                RN.method = RN.method, pooled.RN = FALSE)
        }
        names(out) <- session(object$capthist)   ## 2020-11-04
        out
    }

    ####################################################################
    ## otherwise, this is a non-recursive call for one session...
    else {
        if (ms(object$mask))
            mask <- object$mask[[session]]
        else
            mask <- object$mask
        masktype <- if (inherits(mask, "linearmask")) "linearmask" else "mask"
        
        ########################################################
        ## if necessary, convert vector region to raster
        if (inherits(region, 'mask')) {
            ## includes linearmask
            if (ms(region))
                regionmask <- region[[session]]
            else
                regionmask <- region
        }
        else {
            if (is.null(spacing)) {
                ## use original mask spacing by default
                spacing <- spacing(mask)
            }
            if (masktype == "mask") {
                # generate mask object from polygon region
                region <- boundarytoSF(region)
                bbox <- matrix(st_bbox(region), ncol=2)
                regionmask <- make.mask(bbox, poly = region, buffer = 0,
                    spacing = spacing, type = 'polygon',
                    check.poly = FALSE)
            }
            else {
                stop("linear mask not available")
                # if (!requireNamespace("secrlinear", quietly=TRUE))
                #     stop ("could not load secrlinear")
                # regionmask <- secrlinear::read.linearmask(data = region,
                #     spacing = spacing, spacingfactor = attr(mask, "spacingfactor"))
            }
        }

        if (is.matrix(object$details$userdist)) {
            if (ncol(object$details$userdist) != nrow(regionmask)) {
                warning("userdist matrix incompatible with region mask - ignored")
                object$details$userdist <- NULL
            }
        }
        #######################################################

        cellsize <- getcellsize(regionmask)    ## bug fixed 2022-10-08
        regionsize <- nrow(regionmask) * cellsize

        ngrp <- function(x) sum(getgrpnum(x, object$groups) == group)
        if (ms(object)) {
            if (pooled.RN)
                n <- sum(sapply(object$capthist, ngrp))
            else
                n <- ngrp(object$capthist[[session]])
        }
        else {
            n <- nrow(object$capthist)
        }
        sessnum <- match (session, session(object$capthist))
        #######################################################
        ## conditional likelihood fit without density model
        if (object$CL && is.null(object$model$D)) {
            temp <- derived(object, se.D = se.N) ## inefficient as repeats for each sess
            if (!is.data.frame(temp))
                temp <- temp[[session]]
            D <- temp['D', 'estimate']
            seD <- temp['D', 'SE.estimate']
            EN <- D * regionsize
            if (!se.N) return (EN)    ## and stop here
            seEN <- seD * regionsize
        }

        #######################################################
        ## full likelihood fit or conditional likelihood with relative density model...
        else {
            if (is.null(object$model$D) || is.null(object$link$D))
                stop ("model or link function not found in object")
            if ((object$model$D == ~1) && !userD(object)) {
                predicted <- predict(object)
                if (!is.data.frame(predicted))
                    predicted <- predicted[[1]]
                D <- predicted['D','estimate']
                seD <- predicted['D', 'SE.estimate']

                EN <- D * regionsize
                if (!se.N) return (EN)    ## and stop here
                seEN <- seD * regionsize
            }
            else {
                if (object$CL) {
                    D <- covariates(derivedDsurface(object, regionmask, sessnum))$D.0
                }
                else {
                    D <- predictD (object, regionmask, group, session, parameter = 'D')
                }
                EN <- sum(D) * cellsize
                if (!se.N) return (EN)    ## and stop here
                indx <- object$parindx$D
                
                if (object$CL) {
                    # using approximate conversion of conditional to full model for variances
                    object <- completeDbeta(object, vcv = TRUE)   # vcv unreliable
                }
                dENdphi <- nlme::fdHess (object$fit$par[indx],
                                         betaEN, object = object, region = regionmask, group = group,
                                         session = session)$gradient
                beta.vcv <- object$beta.vcv[indx,indx, drop = FALSE]
                seEN <- (dENdphi %*% beta.vcv %*% dENdphi)^0.5
            }
        }

        #######################################################################
        ## realised N
        ## only makes sense for individual detectors (not unmarked or presence)
        ## assume if we have got this far that SE is required
        if (ms(object)) {
            det <- detector(traps(object$capthist)[[session]])
        }
        else {
            det <- detector(traps(object$capthist))
        }
        if (all(det %in% .localstuff$individualdetectors)) {
            noneuc <- predictD (object, regionmask, group, session, parameter = 'noneuc')
            RN.method <- tolower(RN.method)
            if (RN.method %in% c('mspe', 'poisson')) {
                notdetected <- sumDpdot (object, sessnum, regionmask, D, noneuc,
                    cellsize, constant = FALSE, oneminus = TRUE, pooled = pooled.RN, ncores = ncores)[1]
                RN <- n + notdetected
                if (RN.method == 'mspe') {
                    ## evaluate gradient of RN wrt betas at MLE
                    dNdbeta <- nlme::fdHess (object$fit$par, betaRN, object = object,
                                             region = regionmask)$gradient
                    ## compute variance from gradient & vcv
                    pdotvar <- dNdbeta %*% object$beta.vcv %*% dNdbeta
                    seRN <- (notdetected + pdotvar)^0.5
                }
                else if (RN.method == 'poisson') {
                    seRN <- (seEN^2 - EN)^0.5
                }
            }
            ## RN.method = 'EN'
            else {
                RN <- EN
                seRN <- (seEN^2 - EN)^0.5
            }
        }
        else { 
            RN <- NA
            seRN <- NA 
        }

        #######################################################################
        temp <- data.frame(
            row.names = c('E.N','R.N'),
            estimate = c(EN,RN),
            SE.estimate = c(seEN,seRN))
        ## lower bound added 2011-07-15
        ## vector (0,n) means apply to R.N not E.N
        if (nlowerbound)
            temp <- add.cl (temp, alpha, loginterval, c(0, n))
        else
            temp <- add.cl (temp, alpha, loginterval, c(0, 0))
        
        temp$n <- rep(n, nrow(temp))
        attr(temp, 'regionsize') <- masksize(region) 
        attr(temp, 'region') <- if (keep.region) region else NULL
        temp
    }
}

sumDpdot <- function (object, sessnum = 1, mask, D, noneuc, cellsize, constant = TRUE,
                      oneminus = FALSE, pooled = FALSE, bycluster = FALSE, ncores = NULL)

# Return integral for given model and new mask, D
# 'sessnum' is integer index of session (factor level of the 'session' attribute in capthist)
# object must have at least capthist, detectfn
# D should be scalar or vector of length nrow(mask)

# if 'constant' a much simplified calculation is used, assuming
# constant detection and density, and full detector usage

{
    if (ms(object))
        capthists <- object$capthist[[sessnum]]
    else
        capthists <- object$capthist

    if (ms(mask)) {
        mask <- mask[[sessnum]]
    }
    
    ## Allow for fixed beta parameters 
    beta <- complete.beta(object)
    n       <- max(nrow(capthists), 1)
    s       <- ncol(capthists)
    noccasions <- s

    ##############################################
    ## traps
    if (pooled & ms(object))
        trps <- do.call(rbind, c(traps(object$capthist), list(addusage = TRUE)))
    else
        trps   <- traps(capthists)  ## use session-specific traps
    if (!all(detector(trps) %in% .localstuff$individualdetectors))
        stop ("require individual detector type for sumDpdot")

    ##############################################
    dettype <- detectorcode(trps, noccasions = s)
    nmix    <- getnmix(object$details)
    knownclass <- getknownclass(capthists, nmix, object$hcov)
    groups <- object$groups  
    markocc <- markocc(traps(capthists))
    if (is.null(markocc)) markocc <- rep(1,s)
    MRdata <- list(markocc = markocc, firstocc = -1)
    k <- getk(trps)
    K <- if (length(k)>1) length(k)-1 else k
    binomN <- object$details$binomN
    m      <- length(mask$x)            ## assume session-specific mask...
    ncores <- setNumThreads(ncores)
    grain <- if (ncores==1) 0 else 1
    binomNcode <- recodebinomN(dettype, binomN, telemcode(trps))
    
    ##############################################
    
    if (constant) {
        if (is.null(beta))
            real <- detectpar(object)
        else {

            real <- makerealparameters (object$design0, beta,
                object$parindx, object$link, object$fixed)  # naive
            real <- as.list(real)
            names(real) <- parnames(object$detectfn)
        }
        a <- cellsize * sum(pdot(X = mask, traps = trps, detectfn = object$detectfn,
                             detectpar = real, noccasions = noccasions, ncores = ncores))
        return(a * D)
    }
    else {
        # PIA0 <- object$design0$PIA[sessnum,,1:s,,,drop = FALSE]
        # PIA0 <- object$design0$PIA[sessnum,1:n,1:s,,,drop = FALSE]   ## 2020-11-04
        PIA0 <- object$design0$PIA[sessnum,1:n,1:s,1:K,,drop = FALSE]  ## 2023-02-06
        
        #############################################
        ## trick to allow for changed data 2009 11 20
        ## nmix>1 needs further testing 2010 02 26
        ## NOTE 2010-11-26 THIS LOOKS WEAK
        # if (dim(PIA0)[2] != n) {
        #     PIA0 <- array(rep(PIA0[1,1,,,],n), dim=c(s,K,nmix,n))
        #     PIA0 <- aperm(PIA0, c(4,1,2,3))   ## n,s,K,nmix
        # }
        
        #############################################################
        ## parameter table 
        realparval0 <- makerealparameters (object$design0, beta,
            object$parindx, object$link, object$fixed)  # naive
        if (!is.null(object$fixed$D))
            Dtemp <- object$fixed$D
        else if (object$CL)
            Dtemp <- NA
        else
            Dtemp <- D[1]
        Xrealparval0 <- reparameterize (realparval0, object$detectfn, 
                                        object$details, mask, trps, Dtemp, s)
        #############################################################
        ## usage
        # 2022-01-22
        # if (is.null(usage(trps))) {
        if (is.null(usage(trps)) || object$details$ignoreusage) {
                usage(trps) <- matrix(1, nrow = K, ncol = s)
        }
        used <- (usage(trps) > 1e-10) * 1
        if (any(used==0)) {
            PIA0 <- PIA0 * rep(rep(t(used),rep(n,s*K)),nmix)
        }
        
        usge <- usage(trps)
        
        ############################################################
        pmixn <- getpmix (knownclass, PIA0, Xrealparval0)
        pID <- getpID(PIA0, Xrealparval0, MRdata)
        
        # only works for 1-session fitted object
        miscparm <- getmiscparm(object$details$miscparm, object$detectfn, coef(object)[,1], 
                                object$parindx, object$details$cutval) 
        
        #############################################################
        ## add density as third column of mask
        if (!(length(D) %in% c(1,nrow(mask))))
            stop ("D does not match mask in sumDpdot")

        if (length(D) == 1) {
            useD <- FALSE
        }
        else {
            mask <- cbind (mask, D)   ## loses 'mask' class
            useD <- TRUE
        }
        #############################################################
        ## across all traps, regardless of clusters
        gkhk <- makegk (dettype, object$detectfn, trps, mask, object$details, sessnum, 
                        noneuc, D, miscparm, Xrealparval0, grain, ncores)
            
        haztemp <- gethazard(m, binomNcode, nrow(Xrealparval0), gkhk$hk, PIA0, usge)
        #############################################################
        
        if (bycluster) {
            centres <- cluster.centres(trps)
            nclust <- nrow(centres)
            if (is.null(attr(trps, 'cluster'))) {
                clusterID(trps) <- 1:nclust
            }
        }
        else {
            nclust <- 1
            clusterID(trps) <- rep(1,nrow(trps))
            temptrap <- trps
        }
        pdot <- numeric(nclust)
        for (i in 1:nclust) {
            if (nclust>1) {
                clustok <- as.numeric(clusterID(trps)) == i
                temptrap <- subset(trps, subset = clustok)
                gkhk <- makegk (dettype, object$detectfn, temptrap, mask, object$details, 
                    sessnum, noneuc, D, miscparm, Xrealparval0, grain, ncores)
                usge <- usage(temptrap)
                haztempc <- gethazard(m, binomNcode, nrow(Xrealparval0), gkhk$hk, PIA0[,,,clustok,,drop=FALSE], usge)
            }
            
            ## 2020-05-16
            ## object$design0$individual may be NULL (missing)
            if (any(binomNcode == -2)) {
                CH0 <- nullCH (c(n,s), object$design0$individual)
            }
            else {
                K <- nrow(temptrap)
                CH0 <- nullCH (c(n,s,K), object$design0$individual)
            }
            D <- matrix(D, nrow = m, ncol = 1)
            pd <- integralprw1 (
                cc0 = nrow(Xrealparval0), 
                haztemp = if (nclust>1) haztempc else haztemp, 
                gkhk = gkhk, 
                # pi.density = matrix(1/m, nrow = m, ncol = 1), 
                pi.density = D/sum(D), # bug fix 2024-12-15
                PIA0 = PIA0, 
                CH0 = CH0, 
                binomNcode = binomNcode, 
                MRdata = MRdata,
                grp = rep(1,n),    ## dummy single group
                usge = usge, 
                pmixn = pmixn, 
                pID = pID,
                grain = grain,
                ncores = ncores)
            # else {
            #     pd <- integralprw1poly (detectfn, Xrealparval0, haztemp, gkhk$hk, gkhk$H, pi.density, PIA0, 
            #                             data$CH0, data$binomNcode, data$grp, data$usge, data$mask,
            #                             pmixn, data$maskusage, details$grain, details$minprob, debug = details$debug>3)
            # 
            # }
            ## scale by absolute density (not passed to integralprw1)
            pd <- pd * cellsize * nrow(mask) * mean(D)
            pdot[i] <- length(pd) / sum(1/pd)
        }
        
        if (oneminus) {
            sumD <- ifelse (length(D) == 1, D * nrow(mask), sum(D))
            return(sumD * cellsize - pdot)
        }
        else
          return(pdot)
    }
}
############################################################################################

expected.n <- function (object, session = NULL, group = NULL, bycluster = FALSE,
                        splitmask = FALSE, ncores = NULL) {

    ## Note
    ## splitmask toggles between two methods for clustered detectors:
    ## 1. is to integrate over the whole mask, restricting detection to each cluster in turn
    ## 2. is to split mask into Dirichlet subregions by distance to detector centre
    ## and to integrate over all detectors, assuming those far away will never detect from
    ## within a subregion to which they do not belong
    ## Probably, 1. is more robust

    if (!all(session %in% session(object$capthist)))
        stop ("session incompatible with object")

    if (!is.null(group))
        stop ("not yet working for groups")

    if (!all(group %in% interaction (object$groups)))
        stop ("unrecognised groups")

    if (is.null(session))
        session <- session(object$capthist)

    ####################################################################
    ## if En requested for multiple sessions,
    ## call En recursively for each session
    nsess <- length(session)
    if (nsess > 1) {
        ## predict for each session
        out <- vector('list')
        for (sess in session) {
            out[[sess]] <- expected.n (object, sess, group, bycluster, splitmask)
        }
        out
    }

    ####################################################################
    ## otherwise, this is a non-recursive call for one session...
    else {

        if (ms(object$mask))
            mask <- object$mask[[session]]
        else
            mask <- object$mask
        cellsize <- getcellsize(mask)
        if (ms(object)) {
            n <- nrow(object$capthist[[session]])
            trps <- traps(object$capthist[[session]])
        }
        else {
            n <- nrow(object$capthist)
            trps <- traps(object$capthist)
        }
        sessnum <- match (session, session(object$capthist))

        #######################################################
        ## for conditional likelihood fit,
        if (object$CL && is.null(object$model$D)) {
            temp <- derived(object, se.D = FALSE) ## inefficient as repeats for each sess
            if (!is.data.frame(temp))
                temp <- temp[[session]]
            D <- temp['D', 'estimate']
        }
        #######################################################
        ## for full likelihood fit...
        else {
            if (is.null(object$model$D) | is.null(object$link$D))
                stop ("model or link function not found in object")

            if (object$model$D == ~1 && !object$CL) {
                predicted <- predict(object)
                if (!is.data.frame(predicted))
                    predicted <- predicted[[1]]
                D <- rep(predicted['D','estimate'], nrow(mask))
            }
            else {
                if (object$CL)
                    D <- covariates(derivedDsurface(object, mask, sessnum))$D.0
                else
                    D <- predictD (object, mask, group, session, parameter = 'D')
            }
        }
        #######################################################

        if (is.function(object$details$userdist)) {
          if (object$model$noneuc == ~1) {
            predicted <- predict(object)
            if (!is.data.frame(predicted))
              predicted <- predicted[[1]]
            noneuc <- rep(predicted['noneuc','estimate'], nrow(mask))
          }
          else {
            noneuc <- predictD (object, mask, group, session, parameter = 'noneuc')
          }
        }
        else noneuc <- rep(NA, nrow(mask))

        #################################################################
        if (bycluster) {
            if (splitmask) {
                centres <- cluster.centres(trps)
                nclust <- nrow(centres)
                if (is.null(attr(trps, 'cluster'))) {
                    clusterID(trps) <- 1:nclust
                }
                cluster <- nearesttrap (mask, centres)
                mask <- split (mask, cluster)
                D <- split(D, cluster)
                noneuc <- split(noneuc, cluster)
                out <- numeric (nclust)
                for (i in 1:nclust) {
                    out[i] <- sumDpdot(object = object, sessnum = sessnum,
                                       mask = mask[[i]], D = D[[i]], noneuc = noneuc[[i]], cellsize = cellsize,
                                       constant = FALSE, oneminus = FALSE, ncores = ncores)[1]
                }
            }
            else {
                
                if (any(detector(trps) %in% c('single','multi', 'polygonX', 'transectX'))) {
                    warning("expected.n assumes clusters independent when detectors exclusive (single, multi, polygonX, transectX)")
                }
                out <- sumDpdot(object = object, sessnum = sessnum,
                                mask = mask, D = D, noneuc = noneuc, cellsize = cellsize,
                                constant = FALSE, oneminus = FALSE, bycluster = TRUE, 
                                ncores = ncores)
            }
            out
        }
        else {
            sumDpdot (object, sessnum, mask, D, noneuc, cellsize,
             constant = FALSE, oneminus = FALSE, ncores = ncores)[1]
        }
        #################################################################
    }
}


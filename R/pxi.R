##############################################################################
## package 'secr'
## pxi.R
## 2025-01-01 fork from fxi.R
###############################################################################

sharedData <- function (object, i, sessnum, X, ncores, naive = FALSE) {
    ## temporary fix for lack of fastproximity code
    object$details$fastproximity <- FALSE 
    
    ## data for a single session
    data <- prepareSessionData(object$capthist, object$mask, object$details$maskusage, 
                               object$design, object$design0, object$detectfn, object$groups, 
                               object$fixed, object$hcov, object$details)
    
    sessionlevels <- session(object$capthist)
    beta <- coef(object)$beta
    beta <- fullbeta(beta, object$details$fixedbeta)
    detparindx <- object$parindx[!(names(object$parindx) %in% c('D', 'noneuc'))]
    detlink <- object$link[!(names(object$link) %in% c('D', 'noneuc'))]

    data <- data[[sessnum]]
    reusemask <- is.null(X)
    if (reusemask) {
        X <- data$mask
    }
    else {
        X <- matrix(unlist(X), ncol = 2)
    }
    #----------------------------------------
    # restrict to selected individuals
    xy <- data$xy 
    if (is.null(i)) {
        ok <- 1:nrow(data$CH)
    }
    else {
        ok <- i
        if (!is.null(xy)) {
            ## 2022-02-13 don't want 'no detections on occasion x'
            ch <- suppressWarnings(subset(object$capthist, ok))  
            xy <- getxy(data$dettype, selectCHsession(ch, sessnum))
        }
    }
    if (length(dim(data$CH)) == 2) {
        CH <- data$CH[ok,,drop=FALSE]
    }
    else {
        CH <- data$CH[ok,,,drop=FALSE]
    }
    grp <- data$grp[ok]
    
    ncores <- setNumThreads(ncores)
    grain <- if (ncores==1) 0 else 1;
    
    #----------------------------------------
    # Density
    if (is.null(object$model$D))
        D.modelled <- FALSE
    else {
        if (!is.null(object$fixed$D))
            D.modelled <- FALSE
        else
            D.modelled <- (object$model$D != ~1)
    }
    if (D.modelled) {
        predD <- predictDsurface (object)
        if (ms(object))
            predD <- predD[[sessnum]]
        D <- covariates(predD)$D.0  ## does not apply if groups
        pimask <- D / sum(D)   ## vector of probability mass for each mask cell
    }
    else {
        mm <- nrow(data$mask)
        pimask <- rep(1, mm)  ## could be 1/mm, but as we normalise anyway...
    }
    ## fetch predicted density at each new point X
    ## covariates(session.mask) <- data.frame(pi = pimask)
    if (!is.null(covariates(data$mask)))
        covariates(data$mask) <- cbind(data.frame(pi = pimask), covariates(data$mask))
    else
        covariates(data$mask) <- data.frame(pi = pimask)
    ## does this work for linearmask?
    tmpmask <- suppressWarnings(addCovariates(X, data$mask, strict = TRUE))
    piX <- covariates(tmpmask)$pi
    piX[is.na(piX)] <- 0
    #----------------------------------------
    
    ## TO BE FIXED
    
    # NE <- getD (object$designNE, beta, object$mask, object$parindx, object$link, object$fixed,
    #             levels(data$grp[[1]]), sessionlevels, parameter = 'noneuc')
    # NEX <- getD (object$designNE, beta, X, object$parindx, object$link, object$fixed,
    #             levels(data$grp[[1]]), sessionlevels, parameter = 'noneuc')
    # 
    NE <- NULL
    
    #---------------------------------------------------
    ## allow for scaling of detection
    Dtemp <- if (D.modelled) mean(D) else NA
    if (naive) {
        realparval  <- makerealparameters (object$design0, beta, detparindx,
                                            detlink, object$fixed)
        Xrealparval <- reparameterize (realparval, object$detectfn, object$details,
                                       data$mask, data$traps, Dtemp, data$s)
        PIA <- object$design0$PIA[sessnum, ok, 1:data$s, 1:data$K, ,drop=FALSE]
    }
    else {
        realparval  <- makerealparameters (object$design, beta, detparindx,
                                           detlink, object$fixed)
        Xrealparval <- reparameterize (realparval, object$detectfn, object$details,
                                       data$mask, data$traps, Dtemp, data$s)
        PIA <- object$design$PIA[sessnum, ok, 1:data$s, 1:data$K, ,drop=FALSE]
    }
    
    pmix <- getpmix (data$knownclass[ok], PIA, Xrealparval)  ## membership prob by animal
    
    ## unmodelled beta parameters, if needed
    miscparm <- getmiscparm(object$details$miscparm, object$detectfn, object$beta, 
                            object$parindx, object$details$cutval)
    
    gkhk <- makegk (data$dettype, object$detectfn, data$traps, data$mask, object$details, sessnum, 
                    NE, D, miscparm, Xrealparval, grain, ncores)
    haztemp <- gethazard (data$m, data$binomNcode, nrow(Xrealparval), gkhk$hk, PIA, data$usge)
    
    # return a list
    
    list(
        Xrealparval = Xrealparval,
        haztemp     = haztemp,
        gkhk        = gkhk,
        pimask      = pimask,
        PIA         = PIA,
        CH          = CH,
        grp         = grp,
        pmix        = pmix,
        grain       = grain,
        ncores      = ncores,
        data        = data
    )
    
}

pxi <- function (object, i = NULL, sessnum = 1, X = NULL, ncores = NULL, ...) {
    
    # compute required objects
    sD <- sharedData(object, i, sessnum, X, ncores, naive = TRUE)

    with(sD, {
        # null object
        CH[] <- 0
        pimask <- rep(1,data$m)  
        if (sD$data$dettype[1] %in% c(0,1,2,5,8,13)) {
            prmat <- allhistfxi (data$m, Xrealparval, haztemp, gkhk, pimask, PIA, data$usge,
                                 CH, data$binomNcode, grp, pmix, grain, ncores)
        }
        else {
            prmat <- allhistpolygonfxi (object$detectfn, Xrealparval, haztemp, gkhk$hk, gkhk$H, pimask, PIA, 
                                        CH, xy, data$binomNcode, grp, data$usge, data$mask,
                                        pmix, data$maskusage, grain, ncores, object$details$minprob)
        }
        1-prmat
    })
}

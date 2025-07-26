#############################################################################
## package 'secr'
## derived density from conditional (relativeD) models
## 2025-01-01
## 2025-01-01 pxi forked from fxi.R
## 2025-01-03 derivedDbeta0 groups
## 2025-01-11 se.beta0
## 2025-01-20 replaced derivedDbeta0 by derivedDcoef
## 2025-05-17 addCovariates from mask if needed
## 2025-07-26 pxi included here

#############################################################################

# session-specific

#-------------------------------------------------------------------------------

onek <- function (beta, object, individuals, sessnum)
    ## object is a fitted secr object 
    ## individuals is vector indexing the subset of animals to be used
    # Return k-hat for given beta
    # Only 1 session
{
    out <- sum(1 / esa (object, sessnum, beta, Dweight = TRUE)[individuals])
    # cat('beta ', beta, ' out ', out, '\n')
    out
}
#-------------------------------------------------------------------------------

kgradient <- function (object, individuals, sessnum, ...)
    ## object is a fitted secr object 
    ## individuals is vector indexing the subset of a to be used
{
    nlme::fdHess(object$fit$par, onek, object, individuals, sessnum, ...)$gradient
}
#-------------------------------------------------------------------------------

derivedDcoef <- function (object, sessnum = 1, groups = NULL, se = FALSE) {
    if (is.null(object$model$D) || is.null(object$link$D) || !object$CL) {
        warning ("not relative density model")
        return(NULL)
    }
    else {
        capthist <- object$capthist
        if (ms(capthist)) capthist <- capthist[[sessnum]]
        grp <- secr_group.factor(capthist, groups)
        individuals <- split (1:nrow(capthist), grp)
        ngrp <- length(individuals)   ## number of groups
        
        se.derivedk <- function (selection, object, selected.a, sessnum) {
            A <-  masksize(object$mask)
            s2 <- switch (tolower(object$details$distribution),
                          poisson  = sum (1/selected.a^2),
                          binomial = sum (( 1 - selected.a / A) / selected.a^2))
            kgrad  <- kgradient (object, selection, sessnum)
            vark <- kgrad %*% object$beta.vcv %*% kgrad
            sqrt(vark + s2)
        }
        
        getcoef <- function (selection) {
            selected.a <- esa (object, sessnum, Dweight = TRUE)[selection]
            k <- sum(1 / selected.a)
            if (se) {
                warning ("derivedDcoef() underestimates se(beta0)", call. = FALSE)
                se.k <- se.derivedk (selection, object, selected.a, sessnum)
            }
            else {
                se.k <- NA
            }
            # return on link scale
            beta0 <- transform(k, object$link$D)
            se.beta0 <- se.transform(k, se.k, object$link$D)
            oldcoef <- coef(object)
            alpha <- attr(oldcoef, 'alpha')
            z <- abs(qnorm(1 - alpha/2))
            Dcoef <- c(beta0, se.beta0, beta0 - z * se.beta0, beta0 + z * se.beta0)
            tmp <- rbind(D = Dcoef, oldcoef)
            if (object$link$D == 'identity') {
                tmp[grepl('D.', rownames(tmp)),] <- tmp[grepl('D.', rownames(tmp)),] * tmp[1,1]
            }
            tmp
            
        }
        
        # NB getcoef() gives same (estimate, SE.estimate) as derived(object, Dweight=T), 
        # and takes about the same time
        # May later switch to derived()
        
        if ( ngrp > 1) {
            # multiple groups
            lapply(individuals, getcoef)
        }
        else {    
            # one group
            getcoef(individuals[[1]])
        }
        
        
    }
}
#-------------------------------------------------------------------------------

sharedData <- function (object, i, sessnum, X, ncores, naive = FALSE) {
    ## temporary fix for lack of fastproximity code
    object$details$fastproximity <- FALSE 
    
    ## data for a single session
    data <- secr_prepareSessionData(object$capthist, object$mask, object$details$maskusage, 
                                    object$design, object$design0, object$detectfn, object$groups, 
                                    object$fixed, object$hcov, object$details)
    
    sessionlevels <- session(object$capthist)
    beta <- coef(object)$beta
    beta <- secr_fullbeta(beta, object$details$fixedbeta)
    detparindx <- object$parindx[!(names(object$parindx) %in% c('D', 'noneuc'))]
    detlink <- object$link[!(names(object$link) %in% c('D', 'noneuc'))]
    
    data <- data[[sessnum]]
    reusemask <- is.null(X)
    if (reusemask) {
        X <- data$mask
    }
    else {
        if (!inherits(X, 'mask')) {
            X <- as.data.frame(matrix(unlist(X), ncol = 2))
            names(X) <- c('x','y')
            X <- read.mask(data = X)
            if (!is.null(covariates(data$mask))) {
                X <- addCovariates(X, data$mask)
            }
        }
        data$mask <- X
        data$m <- nrow(X)
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
            xy <- secr_getxy(data$dettype, secr_selectCHsession(ch, sessnum))
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
        predD <- predictDsurface (object, mask = X)
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
    
    # NE <- secr_getD (object$designNE, beta, object$mask, object$parindx, object$link, object$fixed,
    #             levels(data$grp[[1]]), sessionlevels, parameter = 'noneuc')
    # NEX <- secr_getD (object$designNE, beta, X, object$parindx, object$link, object$fixed,
    #             levels(data$grp[[1]]), sessionlevels, parameter = 'noneuc')
    # 
    NE <- NULL
    #---------------------------------------------------
    ## allow for scaling of detection
    Dtemp <- if (D.modelled) mean(D) else NA
    if (naive) {
        realparval  <- secr_makerealparameters (object$design0, beta, detparindx,
                                                detlink, object$fixed)
        Xrealparval <- secr_reparameterize (realparval, object$detectfn, object$details,
                                            data$mask, data$traps, Dtemp, data$s)
        PIA <- object$design0$PIA[sessnum, ok, 1:data$s, 1:data$K, ,drop=FALSE]
    }
    else {
        realparval  <- secr_makerealparameters (object$design, beta, detparindx,
                                                detlink, object$fixed)
        Xrealparval <- secr_reparameterize (realparval, object$detectfn, object$details,
                                            data$mask, data$traps, Dtemp, data$s)
        PIA <- object$design$PIA[sessnum, ok, 1:data$s, 1:data$K, ,drop=FALSE]
    }
    
    pmix <- secr_getpmix (data$knownclass[ok], PIA, Xrealparval)  ## membership prob by animal
    
    ## unmodelled beta parameters, if needed
    miscparm <- secr_getmiscparm(object$details$miscparm, object$detectfn, object$beta, 
                                 object$parindx, object$details$cutval)
    
    gkhk <- secr_makegk (data$dettype, object$detectfn, data$traps, data$mask, object$details, sessnum, 
                         NE, D, miscparm, Xrealparval, grain, ncores)
    haztemp <- secr_gethazard (data$m, data$binomNcode, nrow(Xrealparval), gkhk$hk, PIA, data$usge)
    
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
#-------------------------------------------------------------------------------

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
            prmat <- allhistpolygonfxi (
                object$detectfn, Xrealparval, haztemp, gkhk$hk, gkhk$H, pimask, PIA, 
                CH, xy, data$binomNcode, grp, data$usge, data$mask,
                pmix, data$maskusage, grain, ncores, object$details$minprob)
        }
        1-prmat
    })
}
#-------------------------------------------------------------------------------

derivedDsurface <- function (object, mask = NULL, sessnum = NULL, groups = NULL) {
    Dx <- function(object, mask, sessnum, selection) {
        D <- secr_predictD(object, mask, group = NULL, session = sessnum, parameter = 'D')
        cellsize <- secr_getcellsize(mask)
        px <- pxi(object, sessnum = sessnum, X = mask)   # dim N x m
        D <- matrix(D, ncol = 1)                         # dim m x 1
        intDp <- px %*% D * cellsize                     # dim N x 1
        out <- mask
        covariates(out)$D.0 <- D * sum(1/intDp[selection])
        class(out) <- c("Dsurface", "mask", "data.frame")
        out
    }
    if (ms(object) && is.null(sessnum)) {
        if (is.null(mask)) mask <- object$mask
        if (ms(mask)) mask <- mask[[sessnum]]
      out <- mapply(derivedDsurface, sessnum = 1:length(object$capthist), 
                    MoreArgs = list(object = object, mask = mask), 
                    SIMPLIFY = FALSE)   
      class (out) <- c("Dsurface", "mask", "list") 
      out
    }
    else {
        if (!object$details$relativeD) 
            stop ("derivedDsurface is for relativeD models")
        if (is.null(object$details$fixedbeta) || is.na(object$details$fixedbeta[1])) 
            stop ("derivedDsurface fixedbeta[1] expected for relativeD but not found")
        if (!(object$link$D %in% c("log", "identity")))
            warning ("derivedDsurface relativeD requires log or identity link")
        if (is.null(sessnum)) sessnum <- 1
        if (is.null(mask)) mask <- object$mask
        
        capthist <- object$capthist
        if (ms(capthist)) capthist <- capthist[[sessnum]]
        grp <- secr_group.factor(capthist, groups)
        ind <- 1:nrow(capthist)
        if (length(ind)>0)
            individuals <- split (ind, grp)
        else
            individuals <-  split (numeric(0), grp) ## list of empty grp levels
        ngrp <- length(individuals)   ## number of groups
        if ( ngrp > 1)
            out <- mapply (
                Dx, 
                individuals, 
                MoreArgs = list(object = object, mask = mask, 
                                sessnum = sessnum), 
                SIMPLIFY = FALSE)
        else {
            if (ngrp == 1) {
                out <- Dx(object, mask = mask, sessnum = sessnum)
            }
            else {
                out <- NULL
            }
        }
        out
    }
}

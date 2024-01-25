## package 'secr' 4.5
## makeStart.R
## 2022-04-02, 2023-12-17

makeStart <- function (start = NULL, parindx, capthist, mask, detectfn, link, 
    details = NULL, fixed = NULL, CL = FALSE, anypoly = FALSE, anytrans = FALSE, 
    alltelem = FALSE, sighting = FALSE) {
    
    ############################################
    # Optionally start from previous fit
    ############################################

    if (inherits(start, c('ipsecr', 'secr'))) {
        
        oldbeta <- coef(start)$beta
        fb <- start$details$fixedbeta
        if (is.null(fb)) fb <- rep(NA, length(oldbeta))
        names(fb)[is.na(fb)] <- start$betanames
        oldbeta <- fullbeta(oldbeta, fb) # matches start$parindx
        if (!is.null(details) && !is.null(details$nsim) && details$nsim > 0) {
            start <- oldbeta    ## chat simulations
        }
        else {
            oldnam <- start$betanames
            start <- mapbeta(start$parindx, parindx, oldbeta, NULL)
            
            if (!is.null(details$miscparm)) {
                nb <- length(start)
                start <- c(start, details$miscparm)
                oldnam <- oldnam[oldnam %in% names(details$miscparm)]
                start[oldnam] <- oldbeta[oldnam]
            }
        }
        if (!is.null(details$fixedbeta)) start <- start[is.na(details$fixedbeta)]
        return(start)
    }
    
    ############################################
    # allow for incomplete details from ipsecr
    ############################################
    
    defaultdetails <- list(
        trace = TRUE,
        binomN = 0,                    
        cutval = 0,
        param = 0,
        unmash = FALSE,
        ignoreusage = FALSE,
        autoini = 1,
        nsim = 0,
        fastproximity = FALSE,
        nmix = 1
    )    
    details <- replace (defaultdetails, names(details), details)
    
    NP <- max(unlist(parindx))

    ############################################
    # Start values (model-specific)
    # 'start' is vector of beta values (i.e. transformed) or a list 
    ############################################
    if (is.null(start) | is.list(start)) {
        unmash <- details$unmash
        if (!ms(capthist)) {
            ch <- capthist
            msk <- mask
        }
        else if (tolower(details$autoini)=='all') {
            ch <- uniquerownames(capthist)
            ch <- join(ch)
            attr(ch, 'n.mash') <- rep(1,length(capthist))
            unmash <- TRUE
            msk <- make.mask(traps(ch), nx=32, buffer=5*spacing(traps(ch)))
        }
        else {
            ch <- capthist[[details$autoini]]
            msk <- mask[[details$autoini]]
        }
        rpsv <- fixed$sigma   ## 2021-03-31
        if (is.null(rpsv)) {
            rpsv <- try(RPSV(ch, TRUE), silent = TRUE)
            if (inherits(rpsv, 'try-error')) rpsv <- NA
        }
        start3 <- list(D = NA, g0 = NA, sigma = NA)
        requireautoini <- (is.null(start) | !all(names(parindx) %in% names(start))) & !alltelem
        if (requireautoini) {
            ## not for signal attenuation
            if (!(detectfn %in% c(9,10,11,12,13)) & !anypoly & !anytrans) {
                memo('Finding initial parameter values...', details$trace)
                # specific to session, do not use anytelem
                if (any(detector(traps(ch))=="telemetry")) {
                    if (all(detector(traps(ch))=="telemetry"))
                        stop("cannot compute start from telemetry data; \n",
                            "set manually or select different session with details autoini")
                    ch <- subset(ch, occasions = detector(traps(ch)) != "telemetry")
                }
                if (nrow(ch)<5)
                    stop ("too few values session ", details$autoini, " to determine start; \n",
                        "set manually or select different session with details autoini")
                tempbinomN <- if (details$binomN==1 || details$fastproximity) 
                    max(unlist(usage(traps(capthist)))) else details$binomN
                start3 <- autoini (
                    capthist = ch, 
                    mask = msk, 
                    binomN = tempbinomN,
                    adjustg0 = details$binomN[1]==0 && !details$fastproximity,
                    ignoreusage = details$ignoreusage,
                    ncores = details$ncores)   ## use ncores set previously
                
                if (any(is.na(unlist(start3)))) {
                    warning ("'secr.fit' failed because initial values not found",
                        " (data sparse?); specify transformed values in 'start'", 
                        call. = FALSE)
                    return (NULL)
                }
                if (unmash & !CL) {
                    nmash <- attr(ch, 'n.mash')
                    n.clust <- length(nmash)
                    start3$D <- start3$D / n.clust
                }
                nms <- c('D', 'g0', 'sigma')
                nms <- paste(nms, '=', round(unlist(start3),5))
                memo(paste('Initial values ', paste(nms, collapse=', ')),
                    details$trace)
            }
            else warning ("using default starting values", call. = FALSE)
        }
        #--------------------------------------------------------------
        # assemble start vector
        ## revised 2014-12-04 to avoid sessions with no detections
        n <- nrow(ch)
        
        default <- list(
            D       = ifelse (is.na(start3$D), 1, start3$D),
            g0      = ifelse (is.na(start3$g0), 0.1, start3$g0),
            lambda0 = -log(1-ifelse (is.na(start3$g0), 0.1, start3$g0)),
            sigma   = ifelse (is.na(start3$sigma), rpsv, start3$sigma),
            z       = 5,
            w       = 10,
            pID     = 0.7,
            noneuc  = 50,
            beta0   = details$cutval + 30,
            beta1   = -0.2,
            sdS     = 2,
            b0      = 2,      ## changed from 15 2010-11-01
            b1      = -0.1,
            pmix    = 0,      ## superceded below
            esa     = ifelse (is.na(start3$D), n / 5, n / start3$D),
            a0      = ifelse (is.na(start3$g0), 0.1 * rpsv^2, start3$g0 *
                    start3$sigma^2) / 10000 * 2 * pi
        )
        if (details$param %in% 4:6) {
            default$sigmak <- default$sigma * default$D^0.5
            default$c <- 0 ## but problems if take log(c)
            default$d <- 0.01 
        }
        
        if (detectfn %in% c(6)) {
            default$w <- default$sigma
            default$sigma <- default$sigma/2
        }
        if (detectfn %in% c(7)) {
            default$z <- default$sigma/5
        }
        if (detectfn %in% c(8, 18, 19)) {
            default$z <- 1    ## cumulative gamma HCG, HVP
        }
        if (anypoly | anytrans) {
            default$D <- 2 * nrow(ch) / masksize(msk)
            default$lambda0 <- sum(ch) / nrow(ch) / ncol(ch)
            ## using usage for binomN=1
            if (any(details$binomN == 1) & 
                    all(detector(traps(ch)) %in% c('polygon','transect'))) {
                usge <- usage(traps(ch))
            }
            default$sigma <- rpsv
        }
        
        if (is.na(default$sigma)) default$sigma <- 20
        getdefault <- function (par) {
            transform (default[[par]], link[[par]])
        }
        
        if (is.list(start)) {
            startnames <- names(start)
            default <- replace(default, startnames, start)
        }
        else startnames <- NULL
        
        #########################################
        start <- rep(0, NP)
        for ( i in 1:length(parindx) ) {
            start[parindx[[i]][1]] <- getdefault (names(parindx)[i]) 
        }
        #########################################

        if ((details$nmix>1) && !('pmix' %in% names(fixed)) && !('pmix' %in% startnames))
            start[parindx[['pmix']][1]] <- clean.mlogit((1:details$nmix)-0.5)[2]
        
        if (detectfn %in% c(12,13))
            start <- c(start, 46,3)    ## muN, sdN
        
        # D/ngrp when figure out where to calculate this
        
        if (sighting & is.null(fixed$pID))
            start[parindx$pID[1]] <- getdefault('pID')
        
        # start vector completed
        #--------------------------------------------------------------
    }
    start
}

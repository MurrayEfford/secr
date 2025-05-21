###############################################################################
## package 'secr'
## chatnk.R
## 2022-11-18, 29
## 2023-04-26 to 2023-05-13
## 2025-05-20 adjustVarD moved to separate file
###############################################################################

chat.nk.sess <- function(object, D, capthist, mask, detpar, nsim, 
                         ncores = NULL, seed = NULL, 
                         verbose = TRUE, type = 'Fletcher', multinomial = FALSE) {
    
    ## c-hat for one session
    
    ## potential development:
    ## check PIA0 for variation over occasions, animals, detectors?
    ## substitute detectpar by lookup of Xrealparval0 with PIA0?
    
    noccasions <- dim(capthist)[2]
    traps <- traps(capthist)
    expected.nk <- Enk(
        D          = D, 
        mask       = mask, 
        traps      = traps, 
        detectfn   = object$detectfn, 
        detectpar  = detpar, 
        noccasions = noccasions,
        binomN     = object$binomN, 
        userdist   = object$details$userdist, 
        ncores     = NULL,
        nrepl      = NULL)   # do not simulate expected
    
    np <- length(object$betanames)
    if (np > (nrow(traps)-1)) stop ("c-hat not estimated when np > K-1")
    
    observed.nk <- nk(capthist)
    if (!is.null(nsim) && nsim >= 1) {
        # simulate a list of 'observed' nk vectors
        onesimnk <- function (r) {
            pop <- sim.popn(D, core = mask, model2D = 'IHP')
            ch <- sim.capthist(
                traps      = traps, 
                popn       = pop, 
                detectfn   = detectfn, 
                detectpar  = detpar, 
                noccasions = noccasions, 
                nsessions  = 1, 
                binomN     = binomN)
            nk(ch)  # individuals per detector
        }
        # optionally define cluster for parallel processing
        if (is.null(ncores)) ncores <- setNumThreads()
        detectfn <- object$detectfn
        binomN <- object$binomN
        if (ncores > 1) {
            clustertype <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
            clust <- parallel::makeCluster(ncores, type = clustertype)
            if (clustertype == "PSOCK") {
                clusterExport(clust, c("sim.popn", "sim.capthist",
                                       "onesimnk", "D", "mask", "traps", "detectfn", 
                                       "detpar", "noccasions", "binomN"
                ), environment())
            }
            parallel::clusterSetRNGStream(clust, seed)
            on.exit(parallel::stopCluster(clust))
            simnk <- parallel::parLapply(clust, 1:nsim, onesimnk)
            
        }
        else {
            set.seed(seed)
            simnk <- lapply(1:nsim, onesimnk)
        }
        
        simchat <- unlist(Fletcher.chat(simnk, expected.nk, np, verbose = FALSE, type, multinomial))
        obschat <- Fletcher.chat(observed.nk, expected.nk, np, verbose = FALSE, type, multinomial)
        list(
            type     = type, 
            observed = observed.nk,
            expected = expected.nk,
            chat     = obschat, 
            nsim     = nsim,
            sim.chat = simchat, 
            p        = 1 - rank(c(obschat, simchat))[1] / (nsim+1))
    }
    else {
        Fletcher.chat(observed.nk, expected.nk, np, verbose, type, multinomial)
    }    
}

chat.nk <- function(object, nsim = NULL, ...) {
    det <- unlist(detector(traps(object$capthist)))
    if (!all(det %in% c('multi','proximity','count'))) {
        stop("chat.nk available only for multi, proximity and count detectors")
    }
    warning ("chat.nk was shown by Efford and Fletcher (2025) to be an inadequate",
             " variance inflation factor. Consider using chat.nj")
    
    if (ms(object)) {
        if (object$CL) {
            Dlist <- lapply(derived(object), '[', 'D', 'estimate')
        }
        else if (object$model$D == ~1) {
            Dlist <- lapply(predict(object), '[', 'D', 'estimate')
        }
        else {
            Dlist <- lapply(covariates(predictDsurface(object)), '[[', 'D.0')
        }
        getdet <- function(x) {
            ok <- rownames(x) %in% c('g0','lambda0','sigma','z')
            as.list(setNames(x[ok,'estimate'], rownames(x)[ok]))
        }
        detparlist <- lapply(predict(object), getdet)
        mapply(chat.nk.sess, 
               D = Dlist, 
               capthist = object$capthist,   # 2023-04-23 previously object$capthist[[1]] 
               mask = object$mask, 
               detpar = detparlist, 
               MoreArgs = list(object = object, nsim = nsim, ...), 
               SIMPLIFY = FALSE)
        
    }
    else {
        
        ##-------------------------------------------------------------
        ## Restrict application
        if (!inherits(object, 'secr'))
            stop ("chat.nk expects fitted secr model")
        if (length(object$fixed)>0)
            stop ("chat.nk does not yet accept models with one or more fixed parameters")
        if (!is.null(object$groups))
            stop ("chat.nk does not yet accept models with grouping")
        if (object$details$nmix>1)
            stop ("chat.nk does not yet handle mixtures")
        if (length(table(object$design0$PIA))>1 || length(table(object$design$PIA))>1)
            stop ("chat.nk does not yet handle varying detection probabilities")
        
        ##-------------------------------------------------------------
        
        ## Density (scalar or length = nrow(object$mask))
        
        if (object$CL) {
            D <- derived(object)['D', 'estimate']
        }
        else if (object$model$D == ~1) {
            pred <- predict(object)
            D <- pred['D', 'estimate']
        }
        else {
            D <- covariates(predictDsurface(object))$D.0
        }
        ##-------------------------------------------------------------
        
        capthist <- object$capthist
        mask <- object$mask
        detpar <- detectpar(object)
        chat.nk.sess (object, D, capthist, mask, detpar, nsim = nsim, ...)
        
    }
}

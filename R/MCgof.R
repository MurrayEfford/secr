################################################################################
## package 'secr'
## MCgof.R

## 2024-07-29 started, based on simulate.secr
## 2024-09-03 flexible naming of statfn output
## 2024-09-03 new args usefxi, useMVN
## 2024-09-03 call sim.onepopn from simulate.R if !usefxi
## 2024-09-09 various edits
## 2024-09-11 use mvtnorm instead of MASS

# Murray Efford and Yan Ru Choo
################################################################################

defaultstatfn <- function (capthist) {
    list(
        yik = apply(abs(capthist), c(1,3), sum),
        yi = apply(abs(capthist), 1, sum),
        yk = apply(abs(capthist), 3, sum)
    )
}

# Calculate Freeman-Tukey statistic, from Choo et al. 2024
defaulttestfn <- function(realised, expected){
    sum((sqrt(realised) - sqrt(expected))^2)
}

simfxiAC <- function (object, bytrap) {

    # -----------------------------------
    # sample one location of each _observed_ animal from its pdf
    fxiList <- fxi(object)
    mask <- object$mask
    sp   <- spacing(mask) # cell size
    m <- sapply(fxiList, sample.int, n = nrow(mask), size = 1, replace=TRUE)
    obspop <- mask[m,] + runif(length(m)*2, -sp/2, +sp/2)  # jitter within cells
    class (obspop) <- c('popn', 'data.frame')
    attr(obspop, 'boundingbox') <- attr(mask, 'boundingbox')
    rownames(obspop) <- names(fxiList)
    
    # -----------------------------------
    # sample unobserved AC from their pdf
    
    # density of unobserved AC
    # this code is also in fxTotal
    Dstar <- predictDsurface(object)
    Dstar <- covariates(Dstar)$D.0
    
    # make a null starting popn
    unobspop <- sim.popn(0, attr(obspop, 'boundingbox'), buffer = 0, covariates = NULL) 
    CH <- object$capthist
    n <- nrow(CH)
    
    Nstar <- sum(Dstar) * attr(mask, 'area')
    if (object$details$distribution == 'poisson') {
        Nstar <- rpois(1, Nstar)
    }

    if (Nstar < n) {
        warning("resampled N less than number observed; unobserved not simulated")
    }
    else if (Nstar > n) {
        detpar <- detectpar(object, byclass = TRUE, bytrap = bytrap) 
        ## if bytrap, detpar[[i]] is trap x parameter dataframe for class i
        if (is.na(detpar[[1]]['pmix'])) detpar <- list(c(detpar, pmix = 1))
        for (i in 1:length(detpar)) {
            pd <- pdot(X = mask, traps = traps(CH), detectfn = object$detectfn,
                       detectpar = detpar[[i]], noccasions = ncol(CH))
            # sample locations of unobserved AC for i-th class
            pop <- sim.popn (
                Nbuffer = (Nstar-n) * detpar[[i]]$pmix,  # number
                D       = Dstar * (1 - pd),              # spatial distribution only
                core    = mask, 
                model2D = 'IHP', 
                Ndist   = 'fixed', 
                covariates = NULL
            )
            if (nrow(pop)>0) unobspop <- rbind(unobspop, pop)
        }
        rownames(unobspop) <- paste0('N', 1:nrow(unobspop))
    }

    # -----------------------------------
    # combine observed, unobserved
    popn <- rbind(obspop, unobspop, renumber = FALSE)
    Nobs <- nrow(obspop)
    Nunobs <- nrow(unobspop)  # varies because Dstar varies
    covariates(popn) <- data.frame(obs = rep(c(TRUE,FALSE), c(Nobs, Nunobs)))
    
    # -----------------------------------
    ## add individual covariates that may be needed by sim.detect
    if (!is.null(object$hcov)) {
        ## for unobs, sample with replacement from original hcov field
        oldhcov <- covariates(CH)[,object$hcov]
        covariates(popn)[[object$hcov]] <- 
            sample(oldhcov, size = nrow(popn), replace = TRUE)
        ## then overwrite obs (known) hcov
        covariates(popn)[1:Nobs,object$hcov] <- covariates(object$capthist)[,object$hcov]
    }
    
    popn
}

MCgof.secrlist <- function (object, nsim = 100, statfn = NULL, testfn = NULL,
                        seed = NULL, ncores = 1, clustertype = c("PSOCK","FORK"), 
                        usefxi = TRUE, useMVN = TRUE, Ndist = NULL, quiet = FALSE, ...)
    
{
    out <- lapply(object, MCgof, 
                  nsim = nsim, statfn = statfn, testfn = testfn,
                  seed = seed, ncores = ncores, clustertype = clustertype, 
                  usefxi = usefxi, useMVN = useMVN, Ndist = NULL, quiet = quiet, ...)
    names(out) <- names(object)
}

MCgof.secr <- function (object, nsim = 100, statfn = NULL, testfn = NULL,
                   seed = NULL, ncores = 1, clustertype = c("PSOCK","FORK"), 
                   usefxi = TRUE, useMVN = TRUE, Ndist = NULL, quiet = FALSE, ...)
    
{
    #---------------------------------------------------------------------------
    ## function to run one replicate
    runone <- function(i=1) {
        # -----------------------------------------------
        # randomize parameter values from MVN
        if (useMVN) {
            # MASS::mvrnorm not consistent between platforms
            # possibly due to dependence LAPACK/BLAS
            # object$fit$par <- MASS::mvrnorm(
            #     n     = 1,
            #     mu    = object$fit$par,
            #     Sigma = object$beta.vcv)
            object$fit$par <- mvtnorm::rmvnorm(
                n     = 1,
                mean  = object$fit$par, 
                sigma = object$beta.vcv)[1,]
        }
        # -----------------------------------------------
        # draw randomized AC
        if (usefxi) {
            popn <- simfxiAC(object, bytrap)
        }
        else {
            Darray <- getDensityArray (predictDsurface(object))
            popn <- sim.onepopn(object, Darray)[[1]]
        }
        
        # -----------------------------------------------
        # simulate detections for these parameters and AC
        simCH <- sim.detect(
            object     = object, 
            popnlist   = list(popn), 
            expected   = TRUE, 
            dropzeroCH = FALSE, 
            renumber   = FALSE)
        
        # -----------------------------------------------
        # expected counts for these parameters and AC
        
        # if usefxi = FALSE the 'expected' count is for a random AC
        # and the comparison may not be meaningful
        
        expCH <- attr(simCH, 'expected')
        
        Tsim <- statfn(simCH)
        Texp <- statfn(expCH)
        
        if (!quiet && ncores == 1) setTxtProgressBar(progressbar, i)
        list(Tsim = Tsim, Texp = Texp, par = object$fit$par, popn = popn)
        
    }   # end of runone
    
    #---------------------------------------------------------------------------
    # construct test statistic for simulation 'capti' and marginal counts 'stat'
    ft <- function (capti, stat) {
        N <- nrow(capti$popn)  # varies for each replicate because D varies
        CH <- object$capthist
        n <- nrow(CH)
        # addzeroCH is internal secr function
        paddedCH <- addzeroCH(CH, N - n, prefix = 'N')
        Tobs <- statfn(paddedCH)
        obs <- testfn(Tobs[[stat]],       capti$Texp[[stat]])
        sim <- testfn(capti$Tsim[[stat]], capti$Texp[[stat]])
        c(Tobs = obs, Tsim = sim, p = sim>obs)  
    }   # end of ft
    #---------------------------------------------------------------------------
    
    ## create default statfn if needed
    if (is.null(statfn)) {
        statfn <- defaultstatfn
    }
    ## create default testfn if needed
    if (is.null(testfn)) {
        testfn <- defaulttestfn
    }
    
    #---------------------------------------------------------------------------
    ##  check input
    if (any(c("bn", "bkn", "bkc", "Bkc") %in% tolower(object$vars)))
        stop ("MCgof works only with binary behavioural responses")
    
    if (!inherits(object,'secr'))
        stop ("requires 'secr' object")
    
    if (object$CL)
        stop ("MCgof not implemented for conditional likelihood")
    
    if (ms(object$capthist))
        stop ("MCgof not implemented for multi-session data")
    
    if (!is.null(object$groups)) {
        stop ("MCgof not implemented for groups")
    }
    
    #---------------------------------------------------------------------------
    
    if (any(detector(traps(object$capthist)) == 'single')) {
        # expected values not available for 'single'
        warning("replacing single detector type with multi - results unreliable")
        detector(traps(object$capthist)) <- 'multi'
    }
    #---------------------------------------------------------------------------
    
    # optionally override distribution
    if (!is.null(Ndist)) {
        Ndist <- match.arg(tolower(Ndist), c('poisson', 'fixed'))
        object$details$distribution <- switch (Ndist, 
                         poisson  = "poisson", 
                         fixed    = "binomial")
    }
    
    clustertype <- match.arg(clustertype)
    
    #---------------------------------------------------------------------------
    ## set random seed
    ## copied from simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    #---------------------------------------------------------------------------
    
    ## action starts here
    
    ptm  <- proc.time()
    
    ## Does model have differences in detection parameters between traps?
    ## this is conditional to save time in pdot
    bytrap <- length(unique(apply(object$design0$PIA, 4, mean)))>1
    
    # ----------------------------------------------------------
    # perform nsim simulations
    ## parallel processing optional - often slower
    
    if (is.null(ncores)) ncores <- setNumThreads()  # default is 1
    if (!quiet && ncores == 1) progressbar <- txtProgressBar(0, nsim, style = 3)
    
    if (ncores > 1) {
        if (.Platform$OS.type != "unix") clustertype <- "PSOCK"
        clust <- makeCluster(ncores, type = clustertype)
        if (clustertype == "PSOCK") {
            clusterExport(clust, c("object", "simfxiAC", "statfn", "bytrap"), 
                          environment())
        }
        clusterSetRNGStream(clust, seed)
        on.exit(stopCluster(clust))
        capt <- parLapply(clust, 1:nsim, runone)
        
    }
    else {
        capt <- lapply(1:nsim, runone)
    }
    
    # ----------------------------------------------------------
    # apply function 'ft' for each count vector 'stat'
    # stats correspond to named components of statfn output list
    stats <- names(statfn(object$capthist))
    all <- lapply(stats, function(x) sapply (capt, ft, stat = x))
    names(all) <- stats
    
    # ----------------------------------------------------------
    # wrap up
    out <- list(
        object   = object,
        nsim     = nsim,
        statfn   = statfn,
        testfn   = testfn,
        all      = all,
        proctime = (proc.time() - ptm)[3])
    class(out) <- 'MCgof'
    
    if (!quiet) {
        if (ncores==1) close(progressbar)
        message ("MCgof for ", nsim, " simulations completed in ", 
                 round((proc.time() - ptm)[3],3), " seconds")
        print(summary(out))
    }
    attr(out, 'seed') <- RNGstate      ## save random seed
    invisible(out)
    
}
################################################################################

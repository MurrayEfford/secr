################################################################################
## package 'secr'
## MCgof.R

## 2024-07-29 started, based on simulate.secr

# Murray Efford and Ran Yu Choo
################################################################################

MCgof <- function (object, nsim = 1, statfn = NULL, testfn = NULL,
                   seed = NULL, ncores = NULL, verbose = FALSE, quiet = FALSE)

{
    #---------------------------------------------------------------------------
    ## function to run one replicate
    runone <- function(i=1) {
        
        # randomize parameter values from MVN
        object$fit$par <- MASS::mvrnorm(1,
                                        mu = object$fit$par, 
                                        Sigma = object$beta.vcv)
        
        # sample one location of each observed animal from its pdf
        fxi.list <- fxi.secr(object)
        mask <- object$mask
        sp   <- spacing(mask) # cell size
        m <- sapply(fxi.list, sample.int, n = nrow(mask), size = 1, replace=TRUE)
        obspop <- mask[m,] + runif(length(m)*2, -sp/2, +sp/2)  # jitter within cells
        class (obspop) <- c('popn', 'data.frame')
        attr(obspop, 'boundingbox') <- attr(mask, 'boundingbox')
        
        # generate locations of unobserved animals
        D <- covariates(fx.total(object))$D.nc
        unobspop <- sim.popn (D = D, core = mask, model2D = 'IHP', 
                              Ndist = 'fixed', covariates = NULL)
        
        # combine observed, unobserved
        popn <- rbind(obspop, unobspop)
        Nobs <- nrow(obspop)
        Nunobs <- nrow(unobspop)  # varies because D varies
        covariates(popn) <- data.frame(obs = rep(c(TRUE,FALSE), c(Nobs, Nunobs)))
        
        ## add individual covariates that may be needed by sim.detect
        if (!is.null(object$hcov)) {
            ## for unobs, sample with replacement from original hcov field
            oldhcov <- covariates(object$capthist)[,object$hcov]
            covariates(popn)[[object$hcov]] <- 
                sample(oldhcov, size = nrow(popn), replace = TRUE)
            ## then overwrite obs hcov
            covariates(popn)[1:Nobs,object$hcov] <- covariates(obspop)[,object$hcov]
        }
        
        simCH <- sim.detect(object, popnlist = list(popn), expected = TRUE, dropzero = FALSE)
        expCH <- attr(simCH, 'expected')
        
        Tsim <- statfn(simCH)
        Texp <- statfn(expCH)
        
        if (!quiet) setTxtProgressBar(pb, i)
            
        list(Tsim = Tsim, Texp = Texp, par = object$fit$par, popn = popn)
        
    }   # end of runone
    #---------------------------------------------------------------------------
    
    
    ## create default statfn if needed
    if (is.null(statfn)) {
        statfn <- function (capthist) {
            list(
                yik = apply(abs(capthist), c(1,3), sum),
                yi = apply(abs(capthist), 1, sum),
                yk = apply(abs(capthist),3, sum)
            )
        }
    }
    ## create default testfn if needed
    if (is.null(testfn)) {
        # Calculate Freeman-Tukey statistic, from Choo et al. 2024
        testfn <- function(realised, expected){
            sum((sqrt(realised) - sqrt(expected))^2)
        }
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

    if (detector(traps(object$capthist)) == 'single') {
        # expected values not available for 'single'
        warning("replacing single detector type with multi - results unreliable")
        detector(traps(object$capthist)) <- 'multi'
    }
    
    #---------------------------------------------------------------------------
    ## set random seed
    ## copied from simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    #---------------------------------------------------------------------------
    
    ## action starts here
    
    ptm  <- proc.time()
    
    ## optional parallel processing - often slower
    
    if (is.null(ncores)) ncores <- setNumThreads()
    if (ncores>1) quiet <- TRUE
    if (!quiet) pb <- txtProgressBar(0, nsim, style = 3)

    if (ncores > 1) {
        clustertype <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
        clust <- makeCluster(ncores, type = clustertype)
        if (clustertype == "PSOCK") {
            clusterExport(clust, c("object", "statfn"), environment())
        }
        clusterSetRNGStream(clust, seed)
        on.exit(stopCluster(clust))
        capt <- parLapply(clust, 1:nsim, runone)
        
    }
    else {
        capt <- lapply(1:nsim, runone)
    }

    ft <- function (capti, stat) {
        N <- nrow(capti$popn)  # varies for each replicate because D varies
        paddedCH <- addzeroCH(object$capthist, N - nrow(object$capthist))
        Tobs <- statfn(paddedCH)
        
        # BUT does sim.detect deliver in matching order?
        
        obs <- testfn(Tobs[[stat]],       capti$Texp[[stat]])
        sim <- testfn(capti$Tsim[[stat]], capti$Texp[[stat]])
        c(Tobs = obs, Tsim = sim, simGTobs = sim>obs)
    }
    
    tests <- c('yik','yi','yk')
    all <- lapply(tests, function(x) sapply (capt, ft, stat = x))
    names(all) <- tests
    
    if (!quiet) {
        close(pb)
        message ("completed in ", round((proc.time() - ptm)[3],3), " seconds")
    }
    
    if (verbose) {
        out <- list(
            statfn   = statfn,
            testfn   = testfn,
            all      = all,
            nsim     = nsim,
            means    = sapply(all, apply, 1, mean), 
            proctime = (proc.time() - ptm)[3])
        class(out) <- 'MCgof'
        out
    }
    else {
        sapply(all, apply, 1, mean)
    }
}
################################################################################

plot.MCgof <- function(x, overlay = NULL, ...) {
    main <- c('individual-detector','individual', 'detector')
    onestat <- function (xy, pnum) {
        lim <- range(xy[1:2,]) * c(0.8,1.2)
        MASS::eqscplot(xy['Tsim',], xy['Tobs',], 
                       xlab = 'simulated', ylab = 'observed',
                       xlim = lim, ylim = lim)
        
        abline(0,1, col = 'red')
        mtext(side=3, paste(main[pnum], " p =", round(mean(xy[3,]),3)))
        # optional overlay of points from another MCgof
        if (!is.null(overlay) && inherits(overlay, 'MCgof')) {
            points(overlay[[pnum+1]], ...)
        }
        
    }
    if (!inherits(x, 'MCgof')) stop ("requires verbose output")
    par(pty = 's', mfrow = c(2,2))
    mapply(onestat, x$all, 1:3)
    invisible()
}
################################################################################

summary.MCgof <- function(object, ...) {
    temp <- object[c('nsim', 'means')]
    class(temp) <- 'summary.MCgof'
    temp
}

print.summary.MCgof <- function (x, ...) {
    cat ('nsim =', x$nsim, '\n')
    print(x$means)
}

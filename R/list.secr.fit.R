################################################################################
## package 'secr'
## list.secr.fit.R
## 2024-02-19 supercedes par.secr.fit
## 2024-07-11 robust to partial failure
################################################################################

## keeping it simple
list.secr.fit.old <- function (..., constant = list(), prefix = "fit", names = NULL) {
    fits <- mapply(secr.fit, ..., MoreArgs = constant, SIMPLIFY = FALSE)
    nfits <- length(fits)
    defaultnames <- paste0(prefix, 1:nfits)
    if (is.null(names)) names <- defaultnames
    else if (length(names) != nfits) {
        warning ("number of names does not equal number of fits")
        names <- defaultnames
    }
    if (!is.null(constant$details$LLonly) && constant$details$LLonly) {
        fits  
    }
    else {
        ok <- sapply(fits, inherits, 'secr')
        if (any(!ok)) {
            if (all(!ok)) stop ("no valid fits") 
            else warning (sum(!ok), " fits failed and were dropped")
        }
        secrlist(fits[ok], names = names[ok])
    }
}

list.secr.fit <- function (..., constant = list(), prefix = "fit", names = NULL,
                           ncores = NULL, bymodel = FALSE) {
    ## Remember start time
    ptm       <- proc.time()
    
    if (bymodel && is.null(ncores)) warning ("specify ncores to assign models to cluster (bymodel)")
    # cluster type; override with options(secr.clustertype = "PSOCK")
    clustertype <- getOption("secr.clustertype", 
            default = if (.Platform$OS.type == "unix" && Sys.getenv("RSTUDIO") != "1") 
                "FORK" else "PSOCK")
    calls <- mapply(list, ..., SIMPLIFY = FALSE)
    nfits <- length(calls)
    
    if (is.null(ncores)) bymodel <- FALSE      # nothing to parallelise over
    if (bymodel) {
        if (ncores > parallel::detectCores())
            warning("ncores exceeds the number of available cores")
        if (ncores > nfits) {
            nclust <- nfits
            ncores <- max(1, trunc(ncores/nfits))
            message("list.secr.fit is using a parallel cluster with ", nclust, " workers, ", ncores, " threads per worker")
        }
        else {
            nclust <- ncores
            ncores <- 1
        }
    }
    else {
        nclust <- 0
    }
    if (bymodel && nclust > 1) {
        constant$ncores <- ncores                  # threads per fit
        clust <- parallel::makeCluster(nclust, type = clustertype)
        on.exit(parallel::stopCluster(clust), add = TRUE)
        parallel::clusterEvalQ(clust, library(secr))
        parallel::clusterExport(clust, "constant", envir = environment())
        one <- function(a) try(do.call(secr::secr.fit, c(a, constant)), silent = TRUE)
        environment(one) <- globalenv()
        fits <- parallel::parLapplyLB(clust, calls, one)
    } else {
        if (!is.null(ncores)) constant$ncores <- ncores  # as secr.fit would; NULL leaves its default
        fits <- lapply(calls, function(a) do.call(secr::secr.fit, c(a, constant)))
    }
    
    if (!is.null(constant$details$LLonly) && constant$details$LLonly) {
    }
    else {
        defaultnames <- paste0(prefix, 1:nfits)
        if (is.null(names)) names <- defaultnames
        else if (length(names) != nfits) {
            warning ("number of names does not equal number of fits")
            names <- defaultnames
        }
        ok <- sapply(fits, inherits, 'secr')
        if (any(!ok)) {
            msg <- vapply(fits[!ok], function(f)
                if (inherits(f, "try-error")) conditionMessage(attr(f, "condition")) else "not a secr object",
                character(1))
            
            if (all(!ok)) stop ("no valid fits: ", msg[1])
            else warning (sum(!ok), " fits failed and were dropped: ",
                          paste(unique(msg), collapse = "; "))
            
        }
        fits <- secrlist(fits[ok], names = names[ok])
    }
    attr(fits, 'proctime') <- (proc.time()-ptm)[3]
    fits
}

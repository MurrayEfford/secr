## 2014-02-11
## 2016-06-04 annotation; memo -> message
## 2017-04-04 prefix argument
## 2017-07-22 LB and save.intermediate arguments (code of Mathias Tobler)
## 2017-07-26 seed default changed from 123 to NULL
## 2017-09-20 distinct names for saved intermediate fits

par.secr.fit <- function (arglist, ncores = 1, seed = NULL, 
                          trace = TRUE, logfile = "logfile.txt", prefix = "fit.", 
                          LB = FALSE, save.intermediate = FALSE) {
    
    .Deprecated("list.secr.fit", package="secr", "par.secr.fit deprecated owing to inefficiency of parallel processing; use list.secr.fit",
                old = as.character(sys.call(sys.parent()))[1L])

    ptm  <- proc.time()
    
    ## 'inherits' causes R to search in enclosing frames
    if (is.character(arglist))
        arglist <- mget(arglist, inherits = TRUE)
    
    ## force 'trace' to common value across all components of arglist
    arglist <- lapply(arglist, function (x) {x$trace <- trace; x})
    
    ## ensure args named
    if (is.null(names(arglist)))
        names(arglist) <- paste0("arg", 1:length(arglist))
    
    ## check for capthist, mask, dframe mentioned by name
    ## objects are exported to the worker processes as required
    getnames <- function(obj = 'capthist') {
        tmpnames <- sapply(arglist, function(x) if (is.character(x[[obj]])) x[[obj]] else '')
        unique(tmpnames)
    }
    data <- c(getnames('capthist'), getnames('mask'),getnames('dframe'),getnames('details'))
    data <- data[nchar(data)>0]
    
    ## default details savecall to FALSE across all components of arglist
    arglist <- lapply(arglist, function (x) {
        if (is.null(x$details))
            x$details <- list(savecall = FALSE)
        else if (!('savecall' %in% names(x$details))) {
            x$details[['savecall']] <- FALSE
        }
        x
    })
    
    ###################################################
    ## Based on code of M. Tobler 2017-07
    run.fit <- function(x) {
        fitname <- attr(x,'name')
        fit <- do.call("secr.fit", x)
        if (save.intermediate){
            assign(fitname, fit)
            save(list = fitname, file = paste0(fitname, ".RData"))
        }
        fit
    }
    
    ## capture names
    nameattr <- function (x,n) {attr(x,'name') <- n; x}
    arglist <- mapply(nameattr, arglist, names(arglist), SIMPLIFY = FALSE)
    ###################################################
    
    if (ncores > 1) {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big', 
                             outfile = logfile)
        clusterSetRNGStream(clust, seed)
        clusterExport(clust, c(data, 'secr.fit'), environment())
        
        ## previously
        ## output <- parLapply(clust, arglist, do.call, what = 'secr.fit')
        
        ###################################################
        ## Based on code of M. Tobler 2017-07
        if (LB)
            output <- clusterApplyLB(clust, arglist, run.fit)
        else 
            output <- clusterApply(clust, arglist, run.fit)
        
        ###################################################
        
        stopCluster(clust)
    }
    else {
        set.seed (seed)
        ## output <- lapply(arglist, do.call, what = 'secr.fit')
        output <- lapply(arglist, run.fit)
    }
    
    ## changed from memo() 2016-06-04
    message(paste('Completed in ', round((proc.time() - ptm)[3]/60,3), ' minutes at ',
                  format(Sys.time(), "%H:%M:%S %d %b %Y"), sep=''))
    
    if (inherits(output[[1]], 'secr')) 
        output <- secrlist(output)
    
    ## apply standard naming convention
    names(output) <- paste0(prefix, names(arglist))
    
    output
}

par.derived <- function (secrlist, ncores = 1, ...) {
    .Deprecated("lapply(secrlist, derived)", package="secr", "par.derived deprecated owing to inefficiency of parallel processing; use lapply(secrlist, derived)",
                old = as.character(sys.call(sys.parent()))[1L])

    if (!inherits(secrlist, 'secrlist'))
        stop("requires secrlist input")
    if (ncores > 1) {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
        output <- parLapply(clust, secrlist, derived, ...)
        stopCluster(clust)
    }
    else {
        output <- lapply(secrlist, derived, ...)
    }
    names(output) <- names(secrlist)
    output
}

par.region.N <- function (secrlist, ncores = 1, ...) {
    .Deprecated("lapply(secrlist, region.N)", package=NULL, "par.region.N deprecated owing to inefficiency of parallel processing; use lapply(secrlist, region.N)",
                old = as.character(sys.call(sys.parent()))[1L])

    if (!inherits(secrlist, 'secrlist'))
        stop("requires secrlist input")
    if (ncores > 1) {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
        output <- parLapply(clust, secrlist, region.N, ...)
        stopCluster(clust)
    }
    else {
        output <- lapply(secrlist, region.N, ...)
    }
    names(output) <- names(secrlist)
    output
}

################################################################################


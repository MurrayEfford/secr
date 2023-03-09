################################################################################
## package 'secr'
## setNumThreads.R
## 2023-03-10 moved from utility.R
################################################################################

setNumThreads <- function (ncores, ...) {
    ## environment variable RCPP_PARALLEL_NUM_THREADS is set by 
    ## RcppParallel::setThreadOptions
    
    ## current is NA if variable not previously set
    current <- as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", ""))
    if (missing(ncores)) ncores <- NULL
    if (is.na(current)) {
        if (is.null(ncores)) {
            ncores <-  min(RcppParallel::defaultNumThreads(), 2)
        }
    }
    else if (is.null(ncores) || (ncores == current)) {
        return(current)   ## no need to change 
    }
    if (ncores > RcppParallel::defaultNumThreads()) 
        stop("requested ncores exceeds number available")
    if (ncores<1)
        stop ("specified ncores < 1")
    ncores <- min(ncores, RcppParallel::defaultNumThreads())
    RcppParallel::setThreadOptions(ncores, ...) 
    return(as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "")))
}
##############################################################################

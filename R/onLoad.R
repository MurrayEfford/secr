###############################################################################
## package 'secr'
## onLoad.R
## 2020-02-21, 2021-03-13, 2021-04-13, 2021-05-11, 2021-05-17, 2024-03-02
###############################################################################

.onLoad <- function (libname, pkgname) {
    ## also sets environment variable RCPP_PARALLEL_NUM_THREADS
    defaultncores <- RcppParallel::defaultNumThreads()
    if (defaultncores == 1) {
        RcppParallel::setThreadOptions(1)
    }
    else {
        RcppParallel::setThreadOptions(2)
    }
    
    ## raster and terra are loaded only if needed (loading them is slow); their methods for
    ## mask and Dsurface are registered now if they are already loaded, otherwise when they
    ## are loaded (see raster.R)
    registerWhenLoaded("raster", registerRasterMethods)
    registerWhenLoaded("terra", registerRastMethods)

    ## 2024-03-02
    ## avoid BLAS default as not needed and may trigger valgrind problem
    ## see R-package-devel 2024-03-01
    options(matprod = "internal")
    
    ## TESTING ONLY 2021-05-17
    ## RcppParallel::setThreadOptions(7)
    
    ## following advice of Kevin Ushey 2020-03-18, 2021-05-04
    ## to avoid ASAN/UBSAN errors from parallelFor
    ## use this in CRAN tests see test-initial.R
    ## Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
    ## otherwise
    ## Sys.setenv(RCPP_PARALLEL_BACKEND = "tbb")
    
}

## .onLoad is preferred if actions are required for single functions 
## that may be called without attaching package
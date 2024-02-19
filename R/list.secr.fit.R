################################################################################
## package 'secr'
## list.secr.fit.R
## 2024-02-19 supercedes par.secr.fit
################################################################################

list.secr.fit <- function (arglist, trace = FALSE, prefix = "", seed = NULL, 
                           save.intermediate = FALSE) {
    
    ptm  <- proc.time()
    
    ## ensure args named
    if (is.null(names(arglist))) {
        names(arglist) <- paste0("arg", 1:length(arglist))
    }
    
    ## 'inherits' causes R to search in enclosing frames
    if (is.character(arglist)) {
        arglist <- mget(arglist, inherits = TRUE)
    }
    
    ## force 'trace' to common value across all components of arglist
    arglist <- lapply(arglist, function (x) {x$trace <- trace; x})
    
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
    
    set.seed (seed)
    output <- mapply(run.fit, arglist, SIMPLIFY = FALSE)
    if (inherits(output[[1]], 'secr')) {
        output <- secrlist(output)
    }
    
    message(paste('Completed in ', round((proc.time() - ptm)[3]/60,3), ' minutes at ',
                  format(Sys.time(), "%H:%M:%S %d %b %Y"), sep=''))
    
    
    ## apply standard naming convention
    names(output) <- paste0(prefix, names(arglist))
    
    output
}

# tmp <- list.secr.fit(args, prefix='')
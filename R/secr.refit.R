#############################################################################
## package 'secr'
## refit.R
## 2024-11-23, 2025-01-02, 2025-03-17
#############################################################################

secr.refit <- function (object, ...) {
    cl <- match.call(expand.dots = TRUE)
    newargs <- list(...)
    
    if (is.character(object)) {
        # input from saved progress file
        type <- "saved progress file"
        object <- readRDS(object)
        # use last log entry
        # coefficients start in column 4, after eval, loglik and time
        beta <- unlist(tail(attr(object, 'log'),1)[-(1:3)])
        object$fit <- list(par = beta)
    }
    else {
        # assume secr or list
        type <- "fitted model"
        if ("method" %in% names(newargs) && newargs$method == "none") {
            type <- paste(type, "(variance only")
        }
    }
    
    ## list of input arguments saved in object
    ## omit dframe, verify, biasLimit
    
    argnames <- c("capthist", "model", "mask", "CL", "detectfn", 
     "link", "fixed", "timecov", "sessioncov", "hcov", "groups", "details", 
     "method", "start")
    args <- object[argnames]
    
    # incidental arguments not saved by secr.fit
    args$ncores <- setNumThreads()
    args$verify <- TRUE
    args$trace <- TRUE
    args$biasLimit <- 0.01
    
    # first combine details lists so can supply a partial replacement
    if (!is.null(newargs$details)) {
        args$details <- replace (args$details, names(newargs$details), newargs$details)
        newargs$details <- NULL
    }
        
    # check for change in structure
    structargs <- c("model","CL","detectfn","link", "fixed", "hcov", "groups")
    newstructure <- any(sapply(structargs, function(a) {
        !is.null(newargs[[a]]) && args[[a]] != newargs[[a]]
        }
    ))
    
    # update arguments
    args <- replace (args, names(newargs), newargs)
    
    ## raise binomN from details to full argument
    args$binomN <- object$details$binomN  
    
    ## start at exact beta coefficients
    args$start  <- object$fit$par  
    if (newstructure) {
        
        if(!is.null(object$details$relativeD) && object$details$relativeD) {
            stop ("cannot change structure of relative density model in secr.refit")
        }
            
        if (inherits(object, 'secr')) {
            # rely on makeStart() in secr.fit()
            args$start <- object
            warning("model structure has changed")
        }
        else {
            stop ("cannot change model structure in secr.refit when input is a list")
        }
    }
    
    ## relativeD requires beta0 (overwritten later, but must be present in start vector)
    if (!is.null(object$details$relativeD) && object$details$relativeD) {
        args$start <- c(NA, args$start)
    }
    
    out <- with(args, {
        secr.fit(capthist, model, mask, CL = CL, detectfn = detectfn,
                 binomN = binomN, start = start,link = link, fixed = fixed,
                 timecov = timecov, sessioncov = sessioncov, hcov = hcov,
                 groups = groups, details = details, method = method,
                 verify = verify, biasLimit = biasLimit, trace = trace,
                 ncores = ncores)
    })
    if (args$details$savecall) {
        out$call <- cl
    }
    out
}

# note
# tmp <- mapply(secr.refit, MoreArgs = list(object = secrdemo.0), detectfn = c('HEX','HHN'), SIMPLIFY = FALSE)
# lapply(tmp, predict)

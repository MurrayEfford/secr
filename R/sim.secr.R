############################################################################################
## package 'secr'
## sim.secr.R

## 2019-10-05 converted to cpp
## 2019-10-06 yet to fix: groups?
## 2019-10-06 excludes: knownclass (hcov, h2)
## 2019-10-06 excludes: mixed detector types
## 2019-10-06 excludes: individual covariates
## 2019-10-06 excludes: telemetry
## 2020-11-08 bad rownames for sim.detect when renumber = FALSE
## 2024-07-29 sim.detect.R split out as separate file
## 2024-07-29 simulate.R split out as separate file

############################################################################################

sim.secr <- function (object, nsim = 1, extractfn = function(x)
    c(deviance=deviance(x), df=df.residual(x)), seed = NULL,
    maxperpoly = 100, data = NULL, tracelevel = 1, hessian =
    c('none','auto','fdHess'), start = "true", ncores = NULL, ...) {

## parametric bootstrap simulations based on a fitted secr object
## using simulate.secr (see simulate.R)
    
## extractfn is a required function to extract values from an secr fit
## it should return a vector of named values that does not vary in length
    
## 'hessian' overrides value in object$details
## set hessian = 'auto' if extractfn requires variance-covariance matrix

    hessian <- match.arg(hessian)
    cl <- match.call(expand.dots = TRUE)
    cl <- paste(names(cl)[-1],cl[-1], sep=' = ', collapse=', ' )
    cl <- paste('sim.secr(', cl, ')')

    if (is.null(extractfn)) extractfn <- trim
    test <- extractfn(object, ...)

    if (is.numeric(test)) {
        n.extract <- length(test)
        if (n.extract<=0)
            stop ("invalid 'extractfn'")
    }
    detectnames <- names(object$design0[[1]])   ## names of real detection parameters
    details <- replace(object$details, 'hessian', hessian)
    tracelevel <- as.integer(tracelevel)
    details$trace <- tracelevel > 1
    min.detections <- 1
    .localstuff$iter2 <- 0   

    if (any(detector(traps(object$capthist)) == 'single')) {
        secr_memo ('multi-catch likelihood used for single-catch traps', tracelevel>0)
    }
    
    if (is.null(data)) {
        secr_memo ('sim.secr simulating detections...', tracelevel>0)
        ## use multiple cores only for model fitting 2015-12-02
        data <- simulate(object, nsim = nsim, seed = seed, maxperpoly = maxperpoly)
    }
    else {
        if (any(class(data) != c('secrdata', 'list')))
            stop("invalid data")
    }
    fitmodel <- function (sc) {
        ## i <<- i+1
        ## use second counter so as not to interfere with secrloglik
        .localstuff$iter2 <- .localstuff$iter2 + 1   ## 2016-01-09
        if (tracelevel>0)
            secr_memo (paste('sim.secr fitting replicate', .localstuff$iter2, '...'), TRUE)
        nc <-  sum(counts(sc)$'M(t+1)'[,'Total'])
        if (nc >= min.detections) {
            tempfit <- suppressWarnings( secr.fit(sc, model = object$model, mask = object$mask,
                CL = object$CL, detectfn = object$detectfn, binomN = details$binomN,
                start = start, link = object$link, fixed = object$fixed,
                timecov = object$timecov, sessioncov = object$sessioncov,
                groups = object$groups, dframe = object$dframe, details = details,
                method = object$fit$method, verify = FALSE, biasLimit = NA,
                ncores = ncores) )
            extractfn(tempfit, ...)
        }
        else if (is.list(test)) list() else rep(NA, n.extract)
    }
    
    if (!is.null(start) && !is.numeric(start)) start <- secr_complete.beta(object)
    
    output <-  lapply (data, fitmodel)

    if (is.numeric(test)) {
        output <- do.call(rbind, output)
        output <- data.frame(output)
    }
    else {
        class(output) <- c('secrlist', 'list')
    }

    attr(output,'seed') <- attr(data,'seed')
    attr(output,'call') <- cl
    attr(output,'extractfn') <- extractfn
    output
}
############################################################################################


#############################################################################
## package 'secr'
## summary.traps.R
## 2019-11-10 moved from methods.R
#############################################################################

summary.traps <- function(object, getspacing = TRUE, covariates = FALSE, ...) {
    
    #### NEED TO HANDLE CLUSTER, CLUSTERTRAP 2011-04-12
    if (ms(object)) lapply(object, summary.traps, getspacing = getspacing, ...)
    else {
        if (is.null(object$x) | is.null(object$y))
            stop ("not a valid 'traps' object")
        nd <- nrow(object)
        np <- NA
        sumusage <- NULL
        rangeusage <- NULL
        sumcovar <- NULL
        
        if (all(detector(object) %in% .localstuff$polydetectors)) {
            spacing <- NA
        }
        else {
            spacing <- spacing(object, recalculate = getspacing)
            
            # unclear why following is not for polydetectors?  FAULTY 2012-12-18
            if (is.factor(covariates(object))) {
                susage <- by(usage(object), covariates(object), function(y) apply(y,2,sum))
                sumusage <- matrix(unlist(susage), byrow = T, nrow = length(susage))
                dimnames(sumusage) <- list(levels(covariates(object)), names(susage[[1]]))
            }
            else if (!is.null(usage(object)))  {
                sumusage <- apply(usage(object)>0, 2, sum)
                nocc <- ncol(usage(object))
                rangeusage <- matrix(apply(usage(object), 2, range), ncol=nocc)
                dimnames(rangeusage) <- list(c('min','max'), 1:nocc)
            }
            tempcovar <- covariates(object)
            
            if (!is.null(tempcovar) && covariates) {
                if (ncol(tempcovar)>0 && nrow(tempcovar)>0) 
                    sumcovar <- summary(tempcovar, ...)
            }
        }
        ## defaults
        area <- NA
        totallength <- NA
        np = ndetector(object)
        spacex <- NA
        spacey <- NA
        xrange <- range(object$x)
        yrange <- range(object$y)
        if (all(detector(object) %in% c('polygon', 'polygonX')))
            area <- searcharea(object)
        if (all(detector(object) %in% c('transect', 'transectX')))
            totallength <- transectlength(object)
        tvc <- timevaryingcov(object)
        telemdet <- if (any(detector(object) %in% 'telemetry')) 1 else 0
        temp <- list (
            detector = detector(object),
            ndetector = nd - telemdet,
            npart = np,
            xrange = xrange,
            yrange = yrange,
            spacing = spacing,
            area  = area,
            totallength = totallength,
            markocc = markocc(object),
            telemetrytype = telemetrytype(object),
            rangeusage = rangeusage,
            covar = sumcovar,
            tvc = tvc
        )
        
        class(temp) <- 'summary.traps'
        temp
    }
}
###############################################################################

print.summary.traps <- function (x, terse = FALSE, ...) {
    
    #### NEED TO HANDLE CLUSTER, CLUSTERTRAP 2011-04-12
    
    if (!terse)
        cat ('Object class      ', 'traps', '\n')
    detectorvector <- x$detector
    
    if (length(unique(detectorvector))==1)
        detectorvector <- detectorvector[1]
    cat ('Detector type     ', newstr(x$detector), '\n')
    if (any(x$detector %in% c('telemetry'))) {
        cat ('Telemetry type    ', x$telemetrytype, '\n')
    }
    if (all(x$detector %in% c('polygon', 'polygonX'))) {
        cat ('Number vertices   ', x$ndetector-x$npart, '\n')  ## assume each polygon closed
        cat ('Number polygons   ', x$npart, '\n')
        cat ('Total area        ', sum(x$area), 'ha \n')
    }
    else if (all(x$detector %in% c('transect', 'transectX'))) {
        cat ('Number vertices   ', x$ndetector, '\n')
        cat ('Number transects  ', x$npart, '\n')
        cat ('Total length      ', x$totallength, 'm \n')
    }
    else if (!all(x$detector %in% c('telemetry'))) {
        cat ('Detector number   ', x$ndetector, '\n')
        cat ('Average spacing   ', x$spacing, 'm \n')
    }
    if (!all(x$detector %in% c('telemetry'))) {
        cat ('x-range           ', x$xrange, 'm \n')
        cat ('y-range           ', x$yrange, 'm \n')
    }
    if (!is.null(x$rangeusage)) {
        cat ('\n Usage range by occasion\n')
        print(x$rangeusage, ...)
    }
    if (!is.null(x$markocc)) {
        cat ('\nMarking occasions\n')
        mo <- matrix(as.numeric(x$markocc), byrow=T, nrow = 1,
                     dimnames = list('',1:length(x$markocc)))
        print(mo, ...)
    }
    if (!terse) {
        if (!is.null(x$covar)) {
            cat ('\n')
            cat ('Summary of covariates', '\n')
            print(x$covar, ...)
        }
        if (!is.null(x$usage)) {
            cat ('Usage>0 by occasion', '\n')
            print(x$usage, ...)
        }
        if (!is.null(x$tvc)) {
            cat ('Time-varying covariate(s) name : columns', '\n')
            for (i in 1:length(x$tvc)) cat(names(x$tvc)[[i]], ':', x$tvc[[i]], '\n')
        }
    }
}

###############################################################################

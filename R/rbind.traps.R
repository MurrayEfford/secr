#############################################################################
## package 'secr'
## rbind.traps.R
## Last changed 2017-03-25 checkdet bug: erroneous warning
## 2017-10-25 name suffix option
#############################################################################

rbind.traps <- function (..., renumber = TRUE, addusage, checkdetector = TRUE, suffix = TRUE ) {
    # combine 2 or more traps objects
    # what if multi-session?
    allargs <- list(...)
    check <- function (x) {
        result <- FALSE
        if (!is(x,'traps'))
            stop ("all arguments must be 'traps' objects")
        if (is.null(covariates(x)) != is.null(covariates(allargs[[1]]) ))
            stop ("covariates must be provided for all or none")
        if (is.null(timevaryingcov(x)) != is.null(timevaryingcov(allargs[[1]]) ))
            result <- TRUE
        result
    }
    checkdetect <- function (x) {
        identical(detector(x), detector(allargs[[1]]))
    }
    if (length(allargs) <= 1) {
        if (length(allargs)==1 & ms(allargs[[1]])) {
            allargs <- allargs[[1]]
            temp <- do.call(rbind.data.frame, allargs)
        }
        else {
            if (inherits(allargs[[1]], 'traps'))
                return(allargs[[1]])
            else stop ("bad input to rbind.traps")
        }
    }
    else {
        tvc <- sapply (allargs, check)
        if (any(tvc))
            warning("time-varying covariates differ; using first")
        if (checkdetector) {
            dvar <- sapply(allargs, checkdetect)
            ## if (any(dvar)) 2017-03-27 bug fixed
            if (any(!dvar))
                    warning("inputs differ in detector type; using first")
        }
        temp <- rbind.data.frame(...)
    }
    class(temp) <- c('traps', 'data.frame')
    detector(temp) <- detector(allargs[[1]])
    ## code for originating traps object
    ## uses 'polyID = transectID' behaviour of polyID()
    if (all(detector(temp) %in% .localstuff$polydetectors)) {
        oldtrapsID <- rep(1:length(allargs), sapply(allargs, nrow))
        newpolyID <- unlist(lapply(allargs, function(x) as.character(polyID(x))))
        newpolyID <- factor(paste (oldtrapsID,newpolyID, sep='.'))
        polyID(temp) <- newpolyID
    }
    
    ## covariates 2010 07 05, 2010-08-28
    tempcov <- lapply(allargs, covariates)
    common <- Reduce(intersect, lapply(tempcov, names))
    tempcov <- lapply(tempcov, function(x) x[,common, drop = FALSE])
    covariates(temp) <- do.call(rbind, tempcov)
    timevaryingcov(temp) <- timevaryingcov(allargs[[1]])
    ## clusters  2011-04-12
    tempclus <- lapply(allargs, clusterID)
    if (!is.null(tempclus[[1]])) {
        tempclus <- lapply(tempclus, function(x) as.numeric(as.character(x)) )
        clusterID(temp) <- do.call(c, tempclus)
        clustertrap(temp) <- unlist(lapply(allargs, clustertrap))
    }
    else {
        clusterID(temp) <- NULL
        clustertrap(temp) <- NULL
    }
    
    ## usage
    ## tweaked 2013-01-16
    tempusage <- lapply(allargs, usage)
    if (any(!sapply(tempusage, is.null))) {
        nocc <- unique(unlist(sapply(tempusage, ncol)))
        nr <- unlist(sapply(allargs, nrow))
        if (length(nocc)>1)
            warning ("varying number of occasions; using maximum")
        nocc <- max(nocc)
        fillmissing <- function(x, nr) {
            flush.console()
            if(is.null(x))
                x <- matrix(1, nrow = nr, ncol = nocc)
            else {
                if (ncol(x) < nocc) {
                    x <- cbind(x, matrix(0, nrow = nr, ncol = nocc-ncol(x)))
                }
            }
            ## 2016-01-07 temporary fix
            ## dimnames(x)[[2]] <- 1:ncol(x)
            dimnames(x) <- list(NULL, 1:ncol(x))
            x
        }
        for (i in 1: length(tempusage)) {
            tempusage[[i]] <- fillmissing(tempusage[[i]], nr[i])
        }
        usage(temp) <- do.call(rbind, tempusage)
    }
    ## 2014-03-19
    else if (!missing(addusage)) {
        if (length(addusage) == 1)
            addusage <- rep(addusage, length(allargs))
        if ((length(addusage) != length(allargs)) | any(addusage<1))
            stop ("invalid occasion numbers in addusage")
        makeusage <- function(nr,nc) {
            tmp <- matrix(0, nrow = nr, ncol = nocc)
            tmp[,1:nc] <- 1
            tmp
        }
        nocc <- max(addusage)
        nt <- sapply(allargs, nrow)
        tempusage <- mapply(makeusage, nt, addusage)
        usage(temp) <- do.call(rbind, tempusage)
    }
    
    tempmarkocc <- lapply(allargs, markocc)
    mo <- !sapply(tempmarkocc, is.null)
    if (any(mo)) {
        if (!all(mo))
            warning ("using first non-null markocc")
        wmo <- which(mo)
        mo1 <- tempmarkocc[[wmo[1]]]
        if (length(wmo) > 1) {
            for (i in 2:length(wmo))
                if (any(tempmarkocc[[i]] != mo1))
                    stop ("conflicting markocc")
        }
        markocc(temp) <- mo1
    }
    
    tn <- sapply(allargs, row.names, simplify=F)
    
    
    if (renumber) {
        if (all(detector(temp) %in% .localstuff$polydetectors)) {
            if (!is.null(polyID(temp)))  ## polyID also stores transectID
                polyID(temp) <- factor(as.numeric(polyID(temp)))
            temp <- renamepolyrows(temp)    ## see read.traps
        }
        else
            row.names(temp) <- 1:nrow(temp)
    }
    else {
        if (any(duplicated(unlist(tn)))) {
            if (!suffix) warning("trap ID suffix added to avoid duplication")
            suffix <- TRUE
        }
        if (suffix) {
            for (i in 1:length(tn)) tn[[i]] <- paste(tn[[i]],i,sep='.')
        }
        row.names(temp) <- unlist(tn)
    }
    #    }
    if (!is.null(usage(temp)))
        if (nrow(usage(temp))>0)
            row.names(usage(temp)) <- row.names(temp)
    
    if (!is.null(covariates(temp))) {
        if (nrow(covariates(temp))>0)
            row.names(covariates(temp)) <- row.names(temp)
    }
    
    ## 2012-08-07, 2012-08-24
    if (all(detector(temp) %in% .localstuff$pointdetectors)) {
        spacing(temp) <- spacing(temp, recalculate = TRUE)
    }
    
    temp
}
###############################################################################

trapmerge <- function (trapsM, trapsR, markocc){
    if (!is.null(usage(trapsM)) | !is.null(usage(trapsR)))
        warning("discarding existing usage data")
    KM <- nrow(trapsM)
    KR <- nrow(trapsR)
    S <- length(markocc)
    tmpmarkocc <- as.numeric(markocc>0) ## treat -1 and 0 alike as sighting occasions
    usage(trapsM) <- matrix(tmpmarkocc, byrow = TRUE, nrow = KM, ncol = S)
    usage(trapsR) <- matrix(1-tmpmarkocc, byrow = TRUE, nrow = KR, ncol = S)
    trps <- rbind(trapsM, trapsR, checkdetector = FALSE)
    detector(trps) <- ifelse(markocc>0, detector(trapsM)[1], detector(trapsR)[1])
    ## Note: nrow(trps) == KM + KR
    markocc(trps) <- markocc
    trps
}
# gridM <- make.grid(nx = 3, ny = 3, detector = "multi")
# gridR <- make.grid(nx = 5, ny = 5, detector = "proximity")
# combined <- trapmerge(gridM, gridR, c(1,1,0,0,0))
# summary(combined)
# usage(combined)[1:20,]  ## show usage for first 20 detectors of 34 in the combined layout

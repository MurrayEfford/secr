###############################################################################
## package 'secr'
## addTelemetry.R
## telemetry data manipulation
## 2013-11-18 addTelemetry modified for multi detectors; retains cov if no zerohist
## 2017-01-06 addTelemetry radically revised
## 2017-01-11 read.telemetry adjusted for direct input of telemetry using make.capthist
## moved from telemetry.R 2016-12-28
## 2017-01-27 revised addTelemetry
## 2017-03-25 check for compatible covariates in addTelemetry
## 2017-04-10 bug fixed in preceding check
## 2018-09-30 extended to handle mark-resight data
###############################################################################

## add telemetry occasion(s) to follow existing occasions

addTelemetry <- function (detectionCH, telemetryCH, 
                          type = c('concurrent','dependent','independent'), 
                          collapsetelemetry = TRUE,
                          verify = TRUE) {
    
    ## combine capture histories from telemetry and hair snags etc.
    if (ms(detectionCH) | ms(telemetryCH)) {
        if (!(ms(detectionCH) & ms(telemetryCH)))
            stop ("both detectionCH and telemetryCH should be single- or multi-session")
        if (length(detectionCH) != length(telemetryCH)) {
            if (!all(names(telemetryCH) %in% names(detectionCH)))
                stop ("some telemetry sessions not present in detection data")
            warning ("detectionCH and telemetryCH have differing number of sessions")
            nullTe <- array(dim=c(0,1,1))
            class(nullTe) <- "capthist"
            traps(nullTe) <- traps(telemetryCH[[1]])
            telemetryxy(nullTe) <- NULL
            tmp <- rep(list(nullTe),length(detectionCH))
            class(tmp) <- class(telemetryCH)
            session(tmp) <- session(detectionCH)
            telemetryCH <- replace (tmp, names(telemetryCH), telemetryCH)
        }
        CH <- mapply(addTelemetry, detectionCH, telemetryCH, MoreArgs = 
                         list(type = type, collapsetelemetry = collapsetelemetry, 
                              verify = FALSE))  ## delay verify
        class(CH) <- c('capthist', 'list')
        if (verify) verify(CH)
        CH
    }
    else {
        if (!all(detector(traps(telemetryCH))=='telemetry'))
            stop ("telemetryCH should be of uniform detector type 'telemetry'")
        if (is.null(telemetryCH) || nrow(telemetryCH)==0) {
            return(detectionCH)
        }
        else {
            type <- match.arg(type)
            xylist <- telemetryxy(telemetryCH)
            capdet <- secr_expanddet(detectionCH)
            if (collapsetelemetry)
                teldet <- 'telemetry'
            else
                teldet <- secr_expanddet(telemetryCH)
            
            if (type == "independent") 
                OK <- rep(NA, nrow(telemetryCH))        
            else 
                OK <- match(names(xylist), row.names(detectionCH))
            telemOnly <- sum(is.na(OK))
            
            olddim <- dim(detectionCH)
            
            if (collapsetelemetry) {
                telemCH <- apply(telemetryCH,c(1,3),sum)
                telemCH <- array(telemCH, dim=c(nrow(telemCH),1,1))
            }
            else {
                telemCH <- telemetryCH
            }
            
            #################################
            ## start new CH with all old detections
            newdim <- olddim + c(telemOnly, ncol(telemCH), 1)
            newCH <- array(0, dim = newdim)
            if (olddim[1]>0) {
                newCH[1:olddim[1], 1:olddim[2], 1:olddim[3]] <- detectionCH
            }    
            dimnames(newCH)[[2]] <- 1:ncol(newCH)
            
            class(newCH) <- 'capthist'
            
            #################################
            ## add telemetry detections
            ## if (telemOnly>0)
            if (type == 'independent') {
                names(xylist) <- paste0('T',names(xylist))
                dimnames(newCH)[[1]] <- c(rownames(detectionCH), names(xylist))
            }
            else {
                dimnames(newCH)[[1]] <- c(rownames(detectionCH), names(xylist)[is.na(OK)])
            }
            telemetrd <- match(names(xylist), rownames(newCH))
            newCH[telemetrd,  (olddim[2]+1):newdim[2], newdim[3]] <- telemCH
            
            #################################
            ## consistent order
            if (nrow(newCH)>0) {
                rows <- rownames(newCH)
                xynames <- names(xylist)
                nch <- max(nchar(c(rows, xynames)))
                rows <- str_pad(rows, nch)
                xynames <- str_pad(xynames, nch)
                rownames(newCH) <- rows
                names(xylist) <- xynames
                rows <- sort(rows)
                newCH[,,] <- newCH[rows,,]
                rownames(newCH) <- rows
                xylist <- xylist[rows[rows %in% xynames]]
            }
            #################################
            ## traps attribute of newCH
            oldtraps <- traps(detectionCH)
            nulltraps <- subset(oldtraps,1)  ## duplicate first as dummy location for telemetry
            ## no need for message about differing detector types
            traps(newCH) <- suppressWarnings(rbind(oldtraps, nulltraps))
            
            #################################
            ## detector attributes
            detector(traps(newCH)) <- c(capdet, teldet)
            newusge <- matrix(0, nrow=nrow(traps(newCH)), ncol = ncol(newCH))
            oldusge <- usage(detectionCH)
            if (is.null(oldusge)) 
                oldusge <- matrix(1, nrow=nrow(oldtraps), ncol=ncol(detectionCH))
            ncold <- ncol(oldusge)
            newusge[1:nrow(oldusge), 1:ncold] <- oldusge
            newusge[nrow(newusge), (ncold+1) : ncol(newusge)] <- 1
            usage(traps(newCH)) <- newusge
            telemetrytype(traps(newCH)) <- type
            
            #################################
            ## covariates
            telemcov <- covariates(telemetryCH)
            if (is.null(telemcov) || ncol(telemcov)==0)
                zerohistcov <- matrix(nrow=0,ncol=0)
            else
                zerohistcov <- telemcov[is.na(OK),,drop = FALSE]
            dethistcov <- covariates(detectionCH)
            if ((is.null(zerohistcov) || (nrow(zerohistcov)==0)) &&
                (is.null(dethistcov) || (nrow(dethistcov)==0))) {
                covariates(newCH) <- NULL
            }
            else {
                diffcov <- (is.null(zerohistcov) != is.null(dethistcov)) ||
                    (ncol(dethistcov) != ncol(zerohistcov)) || 
                    !all(names(dethistcov) == names(zerohistcov))
                
                if (diffcov || !all(names(zerohistcov) == names(dethistcov))) {
                    warning ("covariates in telemetryCH do not match detectionCH",
                             " so covariates discarded")
                    covariates(newCH) <- NULL
                }
                else {
                    covariates(newCH) <- rbind(dethistcov, zerohistcov)
                }
            }
            
            #################################
            ## sightings
            if (!is.null(markocc(oldtraps))) {
                markocc(traps(newCH)) <- c(markocc(oldtraps), 1)
                refillT <- function (T) {
                    if (is.matrix(T)) {
                        newT <- matrix(0, nrow=nrow(traps(newCH)), ncol = ncol(newCH))
                        newT[1:nrow(oldtraps), 1:ncol(detectionCH)] <- T
                    }
                    else newT <- T   ## assume summed or NULL
                    newT
                }
                Tu(newCH) <- refillT(Tu(detectionCH))
                Tm(newCH) <- refillT(Tm(detectionCH))
                Tn(newCH) <- refillT(Tn(detectionCH))
            }
            
            #################################
            ## telemetryxy
            telemetryxy(newCH) <- xylist
            
            if (verify) verify(newCH)
            newCH
        }
    }
}
############################################################################################

read.telemetry <- function (file = NULL, data = NULL, covnames = NULL, verify = TRUE, ...) {
    
    fmt <- 'XY'
    inflation <- 1e-8
    nvar <- 5
    
    if (is.null(data)) {
        ## input from text file
        if (is.null(file))
            stop ("must specify either file or data")
        dots <- match.call(expand.dots = FALSE)$...
        
        if (length(file) != 1)
            stop ("requires single 'file'")
        
        filetype <- function(x) {
            nx <- nchar(x)
            tolower(substring(x, nx-3, nx))
        }
        
        countargs <- dots[names(dots) %in% names(formals(count.fields))]
        if (filetype(file) == '.csv')
            countargs$sep <- ','
        countargs$file <- file
        nfield <- max(do.call(count.fields, countargs))
        colcl <- c('character','character',NA,NA,NA, rep(NA,nfield-nvar))
        defaultargs <- list(sep = '', comment.char = '#')
        if (filetype(file)=='.csv') defaultargs$sep <- ','
        captargs <- secr_replacedefaults (defaultargs, list(...))
        captargs <- captargs[names(captargs) %in% names(formals(read.table))]
        capt <- do.call ('read.table', c(list(file = file, as.is = TRUE,
                                              colClasses = colcl), captargs) )
        
    }
    else {
        capt <- data
    }
    
    ## let's be clear about this...
    names(capt)[1:5] <- c('Session','AnimalID','Occ','x','y')
    ## X, Y must be numeric 2016-01-12
    if (any(is.na(capt[,1:nvar])))
        stop ("missing values not allowed")
    
    # cannot have occasions with no detections    
    noccasions <- max(capt$Occ)
    
    readtraps <- function(capt) {
        trps <- t(apply(capt[,4:5],2,mean))
        trps <- as.data.frame(trps)
        dimnames(trps) <- list(1:nrow(trps), c('x','y'))
        class(trps) <- c('traps', 'data.frame')
        attr(trps, 'detector') <- rep('telemetry', noccasions)
        attr(trps, 'telemetrytype') <- 'independent'
        trps
    }
    
    splitcapt <- split(capt, capt[,1])
    trps <- sapply(splitcapt, readtraps, simplify = FALSE)
    if (length(trps)==1)
        trps <- trps[[1]]
    else
        class(trps) <- c("traps","list")
    
    temp <- make.capthist(capt, trps, fmt = fmt,  noccasions = noccasions,
                          covnames = covnames, sortrows = TRUE, cutval = NULL,
                          noncapt = 'NONE')
    
    if (verify)
        verify(temp)
    temp
}
############################################################################################

telemetryxy <- function (object, includeNULL = FALSE) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    
    if (ms(object)) {
        lapply(object, telemetryxy, includeNULL)
    }
    else {
        xylist <- attr(object, 'telemetryxy')
        if (includeNULL) {
            ## expand for untelemetered animals in object
            nullxy <- rep(list(matrix(nrow=0, ncol=2)), nrow(object))
            names(nullxy) <- rownames(object)
            xylist <- c(xylist, nullxy[!(rownames(object) %in% names(xylist))])
        }
        xylist
    }
}
############################################################################################

telemetered <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, telemetered)
    }
    else {
        rownames(object) %in% names(telemetryxy(object))
    }
}
############################################################################################

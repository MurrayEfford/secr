############################################################################################
## package 'secr'
## as.data.frame.R 
## 2017-10-25 changed 
## 2021-05-19 sortorder ksn for polygon detectors etc.
############################################################################################

as.data.frame.capthist <- function (x, row.names = NULL, optional = FALSE, covariates = FALSE, 
                                    fmt = c("trapID", "XY"), ...)
    
{
    if (ms(x)) {
        lapply(x, as.data.frame, row.names = row.names, optional = optional, 
               covariates = covariates, fmt = fmt, ...)
    }
    else {
        oldopt <- options(stringsAsFactors = FALSE)
        on.exit(options(oldopt))
        fmt <- match.arg(fmt)
        n <- nrow(x)
        S <- ncol(x)
        xname <- deparse(substitute(x), control=NULL)
        
        ID <- animalID(x, sortorder = "ksn")
        occ <- occasion(x, sortorder = "ksn")
        session <- rep(session(x),length(ID))
        
        if (is.null(traps(x))) {
            det <- 0
            trap <- rep('1', length(ID))
            trapn <- rep(1, length(ID))
        }
        else {
            det <- detector(traps(x))  # may be vector...
            # trap <- as.character(trap(x), sortorder = "ksn")
            # bug fixed 2021-07-10
            trap <- as.character(trap(x, sortorder = "ksn"))
            trapn <- trap(x, names = FALSE, sortorder = "ksn")
        }
        
        if (det[1] %in% c('polygon','transect','polygonX','transectX')) {
            XY <- xy(x)
            temp <- data.frame (Session=session, ID=ID, Occasion=occ,
                                x = XY$x, y = XY$y)
        }
        else if (det[1] %in% c('signal')) {
            signal <- signal(x)
            if (fmt == "trapID")
                temp <- data.frame (Session = session, ID = ID, Occasion = occ, 
                    TrapID = trap, Signal = signal)
            else
                temp <- data.frame (Session = session, ID = ID, Occasion = occ, 
                    x = traps(x)$x[trapn], y = traps(x)$y[trapn], Signal = signal)
        }
        else if (all(det %in% c('telemetry'))) {
            xyl <- telemetryxy(x)
            xy <- do.call(rbind, xyl)
            ID <- rep(names(xyl), sapply(xyl, nrow))
            temp <- data.frame (Session=session, ID=ID, Occasion=rep(1,nrow(xy)), 
                                x = xy[,1], y = xy[,2])
        }
        else {
            if (fmt == "trapID")
                temp <- data.frame (Session=session, ID=ID, Occasion=occ, TrapID=trap)
            else
                temp <- data.frame (Session=session, ID=ID, Occasion=occ, 
                                    x = traps(x)$x[trapn], y = traps(x)$y[trapn])
        }
        
        if (!is.null(covariates) & !is.null(covariates(x))) {
            if (!is.logical(covariates) | covariates) {
                covs <- covariates(x)
                if (is.character(covariates)) {
                    covlist <- match(covariates, names(covs))
                    covlist <- covlist[!is.na(covlist)]
                }
                else
                    covlist <- names(covs)
                if (length(covlist)>0) {
                    temp <- cbind(temp, covs[match(ID, rownames(x)), covlist, drop = FALSE])
                }
            }
        }
        # ensure nice sort order
        maxl <- max(str_length(temp$ID))
        temp <- temp[order(str_pad(temp$ID,maxl), temp$Occasion),] 
        row.names(temp) <- 1:nrow(temp)
        temp
    }
}
############################################################################################

as.data.frame.traps <- function (x, row.names = NULL, optional = FALSE, usage = FALSE, 
                                 covariates = FALSE, ...) {
    
    if (ms(x)) {
        lapply(x, as.data.frame, row.names = row.names, optional = optional, 
               covariates = covariates, usage = usage, ...)
    }
    else {
        # purge blanks from names
        row.names(x) <- gsub(' ','',row.names(x))
        det <- detector(x)
        poly <- det[1] %in% c('polygon', 'polygonX')
        transect <- det[1] %in% c('transect', 'transectX')
        if (poly) {
            temp <- data.frame (polyID=polyID(x), x=x$x, y = x$y)
        }
        else if (transect) {
            temp <- data.frame (transectID=transectID(x), x=x$x, y = x$y)
        }
        else {
            temp <- data.frame(x = x$x, y = x$y)
            if (!is.null(usage(x)) & usage) {
                xmat <- data.frame(usage(x))
                names(xmat) <- paste0('u', 1:ncol(xmat))
                temp <- cbind(temp,xmat)
            }
        }
        
        covlist <- numeric(0)
        if (!is.null(covariates) & !is.null(covariates(x))) {
            if (!(is.logical(covariates) & !covariates)) {
                covs <- covariates(x)
                if (is.character(covariates)) {
                    covlist <- match(covariates, names(covs))
                    covlist <- covlist[!is.na(covlist)]
                }
                else
                    covlist <- names(covs)
                
                if (length(covlist)>0) {
                    covnames <- paste(covlist, collapse=' ')
                    covs <- covs[, covlist, drop=FALSE]
                    ## assume order of levels of polyID matches order in x
                    if (poly | transect)
                        covs <- covs[as.numeric(polyID(x)), , drop=FALSE]
                    for (i in 1:length(covlist))
                        covs[,i] <- as.numeric(covs[,i])
                    covs <- apply(covs,1,paste, collapse=' ')
                    covs <- paste ('/',covs)
                    temp <- cbind(temp, covs)
                }
            }
        }
        temp
    }
}
###############################################################################

############################################################################################
## package 'secr'
## write.captures.R (was write.capthist - changed 2010-05-02)
## last changed 2009 03 30, 2009 06 11, 2009 07 08 2009 11 17
## revised 2010 04 01 with new argument append and character value allowed for header
## bug fixed 2010-06-07: use row names when single/multi
## 2011-09-09 trivial: ms()
## Write capture histories to text file in DENSITY format
############################################################################################

write.captures <- function (object, file='', deblank = TRUE, header = TRUE,
    append = FALSE, sess = '1', ndec = 2, covariates = FALSE, tonumeric = TRUE, ...)

{
    if (!is(object, 'capthist'))
        stop ("requires a 'capthist' object")

    if (ms(object)) {
        out <- vector('list')
        out[[1]] <- write.captures (object[[1]], file = file, deblank = deblank,
            header = deparse(substitute(object), control=NULL), append = append,
            sess = session(object)[1], ndec = ndec, covariates = covariates,
                        tonumeric = tonumeric, ...)
        for (i in 2:length(object)) {
            out[[i]] <- write.captures (object[[i]], file = file, deblank = deblank,
                header = FALSE, append = TRUE, sess = session(object)[i], ndec = ndec,
                covariates = covariates, tonumeric = tonumeric, ...)
        }
        invisible(out)
    }
    else {
        n <- nrow(object)
        S <- ncol(object)
        objectname <- ifelse (is.character(header),
            header, deparse(substitute(object), control=NULL))
        header <- ifelse (is.character(header), TRUE, header)
        det <- detector(traps(object))

        ## use universal sort order; re-sorted at end
        ID <- animalID(object, sortorder = 'ksn')
        occ <- occasion(object, sortorder = 'ksn')
        session <- rep(sess,length(ID))

        if (det[1] %in% c('polygon','transect','polygonX','transectX')) {
            XY <- xy(object)
            temp <- data.frame (Session=session, ID=ID, Occasion=occ,
                x=round(XY$x,ndec), y=round(XY$y,ndec))
        }
        else if (det[1] %in% c('signal')) {
            signal <- signal(object)
            trap <- trap(object, sortorder = 'ksn')
            temp <- data.frame (Session=session, ID=ID, Occasion=occ, Detector=trap, Signal=signal)
        }
        else if (all(det %in% c('telemetry'))) {
            xyl <- telemetryxy(object)
            xy <- do.call(rbind, xyl)
            ID <- rep(names(xyl), sapply(xyl, nrow))
            trap <- trap(object, sortorder = 'ksn')
            temp <- data.frame (Session=session, ID=ID, Occasion=rep(1,nrow(xy)), 
                                x = round(xy[,1],ndec), y = round(xy[,2],ndec))
        }
        else {
            trap <- trap(object, sortorder = 'ksn')
            temp <- data.frame (Session=session, ID=ID, Occasion=occ, Detector=trap)
        }

        if (!is.null(covariates) & !is.null(covariates(object))) {
            covs <- covariates(object)
            if (is.character(covariates)) {
                covlist <- match(covariates, names(covs))
                covlist <- covlist[!is.na(covlist)]
            }
            else
                covlist <- names(covs)
            if (length(covlist)>0) {
                if (tonumeric) {
                    for (i in 1:length(covlist))
                        covs[,i] <- as.numeric(covs[,i])
                }
                temp <- cbind(temp, covs[match(ID, rownames(object)), covlist, drop=FALSE])
            }
        }

        if (header) {
            cat ("# Capture histories exported from '", objectname, "' \n", sep="", file=file)
            cat ('#', format(Sys.time(), "%a %b %d %X %Y"), '\n', append = TRUE, file=file)
            cat ('#', names(temp), '\n', append = TRUE, file=file)
            append <- TRUE
        }
        if (deblank) temp$Session <- gsub(' ','', temp$Session)
        if (deblank) temp$Session <- gsub(',','', temp$Session)
        if (any(nchar(temp$Session)>17)) {
            warning ("truncating long session names")
            temp$Session <- substring(temp$Session,1,17)
        }
        
        temp <- temp[order(temp$Session, temp$ID, temp$Occasion), ]
        write.table(temp, file = file, row.names = FALSE, col.names = FALSE,
            append = append, quote = FALSE, ...)
    }
}
############################################################################################

#############################################################################
## package 'secr'
## join.R

## join returns single-session object from list of inputs

## 2016-10-10 secr3
## 2017-01-29 telemetry fixes
## 2017-12-23
## 2018-05-12 join uses attr() to avoid copying, and 
## 2020-04-07 join SIMPLIFY = FALSE in call to mapply(condition.usage...)
## 2021-05-27 fixed bug in join|onesession for nonspatial data
#############################################################################

addzerodf <- function (df, oldCH, sess) {
    ## add dummy detection records to dataframe for 'all-zero' case
    ## that arises in sighting-only mark-resight with known marks
    allzero <- apply(oldCH,1,sum)==0
    naz <- sum(allzero)
    if (naz > 0) {
        df0 <- expand.grid(
            newID = rownames(oldCH)[allzero], 
            newocc = NA,
            newtrap = trap(oldCH)[1], 
            alive = TRUE, 
            sess = sess,
            stringsAsFactors = FALSE)
        df$x <- NULL; df$y <- NULL  ## 2021-04-08
        df <- rbind(df,df0)
        if (!is.null(xy(oldCH))) {
            df$x <- c(xy(oldCH)$x, rep(NA, naz))
            df$y <- c(xy(oldCH)$y, rep(NA, naz))
        }
        if (!is.null(signal(oldCH)))  {
            df$signal <- c(signal(oldCH), rep(NA, naz))
        }
    }
    df
}

#-------------------------------------------------------------------------------

join <- function (object, remove.dupl.sites = TRUE, tol = 0.001,
                  sites.by.name = FALSE, drop.sites = FALSE, intervals = NULL, 
                  sessionlabels = NULL, timevaryingcov = NULL) {

    ####################################################################
    onesession <- function (sess) {
        ## form CH as a dataframe
        CH <- object[[sess]]
        if (drop.sites) attr(CH, 'traps') <- NULL   ## 2018-05-12
        ## 2021-05-19 sortorder to match xy
        newID <- animalID(CH, sortorder = 'ksn')
        newocc <- occasion(CH, sortorder = 'ksn') + before[sess]
        newtrap <- if (!is.null(traps(CH))) trap(CH, sortorder = 'ksn')
        else rep(1,length(newID))  ## 2021-05-27
        
        newalive <- alive(CH, sortorder = 'ksn')
        df <- data.frame(
            newID = newID, 
            newocc = newocc, 
            newtrap = newtrap,
            alive = newalive, 
            sess = rep(sess, length(newID)),
            stringsAsFactors = FALSE
        )
        if (!is.null(xy(CH)))
            df[,c('x','y')] <- xy(CH)
        if (!is.null(signal(CH)))
            df[,'signal'] <- signal(CH)

        ## add all-zero (sighting-only) records as required,
        ## otherwise leave unchanged
        addzerodf (df, CH, sess)
    }
    ####################################################################

    condition.usage <- function (trp, i, nocc) {
        if (!is.null(trp)) {
            us <- matrix(0, nrow=nrow(trp), ncol=nnewocc)
            ## don't understand need for this 2019-10-22
            if ('telemetry' %in% detector(trp)) {
                occasions <- outputdetector == 'telemetry'
            }
            else {
                s1 <- c(1, cumsum(nocc)+1)[i]
                s2 <- cumsum(nocc)[i]
                if (any(is.na(c(s1,s2)))) {
                    cat("Houston, we have a problem\n")
                    browser()
                }
                occasions <- s1:s2
            }
            if (is.null(usage(trp))) {
                us[,occasions] <- 1
            }
            else {
                us[,occasions] <- usage(trp)
            }
            usage(trp) <- us
        }
        trp
    }
    ####################################################################
    ## preparing for merge when traps vary... 2015-10-29
    condition.sightings <- function (CH, i, type = 'Tu') {
        T <- attr(CH, type)
        if (!is.null(T)) {
            if (is.matrix(T)) {
                Tnew <- matrix(0, nrow = nrow(traps(CH)), ncol = nnewocc)
                s1 <- c(1, cumsum(nocc)+1)[i]
                s2 <- cumsum(nocc)[i]
                Tnew[,s1:s2] <- T
                attr(CH, type) <- Tnew
            }
        }
        CH
    }
    ####################################################################
    ## mainline
    if (!ms(object) | any(sapply(object, class) != 'capthist'))
        stop("requires multi-session capthist object or list of ",
             "single-session capthist")
    detectorlist <- lapply(object, secr_expanddet)
    outputdetector <- unlist(detectorlist)

    nsession <- length(object)
    nocc <- sapply(object, ncol)
    names(nocc) <- NULL
    nnewocc <- sum(nocc)
    ## cumulative number of preceding occasions
    before <- c(0, cumsum(nocc)[-nsession])
    ##------------------------------------------------------------------
    ## combine capthist as one long dataframe
    df <- lapply(1:nsession, onesession)
    df <- do.call(rbind, df)
    n <- length(unique(df$newID))

    ##------------------------------------------------------------------
    ## resolve traps
    ## first check whether all the same (except usage)
    if (!(drop.sites | is.null(traps(object)))) {
        temptrp <- lapply(traps(object), function(x) {usage(x) <- NULL; x})
        sametrp <- all(sapply(temptrp[-1], identical, temptrp[[1]]))
        telemetrytrap <- function (ch) {
            if ('telemetry' %in% detector(traps(ch))) dim(ch)[3] else 0
        }
        if (sametrp & remove.dupl.sites) {
            newtraps <- temptrp[[1]]
            detector(newtraps) <- outputdetector         ## 2018-05-28
            class(newtraps) <- c("traps", "data.frame")
            if (length(usage(traps(object))) > 0)
                usage(newtraps) <- do.call(cbind, usage(traps(object)))
            ## df$newtrap unchanged
        }
        else {
            temptrp <- traps(object)
            if ('telemetry' %in% outputdetector) {
                # drop all notional 'telemetry' traps and replace at end
                ttraps0 <- sapply(object, telemetrytrap)
                ttraps <- ttraps0[ttraps0>0]
                teltrapno <- ttraps[length(ttraps)] # use last
                df$newtrap[df$newtrap %in% ttraps] <- teltrapno
                dropteltrap <- function (trps, teltrap) {
                    if (teltrap>0) {
                        if (nrow(trps)==1)
                            NULL
                        else
                            subset(trps, (1:nrow(trps)) != teltrap)
                    }
                    else
                        trps
                }
                newteltrap <- subset(temptrp[[1]],1)
                temptrp <- mapply(dropteltrap, temptrp, ttraps0, SIMPLIFY = FALSE)
                rownames(newteltrap) <- teltrapno
                temptrp <- c(temptrp, list(newteltrap))
            }
            else {
                if (!sites.by.name) df$newtrap <- paste(df$newtrap,df$sess, sep=".")
            }
            temptrp <- mapply(condition.usage, temptrp, 1:length(temptrp), 
                              MoreArgs = list(nocc = nocc),
                              SIMPLIFY = FALSE)  ## SIMPLIFY = FALSE added 2020-04-07
            
            ## 2019-10-22 - attempt to cleanup telemetry
            # nt <- length(temptrp)
            # temptrp[1:(nt-1)] <- mapply(condition.usage, temptrp[-nt], 1:(nt-1), MoreArgs=list(nocc=nocc), SIMPLIFY = FALSE)
            # temptrp <- temptrp[!sapply(temptrp, is.null)]
            
            ## 2018-05-11
            ## workaround for large datasets
            if (sites.by.name) {
                temptrp <- do.call(rbind, lapply(temptrp, as.matrix))
                newtraps <- as.data.frame(temptrp[!duplicated(row.names(temptrp)),])
            }
            else {
                newtraps <- do.call(rbind, c(temptrp, renumber = FALSE, checkdetector = FALSE))
            }
            ## end workaround
            detector(newtraps) <- outputdetector
            class(newtraps) <- c("traps", "data.frame")
        }
        if (all(outputdetector %in% .localstuff$polydetectors))
            df$newtrap <- factor(df$newtrap)
        else
            df$newtrap <- factor(df$newtrap, levels=rownames(newtraps))
    }
    else {
        sametrp <- FALSE
    }
    ##------------------------------------------------------------------
    ## ensure retain all occasions
    df$newocc <- factor(df$newocc, levels = 1:nnewocc)
    ## drop any records with missing data
    df <- df[!apply(df,1,function(x) any(is.na(x))),, drop = FALSE]
    ##------------------------------------------------------------------
    ## construct new capthist matrix or array from positive detections
    
    ## 2020-10-20 following can fail (>= 2^31 elements) because 
    ##            many distinct newtrap until these are resolved
    
    nlevels <- sapply(df[,c('newID','newocc','newtrap')], function(x) 
        length(levels(factor(x))))
    if (prod(nlevels) >= 2^31) {
        stop ("More than 2^31 combinations of newID, newocc, newtrap in join(); ",
            "try sites.by.name or drop.sites")
    }
    tempnew <- table(df$newID, df$newocc, df$newtrap, useNA = "no")
    i <- cbind(as.character(df$newID), df$newocc, as.character(df$newtrap))
    tempnew[i] <- tempnew[i] * (df$alive * 2 - 1)

    # alivesign <- tapply(df$alive, list(df$newID,df$newocc,df$newtrap),all)
    # alivesign[is.na(alivesign)] <- TRUE
    # alivesign <- alivesign * 2 - 1
    # tempnew <- tempnew * alivesign
    ##------------------------------------------------------------------
    ## pile on the attributes...
    class(tempnew) <- 'capthist'
    if (!(drop.sites || is.null(traps(object)))) {
        # traps(tempnew) <- newtraps
        attr(tempnew, 'traps') <- newtraps
    }
    attr(tempnew, 'session') <- 1
    neworder <- order (df$newocc, df$newID, df$newtrap)
    ##------------------------------------------------------------------
    ## concatenate marking-and-resighting-occasion vectors
    tempmarkocc <- unlist(markocc(traps(object)))
    if (!is.null(tempmarkocc)) {
        names(tempmarkocc) <- NULL
        markocc(traps(tempnew)) <- tempmarkocc
    }

    ##------------------------------------------------------------------
    ## unmarked and nonID sightings
    ## not yet implemented for varying traps
    if (sametrp & remove.dupl.sites & !drop.sites) {
        ## retain unmarked sightings and nonID sightings if present
        ## ignore if NULL
        Tu <- Tu(object)
        if (!is.null(Tu[[1]])) {
            if (!all(sapply(Tu, is.matrix)))
                Tu(tempnew) <- do.call(sum, Tu)
            else
                Tu(tempnew) <- do.call(cbind, Tu)
        }

        Tm <- Tm(object)
        if (!is.null(Tm[[1]])) {
            if (!all(sapply(Tm, is.matrix)))
                Tm(tempnew) <- do.call(sum, Tm)
            else
                Tm(tempnew) <- do.call(cbind, Tm)
        }
    }
    else {
        ## Tu, Tm not ready yet
        if (!is.null(Tu(object[[1]])) | !is.null(Tm(object[[1]])))
            stop ("join does not yet merge sighting matrices when traps vary")
    }

    ##------------------------------------------------------------------
    ## covariates, xy, signal attributes

    covobj <- covariates(object)
    if (!is.null(covobj)) {
        
        tempcov <- do.call(rbind, covobj)
        if (!is.null(tempcov)) {
            IDcov <- unlist(lapply(object,rownames))
            ## use first match
            tempcov <- tempcov[match(rownames(tempnew), IDcov),,drop = FALSE]
            rownames(tempcov) <- rownames(tempnew)
            if (!is.null(timevaryingcov)) {
                tcv <- vector('list')
                for (i in timevaryingcov) {
                    for (j in 1:nsession) {
                        covname <- paste(i,j,sep=".")
                        tempcov[,covname] <- rep(NA, nrow(tempcov))
                        id <- rownames(object[[j]])
                        tempcov[match(id, rownames(tempcov)),covname] <- covobj[[j]][,i]
                    }
                    tcv[[i]] <- paste(i, 1:nsession, sep='.')
                }
                attr(tempnew, 'timevaryingcov') <- tcv
            }
            attr(tempnew, 'covariates') <- tempcov
            
        }
    }

    ##------------------------------------------------------------------
    ## telemetry fixes

    if ('telemetry' %in% outputdetector) {
        oldtelem <- lapply(object, telemetryxy)
        telnames <- unique(unlist(lapply(oldtelem,names)))
        newtelem <- vector('list', length(telnames))
        names(newtelem) <- telnames
        for (id in telnames) {
            newtelem[[id]] <- do.call(rbind, lapply(oldtelem, '[[', id))
        }
        telemetryxy(tempnew) <- newtelem
    }

    ##------------------------------------------------------------------
    ## negotiate problem that all-zero histories have no xy, signal
    tempdf <- df[neworder,, drop = FALSE]
    if (!is.null(df$x)) {
        xy(tempnew) <- tempdf[!is.na(tempdf$newocc), c('x','y')]
    }
    if (!is.null(df$signal))
        signal(tempnew) <- tempdf[!is.na(tempdf$newocc),'signal']
    ##------------------------------------------------------------------
    ## purge duplicate sites, if requested
    if (remove.dupl.sites & !sametrp & !sites.by.name) {
        tempnew <- reduce(tempnew, span=tol, dropunused = FALSE, verify = FALSE)
    }
    ## remember previous structure, for MARK-style robust design
    tmpintervals <- unlist(sapply(nocc, function(x) c(1,rep(0,x-1))))[-1]
    if (!is.null(intervals)) {
        if (length(intervals) != sum(tmpintervals>0))
            stop("invalid intervals argument")
        tmpintervals[tmpintervals>0] <- intervals
    }
    else if (!is.null(attr(object, 'intervals'))) {
        attrintervals <- attr(object, 'intervals')
        if (length(attrintervals) != sum(tmpintervals>0))
            stop("invalid intervals attribute")
        tmpintervals[tmpintervals>0] <- attrintervals
    }
    if (is.null(sessionlabels)) sessionlabels <- sessionlabels(object)
    if (is.null(sessionlabels)) sessionlabels <- session(object)
    attr(tempnew, 'intervals') <- tmpintervals
    attr(tempnew, 'sessionlabels') <- sessionlabels

    ##------------------------------------------------------------------

    tempnew

}

unjoin <- function (object, intervals, ...) {
    if (missing(intervals) & is.null(attr(object,'intervals')))
        stop ("requires 'intervals' to define sessions")
    if (missing(intervals) )
        intervals <- attr(object,"intervals")
    session <- c(0,cumsum(intervals>0))+1
    nsess <- max(session)
    if (nsess<2) {
        warning ("intervals define only one session")
        return(object)
    }
    newobj <- vector(mode='list', length=nsess)
    for (sess in 1:nsess) {
        newobj[[sess]] <- subset(object, occasions = (session==sess), ...)
    }
    class (newobj) <- c('capthist', 'list')
    session(newobj) <- 1:nsess
    return(newobj)
}


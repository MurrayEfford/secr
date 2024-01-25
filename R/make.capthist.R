############################################################################################
## package 'secr'
## make.capthist.R
## 2010 05 02 (transferred from methods.R) 2010 05 03, 2010-11-21, 2011-01-21
## 2012-02-08 signalnoise
## 2012-02-09 revamped sorting
## 2012-02-12 finished tidy up related to signalframe
## 2012-10-19 telemetry detector type
## 2013-10-29 noccasions derived from usage
## 2015-11-03 Improved handling of 'noncapt' with 'validcapt' - no need for dummy trapID etc
## 2017-01-11 adjusted for direct input of telemetry
## 2019-03-11 transect bug fixed (capttrap order)
## 2020-08-27 coerce third column of captures dataframe to integer
## 2021-09-21 modified for nonspatial detector type
############################################################################################

make.capthist <- function (captures, traps, fmt = c("trapID", "XY"), noccasions = NULL,
    covnames = NULL, bysession = TRUE, sortrows = TRUE, cutval = NULL, tol = 0.01, 
    snapXY = FALSE, noncapt = 'NONE', signalcovariates = NULL)

# captures is a dataframe with the structure:
# fmt = 'trapID'
#   column 1	Session
#   column 2	AnimalID
#   column 3	Occasion
#   column 4	TrapID
#   column 5    Signal   (optional)
#   column 6    Noise    (optional)

# fmt = 'XY'
#   column 1	Session
#   column 2	AnimalID
#   column 3	Occasion
#   column 4	x
#   column 5	y
#   column 6    Signal    (optional)
#   column 7    Noise     (optional)

{
    fmt <- match.arg(fmt)
    session <- captures[,1]
    sessionlevels <- unique(session)  ## retains order
    ## session <- factor(session) ## 2010 04 01 automatically sorts levels
    ## use numeric sort if appropriate 2010 05 02
    if (suppressWarnings( all(!is.na(as.numeric(sessionlevels)))))
        sessionlevels <- sessionlevels[order(as.numeric(sessionlevels))]
    else
        sessionlevels <- sort(sessionlevels)
    session <- factor(session, levels=sessionlevels)
    MS <- bysession & ( length(sessionlevels) > 1)
    if (MS) {  # recursive call of make.capthist
        capturelist <- split (captures, session)
        nsession <- length(capturelist)
        traplist <- inherits(traps, 'list')
        occvector <- length(noccasions)>1
        if (traplist & (length(traps) != nsession))
            stop ("multi-session 'traps' list does not match 'captures'")
        if (occvector & (length(noccasions) != nsession))
            stop ("requires one element in 'noccasions' for each session")

        capthist <- vector('list', nsession)
        for (i in 1:nsession) {
            if (traplist)  trps <- traps[[i]] else trps <- traps
            if (occvector) nocc <- noccasions[i] else nocc <- noccasions
            capthist[[i]]  <- make.capthist (
                captures = capturelist[[i]],
                traps = trps,
                fmt = fmt,
                noccasions = nocc,
                covnames = covnames,
                bysession = FALSE,         ## 2010 04 01
                sortrows = sortrows,
                cutval = cutval,
                tol = tol,
                snapXY = snapXY)
        }
        names(capthist) <- levels(session)
        class(capthist) <- c('capthist', 'list')
        capthist
    }

    else ## single-session call
    {
        ## 2020-08-27 in case has been read as character
        captures[,3] <- as.integer(captures[,3])
        ## 2017-03-28
        ## sort by session, animal & occasion, retaining original order otherwise
        ## by abs(occasion) from 2023-09-21
        captures <- captures[order(captures[,1], captures[,2], abs(captures[,3])),]
        if (missing(traps)) traps <- NULL
        
        if (any(detector(traps) %in% .localstuff$exclusivedetectors)) {
            ## repeated <- duplicated(paste0(captures[,1],captures[,2],captures[,3]))
            ## 2021-01-30 bug fix from Richard Glennie
            repeated <- duplicated(paste(captures[,1], captures[,2], captures[,3], sep = "."))
            if (any(repeated)) {
                if (length(detector(traps)) > 1)
                    occasiondetector <- detector(traps)[abs(captures[,3])]
                else
                    occasiondetector <- detector(traps)
                OK <- (!occasiondetector %in% .localstuff$exclusivedetectors) |
                    !repeated
                if (sum(OK) != nrow(captures))
                    warning("dropping repeat detections within occasions at ",
                            "exclusive detectors (traps)", call. = FALSE)
                captures <- captures[OK, ]
            }
        }
        ## end 2017-03-28

        uniqueID <- unique(captures[,2])
        ## condition inserted 2015-11-03 to avoid need to specify valid trap for noncapt
        validcapt <- uniqueID != noncapt
        
        if (is.null(traps)) {
            captTrap <- rep(1, nrow(captures))
            fmt <- "nonspatial"   # 2023-03-30
        }
        else if (all(validcapt)) {
            if (!(fmt %in% c('trapID','XY')))
                stop ("capture format not recognised")
            if (fmt!='trapID') {
                if (ncol(captures)<5)
                    stop ("too few columns in capture matrix")
                
                if (all(detector(traps) %in% 'telemetry')) {
                    captTrap <- rep(1, nrow(captures))
                }
                else if (all(detector(traps) %in% c('polygon','polygonX'))) {
                    if (nrow(captures)==0)  ## 2016-01-02
                        captTrap <- numeric(0)
                    else
                        captTrap <- xyinpoly(captures[,4:5], traps)
                    if (any(captTrap==0)) {
                        captures <- captures[captTrap>0,]  ## first! 2010-11-17
                        captTrap <- captTrap[captTrap>0]
                        warning ("detections with coordinates outside ",
                                 "polygon(s) were dropped", call. = FALSE)
                    }
                }
                else if (all(detector(traps) %in% c('transect','transectX'))) {
                    if (nrow(captures)==0) {
                        captTrap <- numeric(0)
                    }
                    else {
                        captTrap <- xyontransect(captures[,4:5], traps)
                        ########################################################
                        # implement snapXY for transects 2021-12-11
                        
                        if (any(captTrap==0) && snapXY) {

                            # split by transect
                            bytransect <- strsplit(rownames(traps),'.', fixed=TRUE)
                            ID <- do.call(rbind, bytransect)[,1]
                            vlist <- split(traps, ID)
                            vlist <- lapply(vlist, as.matrix)
                            # each transect as sfg
                            vlist <- lapply(vlist, st_linestring)
                            # combine linestrings in one sfc
                            v <- st_sfc(vlist)     
                            xy <- st_as_sf(captures[,4:5, drop = FALSE], coords = 1:2)
                            xy2 <- snap_points (xy, v, max_dist = tol) # see utility.R
                            distances <- st_distance(xy,v) 
                            distances <- apply(distances,1,min)   # nearest 2022-12-01
                            OK <- distances < tol
                            if (any(!OK)) {
                                cat('points not within tol = ', tol, 'of any transect\n')
                                print(cbind(captures, distances)[!OK,])                       
                            }
                            captures[OK,4:5] <- st_coordinates(xy2)[OK,]
                            
                            warning(call. = FALSE, sum(OK), 
                                " detection(s) snapped to transect(s), maximum distance ", 
                                round(max(distances),2), " m")

                            # repeat transect assignment
                            captTrap <- xyontransect(captures[,4:5], traps)
                        }
                        ########################################################
                    }
                    if (any(captTrap==0)) {
                        warning (call. = FALSE, sum(captTrap==0), 
                            " detection(s) with coordinates not on ",
                            "any transect were dropped")
                        captures <- captures[captTrap>0,]  ## first!! 2019-03-11
                        captTrap <- captTrap[captTrap>0]
                    }
                }
                else {
                    if (!all(detector(traps) %in% .localstuff$pointdetectors))
                        stop("cannot mix point and other detector types")
                    if (snapXY) {
                        ## snap to site
                        dtrap <- distancetotrap(captures[,4:5], traps)
                        if (any(dtrap>tol)) {
                            print(captures[dtrap>tol,1:5])                       
                            stop("capture XY greater than distance tol from nearest trap")
                        }
                        captTrap <- nearesttrap(captures[,4:5], traps)                      
                    }
                    else {
                        trapID    <- interaction(traps$x, traps$y)
                        captTrap  <- match(interaction(captures[,4], captures[,5]), trapID)
                    }
                    if (any(is.na(captTrap))) {
                        print(captures[is.na(captTrap),])
                        stop ("failed to match some capture locations ",
                              "to detector sites")
                    }
                }
            }
            else {
                if (any(detector(traps) %in% .localstuff$polydetectors))
                    stop ("use fmt XY to input detections from polygons or transects or telemetry")
                captTrap <- match(captures[,4], row.names(traps))
                if (any(is.na(captTrap))) {
                    print(captures[is.na(captTrap),])
                    stop ("failed to match some capture locations ",
                          "to detector sites")
                }
            }
        }  ## end of condition for noncapt 2015-11-03
        else captTrap <- numeric(0)

        if (!is.null(usage(traps))) {   ## 2013-10-29
            nocc <- ncol(usage(traps))
            if (nocc <  max(abs(captures[,3])))
                stop("fewer usage fields than max(occasion)")
        }
        else {
            nocc <- max(abs(captures[,3]))
        }
        nocc <- ifelse (is.null(noccasions), nocc, noccasions)
        if (!is.null(traps)) {
            if (is.null(detector(traps)))
                stop ("'traps' must have a detector type e.g. 'multi'")
            if (is.null(cutval) & any(detector(traps)  %in% c('signal','signalnoise')))
                stop ("missing 'cutval' (signal threshold) for signal data")
        }
        wout <- NULL
        ID   <- NULL

        ## optional row sort 2009 09 26, tweaked 2010 05 01
        if (sortrows) {
            if (suppressWarnings( all(!is.na(as.numeric(uniqueID)))))
                rowOrder <- order(as.numeric(uniqueID))
            else
                rowOrder <- order (uniqueID)
            uniqueID <- uniqueID[rowOrder]
        }

        captID <- as.numeric(factor(captures[validcapt,2], levels=uniqueID))
        nID    <- length(uniqueID)
        detectionOrder <- order(captTrap, abs(captures[validcapt,3]), captID)

        ## all dim3 now
        if ((length(detector(traps))>1) & (length(detector(traps)) != nocc))
            stop("detector vector of traps object does not match nocc")
            
        w <- array (0, dim=c(nID, nocc, ndetector(traps)))
        ## drop rows if dummy input row indicates no captures
        if (any(!validcapt)) {
            if (any(validcapt))
                stop ("cannot combine data and noncapt")
            w <- w[FALSE, , , drop = FALSE]
            dimnames(w) <- list(NULL, 1:nocc, 1:ndetector(traps))
        }
        else {
            dimnames(w) <- list(NULL, 1:nocc, 1:ndetector(traps))
            if (nID>0) {
                dimnames(w)[[1]] <- 1:nID   ## 2016-01-02
                temp <- table (captID, abs(captures[,3]), captTrap)
                if ('0' == dimnames(temp)[[2]][1]) {
                    ## drop zero occasion but retain animals for mark-resight
                    temp <- temp[,-1,]
                }
                d <- dimnames(temp)
                d <- lapply(d, as.numeric)
                w[d[[1]], d[[2]], d[[3]]] <- temp
                
                ## fix to retain deads 2010-08-06
                dead <- captures[,3]<0
                deadindices <- cbind(captID[dead], abs(captures[dead,3]), captTrap[dead])
                w[deadindices] <- w[deadindices] * -1
            }
            #################################
        }
        wout <- abind(wout, w, along=1)
        dimnames(wout)[[2]] <- 1:nocc
        attr(wout,'covariates') <- data.frame() ## default 2015-01-06
        
        if (nrow(wout) > 0) {

            dimnames(wout)[[1]] <- uniqueID

            ## check added 2010 05 01
            if (any(is.na(wout))) {
                rw <- row(wout)[is.na(wout)]
                print (wout[rw,,,drop=F])
                stop ("missing values not allowed")
            }
            ## code to input permanent individual covariates if these are present
            zi <- NULL
            startcol <- switch(fmt, nonspatial = 4, trapID = 5, XY = 6, 6)
            if (fmt != "nonspatial" && all(detector(traps) %in% c('signal'))) startcol <- startcol+1
            if (fmt != "nonspatial" && all(detector(traps) %in% c('signalnoise'))) startcol <- startcol+1
            if (ncol(captures) >= startcol) {
                ## zi <- as.data.frame(captures[,startcol:ncol(captures), drop=F])
                zi <- as.data.frame(captures[,startcol:ncol(captures), drop = FALSE], 
                                    stringsAsFactors = TRUE)   ## 2020-05-18
            }
            if (!is.null(zi)) {
                ## find first match of ID with positive value for each covar
                temp <- zi[1:length(uniqueID),,drop=FALSE]
                temp[,] <- NA

                for (j in 1:ncol(zi)) {
                    nonmissing <- function(x) x[!is.na(x)][1]
                    ## levels=uniqueID to retain sort order
                    tempj <- split(zi[,j], factor(captures[,2], levels=uniqueID))
                    tempj2 <- sapply(tempj, nonmissing)
                    ## 2011-01-20 only convert to factor if character
                    if (is.character(tempj2))
                        tempj2 <- factor(tempj2)
                    temp[,j] <- tempj2
                    rownames(temp) <- names(tempj)          ## 2010 02 26
                }
                if (is.null(covnames)) names(temp) <- names(zi)
                else {
                    if (ncol(temp) != length(covnames))
                        stop ("number of covariate names does not match")
                    names(temp) <- covnames
                }
                attr(wout,'covariates') <- temp
            }
        }
        class (wout) <- 'capthist'
        if (!is.null(traps)) {
            traps(wout) <- traps
        }
        session(wout)  <- as.character(captures[1,1])
        if (nrow(wout) > 0 && !is.null(traps)) {
            
            if (all(detector(traps) %in% 'telemetry')) {
                xyl <- split(captures[,4:5], captures[,2], drop = TRUE)
                telemetryxy(wout) <- xyl
            }
            if (all(detector(traps) %in% .localstuff$polydetectors)) {
                xy <- captures[detectionOrder,4:5]
                names(xy) <- c('x','y')
                attr(wout,'detectedXY') <- xy
            }
            if (all(detector(traps) %in% c('signal','signalnoise'))) {
                if (is.null(cutval))
                    stop ("missing value for signal threshold")
                if (fmt=='XY')
                    signl <- captures[,6]
                else
                    signl <- captures[,5]
                signl <- signl[detectionOrder]
                signal(wout) <- signl
                if (detector(traps)[1] %in% 'signalnoise') {
                    if (fmt=='XY')
                        nois <- captures[,7]
                    else
                        nois <- captures[,6]
                    nois <- nois[detectionOrder]
                    noise(wout) <- nois
                }
                if (!is.null(signalcovariates)) {
                    if (!all(signalcovariates %in% names(captures)))
                        stop ("missing signal covariate(s)")
                    ## need to maintain order 2012-09-15
                    attr(wout, 'signalframe') <- cbind(attr(wout, 'signalframe'),
                        captures[detectionOrder,signalcovariates])
                }
                attr(wout, 'cutval')   <- cutval
                ## apply cutval
                wout <- subset(wout, cutval = cutval)
            }
        }
        wout
    }   ## end of single-session call
}
############################################################################################
#############################################################################
## package 'secr'
## Methods for classes traps, capthist and mask
## Last changed
## 2019-10-25 masksize()
## 2019-11-10 summary.traps etc. moved to file summary.traps.R
## 2020-08-27 radical revamp of occasion, trap and animalID functions, 
##            with unified code
## 2020-11-07 alive() rewritten to match preceding
## 2022-01-05 fixed alive() problem with one animal
## 2022-08-23 summary.mask and summary.popn removed to separate files
## 2022-11-15 read.mask separate file
## 2023-04-14 trim.secrlist
## 2023-05-30 shift.mask updates boundingbox; explicitly call secr::shift 
## 2023-08-19 as.popn()  ppp > popn
###############################################################################

# Generic methods for extracting attributes etc

usage      <- function (object, ...) UseMethod("usage")
telemetrytype <- function (object, ...) UseMethod("telemetrytype")
markocc    <- function (object, ...) UseMethod("markocc")
Tu         <- function (object, ...) UseMethod("Tu")
Tm         <- function (object, ...) UseMethod("Tm")
Tn         <- function (object, ...) UseMethod("Tn")
nontarget  <- function (object, ...) UseMethod("nontarget")
clusterID  <- function (object, ...) UseMethod("clusterID")
clustertrap <- function (object, ...) UseMethod("clustertrap")
covariates <- function (object, ...) UseMethod("covariates")
traps      <- function (object, ...) UseMethod("traps")
detector   <- function (object, ...) UseMethod("detector")
spacing    <- function (object, ...) UseMethod("spacing")
session    <- function (object, ...) UseMethod("session")
trim       <- function (object, drop, keep) UseMethod("trim")
timevaryingcov <- function (object, ...) UseMethod("timevaryingcov")

rotate     <- function (object, degrees, centrexy=NULL, ...) UseMethod("rotate")
shift      <- function (object, shiftxy, ...) UseMethod("shift")
flip       <- function (object, lr = FALSE, tb = FALSE, ...) UseMethod("flip")

ms         <- function (object, ...) UseMethod("ms")
detectpar  <- function (object, ...) UseMethod("detectpar")
signal     <- function (object, ...) UseMethod("signal")
noise      <- function (object, ...) UseMethod("noise")

derived    <- function (object, ...) UseMethod("derived")
region.N   <- function (object, ...) UseMethod("region.N")
LLsurface  <- function (object, ...) UseMethod("LLsurface")
intervals  <- function (object, ...) UseMethod("intervals")
sessionlabels <- function (object, ...) UseMethod("sessionlabels")
AICcompatible <- function (object, ...) UseMethod("AICcompatible")

############################################################################################

'intervals<-' <- function (object, value)
    structure (object, intervals = value)

intervals.default <- function (object, ...)       {
    attr(object,'intervals', exact = TRUE)
}

'sessionlabels<-' <- function (object, value) {
    structure (object, sessionlabels = value)
}

sessionlabels.default <- function (object, ...)       {
    attr(object,'sessionlabels', exact = TRUE)
}

############################################################################################


# Default methods for specialised functions

ms.default <- function (object, ...)       {
    inherits(object, 'list')
}
ms.mask <- function (object, ...)       {
    !is.data.frame(object)
}
ms.secr <- function (object, ...)       {
    ms(object$capthist)
}

usage.default <- function (object, ...)       {
    if (ms(object)) lapply(object, usage.default, ...)
    else attr(object,'usage',exact = TRUE)
}

markocc.default <- function (object, ...)       {
    if (ms(object)) lapply(object, markocc.default, ...)
    else attr(object,'markocc',exact = TRUE)
}

telemetrytype.default <- function (object, ...)       {
    if (ms(object)) lapply(object, telemetrytype.default, ...)
    else {
        tmp <- attr(object,'telemetrytype',exact = TRUE)
        if (is.null(tmp)) "none"
        else tmp
    }
}

Tu.default <- function (object, ...)       {
    if (ms(object)) lapply(object, Tu.default, ...)
    else attr(object,'Tu',exact = TRUE)
}

Tm.default <- function (object, ...)       {
    if (ms(object)) lapply(object, Tm.default, ...)
    else attr(object,'Tm',exact = TRUE)
}

Tn.default <- function (object, ...)       {
    if (ms(object)) lapply(object, Tn.default, ...)
    else attr(object,'Tn',exact = TRUE)
}

nontarget.default <- function (object, ...) {
    if (ms(object)) lapply(object, nontarget.default, ...)
    else attr(object,'nontarget',exact = TRUE)
}

sighting <- function(object) {
    mocc <- unlist(markocc(object))
    if (is.null(mocc)) {
        FALSE
    }
    else {
        any(mocc < 1)   ## unlist added 2019-08-13
    }
}


clusterID.default <- function (object, ...)       {
    if (ms(object)) lapply(object, clusterID.default, ...)
    else attr(object,'cluster',exact = TRUE)
}

clustertrap.default <- function (object, ...)       {
    if (ms(object)) lapply(object, clustertrap.default, ...)
    else attr(object,'clustertrap',exact = TRUE)
}

covariates.default <- function (object, ...)  {
    if (ms(object)) lapply(object, covariates.default, ...)
    else attr(object,'covariates',exact = TRUE)
}

timevaryingcov.default <- function (object, ...)  {
    if (ms(object)) {
        lapply(object, timevaryingcov.default, ...)
    }
    else attr(object,'timevaryingcov')
}

traps.default <- function (object, ...)       {
    if (ms(object)) {
        temp <- lapply(object, traps.default, ...)
        if (any(sapply(temp, is.null)))
            NULL
        else {
            class(temp) <- c('traps', 'list')
            temp
        }
    }
    else{
        attr(object,'traps',exact = TRUE)
    }
}

detector.default <- function (object, ...)    {
## assumed constant across MS
    if (ms(object)) {
#        detector.default (object[[1]], ...)
        lapply(object, detector.default, ...)  ## returns list 2017-01-29
    }
    else {
        if (is.null(object)) NULL
        else attr(object,'detector',exact = TRUE)
    }
}

spacing.default <- function (object, ...)    {
    if (is.null(object))
        NULL
    else {
        attr(object,'spacing',exact = TRUE)
    }
}

spacing.traps <- function (object, ..., recalculate = FALSE)    {
    if (ms(object)) {
        sapply(object, spacing.traps, ..., recalculate)
    }
    else {
        if (is.null(object)) {
            NULL
        }
        else {
            temp <- attr(object,'spacing')
            if (is.null(temp) || recalculate) {
                if (nrow(object)>1) {
                    points <- matrix(unlist(object), ncol=2)
                    nearest <- nearestcpp(points, points, non_zero = TRUE)
                    sp <- nearest$distance[nearest$index > -1]
                    median(sp)
                }
                else
                    numeric(0) # NA
            }
            else
                temp
        }
    }
}

spacing.mask <- function (object, ..., recalculate = FALSE)    {
    if (ms(object)) {
        sapply(object, spacing.mask, ..., recalculate = recalculate)
    }
    else {
        if (is.null(object)) NULL
        else {
            temp <- attr(object,'spacing',exact = TRUE)
            if ((is.null(temp) || recalculate) && (nrow(object)>1) ) {
                points <- matrix(unlist(object), ncol=2)
                nearest <- nearestcpp(points, points, non_zero = TRUE)
                sp <- nearest$distance[nearest$index > -1]
                median(sp)
            }
            else
                temp
        }
    }
}

## traps object
polyID <- function (object)    {
    if (ms(object)) {
        polyID (object[[1]])
    }
    else {
        if (inherits(object,'traps')) {
            temp <- attr(object,'polyID')
            if (is.null(temp)) temp <- factor(1:nrow(object))   ## all different
            temp
        }
        else
        if (inherits(object,'capthist')) {
            stop ("use trap() to extract polyID from 'capthist' object")
        }
        else stop ("polyID requires 'traps' object")
    }
}

## traps object
transectID <- function (object)    {
    if (ms(object)) {
        transectID (object[[1]])
    }
    else {
        if (inherits(object,'traps')) {
            if (!all(detector(object) %in% c('transect','transectX')))
                stop ("requires transect detector")
            temp <- attr(object,'polyID',exact = TRUE)
            if (is.null(temp)) temp <- factor(1:nrow(object))
            temp
        }
        else
        if (inherits(object,'capthist')) {
            stop ("use trap() to extract transectID from 'capthist' object")
        }
        else stop ("transectID requires 'traps' object")
    }
}

xy <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, xy)
    }
    else {
        if (any(detector(traps(object)) %in%
            c('polygonX', 'transectX', 'polygon','transect'))) {
            attr(object, 'detectedXY',exact = TRUE)
        }
        else
            NULL
    }
}

alongtransect <- function (object, tol = 0.01) {
    ptalongtransect <- function (i) {
        ## where is point i on its transect k?
        k <- trans[i]
        transectxy <- as.matrix(lxy[[k]])
        nr <- nrow(transectxy)
        alongtransectcpp (
            as.matrix (xyi[i,,drop=FALSE]),
            as.matrix (transectxy),
            as.integer (0),
            as.integer (nr-1),
            as.double (tol))
    }
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    
    if (ms(object)) {
        lapply(object, alongtransect, tol = tol)
    }
    else {
        trps <- traps(object)
        
        if (all(detector(trps) %in% c('transectX', 'transect'))) {
            
            #             trans <- trap(object, names = TRUE)
            #             xyi <- xy(object)
            #             ## 2015-09-02 change to fix occsim: remove 'S' prefix
            #             lxy <- split (trps, levels(transectID(trps)), prefix = "")
            
            trans <- trap(object, names = FALSE, sortorder = 'ksn')
            xyi <- xy(object)
            lxy <- split (trps, transectID(trps))
            
            sapply(1:nrow(xyi), ptalongtransect)
        }
        else
            NULL
    }
}

# under development 2021-12-11, 2022-02-01
distancetotransect <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    
    if (ms(object)) {
        lapply(object, distancetotransect)
    }
    else {
        trps <- traps(object)
        
        if (all(detector(trps) %in% c('transectX', 'transect'))) {
            
            xyi <- xy(object)
            
            # split by transect
            lxy <- split (trps, transectID(trps))
            vlist <- lapply(lxy, as.matrix)
            
            # each transect as sfg
            vlist <- lapply(vlist, st_linestring)
            # combine linestrings in one sfc
            v <- st_sfc(vlist)
            xy <- st_as_sf(as.data.frame(xyi), coords = 1:2)
            neari <- st_nearest_feature(xy, v)
            st_distance(xy,v)[,neari]   # distance to nearest feature

        }
        else
            NULL
    }
}



clusterID <- function (object) {
    if (ms(object)) {
        lapply(object, clusterID)
    }
    else {
        if (inherits(object, 'capthist')) {
            trps <- traps(object)
            clusterID(trps)[trap(object, names = FALSE, sortorder = 'snk')]
        }
        else
            attr(object, 'cluster', exact = TRUE)
    }
}

clustertrap <- function (object) {
    if (ms(object)) {
        lapply(object, clustertrap)
    }
    else {
        if (inherits(object, 'capthist')) {
            trps <- traps(object)
            clustertrap(trps)[trap(object, names = FALSE, sortorder = 'snk')]
        }
        else
        attr(object, 'clustertrap',exact = TRUE)
    }
}

signal.default <- function(object, ...) {
    stop ("only for capthist data")
}

signal.capthist <- function (object, ...) {
    if (ms(object)) {
        lapply(object, signal)
    }
    else {
        if (all(detector(traps(object)) %in% c('signal','signalnoise'))) {
            attr(object, 'signalframe', exact = TRUE)$signal
        }
        else
            NULL
    }
}

noise.default <- function(object, ...) {
    stop ("only for capthist data")
}

noise.capthist <- function (object, ...) {
    if (ms(object)) {
        lapply(object, noise)
    }
    else {
        if (all(detector(traps(object)) %in% c('signalnoise'))) {
            attr(object, 'signalframe',exact = TRUE)$noise
        }
        else
            NULL
    }
}

signalframe <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, signalframe)
    }
    else {
        if (all(detector(traps(object)) %in% c('signal','signalnoise'))) {
            attr(object, 'signalframe',exact = TRUE)
        }
        else
            NULL
    }
}

signalmatrix <- function (object, noise = FALSE, recodezero = FALSE,
    prefix = 'Ch', signalcovariates = NULL, names = NULL) {
    if (ms(object)){
        lapply(object, signalmatrix, noise = noise, recodezero = recodezero,
               prefix=prefix, signalcovariates = signalcovariates, names = names)
    }
    else {
        object <- check3D(object)
        tmpsignal <- object
        if (dim(tmpsignal)[2]>1)
            warning("using only first occasion")
        if (noise)
            sound <- noise(object)
        else
            sound <- signal(object)
        if (!is.null(sound))
            tmpsignal[(tmpsignal>0) | is.na(tmpsignal)] <- sound
        if (recodezero)
            tmpsignal[tmpsignal==0] <- NA
        tmpsignal <- tmpsignal[,1,,drop=FALSE]
        tmpsignal <- as.data.frame(tmpsignal)
        ## added 2013-09-03
        if (!is.null(names)) {
            if (length(names)==ncol(tmpsignal))
                names(tmpsignal) <- names
            else
                names(tmpsignal) <- paste(prefix, 1:ncol(tmpsignal), sep='')
        }
        else {
            names(tmpsignal) <- paste(prefix, rownames(traps(object)), sep='')
        }
        if (!is.null(signalcovariates)) {
            sf <- signalframe(object)
            firstID <- match(rownames(tmpsignal), sf$Selection)
            tmpsignal[,signalcovariates] <- sf[firstID,signalcovariates]
        }
        tmpsignal
    }
}

aliveold <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, alive)
    }
    else {
        if (is.matrix(object)) {
            sign(object)[abs(object)>0] > 0
        }
        else {
            temp <- sign(object)[abs(object)>0] > 0
            rep(temp, abs(object[abs(object)>0]))
        }
    }
}

captmatrix <- function (object, sortorder = c('snk','ksn')) {
  # 2020-08-27 
  # object is 3-D single-session capthist
  # form a matrix with indices of animalID, occasion and detector
    # number of rows is product of these dimensions
    if (length(dim(object)) < 3 ) {
        # stop ("3-D capthist required by animalID(), occasion() and trap()")
        if (length(dim(object))==2)
            object <- array(object, dim=c(dim(object),1))  ## 2020-12-08
    }
    sortorder <- match.arg(sortorder)
    if (nrow(object)<1) {
        stop ("captmatrix expects at least one row in capthist")      
    }
    tmp <- sapply(1:3, slice.index, x = object)
    if (!is.matrix(tmp)) tmp <- matrix(tmp, nrow = 1)
    # add column for number of detections 
    tmp <- cbind(tmp, as.integer(abs(object)))
    # drop unused combinations
    tmp <- tmp[tmp[,4]>0,, drop = FALSE]

    # construct index and sort
    if (sortorder == 'snk') { # occasion, animalID, detector
      ord <- order(tmp[,2], tmp[,1], tmp[,3])
    }
    else if (sortorder == 'ksn') {
      ord <- order(tmp[,3], tmp[,2], tmp[,1])
    }
    tmp[ord,, drop = FALSE]
}

occasion <- function (object, sortorder = c('snk','ksn')) {
  if (!inherits(object, 'capthist'))
    stop ("requires 'capthist' object")
  sortorder <- match.arg(sortorder)
  
  if (ms(object)) {
    lapply(object, occasion, sortorder = sortorder)
  }
  else {
    if (nrow(object) == 0) {
      out <- NULL
    }
    else {
      values <- 1:dim(object)[2]    # safer than following
      # values <- as.numeric(dimnames(object)[[2]])
      tmp <- captmatrix(object, sortorder)
      s <- tmp[,2]    # vector of occasion integer values
      out <- rep(values[s], tmp[,4])
    }
    as.numeric(out)
  }
}

alive <- function (object, sortorder = c('snk','ksn')) {
    # revised 2020-11-07
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
  sortorder <- match.arg(sortorder)
  
    if (ms(object)) {
        lapply(object, alive, sortorder = sortorder)
    }
    else {
        if (nrow(object) == 0) {
            out <- NULL
        }
        else {
            tmp <- captmatrix(object, sortorder)
            ## 2022-01-05 fixed problem with one animal
            ## temp <- sign(object[tmp[,1:3]]) > 0
            temp <- sign(object[tmp[,1:3, drop = FALSE]]) > 0
            out <- rep(temp, tmp[,4])
        }
        as.logical(out)
    }
}

trap <- function (object, names = TRUE, sortorder = c('snk','ksn')) {
  if (!inherits(object, 'capthist'))
    stop ("requires 'capthist' object")
  sortorder <- match.arg(sortorder)
  
  if (ms(object)) {
    lapply(object, trap, names = names, sortorder = sortorder)
  }
  else {
    trps <- traps(object)
    if (is.null(trps)) {
      out <- rep(1, sum(abs(object)))
    }
    else {
        if (nrow(object) == 0) {
            out <- NULL
        }
        else {  
            if (names) {
                values <- rownames(trps)  # but what about polyID?
            }
            else {
                values <- 1:dim(object)[3]
            }
            tmp <- captmatrix(object, sortorder)
            k <- tmp[,3]  # vector of detector ID
            out <- rep(values[k], tmp[,4])
        }
    }    
    if (names) as.character(out) else as.numeric(out)
  }
}
animalID <- function (object, names = TRUE, sortorder = c('snk','ksn')) {
  if (!inherits(object, 'capthist'))
    stop ("requires 'capthist' object")
  sortorder <- match.arg(sortorder)
  if (ms(object)) {
    lapply(object, animalID, names = names, sortorder = sortorder)
  }
  else {
    if (nrow(object) == 0) {
      out <- NULL
    }
    else {
      if (names & !is.null(row.names(object)))  ## 2011-08-18 null check
        values <- row.names(object)
      else
        values <- 1:nrow(object)
      tmp <- captmatrix(object, sortorder)
      n <- tmp[,1]  # vector of animalID
      out <- rep(values[n], tmp[,4])
    }
    if (names) as.character(out) else as.numeric(out)
  }
}

detectionindex <- function (object) {
## detectionindex is non-exported function 2012-02-11
## to which original cell in dim3 capthist object does a detection relate?
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, detectionindex)
    }
    else {
        if (nrow(object) == 0)
           numeric(0)
        else {
            object <- check3D(object)
            x <- object
            x[] <- 1:length(x)
            x[abs(object)==0] <- 0
            x <- aperm(x, c(3,2,1))
            x <- x[x>0]
            count <- aperm(object, c(3,2,1))
            count <- count[abs(count)>0]
            rep(x,count)
        }
    }
}

polyarea <- function (xy, ha = TRUE) {
    if (inherits(xy, c('sf','sfc','sfg'))) {
        temparea <- sum(st_area(xy))
    }
    else if (inherits(xy, 'SpatialPolygons')) {
        temparea <- sum(st_area(st_as_sf(xy)))
    }
    else if (inherits(xy, 'SpatVector')) {
        temparea <- sum(st_area(st_as_sf(xy)))
    }
    else {
        nr <- length(xy$x)
        ## beat problem with integer overflow 2010-11-21
        xy$x <- as.double(xy$x)
        xy$y <- as.double(xy$y)
        sc <- sum (xy$x[-nr] * xy$y[-1] - xy$x[-1] * xy$y[-nr])
        temparea <- abs(sc/2)
    }
    ifelse (ha, temparea/10000, temparea)
}

searcharea <- function (object)    {
    ## requires traps object
    ## discarded some obsolete code 2016-10-08
    if (ms(object)) {
        ## 2011-06-24
        lapply(object, searcharea)
    }
    else {
        if (all(detector(object) %in% c('polygon','polygonX'))) {
            ## assume poly closed
            sapply(split(object, levels(polyID(object))), polyarea)
        }
        else NA
    }
}

transectlength <- function (object)    {
## requires traps object
    if (ms(object)) {
        ## 2011-06-24
        lapply(object, transectlength)
    }
    else {
        if (all(detector(object) %in% c('transect','transectX'))) {
            calclength <- function (xy) {
                nr <- nrow(xy)    ## number of vertices
                segments <- (xy$x[-nr] - xy$x[-1])^2 + (xy$y[-nr] - xy$y[-1])^2
                sum (segments^0.5)
            }
            sapply(split(object, polyID(object)), calclength)
        }
        else NA
    }
}

session.default <- function (object, ...)     {
## bypass session attribute for multi-session objects 2009 12 22
## use names(object) for lists i.e. multi-session objects

    if (ms(object)) {
       temp <- names(object)
    }
    else {
        temp <- attr(object,'session', exact = TRUE)
        if (is.null(temp)) temp <- 1    ## added 2010-02-03
    }
    names(temp) <- NULL
    as.character(temp)      ## 2010 02 25
}

trim.default <- function (object, drop, keep)     {
## drop unwanted named components of a list
## conservative resolution of conflicts between drop & keep
    objnames <- names(object)
    indices <- 1:length(object)
    if (missing(drop)) drop <- indices  # all! but wait...
    if (missing(keep)) keep <- 0
    if (is.character(keep)) keep <- match(keep, objnames)
    if (is.character(drop)) drop <- match(drop, objnames)
    drop <- drop[drop %in% indices[! (indices %in% keep)]]
    ## by index, so have to work from end
    for (i in sort(drop, decreasing = T)) object[[i]] <- NULL
    object
}

###############################################################################

rotate.default <- function (object, degrees, centrexy=NULL, ...) {

    rotatefn <- function (xy) {
        # about centre
        x <- xy[1] - centrexy[1]
        y <- xy[2] - centrexy[2]
        x2 <- x * cos(theta) + y * sin(theta) + centrexy[1]
        y2 <- - x * sin(theta) + y * cos(theta) + centrexy[2]
        c(x2,y2)
      }

##    if (ms(object)) lapply(object, rotate.default, degrees=degrees, centrexy=centrexy, ...)
##    else

    object <- as.matrix(object)
    if (dim(object)[2] <2)
        stop ("requires at least 2 columns")
    if (abs(degrees)>0) {
        if (is.null(centrexy)) centrexy <- c(0,0)
        theta <- 2*pi*degrees/360 # convert to radians
        temp <- t(apply (object,1,rotatefn))
        if (!is.null(dimnames(object)[[2]]))
            dimnames(temp)[[2]] <- dimnames(object)[[2]][1:2]
        temp
    } else object[,1:2]
}
###############################################################################

shift.default <- function (object, shiftxy, ...) {
##    if (ms(object)) lapply(object, shift.default, shiftxy=shiftxy, ...)
##    else

    object <- as.matrix(object[,1:2])
    object[,1] <- object[,1] + shiftxy[1]
    object[,2] <- object[,2] + shiftxy[2]
    object
}
###############################################################################

flip.default <- function (object, lr = FALSE, tb = FALSE, ...) {
##    if (ms(object)) lapply(object, flip.default, lr=lr, tb=tb, ...)
##    else
    object <- as.matrix(object[,1:2])
    if (is.logical(lr)) {
        if (lr) object[,1] <- 2 * mean(object[,1]) - object[,1]  ## flip about mean
    } else
        if (is.numeric(lr)) object[,1] <- 2*lr-object[,1]  ## flip about lr

    if (is.logical(tb)) {
        if (tb) object[,2] <- 2 * mean(object[,2]) - object[,2]  ## flip about mean
    } else
        if (is.numeric(tb)) object[,2] <- 2*tb-object[,2]  ## flip about tb

    object
}
###############################################################################

# Generic methods for replacing values

'usage<-' <- function (object, value) {

    ## extended 2012-12-22
    ## bug : loses session names when called as usage(traps(ovenCH)) <- c(1,10)

    if (ms(object)) {
        ## replicate across sessions as required
        if (!is.list(value)) {
            nonmatrix <- length(dim(value)) != 2
            valuelist <- vector(mode='list', length=length(object))
            if (nonmatrix & (length(value)>2)) {
                if (length(value) != (length(valuelist)+1))
                    stop("invalid value vector - how many occasions?")
                for (i in 1:length(valuelist)) valuelist[[i]] <- value[c(1,i+1)]
            }
            else
                for (i in 1:length(valuelist)) valuelist[[i]] <- value
        }
        else
            valuelist <- value
        temp <- mapply('usage<-', object, valuelist, SIMPLIFY = FALSE)
        class(temp) <- c('traps', 'list')
        temp
    }
    else {
        if (!is.null(value)) {
            nonmatrix <- length(dim(value)) != 2
            if (nonmatrix) {
                value <- unlist(value)
                if (length(value) == 1) {
                    if (is.null(usage(object)))
                        stop("requires number of occasions")
                    else noccasions <- ncol(usage(object))
                }
                else
                    noccasions <- value[2]
                value <- value[1]
            }
            else {
                noccasions <- NULL
                value <- as.matrix(value)  ## coerce dataframe to matrix
            }
            if (is.null(noccasions)) {
                if (is.null(dim(value)))
                    stop("usage value should be matrix or scalar with noccasions")
                else if (nrow(value) != ndetector(object))
                    stop("mismatch between value and number of detectors")
            }
            else {
                ## assume value is scalar
                value <- matrix(value, nrow = ndetector(object), ncol = noccasions)
            }
            if (!is.null(markocc(object))) {
                if (ncol(value) != length(markocc(object)))
                    stop ("replacement value not compatible with existing markocc attribute")
            }
        }
        structure (object, usage = value)
    }
}

'markocc<-' <- function (object, value) {

    if (ms(object)) {
        ## replicate across sessions as required
        if (!is.list(value)) {
            valuelist <- vector(mode='list', length=length(object))
            for (i in 1:length(valuelist)) valuelist[[i]] <- value
        }
        else
            valuelist <- value
        temp <- mapply('markocc<-', object, valuelist, SIMPLIFY = FALSE)
        class(temp) <- c('traps', 'list')
        temp
    }
    else {
        if (!inherits(object, 'traps'))
            stop("markocc attribute can only be assigned for traps objects")
        if (!is.null(value)) {
            if (!is.null(usage(object))) {
                if (length(value) != ncol(usage(object)))
                    stop ("markocc not compatible with existing usage attribute")
            }
            if (!all(value %in% c(-2,-1,0,1)))
                stop ("markocc values must be -2, -1, 0 or 1")
        }
        structure (object, markocc = value)
    }
}

'telemetrytype<-' <- function (object, value) {

    if (ms(object)) {
        temp <- mapply('telemetrytype<-', object, value, SIMPLIFY = FALSE)
        class(temp) <- c('traps', 'list')
        temp
    }
    else {
        if (!inherits(object, 'traps'))
            stop("telemetrytype attribute can only be assigned for traps objects")
        if (!is.null(value)) {
            if (!(value %in% c("none","independent","dependent","concurrent")))
                stop ("telemetrytype value must be none, independent, dependent or concurrent")
        }
        structure (object, telemetrytype = value)
    }
}

'Tu<-' <- function (object, value) {

    if (ms(object)) {
        if (!is.list(value)) {
            stop("replacement of Tu for multisession object requires a list")
        }
        else {
            temp <- mapply('Tu<-', object, value, SIMPLIFY = FALSE)
            class(temp) <- class(object)
            temp
        }
    }
    else {
        if (!is.null(value)) {
            if (is.null(markocc <- markocc(traps(object))))
                stop("cannot assign Tu for capthist in which traps has no markocc attribute")
            if (length(value)>1) {
                if (length(dims <- dim(value)) != 2)
                    stop ("require either total count or traps x occasions matrix")
                if (dims[1] != ndetector(traps(object)))
                    stop ("Tu not compatible with traps attribute")
                if (dims[2] != length(markocc))
                    stop ("Tu not compatible with markocc attribute")
            }
            if (any(value<0))
                stop ("sighting counts cannot be negative")
        }
        structure (object, Tu = value)
    }
}

'Tm<-' <- function (object, value) {

    if (ms(object)) {
        if (!is.list(value)) {
            stop("replacement of Tm for multisession object requires a list")
        }
        else {
            temp <- mapply('Tm<-', object, value, SIMPLIFY = FALSE)
            class(temp) <- class(object)
            temp
        }
    }
    else {
        if (!is.null(value)) {
            if (is.null(markocc <- markocc(traps(object))))
                stop("cannot assign Tm for capthist in which traps has no markocc attribute")
            if (length(value)>1) {
                if (length(dims <- dim(value)) != 2)
                    stop ("require either total count or traps x occasions matrix")
                if (dims[1] != ndetector(traps(object)))
                    stop ("Tm not compatible with traps attribute")
                if (dims[2] != length(markocc))
                    stop ("Tm not compatible with markocc attribute")
            }
            if (any(value<0))
                stop ("sighting counts cannot be negative")
        }
        structure (object, Tm = value)
    }
}

'Tn<-' <- function (object, value) {

    if (ms(object)) {
        if (!is.list(value)) {
            stop("replacement of Tn for multisession object requires a list")
        }
        else {
            temp <- mapply('Tn<-', object, value, SIMPLIFY = FALSE)
            class(temp) <- class(object)
            temp
        }
    }
    else {
        if (!is.null(value)) {
            if (is.null(markocc <- markocc(traps(object))))
                stop("cannot assign Tn for capthist in which traps has no markocc attribute")
            if (length(value)>1) {
                if (length(dims <- dim(value)) != 2)
                    stop ("require either total count or traps x occasions matrix")
                if (dims[1] != ndetector(traps(object)))
                    stop ("Tn not compatible with traps attribute")
                if (dims[2] != length(markocc))
                    stop ("Tn not compatible with markocc attribute")
            }
            if (any(value<0))
                stop ("sighting counts cannot be negative")
        }
        structure (object, Tn = value)
    }
}

'nontarget<-' <- function (object, value) {
    
    if (ms(object)) {
        if (!is.list(value)) {
            stop("replacement of nontarget for multisession object requires a list")
        }
        else {
            if (length(value) != length(object)) {
                stop ("require one non-target matrix for each session")
            }
            temp <- mapply('nontarget<-', object, value, SIMPLIFY = FALSE)
            class(temp) <- class(object)
            temp
        }
    }
    else {
        if (!is.null(value)) {
            if (length(value)>1) {
                if (length(dims <- dim(value)) != 2)
                    stop ("require traps x occasions matrix")
                if (dims[1] != ndetector(traps(object)))
                    stop ("nontarget not compatible with traps attribute")
            }
            if (any(value<0))
                stop ("nontarget values cannot be negative")
            value <- as.matrix(value)
        }
        structure (object, nontarget = value)
    }
}

# 'clusterID<-' <- function (object, value) structure (object, cluster = value)
# 'clustertrap<-' <- function (object, value) structure (object, clustertrap = value)

'covariates<-' <- function (object, value) {
## modified for multi-session data 2010-10-15
    if (is.null(value))
        structure (object, covariates = NULL)
    else {
        if (ms(object)) {
            if (length(object) != length(value))
                stop ("mismatch between multisession object and covariates")
            value <- lapply(value, as.data.frame)
        }
        else {
            value <- as.data.frame(value)
            nrequired <- nrow(object)
            if (inherits(object, 'traps'))
                if (all(detector(object) %in% .localstuff$polydetectors))
                    nrequired <- length(levels(polyID(object)))
            if (nrow(value) != nrequired)
                stop ("length of covariate does not match")
        }
        structure (object, covariates = value)
    }
}

'timevaryingcov<-' <- function (object, value) {
## 2012-10-31, modified 2013-02-07, 2018-05-14, 2019-09-29
    if (is.null(value))
        structure (object, timevaryingcov = NULL)
    else {
        if (ms(object)) {
            temp <- object
            if (!(is.list(value[[1]]) & (length(value) == length(object))))
                value <- list(value)
            for (i in 1:length(object))
                timevaryingcov(temp[[i]]) <- value[[i]]
            temp
        }
        else {
            # if (!inherits(object, 'traps'))
            #     stop("timevaryingcov is for traps objects")
            if (!is.list(value) | is.null(names(value)))
                stop("value should be a list of one or more named vectors")
            if (!is.null(usage(object))) {
                ## traps object with usage
                OK <- sapply(value, function(x) ncol(usage(object)) == length(x))
                if (any(!OK))
                    warning ("mismatch between number of occasions in usage and timevaryingcov")
            }
            else if (inherits(object, 'capthist')) {
                ## capthist object
                nsessions <- length(unique(primarysessions(intervals(object))))
                if (nsessions>1) {    # 2019-09-29
                    OK <- sapply(value, function(x) nsessions == length(x))
                    if (any(!OK))
                        warning ("mismatch between number of primary sessions in object and timevaryingcov")
                }
            }
            if (is.character(value))
                if (!all(value %in% names(covariates(object))))
                    warning ("mismatch between character vector and covariate names")
            structure (object, timevaryingcov = value)
        }
    }
}

'detector<-' <- function (object, value) {
    if (!all(value %in% .localstuff$validdetectors))
        stop ("invalid detector type")
    structure (object, detector = value)
}

'spacing<-' <- function (object, value) {
    ## if (!(is.numeric(value)))
    if (!(is.na(value) || is.numeric(value)))   ## allow NA 2019-01-22
        stop ("non-numeric spacing")
    if (ms(object)) {
        stop ("not sure how to replace spacing of ms object")
    }
    else {
        structure (object, spacing = value)
    }
}

'polyID<-' <- function (object, value) {
    if (!inherits(object,'traps'))
        warning ("polyID requires 'traps' object")
    if (length(value)==1) value <- rep(value, nrow(object))
    value <- factor(value)
    structure (object, polyID = value)
}

'transectID<-' <- function (object, value) {
    if (length(value)==1) value <- rep(value, nrow(object))
    value <- factor(value)
    structure (object, polyID = value)
}

'xy<-' <- function (object, value) {
    if (!is.null(value)) {
        object <- check3D(object)
        polyoccasions <- expanddet(object) %in% .localstuff$polydetectors
        ndetections <- sum(abs(object[,polyoccasions,]))
        if (nrow(value) != ndetections)
            stop ("requires one location per detection")
        if ((sum(polyoccasions)==0) |
            !(inherits(object,'capthist')))
            stop ("requires 'capthist' object with ",
                  "polygon-like detector")
        if (ms(object))
            stop ("requires single-session 'capthist' object")
    }
    structure (object, detectedXY = value)
}

'telemetryxy<-' <- function (object, value) {
    if (!is.null(value)) {
        object <- check3D(object)
        telemoccasions <- expanddet(object) %in% 'telemetry'
        ndetections <- sum(abs(object[,telemoccasions,]))
        if (sum(sapply(value,nrow)) != ndetections)
            stop ("requires one location per detection")
        if ((sum(telemoccasions)==0) |
            !(inherits(object,'capthist')))
            stop ("requires 'capthist' object with ",
                  "telemetry detector")
        if (ms(object))
            stop ("requires single-session 'capthist' object")
    }
    structure (object, telemetryxy = value)
}

'clusterID<-' <- function (object, value) {

    if (ms(object))
        stop ("cluster requires single-session 'traps' object")

    if (length(value)==1) value <- rep(value, nrow(object))

    if (!(inherits(object, 'traps')))
        stop ("requires clustered 'traps' object")

    if (length(value) > 0) {
        if (length(value) != nrow(object))
            stop ("requires one cluster per detector or detector vertex")
        value <- factor(value)   ## moved here 2017-11-15
    }
    else
        value <- NULL

    structure (object, cluster = value)
}

'clustertrap<-' <- function (object, value) {

    if (ms(object))
        stop ("clustertrap requires single-session 'traps' object")

    if (!(inherits(object, 'traps')))
        stop ("requires clustered 'traps' object")

    if (length(value) > 0) {
        if (length(value) != nrow(object))
            stop ("requires one clustertrap per detector or detector vertex")
    }
    else
        value <- NULL

    structure (object, clustertrap = value)
}

'signalframe<-' <- function (object, value) {
    if (is.null(value)) {
        attr(object, 'signalframe') <- NULL
        object
    }
    else {
        value <- as.data.frame(value)
        if (nrow(value) != sum(abs(object)))
            stop ("requires one row per detection")
        if (!('signal' %in% names(value)))
            stop ("value does not contain column 'signal'")
        if (!(all(detector(traps(object)) %in% c('signal','signalnoise'))) |
            !(inherits(object,'capthist')))
            stop ("requires 'capthist' object with 'signal' or 'signalnoise' detector")
        if (ms(object))
            stop ("requires single-session 'capthist' object")
        structure (object, signalframe = value)
    }
}

'signal<-' <- function (object, value) {
    if (length(value) != sum(abs(object)))
        stop ("requires one signal per detection")
    if (!(all(detector(traps(object)) %in% c('signal','signalnoise'))) |
        !(inherits(object,'capthist')))
        stop ("requires 'capthist' object with 'signal' or 'signalnoise' detector")
    if (ms(object))
        stop ("requires single-session 'capthist' object")
    sf <- attr(object, 'signalframe',exact = TRUE)
    if (is.null(sf)) {
        sf <- data.frame(signal = value)
    }
    else
        sf$signal <- value
    structure (object, signalframe = sf)
}

'noise<-' <- function (object, value) {
    if (length(value) != sum(abs(object)))
        stop ("requires one noise per detection")
    if (!((detector(traps(object)) %in% c('signalnoise'))) |
        !(inherits(object,'capthist')))
        stop ("requires 'capthist' object with 'signalnoise' detector")
    if (ms(object))
        stop ("requires single-session 'capthist' object")
    sf <- attr(object, 'signalframe',exact = TRUE)
    if (is.null(sf))
        sf <- data.frame(noise = value)
    else
        sf$noise <- value
    structure (object, signalframe = sf)
}

'traps<-' <- function (object, value) {
    if (!is(value,'traps'))
        stop ("'traps' object required for replacement")
    ## MODIFIED 2010-04-27 2024-02-17
    if (ms(object)) {
        nsess <- length(object)
        temp <- vector(mode='list', nsess)
        if (nsess != length(value))
            stop ("replacement value has wrong length")
        for (i in 1:nsess) temp[[i]] <- `traps<-`(object[[i]], value[[i]])
        class(temp) <- class(object)   # capthist
        temp
    }
    else {
        structure (object, traps = value)
    }
}

'session<-' <- function (object, value) {
    if (ms(object)) {
       if (length(value) != length(object))
           stop ("invalid replacement value")
       for (i in 1:length(object)) session(object[[i]]) <- value[i]   ## 2010 03 26
       structure (object, names = as.character(value))
    }
    else {
        if (length(value) > 1)
            stop ("requires only one session name per session")
        structure (object, session = as.character(value))
    }
}
###############################################################################

######################################
## Class : traps
## defines an array of detectors
## detector may be 'single', 'multi', 'proximity' etc.
######################################

## 2011-10-10 make.grid, makepoly etc. moved to make.grid.R

rotate.traps <- function (object, degrees, centrexy=NULL, ...)
{
##    if (ms(object)) lapply(object, rotate.traps, degrees, centrexy, ...)
##    else

  rotatefn <- function (xy) {
    # about centre
    x <- xy[1] - centrexy[1]
    y <- xy[2] - centrexy[2]
    x2 <- x * cos(theta) + y * sin(theta) + centrexy[1]
    y2 <- - x * sin(theta) + y * cos(theta) + centrexy[2]
    c(x2,y2)
  }

  if (abs(degrees)>0 && nrow(object)>0) {
    if (is.null(centrexy)) centrexy <- c(mean(object$x), mean(object$y))
    theta <- 2*pi*degrees/360 # convert to radians

    traps2 <- data.frame(t(apply (object,1,rotatefn)))
    names(traps2) <- c('x','y')
    attr(traps2,'class')  <- c('traps', 'data.frame')
    detector(traps2)      <- detector(object)
    if (!is.null(usage(object)))
        usage(traps2)         <- usage(object)
    if (!is.null(covariates(object)))
        covariates(traps2) <- covariates(object)
    if (!is.null(timevaryingcov(object)))
        timevaryingcov(traps2) <- timevaryingcov(object)
    if (!is.null(polyID(object)))   ## includes transectID
        polyID(traps2)    <- polyID(object)
  }
  else traps2 <- object
  traps2
}
###############################################################################

rotate.capthist <- function(object, degrees, centrexy=NULL, ...) {
    if (ms(object)) {
        out <- lapply(object, rotate, degrees, centrexy, ...)
        class(out) <- c('capthist', 'list')
        out
    }
    else {
        traps(object) <- rotate(traps(object), degrees, centrexy, ...)
        object
    }
}

shift.traps <- function (object, shiftxy, ...)
{
  object$x <- object$x + shiftxy[1]
  object$y <- object$y + shiftxy[2]
  object
}
###############################################################################

shift.mask <- function (object, shiftxy, ...)
{
  object$x <- object$x + shiftxy[1]
  object$y <- object$y + shiftxy[2]
  bbox <-  attr(object, 'boundingbox',exact = TRUE)
  attr(object, 'boundingbox') <- data.frame(secr::shift(bbox, shiftxy))
  object
}
###############################################################################

flip.traps <- function (object, lr = FALSE, tb = FALSE, ...) {

##    if (ms(object)) lapply(object, flip.traps, lr, tb, ...)
##    else

    if (is.logical(lr)) {
        if (lr) object$x <- 2 * mean(object$x) - object$x  ## flip about centre
    } else
        if (is.numeric(lr)) object$x <- 2*lr - object$x  ## flip about lr

    if (is.logical(tb)) {
        if (tb) object$y <- 2 * mean(object$y) - object$y  ## flip about centre
    } else
        if (is.numeric(tb)) object$y <- 2*tb - object$y  ## flip about tb

    object
}
###############################################################################

rotate.popn <- function (object, degrees, centrexy=NULL, ...)
{
  rotatefn <- function (xy) {
    # about centre
    x <- xy[1] - centrexy[1]
    y <- xy[2] - centrexy[2]
    x2 <- x * cos(theta) + y * sin(theta) + centrexy[1]
    y2 <- - x * sin(theta) + y * cos(theta) + centrexy[2]
    c(x2,y2)
  }
  bbox <-  attr(object, 'boundingbox',exact = TRUE)
  if (abs(degrees)>0) {
    if (is.null(centrexy)) centrexy <- c(mean(bbox$x), mean(bbox$y))
    theta <- 2*pi*degrees/360 # convert to radians

    popn2 <- data.frame(t(apply (object,1,rotatefn)))

    object[,] <- popn2[,]
    attr(object, 'boundingbox') <- data.frame(rotate (bbox, degrees, centrexy))
  }
  object
}
###############################################################################

shift.popn <- function (object, shiftxy, ...)
{
  object$x <- object$x + shiftxy[1]
  object$y <- object$y + shiftxy[2]
  bbox <-  attr(object, 'boundingbox',exact = TRUE)
  attr(object, 'boundingbox') <- data.frame(secr::shift(bbox, shiftxy))
  object
}
###############################################################################

flip.popn <- function (object, lr = FALSE, tb = FALSE, ...) {
    bbox <- attr(object, 'boundingbox',exact = TRUE)
    if (is.logical(lr)) {
        if (lr) object$x <- 2 * mean(bbox$x) - object$x  ## flip about centre
    } else
        if (is.numeric(lr)) object$x <- 2*lr - object$x  ## flip about lr

    if (is.logical(tb)) {
        if (tb) object$y <- 2 * mean(bbox$y) - object$y  ## flip about centre
    } else
        if (is.numeric(tb)) object$y <- 2*tb - object$y  ## flip about tb
    attr(object, 'boundingbox') <- data.frame(flip(bbox, lr=lr, tb=tb, ...))
    object
}
###############################################################################

print.traps <- function(x, ...) {
    if (ms(x)) {
        for (i in 1:length(x)) {
            cat('\n')
            if (!is.null(names(x)))
                cat(names(x)[i], '\n')
            print (x[[i]], ...)
        }
        invisible()
    }
    else {
        temp <- data.frame(row.names=attr(x,'row.names',exact = TRUE), x=x$x, y=x$y)
        print(temp, ...)
    }
}
###############################################################################

subset.traps <- function (x, subset = NULL, occasions = NULL, ...) {
    # subset may be numeric index or logical

    if (ms(x)) {
        temp <- lapply(x, subset.traps, subset=subset, occasions = occasions, ...)
        class(temp) <- c('traps', 'list')
        temp
    }
    else {
        if (is.null(subset)) {
            if (all(detector(x) %in% .localstuff$polydetectors))
                subset <- 1:length(levels(polyID(x)))
            else
                subset <- 1:nrow(x)
        }
        
        if (is.function(subset)) subset <- subset(x, ...)   ## 2020-08-12
        
        ## polygon & transect objects subset by whole polygons or transects
        ## 2011-01-24
        rowsubset <- subset  ## default
        if (all(detector(x) %in% c('polygon', 'polygonX')))
            rowsubset <- as.character(polyID(x)) %in% levels(polyID(x))[subset]
        if (all(detector(x) %in% c('transect', 'transectX')))
            rowsubset <- as.character(transectID(x)) %in% levels(transectID(x))[subset]
        
        ## apply subsetting
        temp <- x[rowsubset,,drop=F]
        class(temp) <- c('traps','data.frame')
        detector(temp) <- detector(x)

        ## 2011-05-09
        if (!is.null(clusterID(x))) {
            clusterID(temp) <- factor(clusterID(x)[rowsubset])
            clustertrap(temp) <- factor(clustertrap(x)[rowsubset])
        }

        ## restore polyiD, transectID, usage, covariates
        if (all(detector(x) %in% c('polygon', 'polygonX')))
            polyID(temp) <- factor(polyID(x)[rowsubset])
        if (all(detector(x) %in% c('transect', 'transectX')))
            transectID(temp) <- factor(transectID(x)[rowsubset])

        if (!is.null(usage(x))) {
            if (is.null(occasions))
                occasions <- 1:ncol(usage(x))
            usage(temp) <- usage(x)[subset,occasions,drop=F]
        }
        if (!is.null(markocc(x))) {
            if (is.null(occasions))
                occasions <- 1:length(markocc(x))
            markocc(temp) <- markocc(x)[occasions]
        }
        if (length(detector(x))>1) {
            if (is.null(occasions))
                occasions <- 1:length(markocc(x))
            detector(temp) <- detector(x)[occasions]
        }
        if (!is.null(covariates(x)))
            covariates(temp) <- covariates(x)[subset,,drop=F]

        if (!is.null(timevaryingcov(x))) {
            timevaryingcov(temp) <- lapply(timevaryingcov(x),
                function(y) y[occasions,drop=FALSE])
        }

        if (all(detector(x) %in% .localstuff$pointdetectors)) {
            spacing(temp) <- spacing(temp, recalculate = TRUE)
        }

        temp
    }
}
###############################################################################


####################################
## Class : capthist
## capture data
####################################

###############################################################################

subset.popn <- function (x, subset = NULL, sessions = NULL, poly = NULL,
    poly.habitat = TRUE, keep.poly = TRUE, renumber = FALSE, ...)
## x - popn object
## subset - vector (character, logical or integer) to subscript rows (dim(1)) of x
## sessions - vector (integer or logical) to subscript sessions
{
    if (ms(x)) {
        if (!is.null(sessions))
            x <- x[sessions]
        out <- vector('list')
        for (i in 1:length(x)) {
            if (is.list(subset))
                sset <- subset[[i]]
            else
                ## sset <- subset
                ## 2018-06-26
                sset <- rep(list(subset), length(x))
            if (is.null(subset)) { ## 2015-03-17
                out[[i]] <- x[[i]]
            }
            else {
                ## 2018-06-26
                ## out[[i]] <- subset(x[[i]], subset[[i]], NULL, renumber, ...)
                out[[i]] <- subset(x[[i]], subset = sset[[i]], sessions = NULL, 
                                   renumber = renumber, ...)
            }
        }
        class(out) <- c('popn', 'list')
        out
    }
    else {
        #-------------------------
        # default subset is all
        if (is.null(subset))
            subset <- 1:nrow(x)
        #-------------------------
        # restrict to a polygon
        # added 2011-10-20

        if (!is.null(poly)) {
            OK <- pointsInPolygon(x, poly)
            if (!poly.habitat)
                OK <- !OK
            subset <- subset[OK]
        }
        
        #-------------------------
        # apply subsetting
        pop <- x[subset,]
        #-------------------------
        if (renumber)
            rownames(pop) <- 1 : nrow(pop)

        class(pop) <- c('popn', 'data.frame')
        attr(pop, 'Ndist') <- NULL     ## no longer known
        attr(pop, 'model2D') <- NULL   ## no longer known
        attr(pop, 'boundingbox') <- attr(x, 'boundingbox')
        if (!is.null(poly) & keep.poly) {
            attr(pop, 'polygon') <- poly
            attr(pop, 'poly.habitat') <- poly.habitat
        }
        if (!is.null(covariates(x))) {
            covariates(pop) <- covariates(x)[subset,,drop=FALSE]
        }
        pop
    }
}

###############################################################################

subset.capthist <- function (x, subset=NULL, occasions=NULL, traps=NULL,
    sessions=NULL, cutval=NULL, dropnullCH=TRUE, dropnullocc=FALSE,
    dropunused = TRUE, droplowsignals = TRUE, dropNAsignals = FALSE,
    cutabssignal = TRUE, renumber=FALSE, ...)  {

## x - capthist object (array with 2 or 3 dimensions)
## subset - vector (character, logical or integer) to subscript rows (dim(1)) of x
## occasions - vector (integer or logical) to subscript occasions
## traps - vector (character, integer or logical) to subscript rows of traps object
## sessions - vector (integer or logical) to subscript sessions
    if (ms(x)) {
        if (is.null(sessions)) sessions <- 1:length(x)
        temp <- lapply (x[sessions], subset.capthist,    ## 2017-11-13 avoid calling subset fn argument
            subset = subset,
            occasions = occasions,
            traps = traps,
            sessions = sessions,   ## inserted 2009 10 01
            cutval = cutval,
            dropnullCH = dropnullCH,
            dropnullocc = dropnullocc,
            dropunused = dropunused,
            droplowsignals = droplowsignals,
            dropNAsignals = dropNAsignals,       ## 2017-01-28
            cutabssignal = cutabssignal,         ## 2017-01-28
            renumber = renumber, ...)
        class(temp) <- c('capthist', 'list')
        if (length(temp) == 1) temp <- temp[[1]]  ## 2009 09 25
        interv <- intervals(x)
        if (!is.null(interv)) {     ## 2018-01-25
            cumi <- cumsum(interv)
            newinterv <- diff(c(0,cumi)[sessions])
            intervals(temp) <- newinterv
        }
        slabels <- sessionlabels(x)
        if (!is.null(slabels)) {    ## 2018-05-10
            sessionlabels(temp) <- slabels[sessions]
        }

        return(temp)
    }
    else {
        rownum <- function (x) {
            if (length(dim(x)) < 1 || dim(x)[1] == 0) NULL
            else 1: (dim(x)[1])
        }
        colnum <- function (x) {
            if (length(dim(x)) < 2 || dim(x)[2] == 0) NULL
            else 1: (dim(x)[2])
        }
        if (is.matrix(x))
            stop("require updated capthist secr >= 3.0")
        trapsx <- secr::traps(x)
        if (is.null(trapsx)) {
           detector <- 'nonspatial'
           nk <- 1
        }
        else {
            detector <- expanddet(x)   # vector 2017-01-28
            nk <- ndetector (trapsx)
        }
        if (is.logical(subset) & (length(subset) != nrow(x)))
            stop ("if 'subset' is logical its length must match number of animals")
        if (is.logical(occasions) & (length(occasions) != ncol(x)))
            stop ("if 'occasions' is logical its length must match number of occasions")
        if (is.logical(traps) & (length(traps) != nk))
            stop ("if 'traps' is logical its length must match number of detectors")
        ## 2022-01-04 if (is.null(occasions)) occasions <- 1:ncol(x)
        if (is.null(occasions)) occasions <- colnum(x)
        if (is.null(traps))  traps <- 1:nk
        if (is.function(subset)) subset <- subset(x, ...)   ## 2017-11-13
        if (is.null(subset)) subset <- 1:nrow(x)

        # if (dropnullCH & (telemetrytype(trapsx) %in% c('independent','concurrent')) &
        #     any(detector[occasions] == 'telemetry'))
        #     warning("dropnullCH = TRUE is probably a mistake when telemetry ",
        #             "independent or concurrent")

        #############################################
        ## coerce subset, traps, occasions to logical
        if (is.character(subset))
            subset <- dimnames(x)[[1]] %in% subset
        else {
            if (is.numeric(subset)) {
                if (all(subset<0))   # negative implies exclusion 2020-09-03
                  ## 2022-01-04 subset <- !(1:nrow(x)) %in% abs(subset)
                  subset <- !(rownum %in% abs(subset))
                else
                  ## 2022-01-04 subset <- (1:nrow(x)) %in% subset
                  subset <- rownum(x) %in% subset
            }
        }
        
        if (is.character(traps)) {
            traps <- rownames(trapsx) %in% traps
        }
        else if (is.function(traps)) {
            traps <- traps(trapsx, ...)
        }            
        else if (!is.logical(traps)) {
            traps <- (1:nk) %in% traps
        }
        if (!is.logical(occasions)) {
          ## 2022-01-04 occasions <- (1:ncol(x)) %in% occasions
          occasions <- colnum(x) %in% occasions
        }
        
        if (any(detector=='telemetry')) {
            if (any(detector[occasions]=='telemetry')) {
                if (!traps[nk])
                    stop ("cannot drop notional telemetry detector from telemetry capthist")
            }
            else traps[nk] <- FALSE  ## forcibly drop notional telemetry detector
        }
        
        #####################################
        ## signaldf is used later...
        if (all(detector %in% c('signal','signalnoise'))) {
          trap <- trap(x, names = F, sortorder = 'ksn')
          occ <- occasion(x, sortorder = 'ksn')
          ID <- animalID(x, names = F, sortorder = 'ksn')
          signaldf <- data.frame(trap = trap, occ = occ, ID = ID,
            attr(x, 'signalframe',exact = TRUE))
          
          #####################################
          ## apply signal threshold if relevant
          if (is.null(cutval)) cutval <- attr(x, 'cutval',exact = TRUE)
          if (cutval < attr(x, 'cutval',exact = TRUE))
            stop ("cannot decrease 'cutval'")
          
          if (cutabssignal) {
            signalOK <- (signal(x) >= cutval)
          }
          else {
            if (is.null(noise(x)))
              stop("could not find noise for relative signal cut")
            signalOK <- ((signal(x)-noise(x)) >= cutval)
          }
          signalOK <- ifelse(is.na(signalOK), !dropNAsignals, signalOK)
          newcount <- table(
            factor(ID, levels = 1:nrow(x))[signalOK],
            factor(occ, levels = 1:ncol(x))[signalOK],
            factor(trap, levels = 1:nk)[signalOK])
          if (nrow(x)>0)  ## 2013-08-17
            x[] <- newcount * sign(x)  ## retain deads, in principle
        }
        
        ###########################
        ## condition missing values
        x[is.na(x)] <- 0
        #################################
        ## optionally drop traps never used on the specified occasions
        if ((detector[1] != 'nonspatial') & dropunused & !is.null(usage(trapsx))) {
            used <- apply(usage(trapsx)[,occasions, drop=FALSE],1,sum) > 0
            traps <- traps & used
        }

        ##############################################
        ## which cells are retained?
        ## for signalframe 2012-02-11
        i <- x
        if (nrow(x) > 0) {
            i[] <- 1:length(i)
            i <- i[subset, occasions, traps, drop = FALSE]
        }

        #################################
        ## prepare to drop null histories
        if (dropnullCH) {
            nonnull <- apply(abs(x[,occasions, traps, drop=FALSE]),1,sum) > 0
        }
        else
            nonnull <- rep(TRUE, nrow(x))
        subset <- subset & !is.na(subset) & nonnull
        if (nrow(x)==0) subset <- 0

        #################################
        ## perform main subset operation
        temp <- x[subset, occasions, traps, drop = FALSE]
        nocc <- ncol(temp)

        #################################
        ## drop null occasions
        OK2 <- rep(T,nocc)
        if (nrow(temp)>0)
        if ((!all(detector %in% c('signal','signalnoise'))) &
            any( apply(abs(temp),2,sum) ==0)) {
            if (dropnullocc) {
                OK2 <- apply(abs(temp),2,sum) > 0
                temp <- temp[,OK2,, drop = FALSE]
                i <- i[,OK2,, drop = FALSE]
            }
            else
                warning ("no detections on occasion(s) ",
                    paste((1:nocc)[apply(abs(temp),2,sum) ==0], collapse=', '), "\n")
            nocc <- dim(temp)[2]  # refresh
        }
        ###################################
        ## attributes
        class(temp) <- 'capthist'
        if (!is.null(trapsx)) { # spatial data
            secr::traps(temp) <- subset (trapsx, traps)
            usage(secr::traps(temp)) <- NULL  ## until we fix markocc
        }
        covariates(temp) <- covariates(x)[subset,,drop = FALSE]
        session(temp) <- session(x)
        attr(temp, 'n.mash') <- attr(x, 'n.mash',exact = TRUE)
        attr(temp, 'centres') <- attr(x, 'centres',exact = TRUE)
        ###################################
        ## mark-resight
        if (!is.null(markocc(trapsx))) {
            markocc(secr::traps(temp)) <- markocc(trapsx)[occasions]
            if (!is.null(Tu(x))) {
                if (is.matrix(Tu(x))) {
                    Tu(temp) <- Tu(x)[traps, occasions, drop = FALSE]
                }
                else {
                    Tu(temp) <- Tu(x)
                    warning ('sightings not separated by occasion; all included')
                }
            }
            if (!is.null(Tm(x))) {
                if (is.matrix(Tm(x))) {
                    Tm(temp) <- Tm(x)[traps, occasions, drop = FALSE]
                }
                else {
                    Tm(temp) <- Tm(x)
                    warning ('nonID sightings not separated by occasion; all included')
                }
            }
            if (!is.null(Tn(x))) {
                if (is.matrix(Tn(x))) {
                    Tn(temp) <- Tn(x)[traps, occasions, drop = FALSE]
                }
                else {
                    Tn(temp) <- Tn(x)
                    warning ('unresolved sightings not separated by occasion; all included')
                }
            }
        }

        if (!is.null(secr::traps(temp))) {
            usage(secr::traps(temp)) <- usage(trapsx)[traps, occasions,
                  drop = FALSE][,OK2, drop = FALSE]  ## drop null occasions
            ## 2022-08-09
            attr(temp, 'nontarget') <- attr(x, 'nontarget')[traps, occasions,
                drop = FALSE][,OK2, drop = FALSE]  ## drop null occasions
        }

        if (length(detector(trapsx))>1)
            detector(secr::traps(temp)) <- detector(trapsx)[occasions]

        ###################################
        ## telemetry
        if (!is.null(secr::traps(temp))) {
            xylist <- telemetryxy(x)
            if (!is.null(xylist) & any (detector(secr::traps(temp)) == 'telemetry')) {
                selectoccasions <- function (xy, id) {
                    originalocc <- rep(1:ncol(x), x[id,,nk])
                    xy[originalocc %in% which(occasions),]
                }
                xylist <- xylist[names(xylist) %in% rownames(x)[subset]]
                if (!all(occasions))
                    xylist <- mapply(selectoccasions, xylist, names(xylist), SIMPLIFY = FALSE)
                telemetryxy(temp) <- xylist
                if (all(detector(secr::traps(temp)) == 'telemetry'))
                    telemetrytype(secr::traps(temp)) <- "independent"
            }
            else telemetrytype(secr::traps(temp)) <- NULL
        }
        ###################################
        ## subset signal of signal capthist
        ## revised 2012-07-28
        if (all(detector %in% c('signal','signalnoise'))) {
            if (nrow(x)>0) {
                if (!droplowsignals) {
                    if (any(x<=0))
                        stop ("droplowsignals = FALSE cannot be applied to CH",
                              " objects with incomplete detection")
                    signalOK <- TRUE
                }
            ## otherwise signalOK remains a logical vector with length equal
            ## to the original number of detections
            ## subset, traps and occasions are logical vectors 2011-01-21, 2011-11-14
            OK <- occasions[signaldf$occ] & subset[signaldf$ID] & traps[signaldf$trap] & signalOK
            signaldf <- signaldf[OK,, drop = FALSE]
            signaldf <- signaldf[order(signaldf$trap, signaldf$occ, signaldf$ID),, drop = FALSE]
            attr(temp, 'signalframe') <- signaldf[,-(1:3), drop = FALSE]
        }
            attr(temp, 'cutval') <- cutval
        }

        ############################################
        ## subset xy of polygon or transect capthist
        if (!is.null(xy(x))) {
            df <- data.frame(trap = trap(x, names = FALSE, sortorder = 'ksn'),
                             occ = occasion(x, sortorder = 'ksn'),
                             ID = animalID(x, names = FALSE, sortorder = 'ksn'),
                             x=xy(x)[,1], y=xy(x)[,2])
            ## subset, traps and occasions are logical vectors 2011-01-21, 2011-11-14
            OK <- occasions[df$occ] & subset[df$ID] & traps[df$trap]
            df <- df[OK,, drop = FALSE]
            df <- df[order(df$trap, df$occ, df$ID),, drop=FALSE]
            attr(temp, 'detectedXY') <- df[,c('x','y')]
        }

        ################################################
        interv <- intervals(x)
        if (!is.null(interv)) {   ## 2018-01-25
            cumi <- cumsum(interv)
            newinterv <- diff(c(0,cumi)[occasions])
            intervals(temp) <- newinterv
        }
        ################################################
        slabels <- sessionlabels(x)
        if (!is.null(slabels)) {    ## 2018-05-10
            if (is.null(intervals(temp))) {
                warning ("sessionlabels but no intervals specified;",
                         " all intervals set to 1.0")
                intervals(temp) <- rep(1, ncol(temp)-1)
            }
            ## bug fixed 2021-04-22
            # sessions <- unique(primarysessions(intervals(temp))[occasions])
            sessions <- unique(primarysessions(interv)[occasions])
            sessionlabels(temp) <- slabels[sessions]
        }
        ################################################

        ## renumber if desired
        if (renumber) {
            if (length(dim(temp))==3)
                dimnames(temp) <- list(1:nrow(temp),1:nocc,NULL)   # renew numbering
            else
                dimnames(temp) <- list(1:nrow(temp),1:nocc)
        }
        temp
    }
}

###############################################################################
##
## MS.capthist and rbind.capthist removed to rbind.capthist.R 13/9/2011
##
###############################################################################

sort.capthist <- function (x, decreasing = FALSE, by = '', byrowname = TRUE, ...) {
    if (ms(x)) {
        newx <- vector(mode='list')
        for (i in 1:length(x)) {
            xi <- subset(x, session=i)
            newx[[i]] <- sort(xi, decreasing=decreasing, by=by)
        }
        temp <- do.call(MS.capthist, newx)
        names(temp) <- names(x)
        session(temp) <- session(x)
        temp
    }
    else {
        if (is.character(by)) {
            if (by == '') {
                by <- vector(mode='list')   ## zero-length
            }
            else {
                if (!all(by %in% names(covariates(x))))
                    stop ("unrecognised sort field(s)")
                by <- as.list(covariates(x)[,by, drop=FALSE])
            }
        }
        else {
            by <- as.list(data.frame(by))
        }
        if (byrowname) by$rownames <- row.names(x)
        by$decreasing <- decreasing
        roworder <- do.call(order, by)
        newrownames <- rownames(x)[roworder]
        tempx <- x   ## for detectionorder
        tempx[tempx!=0] <- 1:length(animalID(x))
        if (length(dim(x))==2) {
            tempx[,] <- tempx[roworder,]
            x[,] <- x[roworder,]
        }
        else {
            tempx[,,] <- tempx[roworder,,]
            x[,,] <- x[roworder,,]
        }
        rownames(x) <- newrownames
        ## covariates
        if (!is.null(covariates(x))) {
            temp <-covariates(x)[roworder,,drop=F]
            rownames(temp) <- newrownames
            covariates(x) <- temp
        }

        detectionorder <- tempx[tempx!=0]
        ## signal
        sf <- attr(x, 'signalframe',exact = TRUE)
        if (!is.null(sf))
            attr(x, 'signalframe') <- sf[detectionorder,,drop=FALSE]
        ## xy
        if (!is.null(xy(x)))
            xy(x) <- xy(x)[detectionorder,,drop=FALSE]
        ## return
        x
    }
}

sort.mask <- function (x, decreasing = FALSE, by = '', byrowname = TRUE, ...) {
    if (ms(x)) {
        oldclass <- class(x)
        if (inherits(by, 'list'))
            temp <- mapply(sort, x, by, MoreArgs = list(decreasing = decreasing, ...), simplify = FALSE)  ## need to evaluate ...!!!
        else
            temp <- lapply(x, sort, decreasing = decreasing, by = by, ...)
        class (temp) <- oldclass
        temp
    }
    else {
        if (is.character(by)) {
            if (by == '') {
                by <- vector(mode='list')   ## zero-length
            }
            else {
                if (!all(by %in% names(covariates(x))))
                    stop ("unrecognised sort field(s)")
                by <- as.list(covariates(x)[,by, drop=FALSE])
            }
        }
        else {
            by <- as.list(data.frame(by))
        }
        if (byrowname) by$rownames <- row.names(x)
        by$decreasing <- decreasing
        roworder <- do.call(order, by)
        newrownames <- rownames(x)[roworder]
        temp <- x[roworder,]
        rownames(temp) <- newrownames
        ## covariates
        if (!is.null(covariates(x))) {
            tempcov <-covariates(x)[roworder,,drop=F]
            rownames(tempcov) <- newrownames
            covariates(temp) <- tempcov
        }

        attr(temp,'type')    <- 'sort'
        attr(temp,'meanSD')  <- attr(x,'meanSD',exact = TRUE)
        attr(temp,'area')    <- attr(x, 'area',exact = TRUE)
        attr(temp,'SLDF')    <- attr(x, 'SLDF',exact = TRUE)
        attr(temp, 'graph')  <- attr(x, 'graph',exact = TRUE)
        attr(temp,'spacing') <- attr(x, 'spacing',exact = TRUE)
        attr(temp,'boundingbox') <- attr(x,'boundingbox',exact = TRUE)
        class(temp) <- class(x)
        temp
    }
}

print.capthist <- function (x,..., condense = FALSE, sortrows = FALSE)
{
    ## recursive if list of capthist
    if (ms(x)) lapply (x, print.capthist, ..., condense = condense,
        sortrows = sortrows)
    else { # strip attributes, but why bother?
        cat('Session = ', session(x), '\n')
        if (is.null(detector(traps(x)))) {
            print.default(x, ...)
            invisible (x)
        }
        else
        if (all(detector(traps(x)) == 'unmarked')) {
            temp <- apply(x>0, 2:3, sum)
            colnames(temp) <- rownames(traps(x))
            print.default(temp, ...)
            invisible (temp)
        }
        else
        if (all(detector(traps(x)) == 'presence')) {
            ## temp <- apply(apply(x>0, 2:3, sum, drop = FALSE) > 0,2,any, drop = FALSE)*1.0
            ## 2015-05-10 retain occasion-specific data
            temp <- (apply(x>0, 2:3, sum, drop = FALSE) > 0)*1.0
            print.default(temp, ...)
            invisible (temp)
        }
        else
        if (condense & all(detector(traps(x)) %in% c('proximity', 'count',
                                                  'signal','signalnoise'))) {
            temp <- apply(x, 3, function(y) y[apply(abs(y),1,sum)>0,, drop=F])
            trapnames <- rownames(traps(x))
            traps <- trapnames[rep(1:length(temp), sapply(temp, nrow))]
            detections <- as.matrix(abind(temp, along=1))
            temp <- data.frame(AnimalID = rownames(detections), Trap = traps, detections,
                stringsAsFactors = FALSE, row.names=NULL)
            names(temp)[3:ncol(temp)] <- 1:(ncol(temp)-2)
            if (sortrows) {
                lab <- temp$AnimalID
                if (suppressWarnings( all(!is.na(as.numeric(lab)))))
                    lab <- as.numeric(lab)
                temp <- temp[order(lab, temp$Trap ), ]
            }
            rownames(temp) <- 1:nrow(temp)
            print(temp, ...)
            invisible (temp)
        }
        else {
            temp <- array (x, dim=dim(x))
            dimnames(temp) <- dimnames(x)
            if (sortrows) {
                lab <- dimnames(temp)[[1]]
                if (suppressWarnings( all(!is.na(as.numeric(lab)))))
                    lab <- as.numeric(lab)
                if (length(dim(x)) == 3) temp <- temp[order(lab),,]
                else temp <- temp[order(lab),]
            }
            print.default(temp,...)
            invisible(temp)
        }
    }
}
############################################################################################

###############################
## Class : mask
## defines a habitat mask
###############################

subset.mask <- function (x, subset, ...) {

    if (ms(x))
        stop ("subset of multi-session mask not implemented")

    # subset may be numeric index or logical
    temp <- x[subset, , drop = FALSE]
    spacing <- attr(x,'spacing', exact = TRUE)
    attr(temp,'type')        <- 'subset'
    attr(temp,'meanSD')      <- getMeanSD(temp)
    attr(temp,'area')        <- attr(x, 'area', exact = TRUE)
    attr(temp,'SLDF')    <- attr(x, 'SLDF', exact = TRUE)
    attr(temp, 'graph')  <- attr(x, 'graph', exact = TRUE)
    attr(temp,'spacing')     <- spacing
    if (!is.null(covariates(x))) covariates(temp) <- covariates(x)[subset,,drop=F]
    xl <- range(temp$x) + spacing/2 * c(-1,1)
    yl <- range(temp$y) + spacing/2 * c(-1,1)
    attr(temp,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
    if (!is.null(attr(temp,'OK',exact = TRUE)))
        attr(temp,'OK') <- attr(temp,'OK',exact = TRUE)[subset]
    ## 2015-10-18
    ## markingpoints identifies a cutoff: the first mp rows are used (for marking)
    ## this must be reset in a subset operation
    if (!is.null(attr(x,'markingpoints',exact = TRUE)))
        attr(temp,'markingpoints') <- sum(rep(1,nrow(x))[subset] <=
                                          attr(x,'markingpoints',exact = TRUE))
    class(temp) <- class(x)
    temp
}
############################################################################################

###############################################################################

####################################################
## Class : secr
## spatially explicit capture-recapture model fit
####################################################

############################################################################################

trim.secr <- function (object, drop = c('call', 'mask', 'designD', 'designNE', 
    'design','design0'), keep = NULL) {
    trim.default(object, drop = drop, keep = keep)
}
############################################################################################

trim.secrlist <- function (object, drop = c('call', 'mask', 'designD', 'designNE', 
    'design','design0'), keep = NULL) {
    out <- lapply(object, trim, drop = drop, keep = keep)
    class(out) <- class(object)
    out
}
############################################################################################

## 2012-11-14, 2017-10-27, 2017-12-13, 2022-04-17, 2024-02-15
secrlist <- function(..., names = NULL) {
    dots <- match.call(expand.dots = FALSE)$...
    allargs <- list(...)
    if (length(allargs)==1 & inherits(allargs[[1]], 'secrlist')) {
        if (!is.null(names)) {
            if (length(allargs[[1]]) != length(names)) stop (length(allargs[[1]]), " names required", call. = FALSE)
            names(allargs[[1]]) <- names
        }
        return (allargs[[1]])
    }
    else {
        if (is.null(names(allargs))) {
            dots2 <- substitute(list(...))[-1]
            names(allargs) <- sapply(dots2, deparse)
        }

        allargs <- lapply(allargs, function(x) if (inherits(x, c('secr','ipsecr'))) list(x) else x)
        temp <- do.call(c, allargs)
        if (!is.null(names)) {
            if (length(temp) != length(names)) stop (length(temp), " names required", call. = FALSE)
            names(temp) <- names
        }
        else if (is.null(names(temp))) {
            names(temp) <- paste("secr", 1:length(temp), sep="")
        }
        if (!all(sapply(temp, function(x) inherits(x, c('secr','ipsecr')))))
            stop ("objects must be of class 'secr' or 'secrlist'")
        class(temp) <- 'secrlist'
        temp
    }
}

############################################################################################
## 2014-02-19
## extract method for secrlist objects
## retains class
'[.secrlist' <- function(x, i) {
 y <- NextMethod("[")
 class(y) <- class(x)
 y
}

############################################################################################

coef.secr <- function (object, alpha=0.05, ...) {
    beta   <- object$fit$par
    if (!is.null(object$beta.vcv))
        sebeta <- suppressWarnings(sqrt(diag(object$beta.vcv)))
    else sebeta <- rep(NA, length(beta))
    z <- abs(qnorm(1-alpha/2))
    temp <- data.frame(
        row.names = object$betanames,
        beta    = beta,
        SE.beta = sebeta,
        lcl = beta - z*sebeta,
        ucl = beta + z*sebeta
    )
    attr(temp, 'alpha') <- alpha
    temp
}
############################################################################################

coef.secrlist <- function (object, alpha=0.05, ...) {
    lapply(object, coef, alpha = alpha, ...)
}
############################################################################################

detectpar.default <- function(object, ...) {
    stop ("only for secr models")
}
## byclass option 2013-11-09
detectpar.secr <- function(object, ..., byclass = FALSE) {
    extractpar <- function (temp) {
        if (!is.data.frame(temp))   ## assume list
        {
            if (byclass)
                lapply(temp, extractpar)
            else
                extractpar(temp[[1]])
        }
        else {
            if (!is.data.frame(temp) |
                (nrow(temp) > length(object$link)))
                stop ("unexpected input to detectpar()")

            temp <- temp[, 'estimate', drop = F]
            temp <- split(temp[,1], rownames(temp))
            temp <- c(temp, object$fixed)
            pnames <- parnames(object$detectfn)
            if (object$details$param == 2)
                pnames[1] <- 'esa'
            if (object$details$param %in% c(4,5)) {
                Dindex <- match ('D', names(temp))
                sigmakindex <- match ('sigmak', names(temp))
                cindex <- match ('c', names(temp))
                temp[[sigmakindex]] <- temp[[sigmakindex]] / temp[[Dindex]]^0.5 + temp[[cindex]]
                names(temp)[sigmakindex] <- 'sigma'
            }
            if (object$details$param %in% c(3,5)) {
                a0index <- match ('a0', names(temp))
                sigmaindex <- match ('sigma', names(temp))
                lambda0 <- temp[[a0index]] / 2 / pi / temp[[sigmaindex]]^2 * 10000
                temp[[a0index]] <- if (object$detectfn %in% 14:19) lambda0 else 1-exp(-lambda0)
                names(temp)[a0index] <- if (object$detectfn %in% 0:8) 'g0'
                else if (object$detectfn %in% 14:19) 'lambda0'
                else stop ('invalid combination of param %in% c(3,5) and detectfn')
            }
            temp <- temp[pnames]
            if ((object$detectfn > 9) & (object$detectfn <14))
                temp <- c(temp, list(cutval = object$details$cutval))
            temp
        }
    }
    if (!inherits(object,'secr'))
        stop ("requires 'secr' object")
    temppred <- predict (object, ...)
    if (ms(object)) {
        temp <- lapply(temppred, extractpar)
        names(temp) <- session(object$capthist)
        temp
    }
    else {
        extractpar(temppred)
    }
}

############################################################################################

print.secr <- function (x, newdata = NULL, alpha = 0.05, deriv = FALSE, call = TRUE, ...) {

    if (!is.null(x$call) & call) {
        cat ('\n')
        print(x$call)
    }

    if (!is.null(x$version)) {
        cat ('secr ', x$version, ', ', x$starttime, '\n', sep='')
    }
    else {   ## for backward compatibility
        cat (x$fitted,'\n')
    }
    cat ('\n')

    print(summary(traps(x$capthist)), terse=TRUE)

    ## 2017-11-02
    if (!is.null(x$details$newdetector)) {
        cat ('New detector type ', x$details$newdetector, '\n')
    }

    cat ('\n')

    ###################
    ## Data description

    if (ms(x$capthist)) {
        print (summary(x$capthist, terse = TRUE))
        Tu <- sapply(x$capthist, function(y) attr(y,'Tu',exact = TRUE))
        det <- detector(traps(x$capthist)[[1]])
        xyl <- NULL
    }
    else {
        det <- detector(traps(x$capthist))
        n  <- nrow(x$capthist)     # number caught
        if (length(dim(x$capthist))>2)
            ncapt <- sum(abs(x$capthist))
        else
            ncapt <- sum(abs(x$capthist)>0)

        defaultmarkresight <- list(Tu='as.is', Tm='as.is', Tn='ignore')
        markresightcontrol <- replacedefaults(defaultmarkresight, x$details$markresight)

        Tu <- Tu(x$capthist)
        Tm <- Tm(x$capthist)
        Tn <- Tn(x$capthist)
        unresolvedocc <- markocc(traps(x$capthist)) == -1
        unmarkedocc <- markocc(traps(x$capthist)) == 0

        if ('g' %in% x$vars) {
            Groups  <- table(group.factor(x$capthist, x$groups))
            temp <- paste (names(Groups), Groups, collapse=', ', sep='=')
            temp <- paste('(',temp,')', sep='')
        }
        else temp <- ''

        cat ('N animals       : ', n, temp, '\n')
        cat ('N detections    : ', ncapt, '\n')

        if (!(is.null(Tu) | markresightcontrol$Tu=='ignore'))
            cat ('N unmarked sght :  ', sum(unlist(Tu)), ' (c-hat ', round(x$details$chat[1],3),  ')\n', sep='')
        if (!(is.null(Tm) | markresightcontrol$Tm=='ignore'))
            cat ('N nonID sghting :  ', sum(unlist(Tm)), ' (c-hat ', round(x$details$chat[2],3), ')\n', sep='')
        if (!(is.null(Tn) | markresightcontrol$Tn=='ignore'))
            cat ('N other sghting : ', sum(unlist(Tn)), '\n')

        cat ('N occasions     : ', ncol(x$capthist), '\n')

        xyl <- telemetryxy(x$capthist)
        if (!is.null(xyl)) {
            ## zeros <- sum(apply(abs(object)>0,1,sum)==0)
            ntelem <- sapply(xyl, nrow)
            nteldet <- if ((nrow(x$capthist) == 0) | (all(detector(traps(x$capthist))=='telemetry')))
                           0 else
                sum(apply(abs(x$capthist)[,-ncol(x$capthist),,drop=FALSE]>0,1,any) [row.names(x$capthist) %in% names(xyl)])
            ## cat ('Known all-zero  : ', zeros, '\n')
            cat ('Telemetry       : ', length(xyl), 'animals,', nteldet, 'detected\n')
            cat ('Telemetry locns : ', paste(range(ntelem), collapse='-'), 'per animal (mean',
                 round(mean(ntelem),2), 'sd', paste(round(sd(ntelem),2), ')',sep=''), '\n')
        }
    }   # end of single-session

    if (any(det %in% .localstuff$countdetectors)) {
        cat ('Count model     :  ')
        if (x$details$binomN == 0) cat ('Poisson \n')
        else if (x$details$binomN == 1) cat ('Binomial, size from usage\n')
        else if (x$details$binomN < 0) cat ('Negative binomial k = ', abs(x$details$binomN), '\n')
        else if (x$details$binomN > 1) cat('Binomial', x$details$binomN, '\n')
    }

    if (!ms(x$capthist)) {
        if (length(maskarea(x$mask))==0)
            cat ('Mask length     : ', masklength(x$mask), 'km \n')
        else
            cat ('Mask area       : ', maskarea(x$mask), 'ha \n')
    }

    ####################
    ## Model description

    ## 2015-03-31
    Npar <- nparameters(x)   ## see utility.R
    AICval <- 2*(x$fit$value + Npar)
    n <- ifelse (ms(x$capthist), sum(sapply(x$capthist, nrow)), nrow(x$capthist))
    AICcval <- ifelse ((n - Npar - 1) > 0,
        2*(x$fit$value + Npar) + 2 * Npar * (Npar+1) / (n - Npar - 1),
        NA)
    cat ('\n')
    cat ('Model           : ', model.string(x$model, x$details$userDfn), '\n')

    ## 2013-06-08
    if (!is.null(x$hcov))
        cat ('Mixture (hcov)  : ', x$hcov, '\n')

    ## 2015-01-16
    if (!is.null(x$details$userdist)) {
        if (is.matrix(x$details$userdist))
            cat ('User distances  :  static (matrix)\n')
        if (is.function(x$details$userdist))
            cat ('User distances  :  dynamic (function)\n')
    }

    cat ('Fixed (real)    : ', fixed.string(x$fixed), '\n')
    cat ('Detection fn    : ', detectionfunctionname(x$detectfn))
    cat ('\n')
    if (!x$CL)
        cat ('Distribution    : ', x$details$distribution, '\n')

    cat ('N parameters    : ', Npar, '\n')
    cat ('Log likelihood  : ', -x$fit$value, '\n')
    cat ('AIC             : ', AICval, '\n')
    cat ('AICc            : ', AICcval, '\n')

    cat ('\n')
    cat ('Beta parameters (coefficients)', '\n')

    print(coef(x), ...)

    if (!is.null(x$fit$hessian)) {
      cat ('\n')
      cat ('Variance-covariance matrix of beta parameters', '\n')
      print (x$beta.vcv, ...)
    }

    # scale newdata covariates... NOT FINISHED 10 05 08
    ## 2018-10-14 NOTE maybe this was never finished!
    ## remove at next maintenance?
    meanSD <- attr(x$mask,'meanSD',exact = TRUE)
    if (!is.null(newdata)) {
         for (i in 1:length(newdata)) {
           ind <- match (names(newdata[i]),names(meanSD))
           if (ind>0 & !is.na(meanSD[1,ind]))
             newdata[[i]] <- (newdata[[i]] - meanSD[1,ind]) / meanSD[2,ind]
         }
     }

    cat ('\n')
    cat ('Fitted (real) parameters evaluated at base levels of covariates', '\n')

    ## 2018-10-14 NOTE realpar is redundant, never used
    ## to be removed at next maintenance
    if (!is.null(x$realpar))
        print( x$realpar )
    else {
        temp <- predict (x, newdata, type = "response", alpha = alpha)
        nd <- length(temp)
        if (is.data.frame(temp)) print(temp, ...)
        else for (new in 1:nd) {
                cat('\n', names(temp)[new],'\n')
                print(temp[[new]], ...)
            }
    }

    #################################
    # Derived parameters
    #################################
    if (deriv) {

        cat ('\n')
        cat ('Derived parameters', '\n')

        temp <- derived(x, alpha=alpha, se.esa = TRUE)
        nd <- length(temp)
        if (is.data.frame(temp)) print(temp, ...)
        else for (new in 1:nd) {
                cat('\n',names(temp)[new],'\n')
                print(temp[[new]], ...)
            }

    }
    cat ('\n')
}
############################################################################################

vcov.secr <- function (object, realnames = NULL, newdata = NULL, byrow = FALSE, ...) {
## return either the beta-parameter variance-covariance matrix
## or vcv each real parameters between points given by newdata (byrow = TRUE)
## or vcv for real parameters at points given by newdata (byrow = TRUE)

    if (is.null(dimnames(object$beta.vcv)))
        dimnames(object$beta.vcv) <- list(object$betanames, object$betanames)

    if (is.null(realnames))
        ## average beta parameters
        return( object$beta.vcv )
    else {
        ## average real parameters
        ## vcv among multiple rows

        if (byrow) {
            ## need delta-method variance of reals given object$beta.vcv & newdata
            if (is.null(newdata)) {
              # newdata <- secr.make.newdata (object)
              newdata <- makeNewData (object)
            }
            nreal <- length(realnames)
            nbeta <- length(object$fit$par)

            rowi <- function (newdatai) {
                reali <- function (beta, rn) {
                    ## real from all beta pars eval at newdata[i,]
                    par.rn <- object$parindx[[rn]]
                    mat <- general.model.matrix(
                        object$model[[rn]], 
                        data = newdatai,
                        gamsmth = object$smoothsetup[[rn]],
                        contrasts = object$details$contrasts)
                    lp <- mat %*% matrix(beta[par.rn], ncol = 1)
                    untransform (lp, object$link[[rn]])
                }
                grad <- matrix(nrow = nreal, ncol = nbeta)
                dimnames(grad) <- list(realnames, object$betanames)
                for (rn in realnames)
                    grad[rn,] <- fdHess (pars = object$fit$par, fun = reali, rn = rn)$gradient
                vcv <- grad %*% object$beta.vcv %*% t(grad)
                vcv
            }

            vcvlist <- list(nrow(newdata))
            for (i in 1:nrow(newdata)) vcvlist[[i]] <- rowi(newdata[i,])
            if (length(vcvlist) == 1) vcvlist <- vcvlist[[1]]
            return(vcvlist)
        }
        else {
            newdata <- as.data.frame(newdata)
            rownames <- apply(newdata, 1, function(x) paste(names(newdata), '=', x, sep='',
                                                            collapse=','))
            vcvlist <- list()
            for (rn in realnames) {
                ## temporary fix 2015-09-30
                if (rn == 'pmix')
                    stop("vcov does not work at present when realname == 'pmix'")
                par.rn <- object$parindx[[rn]]
                mat <- general.model.matrix(
                    object$model[[rn]], 
                    data = newdata,
                    gamsmth = object$smoothsetup[[rn]],
                    contrasts = object$details$contrasts)
                lp <- mat %*% matrix(object$fit$par[par.rn], ncol = 1)
                real <- untransform (lp, object$link[[rn]])
                real <- as.vector(real)
                ## from Jeff Laake's 'compute.real' in RMark...
                deriv.real <- switch(object$link[[rn]],
                    logit = mat * real * (1-real),
                    log = mat * real,
                    identity = mat,
                    sin = mat * cos(asin(2*real-1))/2)
                vcvlist[[rn]] <- deriv.real %*% object$beta.vcv[par.rn, par.rn] %*% t(deriv.real)
                dimnames(vcvlist[[rn]]) <- list(rownames, rownames)
            }
            names (vcvlist) <- realnames
            return (vcvlist)
        }
        ## DIFFERENT VARIANCE TO secr.lpredictor for sigma because there use se.Xuntransfom
    }
}

############################################################################################

## 2017-03-15
## consider case that one of ... is a naked secr object, not a list
## cf secrlist()?
c.secrlist <- function(..., recursive = FALSE) {
    slist <- as.list(...)
    result <- NextMethod('c', slist, recursive = FALSE)
    class(result) <- 'secrlist'
    result
}
############################################################################################

as.mask <- function (x) {
    if (!inherits(x, "traps"))
        stop ("as.mask is defined only for traps objects")
    if (ms(x)) {
        out <- lapply(x, as.mask)
        class(out) <- c('mask', 'list')
        out
    }
    else {
        class(x) <- c('mask', 'data.frame')
        xl <- range(x$x) + spacing(x)/2 * c(-1,1)
        yl <- range(x$y) + spacing(x)/2 * c(-1,1)
        attr(x,'type')        <- "coerced"
        attr(x,'meanSD')      <- getMeanSD (x)
        attr(x,'area')        <- spacing(x)^2 * 0.0001
        attr(x,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
        x
    }
}
############################################################################################

as.popn <- function (x) {
    if (!inherits(x, "ppp")) {        
        stop ("as.popn is defined only for ppp objects")
    }       
    if (requireNamespace("spatstat.geom", quietly = TRUE)) {
        pop <- spatstat.geom::coords(x)
        class(pop) <- c('popn', 'data.frame')
        xl <- x$window$xrange
        yl <- x$window$yrange
        attr(pop, 'boundingbox') <- expand.grid (x=xl, y=yl)[c(1,3,4,2),]
        attr(pop, "Lambda") <- attr(x, "Lambda")   # if present
        pop
    }
    else {
        stop ("as.popn requires package spatstat")
    }
}
############################################################################################

# derived.default <- function (object, ...)    {
#     if (!inherits(object, 'secr'))
#         stop ("no derived method for this class")
#     }
# }

## see derivedMS.R for derived.secr
############################################################################################

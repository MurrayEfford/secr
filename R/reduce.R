###########################################################################################
## package 'secr'
## reduce.R
## 2010-12-01 check for overlapping columns
## 2011-02-08 rewritten to include polygonX and transectX detectors, and simplified
## 2011-03-18 output to unmarked
## 2011-03-21 'by' argument
## 2012-12-17 non-binary usage
## 2012-12-21 re-write reduce.capthist to include spatial lumping
## 2012-12-21 amalgamated 'reduce' methods
## 2013-11-20 telemetry allowed
## 2016-10-06 3D traps
## 2017-03-12 revamp
## 2017-03-25 debug for reduce from transect to count
## 2017-10-20 introducing 'capped' detector
## 2017-11-14 purged unused code for seltrap, selused
## 2018-01-22 made safe for nonspatial data
## 2018-01-25 multisession intervals passed through
## 2018-05-12 faster handling of 'alive'; bug fixed 2018-06-23 3.1.7
## 2019-01-16 bug: clash between usage and markocc in reduce.traps FIXED
## 2021-05-19 sortorder for animalID etc.
## 2022-01-05 fixed alive() problem with one animal
############################################################################################

# From (row) To (column)
#----------------------------------------------------------------------------------------------------
#              single  multi  proximity  count  polygonX transectX  signal  polygon transect capped
# single	&#	&	*	   *      NA       NA         NA      NA      NA       *
# multi		&#	&	*	   *      NA       NA         NA      NA      NA       #
# proximity	&#	&	*	   *      NA       NA         NA      NA      NA       #
# count		&#@	&@	@	   *      NA       NA         NA      NA      NA       #
# polygonX	&#	&	*	   *      &$       NA         NA      NA      NA      NA
# transectX	&#	&	*	   *      NA       &$         NA      NA      NA      NA
# signal        &#	&	*	   @      @        @          $       NA      NA      NA
# polygon     	&#@~	&@~	@~	   *~     @        NA         NA       *      NA      NA
# transect     	&#@~	&@~	@~	   *~     NA       @          NA      NA       *      NA
# capped        &#	&	*	   *      NA       NA         NA      NA      NA       *

#----------------------------------------------------------------------------------------------------
#  * no loss of data
#  # must choose among animals (more than one animal)
#  & must choose among traps (animal caught more than once)
#  @ reduce to binary
#  $ apply rule for combining signals or locations (first, last, random, min, max, mean)
#  ~ form new point detectors from mean of vertices (assumes symmetry)
#  NA not feasible
############################################################################################


reduce     <- function (object, ...) UseMethod("reduce")

reduce.default <- function (object, columns, ...) {
  object <- as.matrix(object)
  if (any(is.na(object)))
      warning ("NAs in input converted to zero")
  firsttrap <- function (y) y[abs(y)>0][1]    # first non-zero
  fnmulti   <- function (occ) apply (object[,occ,drop=F], 1, firsttrap)
  nrow <- nrow(object)
  nnew <- length(columns)
  temp <- sapply (columns, fnmulti)
  temp[is.na(temp)] <- 0
  temp
}

collapseocc <- function (newoccasions, inimatrix) {
    fnused <- function (occ, fn) {
        if (length(occ)>0) {
            temp <- inimatrix[, occ, drop = FALSE]
            apply (temp, 1, fn)
        }
        else NULL
    }
    tempmatrix <- unlist(sapply (newoccasions, fnused, sum))
    matrix(tempmatrix, ncol=length(newoccasions))
}
#############################################################################################################

reduce.traps <- function (object, newtraps = NULL, newoccasions = NULL, span = NULL, rename = FALSE,
    newxy = c("mean", "first"), ...) {
    if (ms(object)) {
        out <- lapply(object, reduce, newtraps = newtraps, newoccasions = newoccasions, span = span, rename = rename,
            newxy = newxy, ...)
        class(out) <- class(object)
        out
    }
    else {
        if (!inherits(object, 'traps'))
            stop ("requires traps object")
        newxy <- match.arg(newxy)
        splitfactor <- 1:ndetector(object)    # default to status quo

        #############################################################
        if (is.null(span) & is.null(newtraps)) {
            newcoord <- object
            newtrapnames <- levels(polyID(object))
            attr(newcoord, 'newtrap') <- newtrapnames
        }
        else {
            if (!all(detector(object) %in% .localstuff$pointdetectors))
                stop ("reduce.traps with 'span' or 'newtraps' is only for point detectors")
            ## allow for vector input, or distance threshold
            if (!is.null(span))
                newtraps <- cutree (hclust(dist(object), ...), h = span)

            if (!is.list(newtraps)) {
                if (is.null(names(newtraps)))
                    names(newtraps) <- 1:length(newtraps)
                newtraps <- split(names(newtraps), newtraps)
            }

            if (any(duplicated(unlist(newtraps))))
                stop("traps should not appear in more than one group")

            if (any(sapply(newtraps, is.character))) {
                ## convert to numeric trap indices
                newtraps <- lapply(newtraps, function(x) match(x, rownames(object)))
            }

            nnew <- length(newtraps)
            if (rename)
                newtrapnames <- 1:nnew
            else if (any(sapply(newtraps, is.character))) {
                newtrapnames <- sapply(newtraps, function(x) paste(x, collapse='+'))
            }
            else {
                namefn <- function(x) paste(rownames(object)[x], collapse='+')
                newtrapnames <- sapply(newtraps, namefn)
            }

            g <- rep(1:nnew, sapply(newtraps,length))
            splitfactor[] <- 0
            splitfactor[unlist(newtraps)] <- g             ## levels are indices of groups 1, 2,...
            splitfactor <- factor(splitfactor)
            grouped <- split.data.frame(object, splitfactor) ## otherwise calls split.traps
            ## or could use for usage and cov below and save repeat splitting...
            if (any (splitfactor == 0)) grouped <- grouped[-1]  ## drop unwanted traps
            names(grouped) <- newtrapnames
            if (newxy == 'mean') {
                grouped <- lapply(grouped, function(df) apply(df,2,mean))
            }
            else {
                grouped <- lapply(grouped, '[', 1,)   ## 2021-04-14
            }
            newcoord <- do.call(rbind, grouped)
            newcoord <- as.data.frame(newcoord)
            class (newcoord)   <- c('traps', 'data.frame')

            sp <- spacing(newcoord, recalculate = TRUE)
            if (!is.null(sp)) spacing(newcoord) <- sp
            temp <- as.numeric(levels(splitfactor)[splitfactor])
            temp[temp==0] <- NA
            attr(newcoord, 'newtrap') <- temp

        }

        #############################################################
        ## what if occasions collapsed?

        newdetector <- detector(object)
        newmarkocc <- markocc(object)
        newusage <- usage(object)
        
        if (!is.null(newoccasions)) {
            ## return one detector type per new occasion
            validdet <- function(occ) {
                if (length(detector(object))==1)
                    return(detector(object))
                else {
                    tempdet <- detector(object)[occ]
                    if (length(unique(tempdet))>1)
                        stop("cannot combine occasions with differing detector type")
                    tempdet[1]
                }
            }
            newdetector <- sapply(newoccasions, validdet)

            ## return one markocc per new occasion
            if (!is.null(markocc(object))) {
                validmocc <- function(occ) {
                    tempmocc <- markocc(object)[occ]
                    if (length(unique(tempmocc))>1)
                        stop("cannot combine occasions with differing markocc code")
                    tempmocc[1]
                }
                newmarkocc <- sapply(newoccasions, validmocc)
            }
        }
        #############################################################

        if (!is.null(usage(object))) {
            # if (is.null(splitfactor)) {
            #     usagelist <- list(usage(object))
            #     newusagenames <- rownames(usage(object))
            # }
            # else {
            usagelist <- split(as.data.frame(usage(object)), splitfactor)
            newusagenames <- newtrapnames
            if (any (splitfactor == 0)) usagelist <- usagelist[-1]  ## drop unwanted traps
            #            }
            newusage <- lapply(usagelist, function(x) apply(x,2,sum))
            newusage <- do.call(rbind, newusage)
            newusage <- as.matrix(newusage)
            dimnames(newusage) <- list(newusagenames, 1:ncol(newusage))
            
            daily <- as.data.frame(newusage)
            names(daily) <- paste('occ', names(daily), sep='')
            
            if (!is.null(newoccasions)) {
                newusage <- collapseocc (newoccasions, newusage)
            }
            
        }
        
        # return to consider covariates when traps combined,
        # using 'daily' from preceding

        if (!is.null(span) || !is.null(newtraps)) {

            if (!is.null(covariates(object))) {
                covlist <- split(covariates(object), splitfactor)
                if (any (splitfactor == 0)) covlist <- covlist[-1]  ## drop unwanted traps
                varying <- sapply(covlist, function(x) apply(x,2,function(y) length(unique(y))>1))
                if (any(varying))
                    warning("covariates vary within groups; using only first")
                temp <- lapply(covlist, head, 1)
                temp <- do.call(rbind, temp)
                covariates(newcoord) <- temp
                if (!is.null(usage(object)))
                    covariates(newcoord) <- cbind(covariates(newcoord), daily)
                else
                    covariates(newcoord)$combined <- sapply(newtraps, length)
            }
            else {  ## new covariate dataframe
                if (!is.null(usage(object)))
                    covariates(newcoord) <- daily
                else
                    covariates(newcoord) <- data.frame(combined = sapply(newtraps, length))
            }
        }
        #############################################################

        markocc(newcoord)  <- NULL   ## to avoid clash in next line; replaced below
        usage(newcoord)    <- newusage
        detector(newcoord) <- newdetector
        markocc(newcoord)  <- newmarkocc
        newcoord
    }
}
#############################################################################################################

poly2point <- function (object, detector = 'count') {
    if (!all(detector(object) %in% c('polygon','polygonX')))
        stop ("requires 'polygon' input")
    if (any(detector %in% .localstuff$polydetectors))
        stop ("requires non-polygon, non-transect output")
    temp <- split(object, polyID(object))
    temp <- lapply(temp, function(df) apply(df,2,mean))
    temp1 <- t(abind(temp, along=2))
    dimnames(temp1) <- list(levels(polyID(object)), c('x','y'))
    temp <- data.frame(temp1, row.names=NULL)
    class (temp)   <- c('traps', 'data.frame')
    detector(temp) <- detector
    usage(temp)    <- usage(object)
    covariates(temp) <- covariates(object)
    attr(temp,'spacex') <- 100 * (searcharea(object)/nrow(temp))^0.5
    attr(temp,'spacey') <- attr(temp,'spacex')
    temp
}
#-----------------------------------------------------------------------------

transect2point <- function (object, detector = 'count') {
    if (!all(detector(object) %in% c('transect','transectX')))
        stop ("requires 'transect' input")
    if (any(detector %in% c('transect', 'transectX')))
        stop ("requires non-transect output")
    temp <- split(object, transectID(object))
    temp <- lapply(temp, function(df) apply(df,2,mean))
    temp1 <- t(abind(temp, along=2))
    dimnames(temp1) <- list(levels(transectID(object)), c('x','y'))
    temp <- data.frame(temp1, row.names=NULL)
    class (temp)   <- c('traps', 'data.frame')
    detector(temp) <- detector
    usage(temp)    <- usage(object)
    covariates(temp) <- covariates(object)
    attr(temp,'spacex') <- mean(transectlength(object))/2   ## arbitrary
    attr(temp,'spacey') <- attr(temp,'spacex')
    temp
}
#-----------------------------------------------------------------------------

## function to make list in which each component is a
## subset of occasions (for use in reduce.capthist)

splitby <- function (x, by) {
    if ((length(x) == 1) & (x[1] > 1))
        x <- 1:x
    if (by < 1)
        stop ("invalid 'by' argument")
    index <- 1:length(x)
    gp <- trunc((index-1)/by) + 1
    split (index, gp)
}
#----------------------------------------------------------------------------------------------

reduce.capthist <- function (object, newtraps = NULL, span = NULL,
    rename = FALSE, newoccasions = NULL, by = 1, outputdetector = NULL,
    select = c("last", "first", "random"), dropunused = TRUE,
    verify = TRUE, sessions = NULL, ...) {

    # newoccasions - list, each component gives occasions to include in new capthist
    # newtraps     - list, each component gives traps to include in new capthist

    #----------------------------------------------------------------------------
    collapse <- function (df) {
        ## reduce data frame to a single row
        if (nrow(df)>1) {
            df$alive <- rep(all(df$alive), nrow(df))
            index <- switch (select, first = 1, last = nrow(df),
                random = sample.int (nrow(df),1) )
            df <- df[index,,drop=FALSE]
        }
        df
    }
    #----------------------------------------------------------------------------

    # main line
    
    if (ms(object)) {
        if (is.null(sessions)) sessions <- 1:length(object)
        temp <- lapply (object[sessions], reduce,
            newoccasions = newoccasions,
            by = by,
            newtraps = newtraps,
            span = span,
            outputdetector = outputdetector,
            select = select,
            rename = rename,
            dropunused = dropunused,
            verify = verify,
            ...)
        class(temp) <- c('capthist', 'list')
        if (length(temp) == 1) temp <- temp[[1]]
        interv <- intervals(object)
        if (!is.null(interv)) {   ## 2018-01-25
            cumi <- cumsum(interv)
            newinterv <- diff(c(0,cumi)[sessions])
            intervals(temp) <- newinterv 
        }
        sessionlabels(temp) <- sessionlabels(object)[sessions]
        return(temp)
    }
    else {
        select <- match.arg(select)
        polygons <- c('polygon','polygonX')
        transects <- c('transect','transectX')
        nrw <- nrow(object)
        if (is.null(newoccasions)) {
            if (tolower(by) == 'all') by <- ncol(object)
            newoccasions <- splitby (1:ncol(object), by)
            if ((ncol(object) %% by) > 0)
                warning ("number of occasions is not a multiple of 'by'")
        }
        if (!is.null(traps(object))) {
            ntrap <- ndetector(traps(object))  ## npoly if 'polygon' or 'transect'
            inputdetector <- detector(traps(object))
        }
        else {
            ntrap <- 1
            inputdetector <- 'nonspatial'
            outputdetector <- 'nonspatial'
        }
        if (length(unique(inputdetector))>1)
            stop("reduce.capthist does not yet handle mixed input detector types")
        inputdetector <- inputdetector[1]
        if (is.null(outputdetector))
            outputdetector <- inputdetector
        else
            if (length(unique(unlist(outputdetector)))>1) {
                warning("reduce.capthist does not yet handle mixed output detector types")
                outputdetector <- unlist(outputdetector[1])
            }
        if (!(outputdetector %in% c(.localstuff$validdetectors, 'nonspatial')))
            stop ("'outputdetector' should be one of ",
                  paste(sapply(.localstuff$validdetectors, dQuote),collapse=','))
        if ((!(inputdetector %in% c('signal','signalnoise'))) & (outputdetector == 'signal'))
                stop ("cannot convert non-signal data to signal data")
        if ((!(inputdetector %in% c('signalnoise'))) & (outputdetector == 'signalnoise'))
                stop ("cannot convert non-signalnoise data to signalnoise data")
        if ((!(inputdetector %in% polygons)) & (outputdetector %in% polygons))
                stop ("cannot convert non-polygon data to 'polygon' data")
        if ((!(inputdetector %in% transects)) & (outputdetector %in% transects))
                stop ("cannot convert non-transect data to 'transect' data")

        ######################################
        ## add usage to count pooled occasions
        ## 2014-10-11
        if (outputdetector %in% .localstuff$countdetectors) {
            if (is.null(usage(traps(object)))) {
                usage(traps(object)) <- matrix(1, nrow = ndetector(traps(object)),
                                               ncol = ncol(object))
            }
        }

        ######################################
        ## check newoccasions
        for (i in length(newoccasions):1) {
            occ <- newoccasions[[i]]
            occ <- occ[occ %in% (1:ncol(object))]  ## discard nonexistent occ
            if (length(occ)==0)
                newoccasions[[i]] <- NULL
            else
                newoccasions[[i]] <- occ
        }
        cumocc <- numeric(0)
        for (i in length(newoccasions):1) {
            if (any (newoccasions[[i]] %in% cumocc))
                warning ("new occasions overlap")
            cumocc <- c(cumocc, newoccasions[[i]])
        }
        nnew <- length(newoccasions)
        newcols <- rep(1:nnew, sapply(newoccasions,length))
        newcols <- factor(newcols)
        ####################################
        ## check and build newtraps
        if (outputdetector != 'nonspatial') {
            trps <- traps(object)
            reducetraps <- !is.null(newtraps) | !is.null(span) | (length(newoccasions) != ncol(object))
            if (reducetraps) {
                trps <- reduce(trps, newtraps = newtraps, newoccasions = newoccasions,
                               span = span, rename = rename, ...)
                newtrapID <- attr(trps, 'newtrap')
                ntrap <- ndetector(trps)
            }
        }
        else {
            trps <- NULL
            reducetraps <- FALSE
        }
        ####################################
        ## build dataframe of observations
        ## 2021-05-19 using ksn as safe universal order regardless of polygon/signal/non
        ## 2022-01-05 fixed bug in alive() when only one animal
        df <- data.frame(
            trap = trap(object, names = FALSE, sortorder = 'ksn'),
            occ = occasion(object, sortorder = 'ksn'),
            ID = animalID(object, names = FALSE, sortorder = 'ksn'),
            alive = alive(object, sortorder = 'ksn'))
        if (reducetraps) {
            df$trap <- newtrapID[df$trap]
        }
        if (any(outputdetector %in% c(polygons, transects))) {
            df$x <- xy(object)[,1]
            df$y <- xy(object)[,2]
        }
        if (!is.null(attr(object,'signalframe'))) {
            df <- cbind(df, attr(object,'signalframe'))
        }
        if (nrow(df) > 0) {
          df$newocc <- newcols[match(df$occ, unlist(newoccasions))]
          if (dropunused) {
            df$newocc <- factor(df$newocc)
            nnew <- length(levels(df$newocc))
          }
          df <- df[!is.na(df$newocc),]                   ## drop null obs
          df$newID <- factor(df$ID)                      ## assign newID
          if (outputdetector %in% .localstuff$exclusivedetectors) {
            ID.occ <- interaction(df$ID, df$newocc)
            dflist <- split(df, ID.occ)
            dflist <- lapply(dflist, collapse)
            df <- do.call(rbind, dflist)
          }
          
          if (outputdetector %in% c('single', 'capped')) {
            occ.trap <- interaction(df$newocc,df$trap)
            dflist <- split(df, occ.trap)
            dflist <- lapply(dflist, collapse)
            df <- do.call(rbind, dflist)
          }
          df$newID <- factor(df$ID)                     ## re-assign newID
          
          ####################################
          ## build new object
          validrows <- (1:nrow(object)) %in% df$ID   ## or newID??? 2012-12-12
          df$trap <- factor(df$trap, levels = 1:ntrap)  ## bug fixed 2017-03-27
          ## drop any records with missing data - just to be sure
          df <- df[!apply(df,1,function(x) any(is.na(x))),, drop = FALSE]
          tempnew <- table(df$newID, df$newocc, df$trap)
          if (! (outputdetector %in% .localstuff$countdetectors)
            & (length(tempnew)>0)) {
            ## convert 'proximity' and 'signal' to binary
            tempnew[tempnew>0] <- 1
          }
          
          # alivesign <- tapply(df$alive, list(df$newID,df$newocc,df$trap),all)
          # alivesign[is.na(alivesign)] <- TRUE
          # alivesign <- alivesign * 2 - 1
          # tempnew <- tempnew * alivesign
          
          ## 2018-05-12 faster...
          ## 2018-06-23 Bug in 3.1.6 because newocc selects non-existent columns
          ## 2018-06-23 fix by enclosing df$newocc in as.character()
          
          ## 2019-11-29 bug fix for count
          if (outputdetector != 'count') {
            i <- cbind(as.character(df$newID), as.character(df$newocc), as.character(df$trap))
            tempnew[i] <- tempnew[i] * (df$alive * 2 - 1)
          }
        }
        else { # if (nrow(tempnew) == 0) {
          # ignores dropunused
          tempnew <- array(dim=c(0, nnew, ntrap))
          validrows <- 0
          df$newocc <- numeric(0)
        }
        ################################
        ## general attributes
        class(tempnew) <- 'capthist'
        session(tempnew) <- session(object)
        attr(tempnew, 'n.mash') <- attr(object, 'n.mash')
        attr(tempnew, 'centres') <- attr(object, 'centres')
        attr(tempnew, 'popn') <- attr(object, 'popn')    # 2023-09-16
        if (outputdetector == 'nonspatial') {
            # NULL   2021-03-31
        }
        else if ((inputdetector %in% polygons) & !(outputdetector %in% polygons))
            traps(tempnew) <- poly2point(trps)
        else
            if ((inputdetector %in% transects) & !(outputdetector %in% transects))
                traps(tempnew) <- transect2point(trps)
        else 
            traps(tempnew) <- trps
        if (outputdetector != 'nonspatial') 
            detector(traps(tempnew)) <- outputdetector

        ################################
        ## covariates and ancillary data

        if (!is.null(covariates(object)))
             covariates(tempnew) <- covariates(object)[validrows,,drop=F]

        detectorder <- order(df$trap, df$newocc,df$ID)  ## CHECK!

        ################################
        ## polygon & transect xy

        if (outputdetector %in% c(polygons, transects))
            xy(tempnew) <- df[detectorder,c('x','y'),drop=FALSE]

        ################################
        ## telemetry xy
        ## NOT FINISHED - NEEDS DISTINCT X,Y COLUMNS IN DF
        if (any(outputdetector == 'telemetry') & (nrow(tempnew>0))) {
            stop("reduce telemetry not ready")
            idlevels <- rowname <- rownames(object)[(1:length(validrows))[validrows]]
            id <- factor(df$ID, levels = idlevels)
            telemetryxy(tempnew) <- split(df[detectorder,c('x','y'),drop=FALSE], id)
        }

        ################################
        ## signalframe & cutval

        if (outputdetector %in% c('signal','signalnoise')) {
            sigcolumns <- names(attr(object,'signalframe'))
            attr(tempnew,'signalframe') <- df[detectorder, sigcolumns, drop=FALSE]
            attr(tempnew, 'cutval') <- attr(object, 'cutval')
        }

        ##################################
        ## sightings
        reduceT <- function (TmTu) {
            if (length(dim(TmTu))==2) {
                ## assume detector x occasion
                collapseocc (newoccasions, TmTu)
                }
            else {
                if (!is.null(TmTu) & ((length(newoccasions)>1) | (length(newoccasions[[1]]) != ncol(object))))
                    warning ("cannot selectively collapse summed sightings, including all")
                TmTu
            }
        }
        Tu(tempnew) <- reduceT(Tu(object))
        Tm(tempnew) <- reduceT(Tm(object))

        ##################################
        ## optionally drop unused detectors
        if (outputdetector != 'nonspatial') {
            if (nrow(tempnew) > 0)
                dimnames(tempnew)[[1]] <- 1:nrow(tempnew)  ## temporary, for animalID in subset
            if (dropunused & !is.null(usage(traps(tempnew)))) {
                OK <- apply(usage(traps(tempnew)), 1, sum) > 0
                tempnew <- subset(tempnew, traps = OK)
            }
            tempnew[is.na(tempnew)] <- 0
        }

        ################################
        interv <- intervals(object) ## 2018-01-25
        if (!is.null(interv)) {
            if (!is.null(newoccasions)) {  ## 2018-05-10
                ## using last-first
                # inti <- sapply(newoccasions, tail, 1)
                # intervals(tempnew) <- interv[inti[-length(inti)]]
                ## using first-first
                cumint <- c(0,cumsum(interv))
                int1 <- sapply(newoccasions, '[', 1)
                intervals(tempnew) <- diff(cumint[int1])
            }
            #warning ("intervals not carried forward")   
            else intervals(tempnew) <- interv
        }
        ################################
        slabels <- sessionlabels(object)  ## 2018-05-10
        if (!is.null(slabels) & !is.null(interv)) {
            sessionlabels(tempnew) <- slabels[unique(primarysessions(interv))]
        }
        
        ################################
        ## dimnames
        if (nrow(tempnew) > 0) {
            indices <- (1:length(validrows))[validrows]
            rowname <- rownames(object)[indices]
        }
        else
            rowname <- NULL

        if (nnew>0) {   ## 2020-11-18 for robustness to zero rows
            dimnames(tempnew) <- list(rowname,1:nnew,NULL)   # renew numbering
        }

        if (verify) verify(tempnew, report=1)
        tempnew
    }
}
############################################################################################

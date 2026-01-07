###############################################################################
## package 'secr'
## rbind.capthist.R
## 2015-01-11 tweak rownames in rbind.capthist (default component names)
## 2015-01-23 replace long object names (>100ch)
## 2017-07-26 revise to make rbind.capthist S3 method
## 2018-05-11 MS.capthist handles intervals, sessionlabels
###############################################################################

flatten <- function(x) {
## x is expected to be a list whose components are non-list and list objects
## create a new list from a concatenation of the non-list objects and the first-order
## components of the lists
    if (!is.list(x))
        stop ("can only flatten a list")
    temp <- lapply(x,
                  function(y)
                      if (is.list(y))
                          return(y)
                      else {
                          return(list(y))
                      }
                 )
    unlist(temp, recursive = FALSE)
}
###############################################################################

MS.capthist <- function (...) {
    # make a list of capthist objects, one for each session
    # modified 7/6/2010 for more general input:
    #    multiple single-session objects (as in old version)
    #    list of single-session objects
    #    combination of existing MS objects
    #    combination of existing MS objects and single-session objects?
    # modified 2/7/2010 so always named
    # modified 12/9/2011 so creates names as needed

    dots <- match.call(expand.dots = FALSE)$...
    ##    oldsess <- unlist(sapply(list(...), session))
    MS <- flatten(list(...))
    
    interv <- lapply(list(...), intervals)
    if (length(interv)>1) {
        ## insert bridging interval
        for (i in 1:(length(interv)-1)) 
            interv[[i]] <- c(interv[[i]], NA)
    }
    class(MS) <- c('capthist', 'list')
    intervals(MS) <- unlist(interv)
    sessionlabels(MS) <- unlist(lapply(list(...), sessionlabels))
    if (is.null(names(MS))) {
        names(MS) <- rep("",length(MS))
        dots2 <- substitute(list(...))[-1]    ## 2017-11-13
        if (length(MS) == length(dots2))
            names(MS) <- sapply(dots2, deparse)
    }
    if (any(duplicated(names(MS))) | any(names(MS)=="")) {
        warning ("session names replaced to avoid duplication or omission")
        names(MS) <- 1:length(MS)
    }
    session(MS) <- names(MS)
    MS
}
###############################################################################

rbind.capthist <- function (..., renumber = TRUE, pool = NULL, verify = TRUE)
    
{
    # dots <- match.call(expand.dots = FALSE)$...
    allargs <- list(...)

    ##############################################################
    # inputnames <- lapply(dots, as.character)
    # if (any(is.na(inputnames) | duplicated(inputnames) | (nchar(inputnames)>100)))
    #     inputnames <- as.character(1:length(allargs))
    # names(allargs) <- inputnames
    ##############################################################
    if (length(allargs)==1) 
        object <- allargs[[1]]
    else {
        object <- allargs
        if (is.null(names(object)) | any(names(object)==""))
            names(object) <- 1:length(object)
    }

    ## Catch singleton - added 2011-09-12
    if ((length(allargs) == 1) & !ms(object) )
        return(object)       ## unchanged!
    if ((length(allargs) == 1) & ms(object) & (length(object) == 1) )
        return(object[[1]])  ## unchanged!
    
    ## Case 1 DEPRECATED 2011-09-12
    ## several lists or a combination of list & elementary capthist objects
    ## concatenate lists, including elementary objects (ignore 'pool')
    ## objects may differ in traps etc.
    
    if ((length(allargs)>1) & any(sapply(allargs, is.list)))
        stop ("invalid input to rbind.capthist; ",
              "use MS.capthist to concatenate sessions")
    
    ## Case 2
    ## a single MS capthist (i.e. a list)
    ## rbind components as identified in 'pool'
    ## recursive call for each component of 'pool'
    
    if((length(allargs)==1) & (is.list(object)) & !is.null(pool)) {
        
        if (!is.list(pool)) {
            pool <- list(combined=1:length(object))
            warning ("list not specified, pooling all components")
        }
        else if (any (sapply(unlist(pool),
                             function(x) length(object[[x]])==0)))
            ## prempted by 'subscript out of bounds'
            stop ("invalid pooling indices")
        
        getpooled <- function (x) {
            temphist <- object[x]
            class(temphist) <- c('capthist', 'list')
            ## recursive call
            rbind(temphist, renumber = renumber, pool=NULL, verify = FALSE)
        }
        
        temp <- lapply (pool, getpooled)
        if (length(temp)==1) {
            temp <- temp[[1]]
            class(temp) <- 'capthist'
        }
        else {
            class (temp) <- c('capthist', 'list')
            if (is.null(names(pool)) | any(names(pool) == ""))
                names(temp) <- sapply(temp,session)
            else {
                session(temp) <- names(pool)
            }
        }
        ## do it once
        if (verify) {
            verify(temp)
        }
        return(temp)
    }
    else {
        
        ## Case 3
        ## 1 to several elementary capthist objects
        ## conventional rbind, given compatible traps, covariates, noccasions
        ## optional renumbering
        check <- function (x) {
            if (!is(x,'capthist'))
                stop ("all arguments must be 'capthist' objects")
            if (is.null(covariates(x)) != is.null(covariates(object[[1]]) ))
                stop ("covariates must be provided for all or none")
            if (is.null(Tu(x)) != is.null(Tu(object[[1]])))
                stop ("unmarked sightings Tu must be provided for all or none")
            if (is.null(Tm(x)) != is.null(Tm(object[[1]])))
                stop ("nonID sightings Tu must be provided for all or none")
            if (any(dim(x)[-1] != dim(object[[1]])[-1]))
                stop ("varying numbers of occasions and/or detectors ",
                      "in rbind.capthist", call. = FALSE)
            notPoolPoly <- !all(detector(traps(object[[1]])) %in% c('polygon','polygonX'))
            if (!identical(traps(x), traps(object[[1]])) & notPoolPoly)
                stop ("cannot pool capthist with different",
                      " detector arrays in rbind.capthist", call. = FALSE)
        }
        sapply (object, check)
        
        ##################################################
        ## form new object
        
        temp <- abind(object, along = 1, hier.names = TRUE) 
        
        ## row names
        if (renumber) {
            row.names(temp) <- 1:nrow(temp)
        }
        else {
            ## use original if unique, otherwise default to hierarchical
            an <- unlist(sapply(object, row.names, simplify = FALSE))
            if (any(duplicated(an))) 
                warning("renaming rows to avoid duplicate names")
            else
                row.names(temp) <- an
        }
        
        class(temp) <- 'capthist'
        
        ##################################################
        ## traps
        trps <- traps(object[[1]])
        mergepoly <- all(detector(trps) %in% c('polygon'))
        if (mergepoly) {
          
          ## 2022-02-13  untested
          matlist <- lapply(traps(object), as.matrix)
          polys <- st_sfc(st_polygon(matlist))
          tmp2 <- st_union(polys) # does not dissolve internal boundary
          trps <- data.frame(st_coordinates(tmp2)[,1:2])
          names(trps) <- c('x','y')
          rownames(trps) <- 1:nrow(trps)
          class(trps) <- c('traps', 'data.frame')
          detector(trps) <- detector( traps(object[[1]]))
          polyID(trps) <- rep(1,nrow(trps))
          ## note any covariates have been abandoned
        }
        traps(temp) <- trps
        
        ##################################################
        ## covariates
        
        tempcov <- covariates(object)
        covnamelist <- lapply (tempcov, names)
        covnames <- Reduce(intersect, covnamelist)
        
        if (length(covnames) > 0) {
            tempcov <- lapply(tempcov, function(x) x[,covnames, drop = FALSE])
            tempcov <- do.call (rbind, tempcov)   ## rbind.data.frame
            row.names(tempcov) <- row.names(temp)
            covariates(temp) <- tempcov
        }
        else
            covariates(temp) <- NULL
        
        ##################################################
        ## sightings
        ## either all-scalar or all-matrix
        
        if (!is.null(Tu(object[[1]])))
            Tu(temp) <- sum(Tu(object))
        if (!is.null(Tm(object[[1]])))
            Tm(temp) <- sum(Tm(object))
        
        ##################################################
        ## polygon or transect coordinates xy
        
        tempxy <-  lapply(object, xy)
        xy(temp) <- do.call(rbind, tempxy)
        
        ##################################################
        ## telemetry coordinates cf join() no merge of identities
        
        if ('telemetry' %in% detector(traps(temp))) {
            newtelem <- lapply(object, telemetryxy)
            ntelem <- sapply(newtelem, length)
            newtelem <- unlist(newtelem, recursive = FALSE)
            names(newtelem) <- paste(rep(names(object), ntelem), names(newtelem), sep='.')
            telemetryxy(temp) <- newtelem
        }
        
        ##################################################
        ## signal
        tempsig <- lapply(object, signalframe)
        signalframe(temp) <- do.call(rbind, tempsig)
        
        if (!is.null(signalframe(temp))) {
            signal(temp) <- do.call(c, tempsig)
            cutvals <- sapply(object, function(x) attr(x,'cutval'))
            attr(temp, 'cutval') <- max(cutvals)
            temp <- subset(temp, cutval = max(cutvals))
        }
        else {
            if (!all(sapply(tempsig, is.null)))
                stop ("signal attribute missing in one or more sessions")
        }
        
        ##################################################
        ## messy problem of correct order of detections
        ## 2021-05-19 NOT RESOLVED for xy, signal
        if (!is.null(xy(temp)) | !is.null(signalframe(temp))) {
            occ <- unlist(lapply(object, occasion))
            ID  <- lapply(object, animalID, names = FALSE)
            maxID <- suppressWarnings(sapply(ID, max))
            nID <- sapply(ID, length)
            ID <- unlist(ID)
            uniqueID <- ID + rep(c(0, cumsum(maxID[-length(maxID)])), nID)
            trp <- unlist(lapply(object, trap))
            neworder <- order (occ, uniqueID, trp)
            if (!is.null(xy(temp))) {
              warning ("secr 4.5 rbind.capthist not tested for polygon data")
                xy(temp) <- xy(temp)[neworder,,drop=F]
            }
            if (!is.null(signalframe(temp))) {
              warning ("secr 4.5 rbind.capthist not tested for signal data")
              signalframe(temp) <- signalframe(temp)[neworder,,drop=F]
            }
        }
        ##################################################
        ## name new sessions
        session (temp) <- paste(names(object), collapse='+')
      
        ##################################################
        ## optionally verify
        if (verify) {
            verify(temp)
        }

        ##################################################
        ## return
        temp
    }
}
###############################################################################

append.capthist <- function (..., synchronous = TRUE)
{
  dots <- match.call(expand.dots = FALSE)$...
  allargs <- list(...)
  names(allargs) <- lapply(dots, as.character)
  if (length(dots)==1) object <- allargs[[1]]
  else object <- allargs
  dims <- sapply(object, dim)
  dimsum <- apply(dims,1,sum)
  if (synchronous) dimsum[2] <- max(dims[2,])
  ch <- array(0, dim=dimsum)
  class (ch) <- 'capthist'
  traplist <- lapply(object,traps)
  startn <- starts <- startk <- 1
  for (i in 1:length(object)) {
    n <- dims[1,i]
    s <- dims[2,i]
    k <- dims[3,i]
    ch[startn:(startn+n-1), starts:(starts+s-1), startk:(startk+k-1)] <- object[[i]][]
    covariates(traplist[[i]])$sub <- factor(rep(i, k), levels=1:length(object))
    if (!is.null(usage(traplist[[i]]))) {
      if (synchronous) {
        # fill extra occasions as needed
        usage(traplist[[i]]) <- cbind(usage(traplist[[i]], 
            matrix(0, nrow=k, ncol = dimsum[2] - dims[2,i])))
      }
      else {
        # construct usage matrix for all occasions
        tmp <- matrix(0, nrow = k, ncol = dimsum[2])
        tmp[,starts:(starts+s-1)] <- usage(traplist[[i]])
        usage(traplist[[i]]) <- tmp
      }
    }
    startn <- startn + n
    startk <- startk + k
    if (!synchronous) starts <- starts + s
    
  }
  cov <- do.call(rbind, lapply(object, covariates))
  if (nrow(cov) > 0)
    covariates(ch) <- cov
  traps(ch) <- do.call(rbind, traplist)
  row.names(ch) <- 1:nrow(ch)
  ch
}
###############################################################################


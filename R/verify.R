############################################################################################
## package 'secr'
## verify.R
## 2009 09 18, 2009 09 19, 2009 09 20, 2009 10 02, 2009 11 05, 2009 11 13
## 2010 05 02 removed erroneous ref to 'areabinary' detector
## 2011 02 15 validpoly
## 2011 03 19 adjustments to allow unmarked; need more complete check of unmarked
## 2012 02 02 signal check applied at level of whole sound
## 2012-10-22 xylist checked
## 2013-05-09 tweak to avoid error in checkcovariatelevels when no covariates
## 2015-01-06 stop if mixture of NULL and non-NULL covariates
## 2015-10-03 resight data
## 2015-10-12 verify.traps error messages
## 2015-11-02 xyinpoly moved to utility.R
## 2016-10-06 secr 3.0 revamped
## 2017-10-18 zerohist check ()
## 2017-10-19 capped detector check

## 2017-01-27 future telemetry checks:
##                cannot have occasions with no detections    
## 2020-08-14 bug in checkcovariatelevels caused by bug in stringsAsFactors
## 2020-08-27 verify.traps did not report markocc checks
## 2020-08-27 verify.capthist returns capture checks invisibly
## 2022-05-07 verify.capthist checks nontarget matrix if present
############################################################################################

verify <- function (object, report, ...) UseMethod("verify")

verify.default <- function (object, report, ...) {
  cat ('no verify method for objects of class', class(object), '\n')
}
############################################################################################

## from is.integer help page
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

overlapcells <- function (xy) {
  vertexinside <- function (a,b) {
    OK <- FALSE
    for (k in 1:4) {
      temp <- insidecpp(unlist(a[k,]), 0, 3, as.matrix(b))
      if (any(temp)) OK <- TRUE
      temp <- insidecpp(unlist(b[k,]), 0, 3, as.matrix(a))
      if (any(temp)) OK <- TRUE
    }
    OK
  }
  spacex <- attr(xy, 'spacex')
  spacey <- attr(xy, 'spacey')
  fuzz <- 1e-10
  spx2 <- spacex/2 - fuzz
  spy2 <- spacey/2 - fuzz
  xy <- as.matrix(xy)
  nr <- nrow(xy)
  if (nr<2)
    FALSE
  else {
    pixel <- matrix(ncol = 2, c(-spx2,-spx2,spx2,spx2,-spx2,-spy2,spy2,spy2,-spy2,-spy2))
    overlap <- matrix(FALSE, nrow = nr, ncol = nr)
    for (i in 1:(nr-1))
      for (j in (i+1):nr)
      {
        verti <- t(apply(pixel, 1, function (x) xy[i,] + x))
        vertj <- t(apply(pixel, 1, function (x) xy[j,] + x))
        if ((length(verti)>0) & (length(vertj)>0))
          overlap[i,j] <- vertexinside(verti,vertj)
      }
    any (overlap, na.rm=T)
  }
}
############################################################################################

overlappoly <- function (xy, polyID) {
  ## overlap of component polygons?
  vertexinside <- function (a,b) {
    OK <- FALSE
    n.a <- nrow(a)
    n.b <- nrow(b)
    a <- as.matrix(a)
    b <- as.matrix(b)
    for (k in 1:n.a) {
      temp <- insidecpp(unlist(a[k,]), 0, n.b-1, as.matrix(b))
      if (any(temp)) OK <- TRUE
    }
    for (k in 1:n.b) {
      temp <- insidecpp(unlist(b[k,]), 0, n.a-1, as.matrix(a))
      if (any(temp)) OK <- TRUE
    }
    OK
  }
  lxy <- split (xy, polyID)
  nr <- length(lxy)
  if (nr<2)
    FALSE
  else {
    overlap <- matrix(FALSE, nrow = nr, ncol = nr)
    for (i in 1:(nr-1))
      for (j in (i+1):nr)
      {
        overlap[i,j] <- vertexinside(lxy[[i]], lxy[[j]])
      }
    any (overlap, na.rm=T)
  }
}
############################################################################################

validpoly <- function (xy, polyID, nx = 500) {
  ## check intersections of perimeters with vertical lines
  OKpoly <- function (xy) {
    cross <- rep(0, nx)
    x <- seq(min(xy[,1]), max(xy[,1]), length=nx)
    for (i in 1:nx) {
      tempx <- xy[,1] - x[i]
      cross[i] <- sum (sign(tempx[-1]) != sign(tempx[-length(tempx)]))
    }
    cross <= 2
  }
  lxy <- split (xy, polyID)
  temp <- lapply(lxy, OKpoly)
  all(unlist(temp))
}
############################################################################################

xyontransect <- function (xy, trps, tol = 0.01) {
  ptontransect <- function (i,k) {
    ## is point i on transect k?
      transectxy <- as.matrix(lxy[[k]])
      nr <- nrow(transectxy)
      ontransectcpp (
          as.matrix (xy[i,,drop=FALSE]),
          as.matrix (transectxy),
          as.integer (0),
          as.integer (nr-1),
          as.double (tol))
  }
  lxy <- split (trps, transectID(trps))
  firsttransect <- function (i) {
    for (k in 1:length(lxy))
      if (ptontransect(i,k)) return(k)
    0
  }
  sapply(1:nrow(xy), firsttransect)
}
############################################################################################

checkcovariatelevels <- function (cov) {
  ## cov is a list of dataframes of covariates
  ## tweak 2013-05-09 to avoid error when no covariates
  xfactor <- function(x) {
    if ((nrow(x)>0) & (ncol(x)>0))
      names(x)[sapply(x,is.factor)]
    else
      NULL
  }
  #############################################
  ## 2020-05-15 convert all character to factor
  cov <- lapply(cov, stringsAsFactors)
  #############################################
  
  factornames <- sapply(cov, xfactor)
  factornames <- unique(unlist(factornames))
  if (length(factornames) > 0) {
    if (any(sapply(cov, function(x) !all(factornames %in% names(x))))) {
      return(FALSE)
    }
    else {
      checklevels <- function(variable) {
        baselevels <- levels(cov[[1]][,variable])
        lev <- lapply(cov, function(x) levels(x[,variable]))
        all(sapply(lev[-1], function(x) identical(x,baselevels)))
      }
      return(all(sapply(factornames, checklevels)))
    }
  }
  else
    return(TRUE)
}
############################################################################################

verify.traps <- function (object, report = 2, ...) {
  
  ## Check internal consistency of 'traps' object
  ##
  ## -- Number of rows in dataframe of detector covariates differs expected
  ## -- Number of detectors in usage matrix differs from expected
  ## -- Occasions with no used detectors
  
  if (!inherits(object, 'traps')) {
    stop ("object must be of class 'traps'")
  }
  
  if (inherits(object, 'list')) {
    temp <- lapply (object, verify, report = min(report,1))
    anyerrors <- any(sapply(temp, function(x) x$errors))
    
    ## check covariate factor levels conform across sessions
    
    if (!all(sapply(covariates(object), is.null))) {
      if (any(sapply(covariates(object), is.null)))
        stop ("mixture of NULL and non-NULL trap covariates in different sessions")
      trapcovariatelevelsOK <- checkcovariatelevels(covariates(object))
      if (!trapcovariatelevelsOK & report>0) {
        warning ('Levels of factor trap covariate(s) differ between sessions - ",
                 "use shareFactorLevels()', call. = FALSE)
      }
    }
    
    if ((report == 2) & !anyerrors)
      cat('No errors found :-)\n')
    invisible(list(errors = anyerrors, bysession = temp))
  }
  else {
    
    poly <- all(detector(object) %in% c('polygon','polygonX'))
    telemonly <- all(detector(object) %in% 'telemetry')
    
    detectorsOK <- TRUE
    usagedetectorsOK <- TRUE
    usagedetectors2OK <- TRUE
    usagenonzeroOK <- TRUE
    markoccOK <- TRUE
    markocc2OK <- TRUE
    areaOK <- TRUE
    polyIDOK <- TRUE
    polyconvexOK <- TRUE
    trapcovariatesOK <- TRUE
    telemOK <- TRUE
    
    if (!is.null(covariates(object)))
      if ((ncol(covariates(object)) == 0 ) |
          (nrow(covariates(object)) == 0 )) covariates(object) <- NULL
    
    ## 1
    trapNAOK <- !any(is.na(object))
    
    ## 2
    if (!is.null(covariates(object)))
      trapcovariatesOK <- nrow(covariates(object)) == ndetector(object)
    
    ## 3
    uniquedetectors <- unique(detector(object))
    nontelemetrydetectors <- uniquedetectors[uniquedetectors != "telemetry"]
    detectorsOK <- (length(uniquedetectors)==1) |
      all(nontelemetrydetectors %in% .localstuff$pointdetectors)
    
    ## 'usage' of traps
    if (!is.null(usage(object))) {
      ## 4
      usagedetectorsOK <- nrow(usage(object)) == ndetector(object)
      
      ## 5
      if (length(detectorcode(object))>1)
        usagedetectors2OK <- length(detectorcode(object)) == ncol(usage(object))
      
      ## 6
      usagecount <- apply(usage(object),2,sum)
      usagenonzeroOK <- !any(usagecount == 0)
    }
    else usagecount <- rep(NA, ncol(object))
    
    ## 7
    if (sighting(object)) {
      if (!is.null(usage(object)))
        markoccOK <- length(markocc(object)) == ncol(usage(object))
      ## 8
      markocc2OK <- all(markocc(object) %in% c(-2,-1,0,1))
    }
    
    ## 9
    if (poly) {
      areaOK <- !overlappoly (object, polyID(object))
    }
    
    ## 10
    if (poly) {
      polyIDOK <- (length(polyID(object)) == nrow(object)) &
        is.factor(polyID(object))
    }
    
    ## 11
    if (poly) {
      polyconvexOK <- validpoly (object, polyID(object))
    }
    
    ## 12
    if (telemonly) {
      telemOK <- nrow(object) == 1
    }
    
    errors <- !all(c(trapNAOK,
      trapcovariatesOK,
      detectorsOK,
      usagedetectorsOK,
      usagedetectors2OK,
      usagenonzeroOK,
      markoccOK,
      markocc2OK,
      areaOK,
      polyIDOK,
      polyconvexOK,
      telemOK))
    
    if (report > 0) {
      if (errors) {
        if (!trapNAOK) {
          cat ('Missing detector coordinates not allowed\n')
        }
        if (!trapcovariatesOK) {
          cat ('Wrong number of rows in dataframe of detector covariates\n')
          cat (uniquedetectors, '\n')
        }
        if (!detectorsOK) {
          cat ('Invalid combination of detector types\n')
          cat ('traps :', ndetector(object), 'detectors\n')
          cat ('usage(traps) :', nrow(usage(object)), 'detectors\n')
        }
        if (!usagedetectorsOK) {
          cat ('Conflicting number of detectors in usage matrix\n')
          cat ('traps :', ndetector(object), 'detectors\n')
          cat ('usage(traps) :', nrow(usage(object)), 'detectors\n')
        }
        if (!usagedetectors2OK) {
          cat ('Conflicting number of occasions in usage matrix\n')
          cat ('detector attribute :', length(detectorcode(object)), 'occasions\n')
          cat ('usage(traps) :', ncol(usage(object)), 'occasions\n')
        }
        if (!usagenonzeroOK) {
          cat ("Occasions when no detectors 'used'\n")
          cat ((1:length(usagecount))[usagecount==0], '\n')
        }
        if (!markoccOK) {
          cat ("Conflicting number of occasions in markocc and usage \n")
          cat ('markocc : ', length(markocc(object)), 'occasions \n')
          cat ('usage : ', ncol(usage(object)), 'occasions \n')
        }
        if (!markocc2OK) {
          cat ("Values in markocc should be -2, -1, 0 or 1 \n")
          cat ('markocc(object) \n')
        }
        if (!areaOK) {
          cat ("Search areas overlap, or no search area specified \n")
        }
        if (!polyIDOK) {
          cat ("Invalid polyID \n")
        }
        if (!polyconvexOK) {
          cat ("The boundary of at least one polygon is concave east-west \n")
        }
        if (!telemOK) {
          cat ("More than one row in notional traps object for telemetry \n")
        }
      }
    }
    
    if ((report == 2) & !errors) message('No errors found :-)')
    
    out <- list(errors = errors,
      trapNAOK = trapNAOK,
      trapcovariatesOK = trapcovariatesOK,
      detectorsOK = detectorsOK,
      usagedetectorsOK = usagedetectorsOK,
      usagedetectors2OK = usagedetectors2OK,
      usagenonzeroOK = usagenonzeroOK,
      markoccOK = markoccOK,
      markocc2OK = markocc2OK,
      areaOK = areaOK,
      polyIDOK = polyIDOK,
      polyconvexOK = polyconvexOK,
      usagecount = usagecount,
      telemOK = telemOK
    )
    
    invisible(out)
    
  }
}
############################################################################################

verify.capthist <- function (object, report = 2, tol = 0.01, ...) {
  
  ## Check internal consistency of 'capthist' object
  ##
  ## -- 'traps' component present
  ## -- verify(traps)
  ## -- No live releases
  ## -- Live detection(s) after reported dead
  ## -- More than one capture in single-catch trap(s)
  
  ## -- Number of rows in 'traps' object not compatible with reported detections
  ## -- Number of rows in dataframe of individual covariates differs from capthist
  ## -- Number of occasions in usage matrix differs from capthist
  ## -- Detections at unused detectors
  
  ## -- resighting attributes Tu, Tm compatible if present
  ## -- no resightings on marking occasions, or new animals on resighting occasions
  if (!inherits(object, 'capthist'))
    stop ("object must be of class 'capthist'")
  if (inherits(object, 'list')) {
    temp <- lapply (object, verify, report = min(report, 1))
    anyerrors <- any(sapply(temp, function(x) x$errors))
    ## check covariate factor levels conform across sessions
    if (!all(sapply(covariates(object), is.null))) {
      if (any(sapply(covariates(object), is.null)))
        stop ("mixture of NULL and non-NULL individual covariates in different sessions")
      covariatelevelsOK <- checkcovariatelevels(covariates(object))
      if (!covariatelevelsOK & report>0) {
        warning ('Levels of factor covariate(s) differ between sessions - ",
                 "use shareFactorLevels()', call. = FALSE)
      }
    }
    if ((report == 2) & !anyerrors)
      cat('No errors found :-)\n')
    invisible(list(errors = anyerrors, bysession = temp))
  }
  else {
    ## preliminaries
    object <- check3D(object)
    detectortype <- expanddet(object)  # vector length noccasions
    telemetrytype <- telemetrytype(traps(object))
    signal <- all(detectortype %in% c('signal','signalnoise'))
    unmarked <- all(detectortype %in% c('unmarked'))
    poly <- all(detectortype %in% c('polygon', 'polygonX'))
    transect <- all(detectortype %in% c('transect', 'transectX'))
    telem <- any(detectortype %in% c('telemetry'))
    
    markocc <- markocc(traps(object))
    sighting <- sighting(traps(object))
    allsighting <- if (sighting) !any(markocc>0) else FALSE
    
    NAOK <- TRUE
    deadOK <- TRUE
    usageOK <- TRUE
    detectorsOK <- TRUE
    usageoccasionsOK <- TRUE
    detectornumberOK <- TRUE
    detectorconflcts <- NULL
    zerohistOK <- TRUE
    singleOK <- TRUE
    multiOK <- TRUE
    cappedOK <- TRUE
    binaryOK <- TRUE
    countOK <- TRUE
    cutvalOK <- TRUE
    signalOK <- TRUE
    xyOK <- TRUE
    xyinpolyOK <- TRUE
    xyontransectOK <- TRUE
    rownamesOK <- TRUE
    sightingsOK <- TRUE
    sightingusageOK <- TRUE
    MOK <- TRUE
    ROK <- TRUE
    telemOK <- TRUE
    nontargetdimOK <- TRUE
    nontargetOK <- TRUE
    
    if (!is.null(covariates(object)))
      if ((ncol(covariates(object)) == 0 ) |
          (nrow(covariates(object)) == 0 ))
        covariates(object) <- NULL
    
    ## 1
    trapspresentOK <- !is.null(traps(object))
    
    ## standalone check of detectors
    ## this is done one session at a time
    ## so does not check between session agreement of covariates
    if (trapspresentOK)
      trapcheck <- verify(traps(object), report = 0)  ## delay reporting
    else
      trapcheck <- list(errors=TRUE)
    
    ## 2
    trapsOK <- !trapcheck$errors
    
    ## 3
    if (length(object)==0) {
      detectionsOK <- unmarked   ## and presence? 2011-10-01
    }
    else  {
      
      detectionsOK <- sum(object[object>0]) > 0
      
      ## 4
      NAOK <- !any(is.na(object))
      
      ## 5
      if (signal) {
        ## must have cutval; deads not allowed
        if (length(attr(object,'cutval')) != 1) cutvalOK <- FALSE
        else {
          maxbyanimal <- tapply(signal(object), animalID(object, sortorder = 'ksn'), max)
          maxbyanimal <- maxbyanimal[!is.na(maxbyanimal)]  ## because NA OK 2012-02-10
          if (any(maxbyanimal < attr(object,'cutval')))
            cutvalOK <- FALSE
        }
        ## 6
        if (length(signal(object)) != sum(abs(object)))
          signalOK <- FALSE
      }
      ## 7
      else {
        fn <- function(x) {   ## tweaked 2017-11-29
          x <- apply(x,1,min)   ## by occasion
          dead <- min(x)<0
          last <- tail(x[x!=0],1)
          if (length(last)==0) last <- 0
          dead & (last>0)
        }
        undead <- apply(object, 1, fn)
        deadOK <- !any(undead)
        if (!deadOK) {
          reincarnated <- object[undead,,, drop=F]
        }
      }
      
      ###################
      ## zero histories
      ## all-zero histories are meaningful only for concurrent telemetry and
      ## sighting-only data with known number marked
      if (telemetrytype != "concurrent" & !allsighting) {
        zerohistOK <- all(apply(abs(object),1,sum)>0)
      }
      
      ###################
      ## binary detectors
      
      detectorcodes <- detectorcode(traps(object), MLonly = FALSE, noccasions = ncol(object))
      
      ## 8
      ## single, multi, capped, proximity, polygonX, transectX, signal and signalnoise must be binary
      multiples <- sum(abs(object[, detectorcodes %in% c(-1,0,1,3,4,5,8,12), , drop = FALSE])>1)
      binaryOK <- multiples == 0
      
      ## 9
      fn <- function (x) duplicated(abs(x)[x!=0])
      ## no more than one obs per animal and per trap on any occasion
      multiplei <- apply(object[, detectorcodes==-1, , drop = FALSE], 1:2, fn)
      multiplek <- apply(object[, detectorcodes==-1, , drop = FALSE], 2:3, fn)
      singleOK <- !any(unlist(multiplei)) & !any(unlist(multiplek))
      
      ## 10
      ## no more than one obs per animal per occasion
      multiplei <- apply(object[, detectorcodes %in% c(0,3,4), , drop = FALSE], 1:2, fn)
      multiOK <- !any(unlist(multiplei))
      
      ## no more than one individual per detector per occasion
      multiplek <- apply(object[, detectorcodes==8, , drop = FALSE], 2:3, fn)
      cappedOK <- !any(unlist(multiplek))
      
      ###################
      ## count detectors
      
      ## 11
      ## count, polygon, transect cannot be dead (?)
      countOK <- all(object[, detectorcodes %in% c(2,6,7), ] >= 0)
      
      ##################
      
    }
    
    ## 12
    if (nrow(object) > 0) {
      if (poly | transect) {
        detectornumberOK <- length(levels(polyID(traps(object)))) == dim(object)[3]
      }
      else {
        detectornumberOK <- dim(object)[3] == nrow(traps(object))
      }
    }
    
    ## 13
    covariatesOK <- ifelse(is.null(covariates(object)),
      TRUE,
      nrow(covariates(object)) == nrow(object))
    
    ## is 'usage' of traps consistent with reported detections?
    if (!is.null(usage(traps(object)))) {
      conflcts <- 0
      
      ## 14
      usageoccasionsOK <- ncol(usage(traps(object))) == ncol(object)
      
      if (detectionsOK) {
        # 2012-12-17
        # notused <- !usage(traps(object))   ## traps x occasions
        notused <- usage(traps(object)) == 0 ## traps x occasions
        if (trapcheck$usagedetectorsOK & usageoccasionsOK) {
          tempobj <- aperm(object, c(2,3,1))   ## occasion, traps, animal sKn
          tempuse <- array(t(usage(traps(object))>0), dim=dim(tempobj)) # repl to fill
          # bug 2017-11-29 && changed to &
          conflcts <- (abs(tempobj)>0) & (tempuse==0)
          tempobjmat <- array(tempobj[,,1], dim= dim(tempobj)[1:2])
          occasion <- rep(row(tempobjmat), dim(tempobj)[3])
          detector <- rep(col(tempobjmat), dim(tempobj)[3])
          ID <- rep(rownames(object), rep(prod(dim(tempobj)[1:2]), nrow(object)))
          detectorconflcts <- as.data.frame(cbind(ID,detector,occasion)[conflcts,])
        }
      }
      
      ## 15
      usageOK <- sum(conflcts)==0
      
    }
    
    if (poly) {
      xy <- xy(object)
      ## 16
      if (all(detectortype %in% c('polygon')))
        xyOK <- nrow(xy) == sum(abs(object))
      else
        xyOK <- nrow(xy) == sum(abs(object)>0)
      # assumes detection sequence is ksn in trap() and xy()
      inpoly <- xyinpoly(xy(object), traps(object))
      inpoly <- inpoly == trap(object, names = FALSE, sortorder = "ksn")  
      ## 17
      xyinpolyOK <- all(inpoly)
    }
    ## 17
    if (transect) {
      xy <- xy(object)
      ## 16 again
      if (all(detectortype == 'transect'))
        xyOK <- nrow(xy) == sum(abs(object))
      else
        xyOK <- nrow(xy) == sum(abs(object)>0)
      ontransect <- xyontransect(xy(object), traps(object), tol = tol)
      ontransect <- ontransect == trap(object, names = FALSE, sortorder = "ksn")
      ## 18
      xyontransectOK <- all(ontransect)
    }
    ## 16 again
    if (telem) {
      telemocc <- detectortype=='telemetry'
      telemdet <- apply(abs(object[,telemocc,,drop=FALSE]),1:2,sum)
      captdet <- apply(abs(object[,!telemocc,,drop=FALSE]),1:2,sum)
      
      telemdet.byoccasion <- apply(telemdet,2,sum)
      telemdet.byanimal <- apply(telemdet,1,sum)
      captdet.byanimal <- apply(captdet,1,sum)
      animaldet <- sapply(telemetryxy(object),nrow)
      telemetrd <- telemetered(object)
      if (any(telemdet.byoccasion == 0))
        # telemetry occasions with no detections
        telemOK <- FALSE
      if (length(animaldet) != sum(telemetrd))
        # mismatch of animals
        telemOK <- FALSE
      else if (any(sort(names(animaldet)) != sort(names(telemdet.byanimal[telemetrd]))))
        # mismatch of names
        telemOK <- FALSE
      else if (any(animaldet != telemdet.byanimal[telemetrd][names(animaldet)]))
        # number of observations does not match
        telemOK <- FALSE
      else if (telemetrytype=='dependent' & any(captdet.byanimal==0))
        # all-zero detection histories for dependent telemetry
        telemOK <- FALSE
      else if (telemetrytype=='independent' & any(captdet.byanimal[telemetrd]!=0))
        # non-zero detection histories for independent telemetry
        telemOK <- FALSE
    }
    ## 19
    ## 2021-03-03 added condition so OK when no rows; fixed 2021-11-16
    rownamesOK <- (!any(duplicated(rownames(object))) &&  
        !is.null(rownames(object))) || (nrow(object)==0) 
    if (!is.null(rownames(object)))
      rownamesOK <- rownamesOK & !any(is.na(rownames(object)))
    
    ## 20
    # superceded by telemOK 2017-01-27
    # if (nrow(object)>0) {
    #     zeros <- apply(abs(object)>0,1,sum)==0
    #     if (!allsighting) {
    #         xyl <- telemetryxy(object)
    #         if (!is.null(xyl) | any(zeros)) {
    #             xylistOK <- all(names(xyl) %in% row.names(object))
    #             if (!all(row.names(object)[zeros] %in% names(xyl)))
    #                 xylistOK <- FALSE
    #         }
    #     }
    # }
    
    ## -- resighting attributes Tu, Tm compatible if present
    ## -- no resightings on marking occasions, or new animals on resighting occasions
    
    ## 21,22
    if (sighting) {
      Tu <- Tu(object)
      Tm <- Tm(object)
      nocc <- ncol(object)
      K <- ndetector(traps(object))
      r <- numeric(nocc)
      usge <- usage(traps(object))
      if (is.null(usge)) usge <- 1
      if (length(markocc) != nocc) sightingsOK <- FALSE
      ## allow scalar summed sighting counts 2015-10-31
      if (!is.null(Tu)) {
        if (length(Tu) > 1) {
          if (any((Tu>0) & (usge==0))) sightingusageOK <- FALSE
          if (ncol(Tu) != nocc) sightingsOK <- FALSE
          if (nrow(Tu) != K) sightingsOK <- FALSE
          r <- r + apply(Tu,2,sum)
        }
        if (any(Tu<0)) sightingsOK <- FALSE
        if (!all(is.wholenumber(Tu))) sightingsOK <- FALSE
      }
      if (!is.null(Tm)) {
        if (length(Tm) > 1) {
          # 2016-10-13 fixed bug Tu>Tm
          if (any((Tm>0) & (usge==0))) sightingusageOK <- FALSE
          if (nrow(Tm) != K) sightingsOK <- FALSE
          if (ncol(Tm) != nocc) sightingsOK <- FALSE
          r <- r + apply(Tm,2,sum)
        }
        if (any(Tm<0)) sightingsOK <- FALSE
        if (!all(is.wholenumber(Tm))) sightingsOK <- FALSE
      }
      ## 23
      if (sightingsOK) {
        if (length(Tu)>1) {
          u <- unlist(counts(object, 'u'))[1:nocc]
          if ( any((u > 0) & (markocc<1) & !allsighting) ) MOK <- FALSE
        }
        ## 24
        if (length(Tm)>1)
          if ( any((r > 0) & (markocc>0))) ROK <- FALSE
      }
    }
    nontarget <- attr(object, 'nontarget', exact = TRUE)
    if (!is.null(nontarget)) {
        nontarget <- as.matrix(nontarget)
        if (nrow(nontarget) != dim(object)[3] ||
                ncol(nontarget) != dim(object)[2])
            nontargetdimOK <- FALSE
        if (nontargetdimOK && detectortype[1] %in% c('single', 'capped')) {
            # capture and nontarget mutually exclusive 
            obs <- t(apply(abs(object), 2:3, sum))
            if (any(obs & nontarget))
                nontargetOK <- FALSE
        }
    }
    
    errors <- !all(c(trapspresentOK,
      trapsOK,
      detectionsOK,
      NAOK,
      cutvalOK,
      signalOK,
      deadOK,
      zerohistOK,
      binaryOK,
      singleOK,
      multiOK,
      cappedOK,
      countOK,
      detectornumberOK,
      covariatesOK,
      usageoccasionsOK,
      usageOK,
      xyOK,
      xyinpolyOK,
      xyontransectOK,
      rownamesOK,
        
      sightingsOK,
      sightingusageOK,
      MOK,
      ROK,
      telemOK,
      nontargetdimOK,
      nontargetOK
        ))
    
    if (report > 0) {
      if (errors) {
        cat ('Session', session(object), '\n')
        
        if (!trapspresentOK) {
          cat ('No valid detectors\n')
        }
        if (!trapsOK) {
          cat ('Errors in traps\n')
          if (!trapcheck$trapNAOK) {
            cat ('Missing detector coordinates not allowed\n')
          }
          if (!trapcheck$trapcovariatesOK) {
            cat ('Wrong number of rows in dataframe of detector covariates\n')
            cat ('traps(capthist) :', ndetector(traps(object)), 'detectors\n')
            cat ('covariates(traps(capthist)) :', nrow(covariates(traps(object))), 'detectors\n')
          }
          if (!trapcheck$usagedetectorsOK) {
            cat ('Conflicting number of detectors in usage matrix\n')
            cat ('traps(capthist) :', ndetector(traps(object)), 'detectors\n')
            cat ('usage(traps(capthist)) :', nrow(usage(traps(object))), 'detectors\n')
          }
          if (!trapcheck$usagenonzeroOK) {
            cat ("Occasions when no detectors 'used'\n")
            cat ((1:length(trapcheck$usagecount))[trapcheck$usagecount==0], '\n')
          }
        }
        
        if (!detectionsOK) {
          cat ('No live releases\n')
        }
        
        if (!NAOK) {
          cat ('Missing values not allowed in capthist\n')
        }
        
        if (!cutvalOK) {
          cat ('Signal less than cutval or invalid cutval\n')
        }
        
        if (!signalOK) {
          cat ('Signal attribute does not match detections\n')
        }
        
        if (!deadOK) {
          cat ('Recorded alive after dead\n')
          print(reincarnated)
        }
        
        if (!zerohistOK) {
          cat ('Empty histories allowed only with concurrent telemetry or sighting-only data\n')
        }
        
        if (!binaryOK) {
          cat ('More than one detection per detector per occasion at binary detector(s)\n')
        }
        
        if (!singleOK) {
          cat ('More than one capture in single-catch trap(s)\n')
        }
        
        if (!multiOK) {
          cat ('Animal trapped at more than one detector\n')
        }
        
        if (!cappedOK) {
          cat ('More than one animal at detector\n')
        }
        
        if (!countOK) {
          cat ('Count(s) less than zero\n')
        }
        
        if (!detectornumberOK) {
          cat ('traps object incompatible with reported detections\n')
          cat ('traps(capthist) :', ndetector(traps(object)), 'detectors\n')
          cat ('capthist :', dim(object)[3], 'detectors\n')
        }
        
        if (!covariatesOK) {
          cat ('Wrong number of rows in dataframe of individual covariates\n')
          cat ('capthist :', nrow(object), 'individuals\n')
          cat ('covariates(capthist) :', nrow(covariates(object)), 'individuals\n')
        }
        if (!usageoccasionsOK) {
          cat ('Conflicting number of occasions in usage matrix\n')
          cat ('capthist :', ncol(object), 'occasions\n')
          cat ('usage(traps(capthist)) :', ncol(usage(traps(object))), 'occasions\n')
        }
        if (!usageOK) {
          cat ("Detections at 'unused' detectors\n")
          print(detectorconflcts)
        }
        if (!xyOK) {
          cat ("Polygon or telemetry detector xy coordinates of detections",
            " do not match counts\n")
        }
        if (!xyinpolyOK) {
          cat ("XY coordinates not in polygon\n")
          print (xy(object)[!inpoly,])
        }
        if (!xyontransectOK) {
          cat ("XY coordinates not on transect\n")
          print (xy(object)[!ontransect,])
        }
        if (!rownamesOK) {
          cat ("Duplicated or missing row names (animal ID)\n")
        }
        if (!sightingsOK) {
          cat("Incompatible dimensions of sighting attributes markocc, Tu or Tm\n")
        }
        if (!sightingusageOK) {
          cat("Sightings at unused detectors\n")
          Tu <- Tu(object)
          Tm <- Tm(object)
          bad <- (Tu>0) & (usage(traps(object))==0) & (length(Tu)>1)
          if (sum(bad)>0) {
            cat("Tu\n")
            print(cbind(Detector = row(bad)[bad], Occasion = col(bad)[bad]))
          }
          bad <- (Tm>0) & (usage(traps(object))==0) & (length(Tm)>1)
          if (sum(bad)>0) {
            cat("Tm\n")
            print(cbind(Detector = row(bad)[bad], Occasion = col(bad)[bad]))
          }
        }
        if (!MOK) {
          cat("New individual(s) on sighting-only occasion\n")
          occ <- split(
            occasion(object, sortorder = 'snk'), 
            animalID(object, names = TRUE, sortorder = 'snk')
          )
          firstocc <- sapply(occ, '[', 1)
          print(firstocc[firstocc %in% (1:ncol(object))[markocc(traps(object))<1]])
        }
        if (!ROK) {
          cat("Sighting(s) on marking-only occasion\n")
        }
          if (!nontargetdimOK) {
              cat ("Nontarget dimensions do not match capthist\n")
          }
          if (!nontargetOK) {
              cat ("Nontarget data conflict with capthist\n")
          }
      }
      
      if ((report == 2) & !errors) message('No errors found :-)')
      
    }
    
    captcheck <- list(
      detectionsOK = detectionsOK,
      NAOK = NAOK,
      cutvalOK = cutvalOK,
      signalOK = signalOK,
      deadOK = deadOK,
      zerohistOK = zerohistOK,
      binaryOK = binaryOK,
      singleOK = singleOK,
      multiOK = multiOK,
      cappedOK = cappedOK,
      countOK = countOK,
      detectornumberOK = detectornumberOK,
      covariatesOK = covariatesOK,
      usageoccasionsOK = usageoccasionsOK,
      usageOK = usageOK,
      xyOK = xyOK,
      xyinpolyOK = xyinpolyOK,
      xyontransectOK = xyontransectOK,
      rownamesOK = rownamesOK,
        
      sightingsOK = sightingsOK,
      sightingusageOK = sightingusageOK,
      MOK = MOK,
      ROK = ROK,
      telemOK = telemOK,
        nontargetdimOK = nontargetdimOK,
        nontargetOK = nontargetOK
    )
    
    out <- list(errors = errors, trapcheck = trapcheck, captcheck = captcheck)
    if (!is.null(detectorconflcts)) out$detections.at.unused.detectors <- detectorconflcts
    invisible(out)
  }
}
############################################################################################

verify.mask <- function (object, report = 2, ...) {
  
  ## Check internal consistency of 'mask' object
  ##
  ## valid x and y coordinates
  ## nrow(covariates) = nrow(object)
  ## ...also look at attributes?
  
  if (!inherits(object, 'mask'))
    stop ("object must be of class 'mask'")
  
  if (inherits(object, 'list')) {
    temp <- lapply (object, verify, report = min(report, 1))
    anyerrors <- any(sapply(temp, function(x) x$errors))
    
    ## check covariate factor levels conform across sessions
    if (!all(sapply(covariates(object), is.null))) {
      covariatelevelsOK <- checkcovariatelevels(covariates(object))
      if (!covariatelevelsOK & report>0) {
        warning ('Levels of factor mask covariate(s) differ between sessions - ",
                 "use shareFactorLevels()', call. = FALSE)
      }
    }
    if ((report == 2) & !anyerrors)
      cat('No errors found :-)\n')
    invisible(list(errors = anyerrors, bysession = temp))
  }
  else {
    
    ## 1
    xyOK <- !(is.null(object$x) | is.null(object$y) | any(is.na(object)))
    xyOK <- xyOK & is.numeric(unlist(object))
    
    ## 2
    
    if (!is.null(covariates(object)))
      covariatesOK <- ifelse (nrow(covariates(object))>0,
        nrow(object) == nrow(covariates(object)), TRUE)
    else
      covariatesOK <- TRUE
    
    errors <- !all(c(xyOK, covariatesOK))
    
    if (report > 0) {
      if (errors) {
        ## cat ('Session', session(object), '\n')
        
        if (!xyOK) {
          cat ('Invalid x or y coordinates in mask\n')
        }
        
        if (!covariatesOK) {
          cat ('Number of rows in covariates(mask) differs from expected\n')
        }
      }
      
      if ((report == 2) & !errors) message('No errors found :-)')
    }
    
    out <- list(errors = errors)
    invisible(out)
  }
}
############################################################################################


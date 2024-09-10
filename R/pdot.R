###############################################################################
## package 'secr'
## pdot.R
## return net detection probability in 'traps' for home ranges centred at X
## 2010 07 01 alpha detection fn
## 2010 10 09 extended for other detection functions
## 2010 10 10 esa.plot added
## 2010 10 11 esa.plot.secr added
## 2010 10 18 etc. pdot.contour, buffer.contour
## 2010 10 24 tweaking poly
## 2010 11 26 usage
## 2010 11 26 check for ms traps
## 2010 12 19 more careful handling of detectpars
## 2011 01 24 debugged pdotpoly
## 2011 02 06 allow polygonX, transectX
## 2011 06 13 moved spatialscale to utility.R
## 2012 12 24 binomN = 'usage'
## 2014-03-26 pdot.contour and buffer.contour extended to multi-session traps
## 2014-10-17 userdist fixes
## 2014-11-17 more userdist fixes
## 2015-05-15 fill argument for contours
## 2016-11-12 pdotpoint, markocc
## 2017-04-04 replaced pdotpoly with hdotpoly
## 2018-06-01 esa.plot type replaces as.density
##            esa.plot conditional argument
##            esa.plot CVpdot
## 2019-01-21 poly.habitat argument in pdot.contour and buffer.contour functions etc.
## 2019-07-29 C++
## 2019-12-28 multithreaded
## 2021-05-19 cv: pmax protects against negative argument to sqrt     
## 2022-11-19 esa.plot in separate file
## 2024-09-07 pdot accepts vector or matrix detectpar for g0/lambda0 and sigma 
##            replicated to fill matrix of dimensions ntraps x noccasions (traps are rows)
###############################################################################

## pdot is used in --

## CVa
## CVpdot
## esa*
## derivedSystematic* (via Fewstervarn)
## fx.total*
## make.mask (pdotmin option)
## reparameterize.esa
## [bias.D  disabled]
## pdot.contour
## MCgof

## * has ncores argument

pdot <- function (X, traps, detectfn = 0, detectpar = list(g0 = 0.2, sigma = 25, z = 1),
                  noccasions = NULL, binomN = NULL, userdist = NULL, ncores = NULL) {

    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coord in col 2

    ncores <- setNumThreads(ncores)
    grain <- if (ncores==1) 0 else 1

    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)
    if ((detectfn > 9) & (detectfn<14) & is.null(detectpar$cutval))
        stop ("requires 'cutval' for detectfn 10:13")
    if (ms(traps))
        stop ("requires single-session traps")

    truncate <- ifelse(is.null(detectpar$truncate), 1e+10, detectpar$truncate)

    ntraps <- ndetector(traps)
    if (!is.null(usage(traps))) {
        usge <- usage(traps)
        if (is.null(noccasions)) {
            noccasions <- ncol(usage(traps))
        }
        else {
            if (noccasions < ncol(usage(traps))) {
                warning ("specified noccasions less than ncol of usage matrix")
            }
            if (noccasions > ncol(usage(traps)))
                stop ("specified noccasions exceeds ncol of usage matrix")
        }
    }
    else {
        if (is.null(noccasions))
            stop("must specify noccasions when traps does not have usage attribute")
        usge <- matrix(1, ntraps, noccasions)
    }

    detectpars <- detectpar[parnames(detectfn)]
    gl0 <- detectpars[[1]]   ## g0 or lambda0
    sig <- detectpars[[2]]   ## sigma
    if (length(detectpars)>2) z <- detectpars[[3]] 
        
    
    otherdetpar <- unlist(detectpar[3:4])  # assume scalar 
    ## g0/lambda0 or sigma is a matrix
    nonscalardetpar <- (is.matrix(gl0) || is.matrix(sig))
    gl0 <- matrix(gl0, nrow = ntraps, ncol = noccasions)
    sig <- matrix(sig, nrow = ntraps, ncol = noccasions)
    
    if ((detectfn>9) & (detectfn<14))  {
        otherdetpar <- c(otherdetpar, detectpar$cutval)
    }
    else {
        
    }
    if (length(detectpars)<1) otherdetpar <- c(otherdetpar,0)
    
    miscparm <- numeric(4);   ## dummy
    
    
    dettype <- detectorcode(traps, noccasions = noccasions)
    binomN <- getbinomN (binomN, detector(traps))
    markocc <- markocc(traps)
    
    if (is.null(markocc)) markocc <- rep(1,noccasions)
    if (!inherits(X, 'mask')) {
        X <- matrix(unlist(X), ncol = 2)
    }
    if (any(detector(traps) %in% c('polygon','polygonX','transect', 'transectX'))) {
        if (!is.null(userdist))
            stop("userdist incompatible with polygon-like detectors")
        if (!(detectfn %in% 14:20))
            stop("pdot requires hazard detectfn for polygon-type detectors")
        k <- table(polyID(traps))   ## also serves transectID
        K <- length(k)              ## number of polygons/transects
        k <-  c(k,0)                ## zero terminate for 
        cumk <- cumsum(c(0,table(polyID(traps))))
        convexpolygon <- TRUE
        dim <- if (any(detector(traps) %in% c('transect', 'transectX'))) 1 else 2
            
        warning("assuming convex polygons in pdot()")
        temp <- hdotpolycpp (
          as.matrix(X),
          as.matrix(traps),
          as.matrix(usge),
          as.integer(markocc),
          as.integer(cumk),
          as.integer(detectfn),
          as.double(detectpars),
          as.logical(convexpolygon),
          as.integer(dim),
          as.integer(grain),
          as.integer(ncores))
        1 - exp(-temp)   ## probability detected at least once, given total hazard
    }
    else {
      distmat2 <- getuserdist (traps, X, userdist, sessnum = NA, NULL, NULL, 
                               miscparm, detectfn == 20)
      #-------------------------------------------------------------
      pdotpointcpp(
        as.matrix(X),
        as.matrix(traps),
        as.matrix(distmat2),
        as.integer(dettype),
        as.matrix(usge),
        as.integer(markocc),
        as.integer(detectfn),
        as.matrix(gl0),             # new
        as.matrix(sig),             # new
        as.double(otherdetpar),     # new    
        as.double(miscparm),
        as.double(truncate^2),
        as.integer(expandbinomN(binomN, dettype)),
        as.integer(grain),
        as.integer(ncores)
      )
    }
}
############################################################################################

pdot.contour <- function (traps, border = NULL, nx = 64, detectfn = 0,
                          detectpar = list(g0 = 0.2, sigma = 25, z = 1),
                          noccasions = NULL, binomN = NULL,
                          levels = seq(0.1, 0.9, 0.1),
                          poly = NULL, poly.habitat = TRUE, plt = TRUE, add = FALSE, fill = NULL, ...) {

    if (ms(traps)) {
        if (length(noccasions) == 1)
            noccasions <- rep(noccasions,length(traps))
        output <- mapply(pdot.contour, traps, detectpar=detectpar, noccasions=noccasions,
                         MoreArgs = list(border = border, nx = nx,
                         detectfn = detectfn, binomN = binomN,
                         levels = levels, poly = poly, poly.habitat = poly.habitat, plt = plt, add = add, ...))
        if (plt)
            invisible(output)
        else
            output
    }
    else {
        if (is.null(border))
            border <- 5 * spatialscale(detectpar, detectfn)
        tempmask <- make.mask (traps, border, nx = nx, type = 'traprect')
        xlevels <- unique(tempmask$x)
        ylevels <- unique(tempmask$y)
        binomN <- getbinomN (binomN, detector(traps))
        z <- pdot(tempmask, traps, detectfn, detectpar, noccasions, binomN)
        if (!is.null(poly)) {
            OK <- pointsInPolygon(tempmask, poly)
            if (poly.habitat)
                z[!OK] <- 0
            else
                z[OK] <- 0
        }
        if (plt) {
            contour (xlevels, ylevels, matrix(z, nrow = nx), add = add, levels = levels, ...)


            ## optional fillin 2015-05-15
            if (!is.null(fill)) {
                z[z < (0.999 * min(levels))] <- NA
                levels <- c(0,levels,1)
                .filled.contour(xlevels, ylevels,  matrix(z, nrow = nx), levels= levels,
                                col = fill)
            }


            invisible(contourLines(xlevels, ylevels, matrix(z, nrow = nx), levels = levels))
        }
        else
            contourLines(xlevels, ylevels, matrix(z, nrow = nx), levels = levels)
    }
}
############################################################################################

buffer.contour <- function (traps, buffer, nx = 64, convex = FALSE, ntheta = 100,
                            plt = TRUE, add = FALSE, poly = NULL, poly.habitat = TRUE,
                            fill = NULL, ...) {
    oneconvexbuffer <- function (buffer) {
        temp  <- data.frame(x = apply(expand.grid(traps$x, buffer * cos(theta)),1,sum),
                       y = apply(expand.grid(traps$y, buffer * sin(theta)),1,sum))
        temp <- temp[chull(temp), ]
        temp <- rbind(temp, temp[1,]) ## ensure closed
        rownames(temp) <- NULL
        if (plt)
            lines(temp,...)
        temp
    }
    if (!(inherits(traps, 'traps') | inherits(traps, 'mask')))
        stop ("requires 'traps' or 'mask' object")

    if (ms(traps)) {
        output <- lapply(traps, buffer.contour, buffer = buffer, nx = nx, convex = convex,
               ntheta = ntheta, plt = plt, add = add, poly = poly, poly.habitat = poly.habitat, ...)
        if (plt)
            invisible(output)
        else
            output
    }
    else {
        if (convex) {
            if (!is.null(poly))
                warning ("'poly' ignored when convex = TRUE")
            ## could use sf etc. to get intersection?
            theta <- (2*pi) * (1:ntheta) / ntheta
            if (!add & plt)
                plot(traps, border = buffer)
            temp <- lapply(buffer, oneconvexbuffer)
            if (plt)
                invisible (temp)
            else
                temp
        }
        else {
            tempmask <- make.mask (traps, max(buffer)*1.2, nx = nx, type = 'traprect')
            xlevels <- unique(tempmask$x)
            ylevels <- unique(tempmask$y)
            z <- distancetotrap(tempmask, traps)
            if (!is.null(poly)) {
                OK <- pointsInPolygon(tempmask, poly)
                if (poly.habitat)
                    z[!OK] <- 1e20
                else
                    z[OK] <- 0
            }
            if (plt) {
                contour (xlevels, ylevels, matrix(z, nrow = nx), add = add,
                         drawlabels = FALSE, levels = buffer,...)
                if (!is.null(fill)) {
                    #z[z < (0.999 * min(levels))] <- NA
                    #levels <- c(0,levels,1)
                    levels <- c(0,buffer)
                    .filled.contour(xlevels, ylevels,  matrix(z, nrow = nx), levels= levels,
                                    col = fill)
                }
                invisible(contourLines(xlevels, ylevels, matrix(z, nrow = nx),
                                       levels = buffer))
            }
            else
                contourLines(xlevels, ylevels, matrix(z, nrow = nx),
                             levels = buffer)
        }
    }
}
################################################################################

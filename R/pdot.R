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

    detectpars <- unlist(detectpar[parnames(detectfn)])
    if ((detectfn>9) & (detectfn<14))  detectpars <- c(detectpars, detectpar$cutval)
    if (length(detectpars)<3) detectpars <- c(detectpars,0)
    miscparm <- numeric(4);   ## dummy

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
        usge <- matrix(1, ndetector(traps), noccasions)
    }
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
      ## distmat2 <- getdistmat2 (traps, X, userdist)
      distmat2 <- getuserdist (traps, X, userdist, sessnum = NA, NULL, NULL, miscparm, detectfn == 20)
      #-------------------------------------------------------------
      pdotpointcpp(
        as.matrix(X),
        as.matrix(traps),
        as.matrix(distmat2),
        as.integer(dettype),
        as.matrix(usge),
        as.integer(markocc),
        as.integer(detectfn),
        as.double(detectpars),
        as.double(miscparm),
        as.double(truncate^2),
        as.integer(expandbinomN(binomN, dettype)),
        as.integer(grain),
        as.integer(ncores)
      )
    }
}
############################################################################################

esa.plot <- function (object, max.buffer = NULL, spacing = NULL, max.mask = NULL, detectfn,
                      detectpar, noccasions, binomN = NULL, thin = 0.1, poly = NULL,
                      poly.habitat = TRUE, session = 1,
                      plt = TRUE, type = c('density', 'esa','meanpdot', 'CVpdot'),
                      n = 1, add = FALSE, overlay = TRUE, conditional = FALSE, ...) {
    type <- match.arg(type)
    if (inherits(object, 'secr')) {
        esa.plot.secr (object, max.buffer, max.mask, thin, poly, poly.habitat, session, plt,
                       type, add, overlay, conditional, ...)
    }
    else {

        if (inherits(object, 'secrlist')) {
            output <- vector('list')
            arg <- list(max.buffer = max.buffer, max.mask = max.mask, thin = thin,
                        poly = poly, poly.habitat = poly.habitat, session = session, plt = plt, type = type,
                        add = add, conditional = conditional)
            extra <- list(...)
            if (!('col' %in% names(extra)))
                extra$col <- c("#000000", rainbow(length(object)))
            arg <- c(arg, extra)
            arg$object <- object[[1]]
            output[[1]] <- do.call( esa.plot.secr, arg)

            if (length(object)>1) {
                for (i in 2:length(object)) {
                    arg$object <- object[[i]]
                    arg$col <- extra$col[i]
                    arg$add <- TRUE
                    output[[i]] <- do.call( esa.plot.secr, arg)
                }
            if (arg$plt) {
                x1 <- par()$usr[1] + (par()$usr[2]-par()$usr[1])*0.65
                y1 <- par()$usr[3] + (par()$usr[4]-par()$usr[3])*0.95
                legend(x1, y1, legend = names(object), lty = 1, col = extra$col)
            }
}
            invisible(output)
        }
        else { if (!inherits(object, 'traps'))
            stop ("requires 'secr' or 'traps' object")
            args <- list(...)
            if(is.null(max.mask)) {
                if (is.null(spacing))
                    spacing <- spacing(object)/3
                max.mask <- make.mask (object, max.buffer, spacing,,, 'trapbuffer', poly, poly.habitat)
            }
            nmask <- nrow(max.mask)
            detectfn <- valid.detectfn(detectfn)
            binomN <- getbinomN (binomN, detector(object))   ## must now be traps object
            a <- pdot (max.mask, object, detectfn, detectpar, noccasions, binomN)
            d <- distancetotrap(max.mask, object)
            ord <- order(d,a)
            cellsize <-  attr(max.mask, 'spacing')^2/10000
            a <- a[ord]

            ## CV 2018-06-01
            mu <- cumsum(a) / (1:nmask)
            ## 2021-05-19 pmax protects against negative argument to sqrt          
            cv <- sqrt(pmax(0,cumsum(a^2)/(1:nmask) - mu^2))/mu
            cumcv <- function(n) {
                an <- a[1:n]
                fx <- an/sum(an)
                mucond <- sum(an * fx)
                cvcond <- sqrt(sum(an^2 * fx) - mucond^2)/mucond
                c(mucond, cvcond)
            }
            ## debug check
            ## tmp <- CVpdot(max.mask, object, detectfn=detectfn, detectpar=detectpar,
            ##   noccasions = noccasions, conditional = TRUE)
            ###################################################
            output <- data.frame(buffer = d[ord], esa =  cumsum(a) * cellsize,
                                 density = n /  cumsum(a) / cellsize,
                                 pdot = a, pdotmin = cummin(a),
                                 meanpdot = mu, CVpdot = cv)

            maxesa <- max(output$esa)
            thinned <- seq(1,  nmask, 1/thin)
            output <- output[thinned,]

            if (conditional) {
                cv <- sapply(thinned, cumcv)
                output$meanpdot <- cv[1,]
                output$CVpdot <- cv[2,]
            }

            if (plt) {
                if (type == 'density') {
                    if (add)
                        lines(output$buffer, n/output$esa, ...)
                    else {
                        xlb <- 'Buffer width  m'
                        ylb <- expression(paste('n / esa(buffer)   ',ha^-1))
                        if ('ylim' %in% names(args))
                            plot(output$buffer, n/output$esa, type = 'l',
                                 xlab = xlb, ylab = ylb, ...)
                        else  ## clunky!
                            plot(output$buffer, n/output$esa, type = 'l',
                                 xlab = xlb, ylab = ylb, ...,
                                 ylim= n / maxesa * c(0.9, 1.2))
                    }
                }
                else if (type == 'esa') {
                    if (add)
                        lines(output$buffer, output$esa, ...)
                    else
                        plot(output$buffer, output$esa, type = 'l',
                             xlab = 'Buffer width  m', ylab = 'esa(buffer)  ha', ...)
                }
                else if (type == 'meanpdot') {
                    if (add)
                        lines(output$buffer, output$meanpdot, ...)
                    else
                        plot(output$buffer, output$meanpdot, type = 'l',
                             xlab = 'Buffer width  m', ylab = 'meanpdot(buffer)', ...)
                }
                else if (type == 'CVpdot') {
                    if (add)
                        lines(output$buffer, output$CVpdot, ...)
                    else
                        plot(output$buffer, output$CVpdot, type = 'l',
                             xlab = 'Buffer width  m', ylab = 'CVpdot(buffer)', ...)
                }

                invisible(output)
            }
            else output
        }
    }
}

###############################################################################

esa.plot.secr <- function (object, max.buffer = NULL, max.mask = NULL,
    thin = 0.1, poly = NULL, poly.habitat = TRUE, session = 1, plt = TRUE, type = 'density',
    add = FALSE, overlay = TRUE, conditional = FALSE, ...) {

    if (!inherits(object,'secr'))
        stop("require secr object")
    MS <- ms(object)
    if (MS) {
        sessnames <- session(object$capthist)
        ## use alphanumeric session ID
        if (is.numeric(session))
            session <- sessnames[session]
    }

    ## recursive call
    if (MS & (length(session) > 1)) {
        esa.plot.outputs <- vector(mode='list')

        for (i in session) {
            addthisone <- ifelse (add | (overlay & (i != session[1])),
                                  TRUE, FALSE)
            esa.plot.outputs[[i]] <- esa.plot.secr (object, max.buffer,
                max.mask, thin, poly, poly.habitat, i, plt, type, addthisone,
                overlay, conditional, ...)
        }
        if (plt)
            invisible(esa.plot.outputs)
        else
            esa.plot.outputs
    }
    ## not recursive
    else {
        if (MS) {
            ## select one session
            trps <- traps(object$capthist[[session]])
            n <- nrow(object$capthist[[session]])
            nocc <- ncol(object$capthist[[session]])
            spacg <- attr(object$mask[[session]], 'spacing')
            detpar <- detectpar(object)[[session]]
            spscale <- spatialscale(object, object$detectfn, session)
        }
        else {
            trps <- traps(object$capthist)
            n <- nrow(object$capthist)
            nocc <- ncol(object$capthist)
            spacg <- attr(object$mask, 'spacing')
            detpar <- detectpar(object)
            spscale <- spatialscale(object, object$detectfn)
        }
        if (is.null(max.mask)) {
            if (is.null(max.buffer)) {
                if (add)
                    max.buffer <- par()$usr[2]  ## width of existing plot
                else {
                    max.buffer <- 5 * spscale
                }
            }
        }
        binomN <- object$details$binomN
        esa.plot (trps, max.buffer, spacg, max.mask, object$detectfn, detpar,
                  nocc, binomN, thin, poly, poly.habitat, session, plt, type, n, add, overlay,
                  conditional, ...)
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

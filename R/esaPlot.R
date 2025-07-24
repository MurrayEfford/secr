###############################################################################
## package 'secr'
## esaPlot.R
## 2022-11-19 separated from pdot.R
## 2024-09-28 renamed
###############################################################################

esa.plot <- function (object, max.buffer = NULL, spacing = NULL, max.mask = NULL, detectfn,
                      detectpar, noccasions, binomN = NULL, thin = 0.1, poly = NULL,
                      poly.habitat = TRUE, session = 1,
                      plt = TRUE, type = c('density', 'esa','meanpdot', 'CVpdot'),
                      n = 1, add = FALSE, overlay = TRUE, conditional = FALSE, ...) {

    .Deprecated("esaPlot", package="secr",
                "esa.plot has been renamed esaPlot; the old version will be removed soon",
                old = as.character(sys.call(sys.parent()))[1L])

    esaPlot (object, max.buffer, spacing, max.mask, detectfn, detectpar,
             noccasions, binomN, thin, poly, poly.habitat, session,plt,
             type, n, add, overlay, conditional, ...)
}

esaPlot <- function (object, max.buffer = NULL, spacing = NULL, max.mask = NULL, detectfn,
    detectpar, noccasions, binomN = NULL, thin = 0.1, poly = NULL,
    poly.habitat = TRUE, session = 1,
    plt = TRUE, type = c('density', 'esa','meanpdot', 'CVpdot'),
    n = 1, add = FALSE, overlay = TRUE, conditional = FALSE, ...) {
    type <- match.arg(type)
    if (inherits(object, 'secr')) {
        esaPlotsecr (object, max.buffer, max.mask, thin, poly, poly.habitat, session, plt,
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
            output[[1]] <- do.call( esaPlotsecr, arg)
            
            if (length(object)>1) {
                for (i in 2:length(object)) {
                    arg$object <- object[[i]]
                    arg$col <- extra$col[i]
                    arg$add <- TRUE
                    output[[i]] <- do.call( esaPlotsecr, arg)
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
            detectfn <- secr_valid.detectfn(detectfn)
            binomN <- secr_getbinomN (binomN, detector(object))   ## must now be traps object
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

esaPlotsecr <- function (object, max.buffer = NULL, max.mask = NULL,
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
        esaPlotOutputs <- vector(mode='list')
        
        for (i in session) {
            addthisone <- ifelse (add | (overlay & (i != session[1])),
                TRUE, FALSE)
            esaPlotOutputs[[i]] <- esaPlotsecr (object, max.buffer,
                max.mask, thin, poly, poly.habitat, i, plt, type, addthisone,
                overlay, conditional, ...)
        }
        if (plt)
            invisible(esaPlotOutputs)
        else
            esaPlotOutputs
    }
    ## not recursive
    else {
        if (MS) {
            ## select one session
            trps <- traps(object$capthist[[session]])
            n <- nrow(object$capthist[[session]])
            nocc <- ncol(object$capthist[[session]])
            spacg <- attr(object$mask[[session]], 'spacing')
            if (!is.numeric(session)) {
                session <- match(session, names(object$capthist))
            }
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
        esaPlot (trps, max.buffer, spacg, max.mask, object$detectfn, detpar,
            nocc, binomN, thin, poly, poly.habitat, session, plt, type, n, add, overlay,
            conditional, ...)
    }
}
#################################################################################

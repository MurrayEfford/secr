#############################################################################
## package 'secr'
## summary.mask.R
## 2022-08-23 moved from methods.R 
#############################################################################

summary.mask <- function(object, ...) {
    
    if (ms(object)) {
        temp <- lapply(object, summary.mask)
        class(temp) <- c('summary.mask', 'list')
        temp
    }
    else {
        if (is.null(object$x) | is.null(object$y))
            stop ("not a valid mask")
        nd <- length(object$x)
        if (length(object$x) != length(object$y))
            stop  ("not a valid mask")
        
        if (!is.null(covariates(object))) {
            sumcovar <- summary(covariates(object), ...)
        } else sumcovar <- NULL
        if (inherits(object, 'linearmask'))
            maskclass <- 'linearmask'
        else if (inherits(object,'Dsurface'))
            maskclass <- 'Dsurface'
        else if (inherits(object,'Rsurface'))
            maskclass <- 'Rsurface'
        else
            maskclass <- 'mask'
        
        ## rearranged 2014-09-06
        temp <- list (
            maskclass = maskclass,
            masktype = attr(object, 'type',exact = TRUE),
            nmaskpoints = nrow(object),
            xrange = range(object$x),
            yrange = range(object$y),
            meanSD = attr(object, 'meanSD',exact = TRUE),
            spacing = attr(object, 'spacing',exact = TRUE),
            cellarea = attr(object, 'area',exact = TRUE),
            boundingbox = attr(object, 'boundingbox',exact = TRUE),
            covar = sumcovar
        )
        class(temp) <- 'summary.mask'
        temp
    }
    
}
############################################################################################


print.summary.mask <- function (x, ...) {
    if (ms(x)) {
        lapply (x, print.summary.mask)
    }
    else {
        cat ('Object class     ', x$maskclass, '\n')
        cat ('Mask type        ', x$masktype, '\n')
        cat ('Number of points ', x$nmaskpoints, '\n')
        cat ('Spacing m        ', x$spacing, '\n')
        if (is.null(x$cellarea))
            cat ('Total length km  ', x$spacing * x$nmaskpoints / 1000, '\n')
        else {
            cat ('Cell area ha     ', x$cellarea, '\n')
            cat ('Total area ha    ', x$cellarea * x$nmaskpoints, '\n')
        }
        cat ('x-range m        ', x$xrange, '\n')
        cat ('y-range m        ', x$yrange, '\n')
        cat ('Bounding box     ','\n')
        print (x$boundingbox, ...)
        cat ('\n')
        if (!is.null(x$covar)) {
            cat ('Summary of covariates', '\n')
            print(x$covar, ...)
        }
    }
}
############################################################################################

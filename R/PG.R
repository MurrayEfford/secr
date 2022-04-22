##############################################################################
## package 'secr'
## PG.R
## 2013-11-23, 2017-12-01, 2022-02-13 (sf)
##############################################################################

PG <- function (CH, poly = NULL, includeNULL = FALSE, plt = FALSE, ...) {
    if (is.null(poly)) {
        poly <- buffer.contour (traps(CH), plt = plt, ...)
        clean <- function(x) {
            ## x may be list or data frame, depending on concave or convex
            if (is.data.frame(x))
                as.matrix(x)
            else  ## assume list
                cbind(x$x, x$y)
        }
        poly <- lapply(poly, clean)
        poly <- st_sfc(st_polygon(poly))
 
    }
    inpoly <- function (xy) {
        if (is.null(xy)) {
            xy <- matrix(0,nrow=1,ncol=2)
            ip <- NA
        }
        else if (nrow(xy) == 0) {
            ip <- NA
        }
        else {
            ip <- pointsInPolygon(xy, poly)
            if (plt) points (xy[,1], xy[,2], pch=c(1,16)[ip+1])
        }
        sum(ip)/nrow(xy)
    }
    if (any(detector(traps(CH)) %in% c('polygon','polygonX', 'transect','transectX'))) {
        xyl <- split(xy(CH), factor(animalID(CH, sortorder = 'ksn'), levels=rownames(CH) )) 
    }
    else {
        xyl <- telemetryxy(CH, includeNULL=includeNULL)
    }
    sapply(xyl, inpoly)
}

##############################################################################

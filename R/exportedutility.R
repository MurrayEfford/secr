################################################################################
## package 'secr'
## exportedutility.R
## 2024-09-16  moved from utility.R
################################################################################

## mean and SD if x numeric
getMeanSD <- function(xy) {
    MeanSD <- function (x) {
        if (is.numeric(x))
            c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
        else
            c(NA,NA)
    }
    as.data.frame (apply(xy, 2, MeanSD))
}

#-------------------------------------------------------------------------------

maskarea <- function (mask, sessnum = 1) {
    if (!ms(mask)) nrow(mask) * attr(mask,'area')
    else nrow(mask[[sessnum]]) * attr(mask[[sessnum]],'area')
}

#-------------------------------------------------------------------------------

masklength <- function (mask, sessnum = 1) {
    if (!ms(mask)) nrow(mask) * attr(mask,'spacing')/1000
    else nrow(mask[[sessnum]]) * attr(mask[[sessnum]],'spacing')/1000
}

#-------------------------------------------------------------------------------

masksize <- function (mask, sessnum = 1) {
    if (inherits(mask, 'linearmask'))
        masklength(mask, sessnum)
    else
        maskarea(mask, sessnum)
}

#-------------------------------------------------------------------------------
edist <- function (xy1, xy2) {
    nr <- nrow(xy1)
    nc <- nrow(xy2)
    x1 <- matrix(xy1[,1], nr, nc)
    x2 <- matrix(xy2[,1], nr, nc, byrow=T)
    y1 <- matrix(xy1[,2], nr, nc)
    y2 <- matrix(xy2[,2], nr, nc, byrow=T)
    sqrt((x1-x2)^2 + (y1-y2)^2)
}

#-------------------------------------------------------------------------------

## least cost paths from mask including barriers to movement
## use edist for equivalent Euclidean distances

nedist <- function (xy1, xy2, mask, inf = Inf, ...) {
    newargs <- list(...)
    if (missing(mask)) mask <- xy2
    noneuc <- covariates(mask)$noneuc
    if (is.null(noneuc)) noneuc <- rep(1, nrow(mask))
    defaultargs <- list(transitionFunction = mean, directions = 16)
    args <- replace(defaultargs, names(newargs), newargs)
    args$x <- raster(mask, values = noneuc)
    if (requireNamespace('gdistance', quietly = TRUE)) {    ## 2015-01-23
        tr <- do.call(gdistance::transition, args)
        tr <- gdistance::geoCorrection(tr, type = "c", multpl = FALSE)
        out <- gdistance::costDistance(tr, as.matrix(xy1), as.matrix(xy2))
    }
    else stop ("package gdistance is required for nedist")
    if (is.finite(inf)) out[!is.finite(out)] <- inf
    out
}

#-------------------------------------------------------------------------------

rlnormCV <- function(n, mean, cv) {
    # n simulated values from log-normal distribution with mean = mean and CV = cv
    sdlog <- log(cv^2 + 1)^0.5
    meanlog <- log(mean) - sdlog^2/2
    rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}
#-------------------------------------------------------------------------------

shareFactorLevels <- function (object, columns = NULL, stringsAsFactors = TRUE) {
    ## stringsAsFactors added 2020-05-16
    if (ms(object)) {
        if (!is.null(covariates(object))) {
            df <- do.call(rbind, covariates(object))
            if (is.null(columns)) {
                columns <- 1:ncol(df)
            }
            if (stringsAsFactors) {
                df[,columns] <- secr_stringsAsFactors(df[,columns, drop = FALSE])
            }
            for (i in columns) {
                if (is.factor(df[,i])) {
                    levelsi <- levels(df[,i])
                    for (sess in 1:length(object)) {
                        covariates(object[[sess]])[,i] <-
                            factor(covariates(object[[sess]])[,i],
                                   levels = levelsi)
                    }
                }
            }
        }
    }
    else {
        # modified 2021-04-27 to apply to covariates, not object itself
        if (!is.null(covariates(object))) {
            if (stringsAsFactors) {
                df <- covariates(object)
                if (is.null(columns)) {
                    columns <- 1:ncol(df)
                }
                df[,columns] <- secr_stringsAsFactors(df[,columns, drop = FALSE])
                covariates(object) <- df
            }
        }
    }
    object
}

#-------------------------------------------------------------------------------

# make.lookup is exported
make.lookup <- function (tempmat) {
    
    ## should add something to protect make.lookup from bad data...
    nrw <- nrow(tempmat)
    ncl <- ncol(tempmat)
    nam <- colnames(tempmat)
    
    df <- is.data.frame(tempmat)
    if (df) {
        lev <- lapply(tempmat, levels)
        tempmat[] <- sapply(tempmat, as.numeric)
        tempmat <- as.matrix(tempmat)
    }
    dimnames(tempmat) <- NULL
    
    temp <- makelookupcpp(tempmat)
    
    lookup <- temp$lookup
    colnames(lookup) <- nam
    if (df) {
        lookup <- as.data.frame(lookup)
        ## restore factors
        for (i in 1: length(lev))
            if (!is.null(lev[[i]]))
                lookup[,i] <- factor(lev[[i]][lookup[,i]], levels = lev[[i]])
    }
    list (lookup=lookup, index=temp$index)
}

#-------------------------------------------------------------------------------

# insertdim is exported
insertdim <- function (x, dimx, dims) {
    ## make vector of values
    ## using x repeated so as to fill array
    ## with dim = dims and the x values occupying dimension(s) dimx
    olddim <- 1:length(dims)
    olddim <- c(olddim[dimx], olddim[-dimx])
    temp <- array (dim=c(dims[dimx], dims[-dimx]))
    tempval <- array(dim=dims[dimx])
    if (length(x) > length(tempval))
        tempval[] <- x[1:length(tempval)]
    else
        tempval[] <- x     ## repeat as needed
    temp[] <- tempval  ## repeat as needed
    if (is.factor(x))
        factor(levels(x), levels=levels(x))[aperm(temp, order(olddim))]   ## 2010 02 25
    else
        as.vector(aperm(temp, order(olddim)))
}
#-------------------------------------------------------------------------------

boundarytoSF <- function (poly) {
    if (is.null(poly)) {
        NULL
    }
    else if(inherits(poly, c('sf','sfc'))) {
        poly <- st_geometry(poly) # extract sfc if not already sfc
        geomtype <- st_geometry_type(poly, by_geometry = FALSE)
        if (geomtype == 'GEOMETRY') {   # 2023-06-02
            geomtype <- st_geometry_type(poly, by_geometry = TRUE)
        }
        if (!all(geomtype %in% c("POLYGON", "MULTIPOLYGON"))) {
            stop ("poly sf/sfc should be of type POLYGON or MULTIPOLYGON")
        }
        poly
    }
    else if (inherits(poly, 'SpatialPolygons')) {   # also SPDF?
        st_as_sfc(poly)
    }
    else if (inherits(poly, 'SpatVector')) {
        st_as_sfc(as(poly,"Spatial"))
    }
    else if (inherits(poly, c('matrix', 'data.frame'))) {
        ## input is 2-column matrix for a single polygon
        poly <- matrix(unlist(poly), ncol = 2)
        poly <- rbind (poly, poly[1,])  ## force closure of polygon
        st_sfc(st_polygon(list(poly)))
    }
    else stop (class(poly), " not valid input to boundarytoSF")
}

#-------------------------------------------------------------------------------

pointsInPolygon <- function (xy, poly, logical = TRUE) {
    # xy is 2-column matrix or data.frame of coordinates
    if (inherits(poly, 'mask')) { 
        if (ms(poly))
            stop ("multi-session masks not supported")
        sp <- spacing(poly)
        minx <- min(poly$x, na.rm = TRUE)
        miny <- min(poly$y, na.rm = TRUE)
        mask <- sweep(poly, MARGIN = 2, FUN = '+', STATS = c(-minx, -miny))
        mask <- round(mask/sp) + 1
        xy <- matrix(unlist(xy), ncol = 2)  ## in case dataframe
        xy <- sweep(xy, MARGIN = 2, FUN = '+', STATS = c(-minx, -miny))
        xy <- round(xy/sp) + 1
        xy[xy<=0] <- NA
        xy[,1][xy[,1]>max(mask$x, na.rm = TRUE)] <- NA
        xy[,2][xy[,2]>max(mask$y, na.rm = TRUE)] <- NA
        
        maskmatrix <- matrix(0, ncol = max(mask$y, na.rm = TRUE), nrow = max(mask$x, na.rm = TRUE))
        maskmatrix[as.matrix(mask)] <- 1:nrow(mask)
        inside <- maskmatrix[as.matrix(xy)]
        inside[is.na(inside)] <- 0
        if (logical)
            inside <- inside > 0
        inside
    }
    else {
        poly <- boundarytoSF(poly)
        if (inherits(poly, c('sf','sfc'))) {
            xy <- st_as_sf(data.frame(xy), coords = 1:2)
            st_crs(xy) <- st_crs(poly)
            apply(st_within(xy, poly, sparse = FALSE), 1, any)
        }
        else {
            stop ("unknown input to pointsInPolygon")
        }
    }
}

#-------------------------------------------------------------------------------

shareFactorLevels <- function (object, columns = NULL, stringsAsFactors = TRUE) {
    ## stringsAsFactors added 2020-05-16
    if (ms(object)) {
        if (!is.null(covariates(object))) {
            df <- do.call(rbind, covariates(object))
            if (is.null(columns)) {
                columns <- 1:ncol(df)
            }
            if (stringsAsFactors) {
                df[,columns] <- secr_stringsAsFactors(df[,columns, drop = FALSE])
            }
            for (i in columns) {
                if (is.factor(df[,i])) {
                    levelsi <- levels(df[,i])
                    for (sess in 1:length(object)) {
                        covariates(object[[sess]])[,i] <-
                            factor(covariates(object[[sess]])[,i],
                                   levels = levelsi)
                    }
                }
            }
        }
    }
    else {
        # modified 2021-04-27 to apply to covariates, not object itself
        if (!is.null(covariates(object))) {
            if (stringsAsFactors) {
                df <- covariates(object)
                if (is.null(columns)) {
                    columns <- 1:ncol(df)
                }
                df[,columns] <- secr_stringsAsFactors(df[,columns, drop = FALSE])
                covariates(object) <- df
            }
        }
    }
    object
}

#-------------------------------------------------------------------------------

## Manually remove some mask points
# simplified 2022-02-03

deleteMaskPoints <- function (mask, onebyone = TRUE, add = FALSE, poly = NULL,
                              poly.habitat = FALSE, ...) {
    ## interface does not work properly in RStudio
    
    if (ms(mask)) {         ## a list of mask objects
        if (inherits(poly, 'list') & (!is.data.frame(poly)))
            stop ("lists of polygons not implemented in 'make.mask'")
        temp <- lapply (mask, deleteMaskPoints, onebyone = onebyone, add = add,
                        poly = poly, poly.habitat = poly.habitat, ...)
        class (temp) <- c('mask', 'list')
        temp
    }
    else {
        plot(mask, add = add, ...)
        if (!is.null(poly)) {
            if (poly.habitat)
                pointstodrop <- (1:nrow(mask))[!pointsInPolygon(mask, poly)]
            else
                pointstodrop <- (1:nrow(mask))[pointsInPolygon(mask, poly)]
        }
        else if (onebyone) {
            cat ('Click to select points; right-click to stop\n')
            flush.console()
            xy <- locator(type = 'p', pch=1, col='red')
            pointstodrop <- if (length(xy$x)==0)
                numeric(0)
            else
                nearesttrap(xy, mask)
        }
        else {
            cat ('Click to select polygon vertices; right-click to stop\n')
            flush.console()
            xy <- locator(type = 'l', col='red')
            xy <- as.data.frame(xy)
            xy <- rbind(xy, xy[1,])
            if (poly.habitat)
                pointstodrop <- (1:nrow(mask))[!pointsInPolygon(mask, xy)]
            else
                pointstodrop <- (1:nrow(mask))[pointsInPolygon(mask, xy)]
        }
        npts <- length(pointstodrop)
        if (npts>0) {
            points(mask[pointstodrop,], pch = 16, col = 'red')
            if(.Platform$OS.type == "windows") {
                pl <- if (npts>1) 's' else ''
                msg <- paste ('Delete ', npts, ' red point',pl, '?', sep='')
                response <-  utils::winDialog(type = "okcancel", msg)
            } else {
                response <- 'OK'
            }
            if (response == 'OK') {
                mask <- subset(mask, -pointstodrop)
                if (npts==1)
                    message("1 point deleted")
                else
                    message(npts, " points deleted")
            }
            else
                message ("point(s) not deleted")
        }
        else
            message ("no points to delete")
        plot(mask, col='green')
        mask
    }
}

#-------------------------------------------------------------------------------
updateCH <- function(object) {
    if (!inherits(object, 'capthist'))
        stop ("requires capthist object")
    if (ms(object)) {
        out <- lapply(object, updateCH)
        class (out) <- c("capthist", "list")
        out
    }
    else {
        if (length(dim(object)) == 3) {
            return(object)
        }
        else {
            K <- secr_ndetector(traps(object))
            ch <- array(0, dim = c(dim(object), K), dimnames = 
                            list(rownames(object), colnames(object), 1:K))
            OK <- as.logical(object!=0)
            animal <- row(object)[OK]
            occ <- col(object)[OK] 
            detn <- object[OK]
            ch[cbind(animal, occ, detn)] <- 1
            traps(ch) <- traps(object)
            class (ch) <- "capthist"
            session(ch) <- session(object)
            ch
        }
    }
}
#-------------------------------------------------------------------------------

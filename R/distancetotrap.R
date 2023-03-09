################################################################################
## package 'secr'
## distancetotrap.R
## 2023-03-10 distancetotrap and nearesttrap moved from utility.R
################################################################################

distancetopoly <- function (X, traps) {
    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coord in col 2
    
    if (is.null(X)) return (NULL) 
    
    detecttype <- detector(traps)
    detecttype <- ifelse (is.null(detecttype), "", detecttype)
    if (!all(detecttype %in% c('polygon', 'polygonX')))
        stop("distancetopoly is for polygon detectors only")
    
    xy <- st_as_sf(data.frame(X), coords=1:2)
    trps <- split(traps, polyID(traps))
    polys <- lapply(trps, boundarytoSF)
    
    dlist <- lapply(polys, st_distance, x = xy)
    matrix(unlist(dlist), ncol = length(dlist))
}

#-------------------------------------------------------------------------------

distancetotrap <- function (X, traps) {
    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coord in col 2
    
    if (is.null(X)) return (NULL)  ## 2022-01-04
    
    X <- as.data.frame(X)
    xy <- st_as_sf(X, coords = 1:2)     # POINTS
    
    nxy <- nrow(X)
    detecttype <- detector(traps)
    detecttype <- ifelse (is.null(detecttype), "", detecttype)
    
    if (all(detecttype %in% c('polygon', 'polygonX'))) {
        trps <- split(traps, polyID(traps))
        polys <- lapply(trps, boundarytoSF)
        dlist <- lapply(polys, st_distance, x = xy)
        dmat <- matrix(unlist(dlist), ncol = length(dlist))
        d <- apply(dmat,1,min)
        return (d)
    }
    
    if (inherits(traps, 'SpatialPolygons')) {
        traps <- st_as_sf(traps)
        d <- st_distance(xy, traps)
        return (d)
    }
    else if (all(detecttype %in% .localstuff$polydetectors)) {
        ## approximate only
        
        traps <- split(traps, polyID(traps))
        trpi <- function (i, n = 100) {
            intrp <- function (j) {
                ## 2020-01-08 dodge issue with polyID in as.data.frame
                ## tmp <- as.data.frame(traps[[i]][j:(j+1),])[,-1]   
                tmp <- data.frame(x = traps[[i]]$x[j:(j+1)], y = traps[[i]]$y[j:(j+1)])
                if (tmp$x[1] == tmp$x[2])
                    data.frame(x=rep(tmp$x[1], n),
                        y=seq(tmp$y[1], tmp$y[2], length=n))
                else {
                    ## 2019-11-30 suppress warnings such as :
                    ## In regularize.values(x, y, ties, missing(ties)) :
                    ## collapsing to unique 'x' values
                    suppressWarnings(data.frame(approx(tmp, n = n)))
                }
            }
            tmp <- lapply(1:(nrow(traps[[i]])-1),intrp)
            do.call(rbind, tmp)
        }
        trps <- do.call(rbind, lapply(1:length(traps), trpi))
        trps <- matrix(unlist(trps), ncol = 2)
    }
    else {
        ## 2015-10-18 added protection
        trps <- matrix(unlist(traps), ncol = 2)
    }
    
    temp <- nearestcpp(as.matrix(X), as.matrix(trps))
    if (all(detecttype %in% c('polygon', 'polygonX'))) {
        inside <- lapply(traps, pointsInPolygon, xy = X)
        inside <- do.call(rbind, inside)
        temp$distance [apply(inside,2,any)] <- 0
    }
    temp$distance
}

#-------------------------------------------------------------------------------

nearesttrap <- function (X, traps) {
    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coord in col 2
    
    if (is.null(X)) return (NULL)  ## 2022-01-04
    
    X <- matrix(unlist(X), ncol = 2)
    nxy <- nrow(X)
    if (inherits(traps, 'SpatialPolygons')) {
        stop ("nearesttrap currently does not accept SpatialPolygons (from 4.5.3)")
        # traps <- sp::coordinates(traps@polygons[[1]]@Polygons[[1]])
        # warning("using only first polygon of SpatialPolygons")
    }
    temp <- nearestcpp(as.matrix(X), as.matrix(traps))
    temp$index
}
#-------------------------------------------------------------------------------

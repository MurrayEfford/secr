###############################################################################
## package 'secr'
## trap.builder.R
## repeat trap layout systematically, by GRTS, or at simple random centres
## across a region
## Also make.systematic, mash(), cluster.counts(), cluster.centres()
## 2011-08-16 (full argument names)
## 2012-01-11 (cutval and signal in mash(); mxy[1:2])
## 2012-02-07 mash noise
## 2012-07-26 mash cleans out trap attributes; names sessions for list input
## 2012-09-17 mash moved to mash.r
## 2013-03-02 exclude, exclmethod arguments added
## 2013-04-20 replace deprecated overlay with over
## 2014-10-25 revamp polygon requirement for grts
## 2014-12-11 set proj4string to NA
## 2014-12-11 method = "GRTS" changed to method == "GRTS" !
## 2018-09-27 coerce CRS of region to CRS()
## 2018-12-18 clipping left to trap.builder()
## 2018-12-19 make.systematic() grid of centres expanded
## 2018-12-29 boundarytoSPDF and boundarytoSP moved to utility.R
## 2018-12-29 make.systematic() argument 'chequerboard'
## 2019-12-20 make.systematic() arguments 'rotate', 'centrexy'
## 2020-01-27 fix bug (saved n = NULL when n missing)
## 2021-10-18 grts temporarily suspended
## 2022-02-05 comprehensive use of sf
## 2022-02-19 grts revived; use sf
###############################################################################

## spsurvey uses sf 2022-01-31

trap.builder <- function (n = 10, cluster, region = NULL, frame =
    NULL, method = c("SRS", "GRTS", "all", "rank"), edgemethod =
    c("clip", "allowoverlap", "allinside", "anyinside", "centreinside"), 
    samplefactor = 2, ranks = NULL, rotation = NULL, detector, 
    exclude = NULL, exclmethod = c("clip", "alloutside", "anyoutside", "centreoutside"), 
    plt = FALSE, add = FALSE, ...) {
    ## region may be -
    ## matrix x,y
    ## sf
    ## sp SpatialPolygonsDataFrame object

    #####################################################
    # 1. form polygon
    # 2. get n random origins SRS, GRTS
    # 3. translate n times; covariate for cluster of origin
    # 4. optionally clip reject edge clusters
    # 5. rbind
    #####################################################
    
    allinside <- function (xy) {
        xy <- st_as_sf(data.frame(xy), coords=1:2)
        st_crs(xy) <- st_crs(region)
        all(st_within(xy, region, sparse = FALSE))
    }

    anyinside <- function (xy) {
        xy <- st_as_sf(data.frame(xy), coords=1:2)
        st_crs(xy) <- st_crs(region)
        any(st_within(xy, region, sparse = FALSE))
    }

    alloutside <- function (xy) {
        xy <- st_as_sf(data.frame(xy), coords=1:2)
        st_crs(xy) <- st_crs(exclude)
        all(!st_within(xy, exclude, sparse = FALSE))
    }

    anyoutside <- function (xy) {
        xy <- st_as_sf(data.frame(xy), coords=1:2)
        st_crs(xy) <- st_crs(exclude)
        any(!st_within(xy, exclude, sparse = FALSE))
    }

    centreinside <- function (xy) {
        xy <- apply(as.matrix(xy),2,mean)
        xy <- matrix(xy, ncol = 2)  # 2022-05-27
        xy <- st_as_sf(data.frame(xy), coords = 1:2)
        st_crs(xy) <- st_crs(region)
        OK <- st_within(xy, region, sparse = FALSE)
        apply(OK,1,any)
    }
    
    centreoutside <- function (xy) {
        xy <- apply(as.matrix(xy),2,mean)
        xy <- matrix(xy, ncol = 2)  # 2022-05-27
        xy <- st_as_sf(data.frame(xy), coords = 1:2)
        st_crs(xy) <- st_crs(region)
        OK <- st_within(xy, region, sparse = FALSE)
        !apply(OK,1,any)
    }

    position <- function (i, cluster) {
        newtraps <- secr::shift(cluster, origins[i,])
        i <- .local$clusteri
        .local$clusteri <- .local$clusteri + 1    ## global update!
        clusterID(newtraps) <- factor(rep(i,nrow(newtraps)), levels=i)
        clustertrap(newtraps) <- as.numeric(polyID(newtraps))
        if (!is.null(rotation)) {
          if (rotation<0)
            rotation <- runif(1) * 360
          newtraps <- secr::rotate(newtraps, rotation, apply(newtraps,2,mean))
        }
        newtraps
    }

    ## mainline
    .local <- new.env()   ## for clusteri
    method <- match.arg(method)
    edgemethod <- match.arg(edgemethod)
    exclmethod <- match.arg(exclmethod)

    if (method == "GRTS" && !is.null(exclude)) {
        stop ("GRTS incompatible with non-null 'exclude'")
    }    
    ## option for single-trap clusters
    if (missing(cluster))
        cluster <- NULL
    if (is.null(cluster)) {
        if (missing(detector))
            detector <- 'multi'
        if (!(detector %in% .localstuff$pointdetectors))
            stop ("solitary detectors must be of a point detector type")
        cluster <- make.grid(nx = 1, ny = 1, detector = detector)
        ## 2019-09-25 block this 
        ## edgemethod <- 'allowoverlap'
    }
    else {
        if ((attr(cluster,'detector') %in% .localstuff$polydetectors) &
            (ndetector(cluster) > 1))
            stop("clusters with multiple polygons or transects not supported")
    }

    if (method == 'all') {
        n <- nrow(frame)
    }
    if (is.null(region)) {
        edgemethod <- 'allowoverlap'
    }

    region <- boundarytoSF(region)   # sfc
    exclude <- boundarytoSF(exclude) # sfc
    
    if (is.null(frame)) {
        if (is.null(region)) {
            stop ("specify at least one of 'region' or 'frame'")
        }

        if (plt & !add) {
            plot(st_geometry(region))   # use plot method for sfc
            if (!is.null(exclude)) {
                plot(st_geometry(exclude), col='lightgrey', add=TRUE, border='lightgrey')
            }
        }
    }
    else {
        if (plt & !add) {
            if (!is.null(region)) {
                plot(st_geometry(region))   # use plot method for sfc
            }
            else {
                MASS::eqscplot (frame, axes = F, xlab = '', ylab = '', pch = 1, cex = 0.5)
            }
            if (!is.null(exclude)) {
                plot(st_geometry(region), border = 'lightgrey', add = TRUE)   # use plot method for sfc
            }
        }
    }

    # define oversample
    ntrial <- max(n * samplefactor, 5)
    
    ####################################
    if (method == 'SRS') {
        
        if (is.null(frame)) {
            pts <- st_sample(region, size = ntrial, type  = "random", exact = TRUE)
            origins <- st_coordinates(pts)
        }
        else {
            if (ntrial > nrow(frame))
                stop ("too few rows in frame for requested sample")
            OK <- sample.int(nrow(frame), ntrial, replace = FALSE)
            origins <- as.matrix(frame[OK, ])
        }
    }
    ####################################
    else if (method == 'GRTS') {
        
        if (!requireNamespace("spsurvey", versionCheck=list(op=NULL, version = ">=5.3.0"), quietly = TRUE)) {
            stop ("package 'spsurvey >= 5.3.0' required for GRTS in trap.builder")
        }
        if (!is.null(region)) {
                sf_frame <- st_sf(region)   # sfc to sf
        }
        else {
            # assume frame of points
            sf_frame <- st_as_sf(data.frame(frame), coords=1:2)
        }

        # ensure valid crs
        crs <- st_crs(sf_frame)
        if (!is.na(crs) && crs$IsGeographic) {   # most likely EPSG 4326
            stop ("region should use projected (Cartesian) coordinates")
        }
        ntrial <- n                   # oversample not allowed

        GRTS.sites <- spsurvey::grts(
            sframe      = sf_frame,
            n_base      = ntrial,
            stratum_var = NULL,        # unstratified
            seltype     = "equal",     # unweighted
            projcrs_check = FALSE,     # override check
            ...
        )
        origins <- st_coordinates(GRTS.sites$sites_base)

    }
    
    ####################################
    else if (method == 'all') {
        if (is.null(frame) || (nrow(frame)==0)) {
            stop ("method = 'all' requires finite frame")
        }
        origins <- as.matrix(frame)
    }
    ####################################
    else if (method == 'rank') {
        if (is.null(frame)) {
            stop ("'rank' requires finite frame")
        }
        if (is.null(ranks)) {
            stop ("'rank' requires ranks")
        }
        nframe <- nrow(frame)
        if (nframe<n)
            stop ("not enough rows in frame for requested n")
        ranks <- ranks + runif(nframe)/(nframe+1)
        frameorder <- rev((1:nframe)[order(ranks)])
        frame <- frame[frameorder,]
        origins <- as.matrix(frame)
    }
    else {
        stop ("method not recognised")
    }
    rownames(origins) <- 1:nrow(origins)
    
    #######################################################
    ## centre cluster on (0,0)
    if (nrow(cluster)>1)
        cxy <- apply(cluster,2,mean)
    else
        cxy <- unlist(cluster)    ## assume one detector
    cluster[,] <- sweep(cluster, MARGIN=2, FUN='-', STATS=cxy)
    #######################################################

    .local$clusteri <- 1    ## updated within position()
    if (method %in% c('all','rank')) {
        ## position all, even if we will later reject some on ranks
        traps <- lapply(1:nrow(frame), position, cluster)
    }
    else {
        traps <- lapply(1:(ntrial), position, cluster)
    }
    ## subset whole clusters
    if ((edgemethod %in% c('alloutside', 'allinside', 'anyinside', 'centreinside')) || 
            (exclmethod %in% c('alloutside', 'anyoutside', 'centreoutside'))) {
        
        if (edgemethod %in% c('allinside', 'anyinside', 'centreinside')) {
            if (is.null(region)) stop ("'edgemethod' requires 'region'")
            # if (method == "GRTS") {
            #     stop ("edge clipping (edgemethod) incompatible with GRTS")
            # }
            OK <- sapply(traps, edgemethod) 
        }
        else {
            OK <- rep(TRUE, length(traps)) ## 2018-12-18
        }
        if (!is.null(exclude) & (exclmethod %in% c('alloutside', 'anyoutside', 'centreoutside'))) {
            OK <- OK & sapply(traps, exclmethod)
        }

        if (method == 'all') {
            n <- sum(OK)
        }
        
        if (sum(OK) < n) {
            stop ("not enough clusters inside polygon")
        }
        
        traps <- traps[OK]  ## subset whole clusters
    }
    traps <- traps[1:n]   ## first n usable clusters
    ## convert list of clusters to flat traps object
    if (n == 1)
        traps <- traps[[1]]
    else {
        traps$renumber <- FALSE
        traps <- do.call(rbind, traps)
    }
    ## drop excluded sites, if requested
    if (edgemethod %in% c('clip', 'centreinside')) {
        xy <- st_as_sf(traps, coords=1:2)
        st_crs(xy) <- st_crs(region)
        # OK matrix, rows are points, cols are polygons
        OK <- st_within(xy, region, sparse = FALSE)  
        OK <- apply(OK,1,any)
        traps <- subset(traps, subset = OK)
    }
    if (!is.null(exclude) && (exclmethod %in% c('clip', 'centreoutside'))) {
        xy <- st_as_sf(traps, coords=1:2)
        st_crs(xy) <- st_crs(exclude)
        OK <- st_within(xy, exclude, sparse = FALSE)  
        OK <- !apply(OK,1,any)
        traps <- subset(traps, subset = OK)
    }

    ## renumber clusters
    oldnames <- unique(clusterID(traps))
    if (attr(cluster,'detector') %in% .localstuff$polydetectors) {
        npoly <- ndetector(traps)
        npercluster <- nrow(cluster)
        polyID(traps) <- factor(rep(1:npoly, rep(npercluster, npoly)))
        clustertrap(traps) <- rep(1, nrow(traps))
        clusterID(traps) <- polyID(traps)
        vertexpart <- rep(rownames(cluster), npoly)
        row.names(traps) <- paste(polyID(traps), vertexpart, sep = '.')
    }
    else {
        clusterID(traps) <- factor(as.numeric(clusterID(traps)))
        if (nrow(cluster) == 1)
            newnames <- clusterID(traps)
        else
            newnames <- paste(clusterID(traps),
                row.names(cluster)[clustertrap(traps)], sep='.')
        row.names(traps) <- newnames
    }
    attr(traps, 'centres') <- origins[oldnames,]

    ####################################
    ## optional plot
    if (plt) {
        plot(traps, add=TRUE)
        invisible(traps)
    }
    else
        traps
    ####################################
}

###############################################################################

make.systematic <- function (n, cluster, region, spacing = NULL,
    origin = NULL, originoffset = c(0,0), chequerboard = c('all','black','white'), 
    order = c('x', 'y', 'xb', 'yb'), 
    rotate = 0, centrexy = NULL, 
    keep.design = TRUE, ...) {
    
    ## 'cluster' is a traps object for one module
    ## 'region' is a survey region
    ## ... arguments passed to trap.builder (rotate, detector)
    temporigin <- origin
    chequerboard <- match.arg(chequerboard)
    order <- match.arg(order)
    region <- boundarytoSF(region)
    if (rotate != 0) {
        ## 2022-02-01 see utility.R for sfrotate
        if (is.null(centrexy)) {
            centrexy <- st_centroid(st_as_sfc(st_bbox(region)))
            centrexy <- st_coordinates(centrexy)
        }
        region <- sfrotate(region, degrees = -rotate, centrexy = centrexy, usecentroid = FALSE)
    }

    bbox <- st_bbox(region)
    wd <- bbox[3]-bbox[1]
    ht <- bbox[4]-bbox[2]
    
    if (missing(cluster)) {
        ## this case is passed to trap builder for single detector placement
        ## if ... does not include detector, detector defaults to 'multi'
        cluster <- NULL
        clwd <- 0
        clht <- 0
    }
    else {
        clwd <- diff(range(cluster$x))
        clht <- diff(range(cluster$y))
    }

    wx <- clwd/2
    wy <- clht/2
    if (!is.null(spacing)) {
        rx <- spacing[1]
        ry <- ifelse(length(spacing)>1, spacing[2], rx)
        nx <- round ((wd-2*wx)/rx + max(spacing)/rx)  ## extra to ensure coverage 2018-10-13
        ny <- round ((ht-2*wy)/ry + max(spacing)/ry)  ## extra to ensure coverage 2018-10-13
    }
    else {
        if (length(n)>1) {
            nx <- n[1]
            ny <- n[2]
        }
        else {
            area <- sum(st_area(region))
            ## 2022-04-27 convert to numeric to avoid units
            cell <- as.numeric(sqrt(area / n))
            nx <- round ((wd - 2*wx) / cell) 
            ny <- round ((ht - 2*wy) / cell)
        }
        rx <- (wd - 2*wx) / nx
        ry <- (ht - 2*wy) / ny
    }
    ## 2022-04-27 convert to numeric to avoid units
    rxy <- as.numeric(c(rx,ry))
    minxy <- bbox[1:2]
    if (is.null(origin)) {
        origin <- runif(2) * rxy + minxy + originoffset
        originbox <- data.frame(
            x = minxy[1] + originoffset[1] + c(0,0,rx,rx),
            y = minxy[2] + originoffset[2] + c(0,ry,ry,0))
    }
    else {
        origin <- unlist(origin)  ## 2018-12-28
        originbox <- NULL
    }

    rowcol <- expand.grid(row = 0:(ny+1), col = 0:(nx+1))
    centres <- data.frame(x = rowcol$col * rx + origin[1],
                          y = rowcol$row * ry + origin[2])
    ##-----------------------------------------------------------------
    nx2 <- nx+2; ny2 <- ny+2
    if (order == 'y') temp <- 1:nrow(centres)
    if (order == 'x') temp <- t(matrix(1:nrow(centres), ncol = nx2))
    if (order == 'yb') {
        temp <- matrix(1:(nx2*ny2), ncol = ny2)
            for (i in seq(2,ny2,2)) temp[,i] <- rev(temp[,i])
    }
    if (order == 'xb') {  
        temp <- t(matrix(1:(nx2*ny2), ncol = nx2))
        for (i in seq(2,nx2,2)) temp[i,] <- rev(temp[i,])
    }
    row.names(centres) <- as.numeric(temp)
    centres <- centres[order(temp),]   
    ##-----------------------------------------------------------------
    
    if (chequerboard != 'all') {
        if (order != 'y') stop ("chequerboard option requires order = 'y'")
        whitesquares <- trunc(rowcol$row + rowcol$col + 0.1) %% 2L == 1L
        if (chequerboard == 'white')
            centres <- centres[whitesquares,]
        else 
            centres <- centres[!whitesquares,]
    }
    
    args <- list(...)
    if (!is.null(args$edgemethod)) {
        if (args$edgemethod %in% c('allinside', 'centreinside','allowoverlap')) {
            sfcentres <- st_as_sf(centres, coords=1:2)
            st_crs(sfcentres) <- st_crs(region)                       
            OK <- st_within(sfcentres, region, sparse = FALSE)
            centres <- st_coordinates(sfcentres)[OK,, drop = FALSE]
        }
    }
    traps <- trap.builder (cluster = cluster, frame = centres, region = region,
        method = 'all', ...)
    
    if (keep.design) {
        design <- list (
            'function' = 'make.systematic',
            n = if (missing(n)) NULL else n, 
            cluster = cluster, 
            region = region, 
            spacing = spacing,
            origin = temporigin, 
            originoffset = originoffset, 
            chequerboard = chequerboard,
            order = order, 
            rotate = rotate, 
            centrexy = centrexy)
        design <- c(design, list(...))
        attr(traps, 'design') <- design
    }
    else {
        attr(traps, 'origin') <- origin
        attr(traps, 'originbox') <- originbox
    }
    
    if (rotate != 0) {
        ## reverse rotation applied to region polygon
        traps <- secr::rotate(traps, rotate, centrexy)
    }
    
    traps
}

###############################################################################
make.lacework <- function (region, spacing = c(100, 20), times = NULL, 
                           origin = NULL, rotate = 0, 
                           radius = NULL, detector = "multi", keep.design = TRUE) {
    spacing <- as.numeric(spacing)
    if (length(spacing) == 1) {
        if (is.null(times)) stop("make.lacework requires times if spacing length 1")
        b <- spacing[1]
        a <- b * times
    }
    else {
        if (length(spacing) != 2) {
            stop("invalid spacing 2-vector in make.lacework")
        }
        a <- spacing[1]
        b <- spacing[2]
    }
    temporigin <- origin
    fraction <- 1.0  ## suspended code
    region <- boundarytoSF(region)
    A <- st_area(region)
    K <- A/a^2 * (2*a/b - 1)
    if (K>5000) stop("Expected number of detectors ", K, " exceeds 5000")
    if (rotate != 0) {
        centrexy <- st_centroid(st_as_sfc(st_bbox(region)))
        centrexy <- st_coordinates(centrexy)
        region <- sfrotate(region, degrees = -rotate, centrexy = centrexy, usecentroid = FALSE)
    }    
    
    bbox <- matrix(st_bbox(region), ncol = 2)                   ## after rotation
    if (is.null(origin)) {
        origin <- bbox[,1]
        origin <- origin - runif(2)*a
    }
    rx <- diff(bbox[1,])
    ry <- diff(bbox[2,])
    n1 <- ceiling((rx + a) / a) + 2
    n2 <- ceiling((ry + a) / b) + 2
    gridx <- make.grid(nx=n1, ny=n2, spacex = a, spacey = b, 
        originxy = origin, detector = detector, ID = 'numxb')
    n1 <- ceiling((ry + a) / a) + 2
    n2 <- ceiling((rx + a) / b) + 2
    gridy <- make.grid(ny=n1, nx=n2, spacey = a, spacex = b, 
        originxy = origin, detector = detector, ID = 'numyb')
    if (fraction < 1) {
        OKx <- ((gridx$y-origin[2]) %% a) < fraction * a
        OKy <- ((gridy$x-origin[1]) %% a) < fraction * a
        gridx <- subset(gridx,OKx)
        gridy <- subset(gridy,OKy)
    }
    grid <- rbind(gridx, gridy, renumber = FALSE, checkdetector = FALSE)
    dupl <- duplicated(round(grid,5))
    crossings <- subset(grid, dupl)
    grid <- subset(grid, !dupl)
    grid <- subset(grid, pointsInPolygon(grid, region))
    crossings <- subset(crossings, pointsInPolygon(crossings, region))
    if (!is.null(radius)) {
        grid <- subset(grid, distancetotrap(grid, crossings)<=radius)
    }
    rownames(grid) <- sapply(lapply(strsplit(rownames(grid), ".", TRUE), rev), paste, collapse='.')
    if (rotate != 0) {
        grid <- secr::rotate(grid, degrees = rotate, centrexy = centrexy)
        crossings <- secr::rotate(crossings, degrees = rotate, centrexy = centrexy)
    }
    attr(grid, 'crossings') <- as.matrix(crossings)
    if (keep.design) {
        design <- list (
            'function' = 'make.lacework',
            region = region, 
            spacing = spacing,
            origin = temporigin, 
            rotate = rotate, 
            radius = radius,
            detector = 'multi')
        attr(grid, 'design') <- design
    }
    
    grid
}
###############################################################################

cluster.counts <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires capthist object")
    clust <- clusterID(traps(object))
    if (is.null(clust) | (length(clust) ==0) ) {
        clust <- factor(1: nrow(traps(object)))
        warning ("clusters not defined, so treating each detector as a cluster")
    }
    cl <- clust[trap(object, names = FALSE, sortorder = 'snk')]
    tmp <- data.frame(ID=animalID(object, sortorder = 'snk'), cluster = cl)
    sapply(split(tmp,tmp$cluster), function(x) length(unique(x$ID)))
}
###############################################################################

cluster.centres <- function (object) {
    if (!inherits(object, 'traps'))
        stop ("requires traps object")
    clust <- clusterID(object)
    if (is.null(clust) | (length(clust) ==0) ) {
        clust <- factor(1: nrow(object))
        warning ("clusters not defined, so treating each detector as a cluster")
    }
    data.frame(x = tapply(object$x,clust,mean),
        y = tapply(object$y,clust,mean))
}
###############################################################################

# lacework test
# library(secr)
# bx <- data.frame(x=c(0,0,270,270,0), y = c(0,270,270,0,0))
# plot(make.lacework(bx, c(20,5)), border=10)
# lines(bx)
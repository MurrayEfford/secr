###############################################################################
## package 'secr'
## make.mask.R
## 2011 10 10 transferred from methods.R
## 2012 04 10 fixed bug in ymax of bounding box
## 2012 04 11 added 'rectangular' mask type
## 2014-03-22 all polygon type with no 'traps'
## 2016-02-17 allow SpatialPolygons
## 2016-03-02 allow poly to be mask object
## 2016-03-21 random.origin
## 2016-03-23 cell.overlap, using cellPointsInPolygon
## 2017-11-13 stop if provided with capthist
## 2022-02-01 use sf functions
###############################################################################

getCentres <- function (xy) {
  nrxy <- nrow(xy)
  if (nrxy > 1)
    (xy[-1,] + xy[-nrxy,]) / 2
  else
    xy
}
cellPointsInPolygon <- function (object, poly, cell.overlap = c("centre","any","all"), spacing) {
  if (inherits(object, "mask") & missing(spacing))
    spacing <- spacing(object)
  cell.overlap <- match.arg(cell.overlap)
  if (cell.overlap == "centre") {
    inside <- pointsInPolygon(object, poly)  # vector
  }
  else {
    sp2 <- spacing/2 * c(-1, +1)
    dxy <- expand.grid(dx=sp2, dy=sp2)
    vertices <- lapply(1:4, function(i) sweep(object, MARGIN=2, STATS=unlist(dxy[i,]), FUN="+"))
    vertices <- array(unlist(vertices), dim = c(nrow(object), 2, 4))
    inside <- apply(vertices, 3, pointsInPolygon, poly)
  }
  if (cell.overlap=="any")
    inside <- apply(inside,1,any)
  if (cell.overlap=="all")
    inside <- apply(inside,1,all)
  inside
}

make.mask <- function (traps, buffer = 100, spacing = NULL, nx = 64, ny = 64,
  type = c("traprect", "trapbuffer", "pdot", "polygon", "clusterrect",
    "clusterbuffer", "rectangular", "polybuffer"), 
  poly = NULL, poly.habitat = TRUE, cell.overlap = c("centre","any","all"), 
  keep.poly = TRUE, check.poly = TRUE, pdotmin = 0.001, random.origin = FALSE, ...)
{
  type <- match.arg(type)
  if (missing(traps)) traps <- NULL
  if (ms(traps)) {         ## a list of traps objects
    if (inherits(poly, 'list') & (!is.data.frame(poly)))
      stop ("lists of polygons not implemented in 'make.mask'")
    ## 2014-09-20 2018-10-11 now passes keep.poly and check.poly
    temp <- lapply (traps, make.mask, buffer = buffer, spacing = spacing, nx = nx, ny = ny,
      type = type, poly = poly, poly.habitat = poly.habitat, keep.poly = keep.poly,
      check.poly = check.poly,  pdotmin = pdotmin, random.origin = random.origin, ...)
    class (temp) <- c('mask', 'list')
    temp
  }
  else {
    dots <- match.call(expand.dots = FALSE)$...
    if ((length(dots)==0) & (type == 'pdot'))
      warning ("no detection parameters supplied; using defaults", call. = FALSE)
    
    if (inherits(traps, 'capthist'))
      stop("argument should be a traps object or equivalent, not capthist")
    ## extend buff to allow for random origin 2016-03-21
    if (random.origin) {
      if (is.null(spacing))                     
        stop ("random.origin requires that spacing to be specified")
      if (!(type %in% c("traprect","trapbuffer","polygon", "rectangular")))
        stop ("random.origin not implemented for this mask type")
      if (type == "polygon") {
        buffer <- 0    ## to ensure no change in behaviour secr 2.10.3
      }
      offset <- runif(2)
      buffx <- c(-buffer-offset[1]*spacing,+buffer+(1-offset[1])*spacing)
      buffy <- c(-buffer-offset[2]*spacing, +buffer+(1-offset[2])*spacing)
    }
    else {
      buffx <- buffy <- c(-buffer,+buffer)
    }  
    
    if (!inherits(poly, "mask")) {    
      poly <- boundarytoSF(poly)  # standardise on sfc and mask
    }
    
    if (is.null(traps)) check.poly <- FALSE
    
    if (is.null(traps) & !(type %in% c('polygon','polybuffer')))
      type <- 'rectangular'
    if (type == 'rectangular') {
      if (is.null(spacing))
        stop ("require spacing for rectangular mask")
      xl <- c(0, spacing * nx)
      yl <- c(0, spacing * ny)
    }
    else if (type %in% c('polygon', 'polybuffer')) {
      if (is.null(poly))
        stop ("mask polygon must be supplied")
      if (!poly.habitat)
        stop ("types 'polygon' and 'polybuffer' not compatible with nonhabitat")
      if (inherits(poly, "mask")) {     
        xl <- range(poly[,1]) + c(-1,1) * spacing(poly)/2
        yl <- range(poly[,2]) + c(-1,1) * spacing(poly)/2
      }
      else  {    # sfc
        xl <- st_bbox(poly)[c(1,3)]
        yl <- st_bbox(poly)[c(2,4)]
      }
    }
    else {
      xl <- range(traps$x)
      yl <- range(traps$y)
    }
    ## 2016-03-21 for random
    xl <- xl + buffx
    yl <- yl + buffy
    
    if (is.null(spacing)) spacing <- diff(xl) / nx
    
    if (type %in% c('clusterrect', 'clusterbuffer')) {
      ID <- clusterID(traps)
      meanx <- unique(tapply(traps$x, ID, mean))
      meany <- unique(tapply(traps$y, ID, mean))
      cluster <- subset(traps, subset = clusterID(traps)==1) ## extract a single cluster
      ## assume identical wx, wy are half-width and half-height of a box
      ## including the cluster and the rectangular buffer
      wx <- diff(range(cluster$x)) / 2 + buffer
      wy <- diff(range(cluster$y)) / 2 + buffer
      wx <- round(wx/spacing) * spacing   ## to make symmetrical
      wy <- round(wy/spacing) * spacing   ## to make symmetrical
      dx <- seq(-wx,wx,spacing)
      dy <- seq(-wy,wy,spacing)
      x <- as.numeric(outer(FUN='+', dx, meanx))
      y <- as.numeric(outer(FUN='+', dy, meany))
    }
    else {
      x <- seq(xl[1] + spacing/2, xl[2], spacing)
      y <- seq(yl[1] + spacing/2, yl[2], spacing)
      
    }
    
    mask   <- expand.grid (x=x, y=y)
    attr(mask,'out.attrs') <- NULL   ## added 2009 07 03
    
    if (type=='trapbuffer') {
      ## appropriate convex buffer 2011-01-22
      ## (this re-use of nx may not be appropriate)
      
      if (!is.null(detector(traps)) &    ## 2017-01-27
          all(detector(traps) %in% c('polygon','polygonX'))) {
        temp <- bufferContour(traps, buffer = buffer, nx = nx,
          convex = T, plt = F)
        OK <- array(dim=c(length(x), length(y), length(temp)))
        for (i in 1:length(temp)) {
          OK[,,i] <- pointsInPolygon(mask, temp[[i]][,c('x','y')])
        }
        OK <- apply(OK, 1:2, any)
        mask <- mask[OK,,drop=F]
      }
      else {
        mask <- mask[distancetotrap(mask, traps) <= buffer,]
      }
    }
    
    if (type=='polybuffer') {
      if (inherits(poly, 'mask')) {
        polymask <- make.mask(type = 'polygon', poly = poly, keep.poly = FALSE, check.poly = FALSE,
          spacing = spacing/2)   # arbitrary
        mask <- mask[distancetotrap(mask, polymask) <= buffer,]
      }
      else {  # sfc
        bufferedpoly <- st_buffer(poly, dist = buffer)
        mask <- mask[pointsInPolygon(mask, bufferedpoly),]
      }
    }
    
    if (type=='clusterbuffer') {
      mask <- mask[distancetotrap(mask, traps) <= buffer,]
    }
    
    if (type=='pdot') {
      dettype <- detectorcode(traps)
      if (!all(dettype %in%  c(0,1,2,5,8,13))) {
          stop ("type pdot is available only for point detectors")
      }
      OK <- pdot(mask, traps = traps, ...) > pdotmin
      edge <- function (a,b) any (abs(a-b) < (spacing))
      mask <- mask[OK,]
      attr(mask,'pdotmin') <- pdotmin   # save nominal threshold
      if (edge(mask[,1],xl[1]) |
          edge(mask[,1],xl[2]) |
          edge(mask[,2],yl[1]) |
          edge(mask[,2],yl[2]))
        warning ("'pdot' mask may have been truncated; ",
          "possibly increase buffer", call. = FALSE)
    }
    
    if ((!is.null(poly)) & (type != 'polybuffer' )) {
      inpoly <- cellPointsInPolygon(mask, poly, cell.overlap, spacing)
      if (poly.habitat) {
        mask <- mask[inpoly,]
        if (check.poly) {
          if (any (!pointsInPolygon(traps, poly)))
            warning ("some traps are not inside habitat polygon(s)", call. = FALSE)
        }
      }
      else {
        mask <- mask[!inpoly,]
        if (check.poly)
          if (any (pointsInPolygon(traps, poly)))
            warning ("some traps are inside non-habitat polygon(s)", call. = FALSE)
      }
      if (keep.poly) {
        attr(mask, 'polygon') <- poly   # save
        attr(mask, 'poly.habitat') <- poly.habitat   # save
      }
    }
    
    if (nrow(mask)>0) {
      xl <- range(mask$x) + spacing/2 * c(-1,1)
      yl <- range(mask$y) + spacing/2 * c(-1,1)
    }
    else {
      warning("no points in mask", call. = FALSE)
      xl <- c(-Inf,Inf)
      yl <- c(-Inf,Inf)
    }
    
    attr(mask,'type')        <- type
    attr(mask,'meanSD')      <- getMeanSD (mask)
    attr(mask,'area')        <- spacing^2 * 0.0001
    attr(mask,'spacing')     <- spacing
    attr(mask,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
    class(mask)  <- c('mask', 'data.frame')
    
    mask
  }
}
###############################################################################

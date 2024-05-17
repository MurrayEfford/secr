# 2024-05-18

make.spcosa <- function (n, cluster, region, rotate = 0, randomize = FALSE, maxtries = 100, keep.mask = FALSE, ...) {
  allinside <- function (xy, subregion) {
    xy <- st_as_sf(data.frame(xy), coords=1:2)
    st_crs(xy) <- st_crs(subregion)
    all( apply(st_within(xy, subregion, sparse = FALSE),1,any))
  }
  
  if (!requireNamespace("spcosa", quietly = TRUE)) {
    stop ("package 'spcosa' is required for make.spcosa")
  }
  
  msk <- make.mask(type = 'polygon', buffer = 0, poly = region, keep.poly = FALSE, ...)
  pix <- sp::SpatialPixels(sp::SpatialPoints(as.matrix(msk)))
  aa <- spcosa::stratify(pix, n, equal = TRUE)  ## SLOW
  covariates(msk)$stratum <- aa@stratumId + 1
  if (randomize) {
    onestra <- function (stra) {
      msk <- subset(msk, covariates(msk)$stratum == stra)
      subregion <- gridCells(msk)
      nm <- nrow(msk)
      OK <- FALSE
      for (i in 1:maxtries) {
        tmp <- trap.builder(cluster = cluster, frame = msk[sample.int(nm,1),], 
                            method = "all", region = subregion, rotation = rotate, 
                            edgemethod = 'allowoverlap')
        OK <- allinside(tmp, subregion)
        if (OK) break
      }
      if (!OK) stop("random location within subregion not found in ", maxtries, ' attempts')
      tmp
    }
    clusterlist <- lapply(1:n, onestra)
    cosa <- do.call(rbind, clusterlist)
    # stop ("option not yet available")
  }
  else {
    cent <- sp::coordinates(spcosa::getCentroid(aa))
    cosa <- trap.builder(cluster = cluster, frame = cent, method = 'all', rotation = rotate)
  }
  if (keep.mask) attr(cosa, 'mask') <- msk
  cosa

}



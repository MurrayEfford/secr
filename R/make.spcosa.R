# 2024-05-18

make.spcosa <- function (n, cluster, poly, rotate = 0, randomize = FALSE, keep.mask = FALSE, ...) {
  
  if (!requireNamespace("spcosa", quietly = TRUE)) {
    stop ("package 'spcosa' is required for make.spcosa")
  }
  
  msk <- make.mask(type = 'polygon', buffer = 0, poly = poly, keep.poly = FALSE, ...)
  pix <- sp::SpatialPixels(sp::SpatialPoints(as.matrix(msk)))
  aa <- spcosa::stratify(pix, n, equal = TRUE)  ## SLOW
  covariates(msk)$stratum <- aa@stratumId + 1
  if (randomize) {
    stop ("option not yet available")
  }
  else {
    cent <- sp::coordinates(spcosa::getCentroid(aa))
  }
  cosa <- trap.builder(cluster = cluster, frame = cent, method = 'all', rotation = rotate)
  if (keep.mask) attr(cosa, 'mask') <- msk
  cosa

}



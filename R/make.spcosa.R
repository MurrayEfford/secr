# 2024-05-17

make.spcosa <- function (n, cluster, poly, rotate = 0, keep.mask = FALSE, ...) {
    msk <- make.mask(type = 'polygon', buffer = 0, poly = poly, keep.poly = FALSE, ...)
    pix <- sp::SpatialPixels(sp::SpatialPoints(as.matrix(msk)))
    aa <- spcosa::stratify(pix, n, equal = TRUE)  ## SLOW
    covariates(msk)$stratum <- aa@stratumId + 1
    cent <- sp::coordinates(spcosa::getCentroid(aa))
    cosa <- trap.builder(cluster = cluster, frame = cent, method = 'all', rotation = rotate)
    if (keep.mask) attr(cosa, 'mask') <- msk
    cosa
}



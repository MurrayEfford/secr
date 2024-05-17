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

source('d:/density communication/secrbook/figures/setup.R')
polyexample1 <- read.traps(file = 'data/polygonexample.txt', detector = 'polygon')
grid <- make.grid(5,2, spacing = 10, detector = 'proximity')

set.seed(1234)
test <- make.spcosa (5, grid, polyexample1, keep.mask = T, rotate = 0)
set.seed(1234)
test1 <- make.spcosa (5, grid, polyexample1, keep.mask = T, rotate = -1)

par(mfrow=c(1,2))

plot(attr(test,'mask'), cov='stratum', col=qual8, legend = FALSE)
plot(test,add=T)
plot(attr(test1,'mask'), cov='stratum', col=qual8, legend = FALSE)
plot(test1,add=T)

polyID(test)
clusterID(test)
clustertrap(test)

nrow(test)

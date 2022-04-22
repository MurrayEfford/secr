## 2022-01-05 Started 
## 2022-03-15 tests using 'terra' suppressed 

library(secr)
library(sf)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

####################################################################

set.seed(123)
detectors <- make.grid (nx = 6, ny = 6, detector = "count")
CH1 <- sim.capthist (detectors, popn = list(D = 0.5, buffer = 100), 
    detectfn = 'HHN', detectpar = list(lambda0 = 0.2, sigma = 25), 
    noccasions = 4, nsessions = 3, renumber = FALSE)

test_that("reduce.capthist accepts zero, one, or more animals per session", {
    sum1 <- summary(CH1, terse = TRUE)
    CH2 <- reduce(CH1, outputdetector = "count", by = "all", 
        verify = FALSE, dropunused = FALSE)
    sum2 <- summary(CH2, terse = TRUE)
    expect_equal(sum1[2:4,], sum2[2:4,], tolerance = 1e-4, check.attributes = FALSE)
})

####################################################################

## 2022-01-21 4.5.2

test_that("home range methods work with empty capthist", {
    ch0 <- subset(captdata, 0)
    expect_equal(dbar(ch0), NA)
    expect_equal(MMDM(ch0), NA)
    expect_equal(ARL(ch0), NA)
    expect_equal(RPSV(ch0), NA)
    expect_equal(centroids(ch0), NA)  # avoid conflict terra::centroids
    expect_equal(length(moves(ch0)), 0)
})
####################################################################

test_that("correct RPSV", {
    expect_equal(RPSV(captdata, CC = TRUE), 25.62888, tolerance = 1e-4)
})
####################################################################

# implicitly tests boundarytoSF()

test_that("points in polygon", {
    set.seed(123)
    xy <- matrix(runif(20)*10, ncol = 2)
    region1 <- cbind(x = c(0,6,6,0,0), y = c(0,0,6,6,0))  
    inside <- c(F,F,F,F,F,F,T,F,T,F)
    inside1 <- pointsInPolygon(xy, region1)            # matrix
    region2 <- st_sfc(st_polygon(list(region1)))
    inside2 <- pointsInPolygon(xy, region2)            # sfc
    region3 <- st_sf(geometry=region2)
    inside3 <- pointsInPolygon(xy, region3)            # sf
    region4 <- as(region3, "Spatial")
    inside4 <- pointsInPolygon(xy, region4)            # SpatialPolygonsDataFrame
    # region5 <- terra::vect(region4)
    # inside5 <- pointsInPolygon(xy, region5)            # SpatVector
    expect_equal(inside, inside1) 
    expect_equal(inside, inside2)
    expect_equal(inside, inside3)
    expect_equal(inside, inside4)
    # expect_equal(inside, inside5)
})
####################################################################

## 2022-02-13 4.5.3

test_that("addCovariates raster datasource", {
    
    spatialdata0 <- possummask                                  # mask 
    spatialdata1 <- raster(possummask, cov = 'd.to.shore')      # RasterLayer
    spatialdata2 <- as(spatialdata1, "SpatialGridDataFrame")    # SpatialGridDataFrame
    # spatialdata3 <- rast(possummask, cov = 'd.to.shore')        # SpatRaster
    
    xy <- data.frame(x = c(2697649, 2698074, 2698438), y = c(6078402, 6078019, 6077678))
    trp <- read.traps(data = xy)
    
    # mask
    expect_warning(cov0 <- covariates(addCovariates(trp, spatialdata0, strict = TRUE))[,1])
    
    # RasterLayer
    expect_warning(cov1 <- covariates(addCovariates(trp, spatialdata1))[,1])
    
    # SpatialGridDataFrame
    expect_warning(cov2 <- covariates(addCovariates(trp, spatialdata2))[,1])
    
    # SpatRaster
    # expect_warning(cov3 <- covariates(addCovariates(trp, spatialdata3))[,1])
    
    target <- c(NA, 416.74212650, 815.38825108)
    expect_equal(cov0, target) 
    expect_equal(cov1, target) 
    expect_equal(cov2, target) 
    # expect_equal(cov3, target) 
})
####################################################################

test_that("addCovariates vector datasource", {
    
    ## sf
    spatialdata4 <- st_read(system.file('extdata/OVforest.shp', package = 'secr'), quiet = TRUE)
    ## shapefile name
    spatialdata5 <- system.file('extdata/OVforest.shp', package = 'secr')
    ## SpatialPolygonsDataFrame
    # warning from discarded proj4 
    spatialdata6 <- suppressWarnings(as(spatialdata4, "Spatial"))
    
    xy <- data.frame(
        x = c(2674961, 2674864, 2674890, 2675075, 2674563),
        y = c(5982263, 5982470, 5982721, 5982919, 5982797))
    trp <- read.traps(data = xy)

    # warnings from NA returned for point outside polygons
    expect_warning(cov4 <- covariates(addCovariates(trp, spatialdata4))$forest)
    expect_warning(cov5 <- covariates(addCovariates(trp, spatialdata5))$forest)
    expect_warning(cov6 <- covariates(addCovariates(trp, spatialdata6))$forest)
    
    target <- c("nonbeech", "beech", "nonbeech", NA, "nonbeech")
    
    expect_equal(cov4, target) 
    expect_equal(cov5, target) 
    expect_equal(cov6, target) 
    
})
####################################################################

## future tests

# subset
# reduce          # tested implicitly with fastproximity
# join
# rbind

## for data types
# point  
# polygon
# transect
# signal

## using RPSV as test criterion?

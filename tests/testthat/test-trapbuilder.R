# test GRTS

library(secr)
library(sf)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

shapefilename <- system.file("extdata/OVforest.shp", package = "secr")
OVforest <- sf::st_read(shapefilename, quiet = TRUE)[1:2,1]

test_that("method = 'SRS' works", {
    set.seed(1235)
    # simple random sample
    expect_silent(sample0 <- trap.builder(n = 32, region = OVforest, method = 'SRS'))
    expect_equal(spacing(sample0), 76.3612953003, tolerance = 1e-4, check.attributes = FALSE)
})

# conditional from 2025-03-27 (see )
if (requireNamespace("spsurvey", versionCheck=list(op=NULL, version = ">=5.3.0"), quietly = TRUE)) {
    
    test_that("method = 'GRTS' works", {
        set.seed(1235)
        # GRTS sample
        expect_silent(sample1 <- trap.builder(n = 32, region = OVforest, method = 'GRTS'))
        expect_equal(spacing(sample1), 94.7225747307, tolerance = 1e-4, check.attributes = FALSE)
    })
    
}

## Started 2021-12-14

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

####################################################################
## Prompted by problems with multi-part transects, fixed in 4.5.0
## John Nettles post phidot 7/12/21

## cut-down version of transect example in secr-polygondetectors.pdf

x <- seq(0, 2*pi, length = 10)
transect <- make.transect(x = x*100, y = sin(x)*300, exclusive = FALSE)
transectCH <- sim.capthist(transect, popn = list(D = 2, buffer = 300),
    detectfn = 'HHN', detectpar = list(lambda0 = 1.0, sigma = 50), 
    binomN = 0, seed = 123)
msk <- make.mask(transect, buffer = 200, spacing = 30, type = 'trapbuffer')
st <- list(D = 2, lambda0 = 1, sigma = 50)

test_that("correct likelihood (Poisson transect data)", {
    
    # 3 parts, each 450 m
    expect_warning(ch3 <- snip(transectCH, by = 450, keep.incomplete = FALSE))
    
    LL3 <- secr.fit(ch3, mask = msk, detectfn = 'HHN', binomN = 0, 
        details = list(LLonly = TRUE), start = st)
    expect_equal(LL3, -1674.425557, tolerance = 1e-4, check.attributes = FALSE)
    
    # discretized, 48 x 30-m intervals
    chd <- discretize(ch3, spacing = 30,  outputdetector = 'count')
    # If omit detectfn = 'HHN' then detectfn = 'HN'; start lambda0 = 1.0
    # is ignored and generates varying default start value for g0, 
    # hence erratic LL
    LLd <- secr.fit(chd, mask = msk, detectfn = 'HHN', binomN = 0, 
        details = list(LLonly = TRUE), start = st)
    expect_equal(LLd, -450.8793538, tolerance = 1e-4, check.attributes = FALSE)
    
})

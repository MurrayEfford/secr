# test CL density model

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

msk <- make.mask(traps(captdata), buffer=100, type='trapbuffer', nx = 32)
fitrD  <- secr.fit(captdata, CL = TRUE, mask = msk, model= D~x, trace = FALSE)
fitrDi <- secr.fit(captdata, CL = TRUE, mask = msk, model= D~x, trace = FALSE, 
                   link=list(D='identity'))

test_that("relativeD estimates correct", {
    expect_equal(coef(fitrD)[1,1], 0.01261994, tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(derivedDbeta0(fitrD), 1.700941, tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(derivedDbeta0(fitrDi), 5.479395, tolerance = 1e-4, check.attributes = FALSE)
})

test_that("region.N correct with relativeD", {
    expect_equal(region.N(fitrD)[1,1:4], c(116.3148,1.874043, 112.6993, 120.0462), 
                 tolerance = 1e-4, check.attributes = FALSE)
})

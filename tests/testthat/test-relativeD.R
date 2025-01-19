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
    expect_equal(derivedDcoef(fitrD)[1,1], 1.700941, tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(derivedDcoef(fitrDi)[1,1], 5.479395, tolerance = 1e-4, check.attributes = FALSE)
})

test_that("region.N correct with relativeD", {
    expect_warning(Nhat <- region.N(fitrD)[1,1:2])
    expect_equal(Nhat, c(116.314763,NA), tolerance = 1e-4, check.attributes = FALSE)
})

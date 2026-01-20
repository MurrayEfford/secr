# test Dlambda model
# 2026-01-20

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

msk <- make.mask(traps(ovenCH[[1]]), buffer = 300, type='trapbuffer', nx = 32)

test_that("Dlambda likelihood correct", {
    LL.0 <- secr.fit(ovenCH, model = D~1, mask = msk, trace = TRUE, 
                    details = list(Dlambda = TRUE, LLonly = TRUE), 
                    start = c(0,0,-3.5,4), ncores = 2)
    LL.CL <- secr.fit(ovenCH, model = D~1, mask = msk, trace = TRUE, 
                    details = list(Dlambda = TRUE, LLonly = TRUE), 
                    start = c(0,0,-3.5,4), ncores = 2, CL = TRUE)
    expect_equal(LL.0,  -956.66706, tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(LL.CL, -946.41862, tolerance = 1e-4, check.attributes = FALSE)
})

# now with fitted models included in package
# ovenbird.model.D, ovenbird.model.Dl
test_that("Ovenbird Dlambda correct", {
    expect_equal(
        predictDlambda(ovenbird.model.Dl)[,1],  
        c(1.0322172, 0.9381376, 0.9381376, 0.9381376, 0.9381376), 
        tolerance = 1e-5, check.attributes = FALSE)
    expect_equal(
        predictDlambda(ovenbird.model.Dl)[2,1],  
        exp(coef(ovenbird.model.D)[2,1]), 
        tolerance = 1e-5, check.attributes = FALSE)
})

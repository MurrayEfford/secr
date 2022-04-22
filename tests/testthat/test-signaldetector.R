## Started 2022-01-05

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

####################################################################

test_that("signal likelihoods correct", {
    
    signalCH.525 <- subset(signalCH, cutval = 52.5)
    omask <- make.mask(traps(signalCH), buffer = 200)
    ostart <- c(log(20), 80, log(0.1), log(2))
    ovensong.LL.1 <- secr.fit( signalCH.525, mask = omask, start = ostart, 
        detectfn = 11, details = list(LLonly = TRUE))
    ovensong.LL.2 <- secr.fit( signalCH.525, mask = omask, start = ostart, 
        detectfn = 10, details = list(LLonly = TRUE))
    
    expect_equal(ovensong.LL.1, -1637.240249, tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(ovensong.LL.2, -846.6429107, tolerance = 1e-4, check.attributes = FALSE)
    
})

    
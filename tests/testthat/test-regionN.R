# test region.N

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

msk1 <- make.mask(traps(captdata), buffer=100, spacing = 5)
msk2 <- make.mask(traps(captdata), buffer=100, spacing = 10)
N1 <- region.N(secrdemo.0, region = msk1)[1,]
N2 <- region.N(secrdemo.0, region = msk2)[1,]

test_that("region.N independent of mask spacing", {
    expect_equal(N1, N2, tolerance = 1e-4, check.attributes = FALSE)
})

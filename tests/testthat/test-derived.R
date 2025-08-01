## Started 2022-11-20
## tests of derived quantities

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

test_that("correct chat.nk ovenbird.model.1", {
    expect_warning(chatnk <- sapply(chat.nk(ovenbird.model.1), '[[', 'chat'))
    expect_equal(chatnk, c(0.8312778, 1.7388806, 1.3127876, 1.2391087, 1.1826071), 
        tolerance = 1e-5, check.attributes = FALSE)
})

test_that("correct derived density secrdemo.CL", {
    expect_equal(derived(secrdemo.CL)['D',], 
        c(5.4798074, 0.6445541, 4.3549957, 6.8951365, 0.11470787, 0.02602675, 0.11762349),
        tolerance = 1e-6, check.attributes = FALSE)
})

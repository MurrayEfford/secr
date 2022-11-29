## Started 2022-11-20
## tests of derived quantities

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

test_that("correct chat.nk secrdemo.0", {
    expect_equal(chat.nk(secrdemo.0)$chat, 1.33591167, 
        tolerance = 1e-6, check.attributes = FALSE)
})

test_that("correct derived density secrdemo.CL", {
    expect_equal(derived(secrdemo.CL)['D',], 
        c(5.4798074, 0.64455399, 4.3549958, 6.8951362, 0.11470787, 0.026026665, 0.11762348), 
        tolerance = 1e-6, check.attributes = FALSE)
})

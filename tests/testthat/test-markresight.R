## 2022-03-06

library(secr)
set.seed(123)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

grid <- make.grid(detector = 'proximity')

# variable usage: only alternate detectors used on second occasion
usage(grid) <- matrix(1, nrow = 36, ncol = 4)
usage(grid)[,2] <- rep(0:1,18)

msk <- make.mask(grid, buffer = 100, nx = 20, type = 'trapbuffer')

# all sighting
markocc(grid) <- c(0, 0, 0, 0)

MRCH <- sim.resight(
    traps     = grid, 
    popn      = list(D = 5, pID = 1.0),
    detectpar = list(g0 = 0.3, sigma = 25), 
    nonID     = FALSE,
    unsighted = TRUE
)
    
argssecr <- list(
    capthist = MRCH, 
    mask     = msk, 
    detectfn = 'HN',
    start    = list(D = 5, g0 = 0.3, sigma = 25),
    details  = list(LLonly = TRUE)
)

test_that("correct likelihood (all sighting mark-resight)", {
    LL <- do.call(secr.fit, argssecr)[1]
    expect_equal(LL, -182.133065, tolerance = 1e-4, check.attributes = FALSE)
})

## Started 2021-11-09
## revised 2023-03-09

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

####################################################################
## Ian Durbach Boost bug with Poisson count data 2021-11-03
## problem arose from switch to Boost distributions in C++ functions 
## gpois, pski: when lambda zero Boost result is NAN instead of zero
## previously returned by dpois

# Poisson counts

set.seed(123)
detectors <- make.grid (nx = 6, ny = 8, detector = "count")
CHpois <- sim.capthist (detectors, popn = list(D = 10, buffer = 100), 
    detectpar = list(g0 = 0.2, sigma = 25), noccasions = 1)

test_that("correct likelihood (Poisson count data)", {
    args <- list(capthist = CHpois, detectfn = 'HN', buffer = 100, 
        start = c(2.121835788, -1.040097594,  3.201521728), verify = FALSE)
    
    args$details <- list(LLonly = TRUE, fastproximity = TRUE)
    LL1 <- do.call(secr.fit, args)[1]
    expect_equal(LL1, -122.138538, tolerance = 1e-4, check.attributes = FALSE)
    
    args$details <- list(LLonly = TRUE, fastproximity = FALSE)
    LL2 <- do.call(secr.fit, args)[1]
    # expect_equal(LL2, -122.1523842, tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(LL2, -122.1385378, tolerance = 1e-4, check.attributes = FALSE)  # 4.5.9
    
    args$detectfn <- 'HHN'  
    args$start <- c(2.121836235, -1.342742898, 3.201525519)
    
    args$details <- list(LLonly = TRUE, fastproximity = TRUE)
    LL3 <- do.call(secr.fit, args)[1]
    expect_equal(LL3, -122.1523842, tolerance = 1e-4, check.attributes = FALSE)

    args$details <- list(LLonly = TRUE, fastproximity = FALSE)
    LL4 <- do.call(secr.fit, args)[1]
    # expect_equal(LL4, -122.3167724, tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(LL4, -122.1523842, tolerance = 1e-4, check.attributes = FALSE) # 4.5.9
})

###############################################################################

# binomial counts

set.seed(123)
CHbinom <- sim.capthist (detectors, popn = list(D = 10, buffer = 100), 
    detectpar = list(g0 = 0.1, sigma = 25), binomN = 5, noccasions = 1)

test_that("correct likelihood (binomial count data)", {
    args <- list(capthist = CHbinom, detectfn = 'HN', binomN=5, buffer = 100, 
        start = c(2.3648, -2.4498, 3.3361 ), verify = FALSE,
        details = list(LLonly = TRUE, fastproximity = FALSE))
    LL1 <- do.call(secr.fit, args)[1]
    expect_equal(LL1, -250.716399 , tolerance = 1e-4, check.attributes = FALSE)
    usage(traps(CHbinom)) <- matrix(5, nrow = nrow(traps(CHbinom)), ncol = 1)
    args$capthist <- CHbinom
    args$binomN <- 1  # binomial size from usage
    LL2 <- do.call(secr.fit, args)[1]
    expect_equal(LL2, -250.716399 , tolerance = 1e-4, check.attributes = FALSE)
})    



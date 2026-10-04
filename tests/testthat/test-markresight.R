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

## analytic overdispersion (details$chatmethod = "analytic") for sighting-only data,
## binary detectors with variable usage; the simulated c-hat (nsim = 40000) was 7.59
## (Monte Carlo error about 1%)
test_that("analytic c-hat (all sighting mark-resight)", {
    args <- argssecr
    args$details <- list(chatmethod = "analytic", chatonly = TRUE)
    chat <- do.call(secr.fit, args)
    expect_equal(chat[1, "Tu"], 7.6817, tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(chat[1, "Tu"], 7.6, tolerance = 0.03, check.attributes = FALSE)
    expect_equal(chat[1, "Tm"], 1, check.attributes = FALSE)
    ## not available when there are marking occasions
    MRCHm <- MRCH
    markocc(traps(MRCHm)) <- c(1, 0, 0, 0)
    args$capthist <- MRCHm
    args$verify <- FALSE      # the altered data are not consistent
    expect_error(do.call(secr.fit, args), "sighting-only")
})

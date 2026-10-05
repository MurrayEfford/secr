## 2026-07-04

library(secr)

# needed for consistent randoms after test.trapbuilder.R
suppressWarnings(RNGkind(sample.kind = "Rounding"))

set.seed(123)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

################################################
# Code to simulate capthist objects (trCH, teCH)

# detectors
te <- make.telemetry()
tr <- make.grid(detector = "proximity", nx = 8, ny = 8)

# spatial population
# block gratuitous covariate 'sex'
pop4 <- sim.popn(tr, D = 5, buffer = 100, seed = 567, covariates = NULL)

# select 12 telemetered individuals from larger population
pop4C <- subset(pop4, sample.int(nrow(pop4), 12))

# renumber = FALSE (keep original animalID) needed for matching
trCH <- sim.capthist(tr,  popn = pop4, renumber = FALSE, 
                     detectfn = "HHN", detectpar = list(lambda0 = 0.1, 
                                                        sigma = 25), seed = 123)
teCH <- sim.capthist(te, popn = pop4C, renumber = FALSE, 
                     detectfn = "HHN", detectpar = list(lambda0 = 1, sigma = 25),
                     noccasions = 10, seed = 345)

################################################

msk <- make.mask(traps(trCH), buffer = 100, type = 'trapbuffer', nx = 32)
argssecr <- list(mask = msk, detectfn = 'HHN', CL = TRUE,
                 start = list(lambda0 = 0.1, sigma = 25),
                 details = list(LLonly = TRUE, safeLL = TRUE, uselog = TRUE, 
                                debug = 0, fastproximity = FALSE))

test_that("correct standalone likelihood", {
    argssecr$capthist <- teCH
    LL <- do.call(secr.fit, argssecr)[1]
    expect_equal(LL, -1159.6049, tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct combined likelihood, independent telemetry", {
    combinedCHI <- addTelemetry(trCH, teCH, type = "independent")
    argssecr$capthist <- combinedCHI
    LL <- do.call(secr.fit, argssecr)[1]
    expect_equal(LL, -1422.7557, tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct combined likelihood, independent telemetry in separate session", {
    combinedCHS <- MS.capthist(trCH, teCH)   # session independence
    argssecr$capthist <- combinedCHS
    LL <- do.call(secr.fit, argssecr)[1]
    expect_equal(LL, -1422.7557, tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct combined likelihood, concurrent telemetry", {
    combinedCHC <- addTelemetry(trCH, teCH, type = "concurrent")
    argssecr$capthist <- combinedCHC
    LL <- do.call(secr.fit, argssecr)[1]
    expect_equal(LL, -1406.7027, tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct combined likelihood, dependent telemetry", {
    expect_warning(combinedCHD <- addTelemetry(trCH, teCH, type = "dependent") )
    argssecr$capthist <- combinedCHD
    LL <- do.call(secr.fit, argssecr)[1]
    expect_equal(LL, -824.33822, tolerance = 1e-4, check.attributes = FALSE)
})


## With telemetry, details safeLL and uselog default to TRUE: the product of many fixes at metre
## scale underflows otherwise (secr 5.5.1)
test_that("log-sum likelihood is the default with telemetry", {
    set.seed(2)
    tr  <- make.grid(nx = 6, ny = 6, spacing = 12000, detector = "proximity")
    msk <- make.mask(tr, buffer = 40000, spacing = 4000)
    pop <- sim.popn(D = 2e-4, core = tr, buffer = 40000, seed = 3)
    ch  <- sim.capthist(tr, popn = pop, detectfn = "HHN", renumber = FALSE, noccasions = 6, seed = 4,
                        detectpar = list(lambda0 = 0.2, sigma = 10000))
    tepop <- subset(pop, row.names(pop) %in% row.names(ch)[1:8])
    teCH  <- sim.capthist(make.telemetry(), popn = tepop, detectfn = "HHN", renumber = FALSE, 
                          noccasions = 150, seed = 5, detectpar = list(lambda0 = 1, sigma = 10000))
    comb <- suppressWarnings(addTelemetry(ch, teCH, type = "concurrent"))
    LL <- function (...) as.numeric(secr.fit(comb, mask = msk, detectfn = "HHN", CL = TRUE, trace = FALSE,
                                             start = list(lambda0 = 0.2, sigma = 10000),
                                             details = list(LLonly = TRUE, ...)))
    expect_equal(LL(), -27028.239, tolerance = 1e-6)
    expect_equal(LL(), LL(safeLL = TRUE, uselog = TRUE))
    expect_lt(LL(safeLL = FALSE, uselog = FALSE), -1e9)    # explicit FALSE is respected (underflow)
})

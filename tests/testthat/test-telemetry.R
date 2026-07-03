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


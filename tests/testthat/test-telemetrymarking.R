# test telemetry type 'marking': collared animals (telemetry plus camera cues) and
# unmarked camera cues Tu; secr 5.5.1

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

## simulated data: 10 collared animals x 100 fixes, 81 count cameras, 5 occasions,
## unmarked cue total Tu = 58
ch <- readRDS(test_path("telemetry_marking_ch.RDS"))
trps <- make.grid(detector = "count", nx = 9, ny = 9, spacing = 90/8, origin = c(30, 30))
msk <- make.mask(trps, buffer = 30, type = "traprect", spacing = 5)

## c-hat is supplied (analytic c-hat for 'marking' is not yet available);
## value from the RTMB reference model
fit <- secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE,
                details = list(safeLL = TRUE, uselog = TRUE, chat = 4.41))
est <- predict(fit)

## Reference: RTMB model (telemetry_TMB5.R) with Tu likelihood Poisson/chat, chat = 4.41,
## mask-free integration over telemetry: D 32.26 (SE 7.69), lambda0 0.6505, sigma 3.943.
## Coarse mask spacing biases sigma and lambda0 (telemetry density discretised), so
## the tolerances here are loose; D is robust.
test_that("marking model agrees with RTMB reference", {
    expect_equal(est["D", "estimate"], 32.26, tolerance = 0.05, check.attributes = FALSE)
    expect_equal(est["D", "SE.estimate"], 7.69, tolerance = 0.05, check.attributes = FALSE)
    expect_equal(est["sigma", "estimate"], 3.943, tolerance = 0.10, check.attributes = FALSE)
})

## secr self-regression at this mask
test_that("marking model unchanged", {
    expect_equal(est["D", "estimate"], 31.7574, tolerance = 1e-3, check.attributes = FALSE)
    expect_equal(est["lambda0", "estimate"], 0.56940, tolerance = 1e-3, check.attributes = FALSE)
    expect_equal(est["sigma", "estimate"], 4.24709, tolerance = 1e-3, check.attributes = FALSE)
})

## Gauss-Hermite integration over the activity centres of telemetered animals
## does not depend on the mask spacing and matches the RTMB reference closely
fitGH <- secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE,
                  details = list(safeLL = TRUE, uselog = TRUE, chat = 4.41,
                                 telemetryint = "GH"))
estGH <- predict(fitGH)

test_that("telemetryint = 'GH' agrees with RTMB reference", {
    expect_equal(estGH["D", "estimate"], 32.26, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(estGH["lambda0", "estimate"], 0.6505, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(estGH["sigma", "estimate"], 3.943, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(estGH["D", "SE.estimate"], 7.69, tolerance = 0.05, check.attributes = FALSE)
})

test_that("telemetryint input checks", {
    expect_error(secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE,
                          details = list(telemetryint = "none")),
                 "should be 'mask' or 'GH'")
    ch1 <- ch
    telemetrytype(traps(ch1)) <- "concurrent"
    expect_error(secr.fit(ch1, detectfn = "HHN", mask = msk, trace = FALSE,
                          details = list(telemetryint = "GH")),
                 "only for telemetrytype 'marking'")
})

## Unidentified marked cues Tm (pID estimated) and collar windows (atrisk)
## simulated: 15 collars, 5 occasions, q = 0.6, collar active on 60% of occasions,
## Tu = 188, Tm = 8, 14 identified detections.
chq <- readRDS(test_path("telemetry_marking_q.RDS"))
chatq <- c(2.498894, 1, 1)          # RTMB value
detq  <- list(safeLL = TRUE, uselog = TRUE, chat = chatq, telemetryint = "GH")
fitq  <- secr.fit(chq, detectfn = "HHN", mask = msk, trace = FALSE, details = detq)
estq  <- predict(fitq)

## Reference: RTMB model (telemetry_TMB8.R, negative binomial Tu) D 79.25 (SE 17.23),
## lambda0 0.6656, sigma 3.9629, q 0.6369
test_that("Tm, pID and atrisk agree with RTMB reference", {
    expect_equal(estq["D", "estimate"], 79.25, tolerance = 0.01, check.attributes = FALSE)
    expect_equal(estq["lambda0", "estimate"], 0.6656, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(estq["sigma", "estimate"], 3.9629, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(estq["pID", "estimate"], 0.6369, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(estq["pID", "SE.estimate"], 0.1024, tolerance = 0.02, check.attributes = FALSE)
})

test_that("atrisk has an effect, and all-ones atrisk is the same as none", {
    ch1 <- chq
    atrisk(ch1) <- matrix(1, nrow(chq), ncol(chq))
    ch0 <- chq
    atrisk(ch0) <- NULL
    fit1 <- secr.fit(ch1, detectfn = "HHN", mask = msk, trace = FALSE, details = detq)
    fit0 <- secr.fit(ch0, detectfn = "HHN", mask = msk, trace = FALSE, details = detq)
    expect_equal(logLik(fit1), logLik(fit0), tolerance = 1e-6)
    ## ignoring the windows over-states density (here by about two thirds)
    expect_gt(predict(fit0)["D", "estimate"], 1.3 * estq["D", "estimate"])
})

test_that("atrisk attribute checks and subsetting", {
    expect_error(atrisk(chq) <- matrix(1, 3, 3), "same number of animals and occasions")
    expect_error(atrisk(chq) <- matrix(-1, nrow(chq), ncol(chq)), "non-negative")
    sub <- suppressWarnings(subset(chq, 1:5))     # warns of occasions without detections
    expect_equal(dim(atrisk(sub)), c(5, ncol(chq)))
    expect_equal(atrisk(sub), atrisk(chq)[1:5, , drop = FALSE], check.attributes = FALSE)
})

test_that("marking input checks", {
    ch0 <- ch
    Tu(ch0) <- NULL
    expect_error(secr.fit(ch0, detectfn = "HHN", mask = msk, trace = FALSE),
                 "requires counts Tu")
    expect_error(secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE,
                          details = list(nsim = 10)),
                 "not available for telemetrytype 'marking'")
})

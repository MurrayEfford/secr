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

test_that("marking input checks", {
    ch0 <- ch
    Tu(ch0) <- NULL
    expect_error(secr.fit(ch0, detectfn = "HHN", mask = msk, trace = FALSE),
                 "requires counts Tu")
    expect_error(secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE,
                          details = list(nsim = 10)),
                 "not available for telemetrytype 'marking'")
})

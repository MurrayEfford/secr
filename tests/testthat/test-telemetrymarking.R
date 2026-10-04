# test telemetry type 'marking': collared animals (telemetry plus camera cues) and
# unmarked camera cues Tu; secr 5.5.1
#
# Most checks evaluate a single log-likelihood at fixed parameter values (fast).
# One model is fitted in the default tests; further fits run when not on CRAN.

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

## ch:  10 collared animals x 100 fixes, 81 count cameras, 5 occasions, unmarked cues Tu = 58
## chq: 15 collars x 50 fixes, 5 occasions, unidentified marked cues Tm = 8 (pID = 0.6),
##      collars working on 60% of occasions (attribute 'marked'), Tu = 188
ch  <- readRDS(test_path("telemetry_marking_ch.RDS"))
chq <- readRDS(test_path("telemetry_marking_q.RDS"))
trps <- make.grid(detector = "count", nx = 9, ny = 9, spacing = 90/8, origin = c(30, 30))
msk <- make.mask(trps, buffer = 30, type = "traprect", spacing = 8)

## c-hat is supplied (analytic c-hat for 'marking' is not yet available); values from RTMB
chat1 <- 4.41
chatq <- c(2.498894, 1, 1)

## log-likelihood at given parameter values (transformed scale: log D, log lambda0,
## log sigma, logit pID)
LL <- function (x, beta, chat, int = "GH") {
    as.numeric(secr.fit(x, detectfn = "HHN", mask = msk, trace = FALSE, start = beta,
                        details = list(safeLL = TRUE, uselog = TRUE, chat = chat,
                                       telemetryint = int, LLonly = TRUE)))
}

## Estimates from the RTMB reference model (telemetry_TMB5.R with Tu likelihood Poisson/chat
## for ch; telemetry_TMB8.R for chq), which integrates telemetered ACs exactly
beta1 <- log(c(32.261, 0.65045, 3.94255))
betaq <- c(log(c(79.25, 0.6656, 3.9629)), qlogis(0.6369))

test_that("log-likelihood unchanged", {
    expect_equal(LL(ch,  beta1, chat1),         -5717.931939, tolerance = 1e-6)
    expect_equal(LL(ch,  beta1, chat1, "mask"), -5891.899452, tolerance = 1e-6)
    expect_equal(LL(chq, betaq, chatq),         -4364.148674, tolerance = 1e-6)
})

## No 1% change in any parameter increases the log-likelihood at the RTMB estimates,
## i.e. the maximum is within 0.5% of the RTMB values (D and lambda0 are weakly determined)
test_that("log-likelihood maximal at RTMB estimates", {
    check <- function (x, beta, chat) {
        l0 <- LL(x, beta, chat)
        for (j in seq_along(beta)) {
            for (h in c(-0.01, 0.01)) {
                b <- beta
                b[j] <- if (j == 4) qlogis(plogis(beta[j]) * (1 + h)) else beta[j] + h
                expect_gt(l0, LL(x, b, chat))
            }
        }
    }
    check(ch, beta1, chat1)
    check(chq, betaq, chatq)
})

test_that("marked: exposure windows matter, all-ones marked is the same as none", {
    ch1 <- chq
    marked(ch1) <- matrix(1, nrow(chq), ncol(chq))
    ch0 <- chq
    marked(ch0) <- NULL
    expect_equal(LL(ch1, betaq, chatq), LL(ch0, betaq, chatq), tolerance = 1e-8)
    expect_lt(abs(LL(chq, betaq, chatq) - LL(ch0, betaq, chatq)), 100)   # finite
    expect_false(isTRUE(all.equal(LL(chq, betaq, chatq), LL(ch0, betaq, chatq))))
})

test_that("marked attribute checks and subsetting", {
    expect_error(marked(chq) <- matrix(1, 3, 3), "same number of animals and occasions")
    expect_error(marked(chq) <- matrix(-1, nrow(chq), ncol(chq)), "non-negative")
    sub <- suppressWarnings(subset(chq, 1:5))     # warns of occasions without detections
    expect_equal(dim(marked(sub)), c(5, ncol(chq)))
    expect_equal(marked(sub), marked(chq)[1:5, , drop = FALSE], check.attributes = FALSE)
})

test_that("starting values from telemetry and counts", {
    st <- secr:::secr_startmarking(chq, msk)        # true sigma 4, D 53
    expect_equal(st$sigma, 4, tolerance = 0.05)
    expect_equal(st$D, 79, tolerance = 0.2)
})

test_that("marking input checks", {
    ch0 <- ch
    Tu(ch0) <- NULL
    expect_error(secr.fit(ch0, detectfn = "HHN", mask = msk, trace = FALSE),
                 "requires counts Tu")
    expect_error(secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE,
                          details = list(nsim = 10)),
                 "not available for telemetrytype 'marking'")
    expect_error(secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE,
                          details = list(telemetryint = "none")),
                 "should be 'mask' or 'GH'")
    ch1 <- ch
    telemetrytype(traps(ch1)) <- "concurrent"
    expect_error(secr.fit(ch1, detectfn = "HHN", mask = msk, trace = FALSE,
                          details = list(telemetryint = "GH")),
                 "only for telemetrytype 'marking'")
})

## One fit with unidentified marked cues, estimated pID, collar windows and
## Gauss-Hermite integration (about 7 s). Reference: RTMB TMB8 (negative binomial Tu,
## so D differs slightly): D 79.25 (SE 17.23), lambda0 0.6656, sigma 3.9629, pID 0.6369 (SE 0.1024)
test_that("fit with Tm, pID, marked and GH agrees with RTMB reference", {
    fitq <- secr.fit(chq, detectfn = "HHN", mask = msk, trace = FALSE,
                     details = list(safeLL = TRUE, uselog = TRUE, chat = chatq,
                                    telemetryint = "GH"))
    est <- predict(fitq)
    expect_equal(est["D", "estimate"], 79.25, tolerance = 0.01, check.attributes = FALSE)
    expect_equal(est["lambda0", "estimate"], 0.6656, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(est["sigma", "estimate"], 3.9629, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(est["pID", "estimate"], 0.6369, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(est["pID", "SE.estimate"], 0.1024, tolerance = 0.02, check.attributes = FALSE)
})

## Slower checks
test_that("fits of ch: mask summation and GH integration", {
    skip_on_cran()
    ## RTMB (Poisson/chat, exact telemetry integration): D 32.26 (SE 7.69), lambda0 0.6505, sigma 3.943
    ## GH is close to this at any mask spacing; mask summation is not (coarse mask biases sigma, lambda0)
    msk5 <- make.mask(trps, buffer = 30, type = "traprect", spacing = 5)
    fitGH <- secr.fit(ch, detectfn = "HHN", mask = msk5, trace = FALSE,
                      details = list(safeLL = TRUE, uselog = TRUE, chat = chat1, telemetryint = "GH"))
    est <- predict(fitGH)
    expect_equal(est["D", "estimate"], 32.26, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(est["lambda0", "estimate"], 0.6505, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(est["sigma", "estimate"], 3.943, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(est["D", "SE.estimate"], 7.69, tolerance = 0.05, check.attributes = FALSE)
    fit <- secr.fit(ch, detectfn = "HHN", mask = msk5, trace = FALSE,
                    details = list(safeLL = TRUE, uselog = TRUE, chat = chat1))
    est <- predict(fit)
    expect_equal(est["D", "estimate"], 32.26, tolerance = 0.05, check.attributes = FALSE)
    expect_equal(est["sigma", "estimate"], 3.943, tolerance = 0.10, check.attributes = FALSE)
})

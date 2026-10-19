# details$userDfn: session-specific density, with and without session covariates
# (Borchers 2026-10-15: userDfn was always called with nsession = 1)
# Log-likelihoods are evaluated at fixed parameter values (LLonly), so the tests are fast

library(secr)

test_that("userDfn reproduces D ~ year across sessions", {
    msk <- make.mask(traps(ovenCH), buffer = 300, nx = 16, type = "trapbuffer")
    sesscov <- data.frame(year = 0:4)
    st <- c(-1, 0.1, -0.85, 3.4)    # D intercept, year, logit g0, log sigma
    LL <- function(...) {
        as.numeric(secr.fit(ovenCH, model = list(D ~ year, g0 ~ 1, sigma ~ 1), mask = msk,
                            sessioncov = sesscov, start = st, trace = FALSE,
                            details = list(LLonly = TRUE, ...)))
    }
    ref <- LL()

    ## documented form: array [nmask, ngroup, nsession]; no sessioncov argument
    make <- function(sc) function(Dbeta, mask, ngroup, nsession) {
        if (is.character(Dbeta)) {
            if (Dbeta == "name") return("loglinear")
            if (Dbeta == "parameters") return(c("b0", "b1"))
        }
        Dvec <- exp(Dbeta[1] + Dbeta[2] * sc$year)
        array(rep(Dvec, each = nrow(mask) * ngroup), dim = c(nrow(mask), ngroup, nsession))
    }
    expect_equal(LL(userDfn = make(sesscov)), ref, tolerance = 1e-8)

    ## sessioncov supplied by secr.fit, including automatic 'Session'
    f2 <- function(Dbeta, mask, ngroup, nsession, sessioncov) {
        if (is.character(Dbeta)) {
            if (Dbeta == "name") return("loglinear")
            if (Dbeta == "parameters") return(c("b0", "b1"))
        }
        ## one row per session in this call; 'Session' agrees with 'year' (= Session here)
        stopifnot(nrow(sessioncov) == nsession, all(sessioncov$Session == sessioncov$year),
                  is.factor(sessioncov$session), nlevels(sessioncov$session) == 5,
                  as.integer(sessioncov$session) == sessioncov$year + 1)
        Dvec <- exp(Dbeta[1] + Dbeta[2] * sessioncov$year)
        array(rep(Dvec, each = nrow(mask) * ngroup), dim = c(nrow(mask), ngroup, nsession))
    }
    expect_equal(LL(userDfn = f2), ref, tolerance = 1e-8)

    ## session-specific masks that differ: one call per session
    msk2 <- msk
    msk2[[3]] <- make.mask(traps(ovenCH)[[3]], buffer = 250, nx = 14, type = "trapbuffer")
    st2 <- st
    ref2 <- as.numeric(secr.fit(ovenCH, model = list(D ~ year, g0 ~ 1, sigma ~ 1), mask = msk2,
                                sessioncov = sesscov, start = st, trace = FALSE,
                                details = list(LLonly = TRUE)))
    LL2 <- as.numeric(secr.fit(ovenCH, model = list(D ~ year, g0 ~ 1, sigma ~ 1), mask = msk2,
                               sessioncov = sesscov, start = st, trace = FALSE,
                               details = list(LLonly = TRUE, userDfn = f2)))
    expect_equal(LL2, ref2, tolerance = 1e-8)
})

test_that("legacy userDfn returning a single session still works; bad dimensions fail", {
    msk <- make.mask(traps(ovenCH), buffer = 300, nx = 16, type = "trapbuffer")
    fit <- function(f) secr.fit(ovenCH, model = list(D ~ 1, g0 ~ 1, sigma ~ 1), mask = msk,
                                start = c(-1, -0.85, 3.4), trace = FALSE,
                                details = list(LLonly = TRUE, userDfn = f))
    ref <- as.numeric(secr.fit(ovenCH, model = list(D ~ 1, g0 ~ 1, sigma ~ 1), mask = msk,
                               start = c(-1, -0.85, 3.4), trace = FALSE,
                               details = list(LLonly = TRUE)))
    mk <- function(dimfn) function(Dbeta, mask, ngroup, nsession) {
        if (is.character(Dbeta)) {
            if (Dbeta == "name") return("constant")
            if (Dbeta == "parameters") return("b0")
        }
        dimfn(exp(Dbeta[1]), mask, ngroup, nsession)
    }
    old <- mk(function(D, mask, ngroup, nsession) array(D, dim = c(nrow(mask), ngroup, 1)))
    expect_equal(as.numeric(fit(old)), ref, tolerance = 1e-8)
    bad <- mk(function(D, mask, ngroup, nsession) rep(D, 7))
    expect_error(fit(bad), "expected dim")
})

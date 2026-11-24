# start = fitted model must carry the pmix coefficient(s) of a finite mixture
# (bug: secr_mapbeta matched real parameters by names(beta), and the first pmix coefficient
#  is named 'pmix.h22', so pmix was reset to 0)

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

test_that("pmix is transferred when start is a fitted model", {
    skip_on_cran()
    msk <- make.mask(traps(captdata), buffer = 100, nx = 32)
    fit <- suppressWarnings(
        secr.fit(captdata, mask = msk, model = list(g0 ~ h2), trace = FALSE))
    LLstart <- function (...) {
        as.numeric(suppressWarnings(secr.fit(captdata, mask = msk, trace = FALSE,
                                  details = list(LLonly = TRUE), ...)))
    }
    ## same model: the log-likelihood at the start is the maximised log-likelihood
    expect_equal(LLstart(model = list(g0 ~ h2), start = fit), as.numeric(logLik(fit)),
                 tolerance = 1e-6)
    ## more complex model: the old coefficients (including pmix) are kept, new ones are zero
    expect_equal(LLstart(model = list(g0 ~ h2, sigma ~ h2), start = fit),
                 as.numeric(logLik(fit)), tolerance = 1e-6)
    ## fixed pmix fit as the start of a free-pmix model
    fitfx <- suppressWarnings(
        secr.fit(captdata, mask = msk, model = list(g0 ~ h2), fixed = list(pmix = 0.3),
                 trace = FALSE))
    expect_equal(LLstart(model = list(g0 ~ h2), start = fitfx), as.numeric(logLik(fitfx)),
                 tolerance = 1e-6)
})

test_that("mixture model started from a fit without a mixture leaves the stationary point", {
    skip_on_cran()
    ## at the maximum of the one-class model the likelihood is stationary in the class effects;
    ## makeStart supplies default proportions and small unequal class effects
    msk <- make.mask(traps(captdata), buffer = 100, nx = 32)
    fit0  <- suppressWarnings(secr.fit(captdata, mask = msk, trace = FALSE))
    fresh <- suppressWarnings(secr.fit(captdata, mask = msk, model = list(g0 ~ h2), trace = FALSE))
    fit   <- suppressWarnings(secr.fit(captdata, mask = msk, model = list(g0 ~ h2), start = fit0,
                                       trace = FALSE))
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(fresh)), tolerance = 1e-4)
    expect_true(as.numeric(logLik(fit)) > as.numeric(logLik(fit0)))
    expect_false(anyNA(coef(fit)$SE.beta))
})

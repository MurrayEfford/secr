# test fixed pmix in finite mixture models: the value is the proportion in the
# second latent class, so that the first class has 1 - pmix (secr 5.5.1)

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

detector(traps(captdata)) <- "multi"  # to dodge "single" warning
msk <- make.mask(traps(captdata), buffer = 100, type = "trapbuffer", nx = 24)
mod <- list(g0 ~ h2, sigma ~ h2)

## log-likelihood at given beta (free pmix: D, g0, g0.h22, sigma, sigma.h22, pmix.h22)
LL <- function (beta, fixed = NULL) {
    ## start has an element for a fixed pmix; its value is ignored (here deliberately wrong).
    ## method = "none" evaluates the model at start without maximisation
    if (!is.null(fixed$pmix)) beta <- c(beta, 5)
    as.numeric(logLik(secr.fit(captdata, mask = msk, model = mod, fixed = fixed, trace = FALSE, 
                               start = beta, method = "none", details = list(hessian = FALSE))))
}
beta <- c(log(5.5), qlogis(0.3), qlogis(0.2) - qlogis(0.3), log(30), log(40) - log(30))

test_that("fixed pmix equals a free pmix with the same value", {
    expect_equal(LL(c(beta, qlogis(0.3))), LL(beta, fixed = list(pmix = 0.3)), tolerance = 1e-10)
    expect_equal(LL(c(beta, qlogis(0.8))), LL(beta, fixed = list(pmix = 0.8)), tolerance = 1e-10)
    ## the value refers to the second class, so 0.3 and 0.7 differ
    expect_false(isTRUE(all.equal(LL(beta, fixed = list(pmix = 0.3)),
                                  LL(beta, fixed = list(pmix = 0.7)))))
})

test_that("fixed pmix: input checks", {
    expect_error(secr.fit(captdata, mask = msk, model = mod, fixed = list(pmix = 1.2),
                          trace = FALSE), "between 0 and 1")
    expect_error(secr.fit(captdata, mask = msk, model = list(g0 ~ h3), fixed = list(pmix = 0.4),
                          trace = FALSE), "must equal 0.3333")
})

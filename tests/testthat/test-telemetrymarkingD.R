# test telemetry type 'marking' with a density model D ~ covariate: conditional (shape)
# likelihood for the density coefficients, derived() Dcw, sandwich variance; secr 5.5.1

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

## chd: 20 collared animals x 100 fixes, 81 count cameras (about 10% outages), 10 occasions,
##      density D(x) = D0 exp(0.7 z(x)) with z the standardised x coordinate, Tu = 272
chd  <- readRDS(test_path("telemetry_marking_dx.RDS"))
mskd <- make.mask(traps(chd), buffer = 30, type = "traprect", spacing = 8)
covariates(mskd) <- data.frame(z = as.numeric(scale(mskd$x)))

## overdispersion: total count, Tm (unused here), spread across detectors (RTMB values)
chatd <- c(6.431146, 1, 7.634429)

LLd <- function (beta, chat = chatd) {
    as.numeric(secr.fit(chd, detectfn = "HHN", mask = mskd, model = list(D ~ z), trace = FALSE,
                        start = beta,
                        details = list(safeLL = TRUE, uselog = TRUE, chat = chat,
                                       telemetryint = "GH", LLonly = TRUE)))
}
## RTMB TMB8 (likelihood = "cond"): D0 50.31, gamma 0.5914, lambda0 0.9523, sigma 3.9964
## beta: log D0, gamma, log lambda0, log sigma
betad <- c(log(50.3139), 0.5914, log(0.9523), log(3.9964))

test_that("density model: log-likelihood unchanged", {
    expect_equal(LLd(betad), -11759.338687, tolerance = 1e-6)
    ## the spread of unmarked cues across detectors is what identifies the density coefficient:
    ## scaling that component out (chat[3] very large) changes the log-likelihood
    expect_equal(LLd(betad, chat = c(6.431146, 1, 1e9)), -11606.131226, tolerance = 1e-6)
})

## no 5% change in any parameter increases the log-likelihood, i.e. the maximum is
## within about 2.5% of the RTMB estimates (the mask is coarse here, and RTMB uses a
## negative binomial for the total Tu)
test_that("density model: log-likelihood maximal near RTMB estimates", {
    l0 <- LLd(betad)
    for (j in seq_along(betad)) {
        for (h in c(-0.05, 0.05)) {
            b <- betad
            b[j] <- b[j] + h
            expect_gt(l0, LLd(b))
        }
    }
})

## Analytic overdispersion, details$chatmethod = "analytic", at the RTMB estimates.
## RTMB (fixed number of animals, mask spacing 2.5) at the same parameter values:
## chat for the total 6.471; Pearson dispersion across detectors 7.63
test_that("analytic chat agrees with RTMB", {
    detc <- list(safeLL = TRUE, uselog = TRUE, telemetryint = "GH", distribution = "binomial",
                 chatmethod = "analytic", chatonly = TRUE)
    cs <- secr.fit(chd, detectfn = "HHN", mask = mskd, model = list(D ~ z), trace = FALSE,
                   start = betad, details = detc)
    expect_equal(colnames(cs), c("Tu", "Tm", "shape"))
    expect_equal(cs[1, "Tu"], 6.471, tolerance = 0.03)
    expect_equal(cs[1, "Tm"], 1)
    expect_equal(cs[1, "shape"], 7.63, tolerance = 0.03)
    ## a Poisson number of animals adds variance relative to a fixed number
    detc$distribution <- "poisson"
    csp <- secr.fit(chd, detectfn = "HHN", mask = mskd, model = list(D ~ z), trace = FALSE,
                    start = betad, details = detc)
    expect_gt(csp[1, "Tu"], cs[1, "Tu"])
})

test_that("chatmethod input checks", {
    expect_error(secr.fit(chd, detectfn = "HHN", mask = mskd, model = list(D ~ z), trace = FALSE,
                          details = list(chatmethod = "none")), "should be 'simulate' or 'analytic'")
    ch1 <- chd
    telemetrytype(traps(ch1)) <- "concurrent"
    expect_error(secr.fit(ch1, detectfn = "HHN", mask = mskd, trace = FALSE,
                          details = list(chatmethod = "analytic")),
                 "only for telemetrytype 'marking'")
    expect_error(secr.fit(chd, detectfn = "HHN", mask = mskd, model = list(D ~ z), trace = FALSE,
                          details = list(nsim = 10)), "chatmethod = 'analytic'")
})

test_that("density model: fit, Dcw from derived(), sandwich variance", {
    skip_on_cran()
    fit <- secr.fit(chd, detectfn = "HHN", mask = mskd, model = list(D ~ z), trace = FALSE,
                    details = list(safeLL = TRUE, uselog = TRUE, chat = chatd, telemetryint = "GH"))
    cf <- coef(fit)
    ## RTMB: gamma 0.5914 (SE 0.2629); coarse mask here
    expect_equal(cf["D.z", "beta"], 0.5914, tolerance = 0.05)
    expect_equal(cf["D.z", "SE.beta"], 0.2629, tolerance = 0.05)

    ## Dcw: independent calculation from the fitted surface and the cue-rate weights
    p  <- predict(fit, newdata = data.frame(z = 0))
    lam0 <- p["lambda0", "estimate"]; sg <- p["sigma", "estimate"]
    trp <- traps(chd)
    K1  <- nrow(trp) - 1                       # omit the notional telemetry detector
    S1  <- ncol(usage(trp)) - 1                # and the telemetry occasion
    d2 <- outer(mskd$x, trp$x[1:K1], "-")^2 + outer(mskd$y, trp$y[1:K1], "-")^2
    L  <- as.vector(exp(-d2 / (2 * sg^2)) %*% rowSums(usage(trp)[1:K1, 1:S1])) * lam0
    Dx <- exp(cf["D", "beta"] + cf["D.z", "beta"] * covariates(mskd)$z)
    dcw <- derived(fit)
    expect_equal(rownames(dcw), "Dcw")
    expect_equal(dcw["Dcw", "estimate"], sum(Dx * L) / sum(L), tolerance = 1e-4)
    expect_true(is.finite(dcw["Dcw", "SE.estimate"]) && dcw["Dcw", "SE.estimate"] > 0)
    expect_true(dcw["Dcw", "lcl"] < dcw["Dcw", "estimate"] && dcw["Dcw", "ucl"] > dcw["Dcw", "estimate"])

    ## sandwich: larger than the model-based SE because cues from one animal are correlated
    ## across detectors; RTMB 0.2736
    sw <- secr:::secr_shapeSandwich(fit)
    expect_equal(rownames(sw), "D.z")
    expect_gt(sw$SE.sandwich, sw$SE.model)
    expect_lt(sw$SE.sandwich, 2 * sw$SE.model)
    expect_equal(sw$SE.sandwich, 0.2736, tolerance = 0.1)

    ## second pass with the analytic overdispersion (secr.refit keeps the model and starts at the fit)
    fit2 <- secr.refit(fit, details = list(chatmethod = "analytic"), trace = FALSE)
    expect_equal(unname(fit2$details$chat[1]), 10.4, tolerance = 0.05)   # Poisson N (default); 6.5 for fixed N
    expect_equal(unname(fit2$details$chat[3]), 7.6, tolerance = 0.05)
    expect_equal(coef(fit2)["D.z", "beta"], cf["D.z", "beta"], tolerance = 0.02)
})

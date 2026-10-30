# test finite mixtures with three classes: mlogit standard errors of pmix, fixed pmix as a
# vector of class proportions, and starting values (secr 5.5.2)

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

softmax <- function (eta) exp(eta) / sum(exp(eta))

test_that("mlogit delta method agrees with finite differences", {
    p <- c(0.5, 0.3, 0.2)
    eta <- log(p / p[1])
    V <- matrix(0, 3, 3)
    V[2:3, 2:3] <- matrix(c(0.04, 0.01, 0.01, 0.09), 2)
    ## d logit(p[j]) / d eta[k] by central differences
    J <- sapply(1:3, function(k) {
        h <- 1e-6
        e1 <- eta; e1[k] <- e1[k] + h
        e2 <- eta; e2[k] <- e2[k] - h
        (qlogis(softmax(e1)) - qlogis(softmax(e2))) / (2 * h)
    })
    expect_equal(secr:::secr_mlogitvcv(V, p, factor(c(1, 1, 1)), c(TRUE, FALSE, FALSE)),
                 J %*% V %*% t(J), tolerance = 1e-6)
    ## two latent models (e.g. two sessions)
    p6 <- c(p, 0.6, 0.3, 0.1)
    V6 <- matrix(0, 6, 6)
    V6[2:3, 2:3] <- V[2:3, 2:3]
    V6[5:6, 5:6] <- V[2:3, 2:3]
    V6[2:3, 5:6] <- 0.005
    V6[5:6, 2:3] <- 0.005
    got <- secr:::secr_mlogitvcv(V6, p6, factor(rep(1:2, each = 3)), rep(c(TRUE, FALSE, FALSE), 2))
    expect_equal(got[1:3, 1:3], J %*% V %*% t(J), tolerance = 1e-6)
    expect_equal(got[1:3, 4:6], t(got[4:6, 1:3]), tolerance = 1e-12)
    ## two classes: the reference class has the same standard error as the other class
    got2 <- secr:::secr_mlogitvcv(matrix(c(0, 0, 0, 0.04), 2), c(0.4, 0.6), factor(c(1, 1)), c(TRUE, FALSE))
    expect_equal(sqrt(diag(got2)), c(0.2, 0.2))
})

## known classes (hcov with three levels) so that pmix is well determined
set.seed(1)
ch3 <- captdata
covariates(ch3) <- data.frame(g = factor(sample(c("a", "b", "c"), nrow(ch3), TRUE, prob = c(.5, .3, .2))))
detector(traps(ch3)) <- "multi"
msk3 <- make.mask(traps(ch3), buffer = 100, type = "trapbuffer", nx = 16)
mod3 <- list(sigma ~ h3)

## log-likelihood at beta (D, g0, sigma, sigma.h32, sigma.h33, pmix.h32, pmix.h33)
LL3 <- function (beta, fixed = NULL) {
    as.numeric(logLik(secr.fit(ch3, mask = msk3, hcov = "g", model = mod3, fixed = fixed, trace = FALSE,
                               start = beta, method = "none", details = list(hessian = FALSE))))
}
b5 <- c(log(5), qlogis(0.3), log(30), 0.1, -0.1)

test_that("fixed pmix for three classes: proportions in classes 2 and 3, or in all classes", {
    skip_on_cran()
    pf <- c(0.5, 0.3, 0.2)
    free <- LL3(c(b5, log(pf[-1] / pf[1])))
    expect_equal(LL3(c(b5, 5, 5), fixed = list(pmix = pf[-1])), free, tolerance = 1e-10)
    expect_equal(LL3(c(b5, 5, 5), fixed = list(pmix = pf)), free, tolerance = 1e-10)
    ## a single value equal to 1/3 still applies to every class
    expect_equal(LL3(c(b5, 0, 0)), LL3(b5, fixed = list(pmix = 1/3)), tolerance = 1e-10)
    expect_false(isTRUE(all.equal(free, LL3(c(b5, 5, 5), fixed = list(pmix = c(0.2, 0.3))))))
})

test_that("fixed pmix for three classes: input checks", {
    f <- function (p) secr.fit(ch3, mask = msk3, hcov = "g", model = mod3, fixed = list(pmix = p), trace = FALSE)
    expect_error(f(0.4), "must equal 0.3333")
    expect_error(f(c(0.7, 0.6)), "sum to 1")
    expect_error(f(c(0.1, 0.2, 0.3, 0.4)), "1, 2 or 3 values")
    expect_error(f(c(0.5, 1.2)), "between 0 and 1")
})

test_that("three-class starting proportions are unequal", {
    skip_on_cran()
    f0 <- secr.fit(ch3, mask = msk3, hcov = "g", model = mod3, trace = FALSE, method = "none",
                   details = list(hessian = FALSE))
    b <- coef(f0)[f0$parindx$pmix, "beta"]
    expect_equal(unname(b), log(c(2, 1) / 3))
})

test_that("predicted pmix for three classes has standard errors from the delta method", {
    skip_on_cran()
    fit <- secr.fit(ch3, mask = msk3, hcov = "g", model = mod3, trace = FALSE)
    pr <- predict(fit)
    se <- sapply(1:3, function(j) pr[[j]]["pmix", "SE.estimate"])
    expect_false(anyNA(se))
    ## independent calculation from the coefficients and their variance matrix
    ix <- fit$parindx$pmix
    b <- coef(fit)[ix, "beta"]
    vb <- vcov(fit)[ix, ix]
    eta <- c(0, b)
    V <- matrix(0, 3, 3)
    V[2:3, 2:3] <- vb
    J <- sapply(1:3, function(k) {
        h <- 1e-6
        e1 <- eta; e1[k] <- e1[k] + h
        e2 <- eta; e2[k] <- e2[k] - h
        (qlogis(softmax(e1)) - qlogis(softmax(e2))) / (2 * h)
    })
    ## predict() gives the standard error on the probability scale: p (1 - p) times the logit-scale SE
    p <- softmax(eta)
    expect_equal(se, p * (1 - p) * sqrt(diag(J %*% V %*% t(J))), tolerance = 1e-4)
    expect_equal(sapply(1:3, function(j) pr[[j]]["pmix", "estimate"]), softmax(eta), tolerance = 1e-6)
    ## proportions of the observed class composition
    expect_equal(sum(softmax(eta)), 1)
})

## Three classes with telemetry type 'marking' and the class of every collared animal known.
## Splitting class M of the sex data set into two classes with identical parameters leaves the
## likelihood unchanged apart from the multinomial term for the known classes:
## log L3 - log L2 = n2 log(p2 / (p2 + p3)) + n3 log(p3 / (p2 + p3))
test_that("marking with three known classes reduces to two classes", {
    skip_on_cran()
    chs <- readRDS(test_path("telemetry_marking_sex.RDS"))
    trps <- make.grid(detector = "count", nx = 9, ny = 9, spacing = 90/8, origin = c(30, 30))
    msk <- make.mask(trps, buffer = 30, type = "traprect", spacing = 8)
    sx <- as.character(covariates(chs)$sex)
    iM <- which(sx == "M")
    sx3 <- sx
    sx3[iM[seq(2, length(iM), 2)]] <- "N"
    ch3 <- chs
    covariates(ch3)$sex <- factor(sx3)
    det <- list(safeLL = TRUE, uselog = TRUE, chat = c(4.758513, 1, 1), telemetryint = "GH", LLonly = TRUE)
    ## beta: log D, log lambda0 (F), lambda0 M - F, log sigma (F), sigma M - F, logit pID
    b2 <- c(log(51.33), log(0.9929), log(0.9542/0.9929), log(3.5025), log(4.5877/3.5025), qlogis(0.6521))
    p <- c(0.5, 0.3, 0.2)
    L2 <- as.numeric(secr.fit(chs, detectfn = "HHN", mask = msk, trace = FALSE, hcov = "sex",
                              model = list(sigma ~ h2, lambda0 ~ h2),
                              start = c(b2, qlogis(p[2] + p[3])), details = det))
    L3 <- as.numeric(secr.fit(ch3, detectfn = "HHN", mask = msk, trace = FALSE, hcov = "sex",
                              model = list(sigma ~ h3, lambda0 ~ h3),
                              start = c(b2[1:3], b2[3], b2[4:5], b2[5], b2[6], log(p[-1] / p[1])),
                              details = det))
    n <- table(factor(sx3, levels = c("F", "M", "N")))
    expect_equal(L3 - L2, n[["M"]] * log(p[2] / (p[2] + p[3])) + n[["N"]] * log(p[3] / (p[2] + p[3])),
                 tolerance = 1e-8)
})

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
LL <- function (x, beta, chat, int = "GH", CL = FALSE) {
    as.numeric(secr.fit(x, detectfn = "HHN", mask = msk, trace = FALSE, start = beta, CL = CL,
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
    marked(ch1) <- matrix(1, nrow(chq), ncol(chq) - 1)   # non-telemetry occasions
    ch0 <- chq
    marked(ch0) <- NULL
    expect_equal(LL(ch1, betaq, chatq), LL(ch0, betaq, chatq), tolerance = 1e-8)
    expect_lt(abs(LL(chq, betaq, chatq) - LL(ch0, betaq, chatq)), 100)   # finite
    expect_false(isTRUE(all.equal(LL(chq, betaq, chatq), LL(ch0, betaq, chatq))))
})

test_that("marked attribute checks and subsetting", {
    expect_error(marked(chq) <- matrix(1, 3, 3), "one column per non-telemetry occasion")
    expect_error(marked(chq) <- matrix(-1, nrow(chq), ncol(chq) - 1), "non-negative")
    sub <- suppressWarnings(subset(chq, 1:5))     # warns of occasions without detections
    expect_equal(dim(marked(sub)), c(5, ncol(chq) - 1))
    expect_equal(marked(sub), marked(chq)[1:5, , drop = FALSE], check.attributes = FALSE)
})

test_that("starting values from telemetry and counts", {
    st <- secr:::secr_startmarking(chq, msk)        # true sigma 4, D 53
    expect_equal(st$sigma, 4, tolerance = 0.05)
    expect_equal(st$D, 79, tolerance = 0.2)
})

test_that("addTelemetry sets markocc for type 'marking' unless already set", {
    set.seed(1)
    tr  <- make.grid(nx = 5, ny = 5, spacing = 10, detector = "count")
    pop <- sim.popn(D = 50, core = make.mask(tr, buffer = 20, spacing = 5), buffer = 0)
    dp  <- list(lambda0 = 1, sigma = 5)
    sightCH <- sim.capthist(tr, popn = pop, detectfn = "HHN", detectpar = dp,
                            noccasions = 3, renumber = FALSE)
    teCH <- sim.capthist(make.telemetry(), popn = pop, detectfn = "HHN", detectpar = dp,
                         noccasions = 10, renumber = FALSE)
    expect_null(markocc(traps(sightCH)))
    expect_equal(markocc(traps(addTelemetry(sightCH, teCH, type = "marking"))), c(0, 0, 0, 1))
    markocc(traps(sightCH)) <- c(0, 1, 0)     # an existing markocc is extended, not replaced
    expect_equal(markocc(traps(addTelemetry(sightCH, teCH, type = "marking"))), c(0, 1, 0, 1))
})

test_that("make.marking ignores Occasion 0 rows (no detection) and checks the IDs", {
    tr <- make.grid(nx = 5, ny = 5, spacing = 100, detector = "count")
    usage(tr) <- matrix(1, nrow(tr), 4)
    set.seed(123)
    tel <- data.frame(ID = rep(c("A","B","C"), each = 10), x = rnorm(30, 200, 40),
                      y = rnorm(30, 200, 40))
    det <- data.frame(Session = 1, ID = c("A","A","C"), Occasion = c(1,3,2),
                      Detector = rownames(tr)[c(1,7,12)])
    ch0 <- suppressMessages(make.marking(tr, det, tel))
    expect_equal(rownames(ch0), c("A","B","C"))          # B was never detected
    ## placeholder row for B, as make.capthist would have it
    place <- function (id) data.frame(Session = 1, ID = id, Occasion = 0, Detector = rownames(tr)[1])
    ch1 <- suppressMessages(make.marking(tr, rbind(det, place("B")), tel))
    expect_equal(ch1, ch0)
    ## not allowed for an animal that also has detections, or that has no telemetry
    expect_error(make.marking(tr, rbind(det, place("A")), tel), "also have detections")
    expect_error(make.marking(tr, rbind(det, place("D")), tel), "without telemetry")
    ## columns 5 onwards of detections are individual covariates (as in make.capthist)
    detc <- cbind(rbind(det, place("B")), Sex = c("F","F","M","M"))
    ch2  <- suppressMessages(make.marking(tr, detc, tel))
    expect_equal(as.character(covariates(ch2)$Sex), c("F","M","M"))
    expect_equal(rownames(covariates(ch2)), c("A","B","C"))
})

test_that("addTelemetry warns when type 'marking' has no telemetered animals", {
    set.seed(1)
    tr  <- make.grid(nx = 5, ny = 5, spacing = 10, detector = "count")
    pop <- sim.popn(D = 50, core = make.mask(tr, buffer = 20, spacing = 5), buffer = 0)
    dp  <- list(lambda0 = 1, sigma = 5)
    sightCH <- sim.capthist(tr, popn = pop, detectfn = "HHN", detectpar = dp,
                            noccasions = 3, renumber = FALSE)
    teCH <- sim.capthist(make.telemetry(), popn = pop, detectfn = "HHN", detectpar = dp,
                         noccasions = 10, renumber = FALSE)
    noTe <- suppressWarnings(subset(teCH, rep(FALSE, nrow(teCH)), dropnullCH = FALSE))
    expect_equal(nrow(noTe), 0)
    expect_warning(out <- addTelemetry(sightCH, noTe, type = "marking"), "make.marking")
    expect_identical(out, sightCH)                       # returned unchanged
    expect_no_warning(addTelemetry(sightCH, noTe, type = "concurrent"))
})

test_that("Gauss-Hermite integration is accurate for HEX (few fixes per collar)", {
    ## simulated HEX data, 5 fixes per collar (the hardest case); the default 7 x 7 rule is
    ## compared with a 25 x 25 rule by forcing the number of nodes
    set.seed(1)
    tr  <- make.grid(nx = 9, ny = 9, spacing = 11, detector = "count")
    msk <- make.mask(tr, buffer = 40, spacing = 5, type = "traprect")
    pop <- sim.popn(D = 50, core = msk, buffer = 0)
    collared <- 1:nrow(pop) %in% sample(nrow(pop), 15)
    dp <- list(lambda0 = 1, sigma = 4)
    sightCH <- sim.capthist(tr, popn = subset(pop, collared), detectfn = "HEX", detectpar = dp,
                            noccasions = 5, renumber = FALSE)
    teCH <- sim.capthist(make.telemetry(), popn = subset(pop, collared), detectfn = "HEX",
                         detectpar = dp, noccasions = 5, renumber = FALSE)
    uCH <- sim.capthist(tr, popn = subset(pop, !collared), detectfn = "HEX", detectpar = dp,
                        noccasions = 5, renumber = FALSE)
    Tu <- t(apply(uCH, 2:3, sum))
    ch <- addTelemetry(sightCH, teCH, type = "marking", verify = FALSE)
    ch <- addSightings(ch, unmarked = Tu, verify = FALSE)
    LLgh <- function () {
        as.numeric(secr.fit(ch, mask = msk, detectfn = "HEX", trace = FALSE, start = log(c(50, 1, 4)),
                            details = list(safeLL = TRUE, uselog = TRUE, telemetryint = "GH", LLonly = TRUE)))
    }
    ll7 <- LLgh()
    orig <- secr:::secr_telemGH
    testthat::local_mocked_bindings(secr_telemGH = function(...) orig(..., ngh = 25), .package = "secr")
    ll25 <- LLgh()
    expect_lt(abs(ll7 - ll25), 0.1)
    ## the HEX telemetry density at the nodes is in closed form (sum of distances to the animal's
    ## own fixes); mask summation (telemetryint = "mask") evaluates every fix directly, and
    ## agrees closely when the mask is fine and wide enough (2026-11-03). The absolute
    ## log-likelihood of simulated data is not used because it differs between platforms
    fine <- make.mask(tr, buffer = 90, spacing = 1.5, type = "traprect")
    llfine <- function (int) {
        as.numeric(secr.fit(ch, mask = fine, detectfn = "HEX", trace = FALSE, start = log(c(50, 1, 4)),
                            details = list(safeLL = TRUE, uselog = TRUE, telemetryint = int, LLonly = TRUE)))
    }
    expect_lt(abs(llfine("GH") - llfine("mask")), 0.05)
})

test_that("per-animal mask limit (maxdistance) is keyed by animal ID, not by row number", {
    ## IDs that are not 1..n, and animals with telemetry only (prefix 'T')
    set.seed(1)
    tr   <- make.grid(nx = 5, ny = 5, spacing = 10, detector = "count")
    msk0 <- make.mask(tr, buffer = 20, spacing = 5)
    pop  <- sim.popn(D = 50, core = msk0, buffer = 0)
    dp   <- list(lambda0 = 1, sigma = 5)
    sightCH <- sim.capthist(tr, popn = pop, detectfn = "HHN", detectpar = dp,
                            noccasions = 3, renumber = FALSE)
    teCH <- sim.capthist(make.telemetry(), popn = pop, detectfn = "HHN", detectpar = dp,
                         noccasions = 10, renumber = FALSE)
    ch <- addTelemetry(sightCH, teCH, type = "marking")
    expect_true(any(grepl("^ *T", rownames(ch))))
    expect_false(identical(trimws(rownames(ch)), as.character(seq_len(nrow(ch)))))

    thr <- 12
    mb  <- secr:::secr_maskboolean2(ch, msk0, thr)
    n   <- nrow(ch)
    expect_equal(mb$mask_id, 0:(n - 1))
    expect_length(mb$mask_offsets, n + 1)

    ## expected: mask points within thr of the centroid of the animal's detections and fixes
    trxy <- as.matrix(traps(ch))[seq_len(nrow(tr)), ]
    mxy  <- as.matrix(msk0)
    xyl  <- telemetryxy(ch)
    for (i in seq_len(n)) {
        id  <- rownames(ch)[i]
        cnt <- colSums(unclass(ch)[i, 1:3, seq_len(nrow(tr)), drop = FALSE], dims = 2)  # events per camera
        pts <- rbind(trxy[rep(seq_along(cnt), cnt), , drop = FALSE], as.matrix(xyl[[id]]))
        ctr <- colMeans(pts)
        expected <- as.vector(which((mxy[, 1] - ctr[1])^2 + (mxy[, 2] - ctr[2])^2 <= thr^2))
        got <- mb$mask_indices[(mb$mask_offsets[i] + 1):mb$mask_offsets[i + 1]] + 1
        expect_equal(as.vector(got), expected)
    }
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
    ## binary detectors are not supported
    chp <- ch
    dd <- detector(traps(chp))
    detector(traps(chp)) <- ifelse(dd == "count", "proximity", dd)
    expect_error(secr.fit(chp, detectfn = "HHN", mask = msk, trace = FALSE),
                 "requires count detectors")
})

## Sex: hcov with latent classes known for every marked animal.
## chs: 24 collars (11 F, 13 M) x 60 fixes, 6 occasions, sigma F 3.5 M 4.5, lambda0 F 0.8 M 1.0,
##      pID 0.6, collars working on 70% of occasions, Tu = 182, Tm = 24.
chs <- readRDS(test_path("telemetry_marking_sex.RDS"))
chats <- c(4.758513, 1, 1)
## log-likelihood at beta (log D, log lambda0 F, lambda0 M - F, log sigma F, sigma M - F, logit pID, 
## and, if pmix is free, logit of the proportion of class M)
LLs <- function (beta, fixed = NULL, x = chs, method = NULL) {
    det <- list(safeLL = TRUE, uselog = TRUE, chat = chats, telemetryint = "GH")
    if (is.null(method)) {
        ## single evaluation (fast)
        det$LLonly <- TRUE
        as.numeric(secr.fit(x, detectfn = "HHN", mask = msk, trace = FALSE, hcov = "sex", 
                            fixed = fixed, model = list(sigma ~ h2, lambda0 ~ h2), start = beta, details = det))
    }
    else {
        ## as secr.fit evaluates at start without maximising; start includes a (here ignored) 
        ## element for a fixed pmix; slower than LLonly but applies fixedbeta
        det$hessian <- FALSE
        as.numeric(logLik(secr.fit(x, detectfn = "HHN", mask = msk, trace = FALSE, hcov = "sex", 
                                   fixed = fixed, model = list(sigma ~ h2, lambda0 ~ h2), start = beta, 
                                   method = "none", details = det)))
    }
}
## RTMB TMB8 estimates with population sex ratio fixed at 0.5: D 51.33,
## lambda0 F 0.9929 M 0.9542, sigma F 3.5025 M 4.5877, q 0.6521
## beta: log D, log lambda0 (F), lambda0 M - F, log sigma (F), sigma M - F, logit pID
betas <- c(log(51.33), log(0.9929), log(0.9542/0.9929), log(3.5025), log(4.5877/3.5025), qlogis(0.6521))
betas5 <- c(betas, qlogis(0.5))      # pmix free, evaluated at equal classes

test_that("sex classes: log-likelihood unchanged and maximal at RTMB estimates", {
    ## pmix free but evaluated at equal classes (RTMB comparison uses a population ratio of 0.5)
    expect_equal(LLs(betas5), -8456.278150, tolerance = 1e-6)
    ## pmix at the proportion of class M among the marked animals
    expect_equal(LLs(c(betas, qlogis(13/24))), -8456.216831, tolerance = 1e-6)
    ## no 2% change in any parameter increases the log-likelihood (D differs from RTMB by 0.7%
    ## because RTMB uses a negative binomial for Tu)
    l0 <- LLs(betas5)
    for (j in seq_along(betas)) {
        for (h in c(-0.02, 0.02)) {
            b <- betas5
            b[j] <- if (j == 6) qlogis(plogis(betas[j]) * (1 + h)) else betas[j] + h
            expect_gt(l0, LLs(b))
        }
    }
})

test_that("sex classes: fixed pmix is the proportion in the second class", {
    skip_on_cran()
    ## fixed pmix = 0.4 is a free pmix at logit(0.4): class proportions 0.6 and 0.4
    fixed4 <- LLs(c(betas, 5), fixed = list(pmix = 0.4), method = "none")   # 5 is ignored
    expect_equal(fixed4, LLs(c(betas, qlogis(0.4))), tolerance = 1e-10)
    fixed6 <- LLs(c(betas, 5), fixed = list(pmix = 0.6), method = "none")
    expect_false(isTRUE(all.equal(fixed4, fixed6)))
})

test_that("sex classes: derived() reports Dcw (known classes of the collared animals)", {
    skip_on_cran()
    ## failed with 'number of items to replace is not a multiple of replacement length'
    ## when the number of collared animals exceeded one
    fits <- secr.fit(chs, detectfn = "HHN", mask = msk, trace = FALSE, hcov = "sex",
                     model = list(sigma ~ h2, lambda0 ~ h2), start = betas5,
                     details = list(safeLL = TRUE, uselog = TRUE, chat = chats, telemetryint = "GH"))
    dd <- derived(fits)
    expect_equal(dd["D", "estimate"], predict(fits)[[1]]["D", "estimate"], tolerance = 1e-6)
    expect_true(is.finite(dd["D", "SE.estimate"]))
})

test_that("sex classes: input checks", {
    ch1 <- chs
    covariates(ch1)$sex[1] <- NA
    expect_error(LLs(betas5, x = ch1), "requires the class")
})

## One fit with unidentified marked cues, estimated pID, collar windows and
## Gauss-Hermite integration (about 7 s). Reference: RTMB TMB8 (negative binomial Tu,
## so D differs slightly): D 79.25 (SE 17.23), lambda0 0.6656, sigma 3.9629, pID 0.6369 (SE 0.1024)
test_that("fit with Tm, pID, marked and GH agrees with RTMB reference", {
    skip_on_cran()      # about 4 s; the likelihood anchors above run on CRAN
    fitq <-secr.fit(chq, detectfn = "HHN", mask = msk, trace = FALSE,
                     details = list(safeLL = TRUE, uselog = TRUE, chat = chatq,
                                    telemetryint = "GH"))
    est <- predict(fitq)
    expect_equal(est["D", "estimate"], 79.25, tolerance = 0.01, check.attributes = FALSE)
    expect_equal(est["lambda0", "estimate"], 0.6656, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(est["sigma", "estimate"], 3.9629, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(est["pID", "estimate"], 0.6369, tolerance = 0.005, check.attributes = FALSE)
    expect_equal(est["pID", "SE.estimate"], 0.1024, tolerance = 0.02, check.attributes = FALSE)
    ## derived() reports Dcw, which equals D when density is uniform
    dd <- derived(fitq)
    expect_equal(rownames(dd), c("Lambda", "D"))
    ## variance components: CVD^2 = CVcues^2 + CVLambda^2
    expect_equal(dd["D", "CVcues"]^2 + dd["D", "CVLambda"]^2, dd["D", "CVD"]^2, tolerance = 1e-6)
    expect_true(dd["D", "CVcues"] > 0 && dd["D", "CVLambda"] > 0 && is.finite(dd["Lambda", "SE.estimate"]))
    expect_equal(dd["D", "estimate"], est["D", "estimate"], tolerance = 1e-5)
})

## Conditional likelihood (CL = TRUE): Tu is ignored (it informs only density) but the unidentified
## marked cues Tm, which do not involve density, are retained so that lambda0 and pID can be
## separated (previously the CL variance was NaN). Swapping the Tm counts of two detectors changes
## the full and CL log-likelihoods by the same amount, as the pooled Tu term is untouched.
test_that("CL keeps the unidentified marked cues Tm", {
    tm <- Tm(chq); K <- nrow(tm); S <- ncol(tm)
    full <- which(rowSums(usage(traps(chq))[1:K, 1:S]) == S)       # detectors used on every occasion
    i <- full[1]
    j <- full[which(rowSums(tm)[full] != rowSums(tm)[i])[1]]
    tm[c(i, j), ] <- tm[c(j, i), ]
    chq2 <- chq
    Tm(chq2) <- tm
    bcl <- betaq[-1]                                                 # no density intercept
    dfull <- LL(chq2, betaq, chatq) - LL(chq, betaq, chatq)
    dcl   <- LL(chq2, bcl, chatq, CL = TRUE) - LL(chq, bcl, chatq, CL = TRUE)
    expect_gt(abs(dcl), 0.01)
    expect_equal(dcl, dfull, tolerance = 1e-8)
    ## nothing to scale: overdispersion is 1
    chat <- secr.fit(chq, detectfn = "HHN", mask = msk, trace = FALSE, start = bcl, CL = TRUE,
                     details = list(safeLL = TRUE, uselog = TRUE, telemetryint = "GH",
                                    chatmethod = "analytic", chatonly = TRUE))
    expect_equal(as.numeric(chat), rep(1, 3))
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

## region.N: telemetry 'marking' has no marking occasions (the collared animals are marked from the
## start); this used to fail in subset.capthist. Realised N is not estimable (NA), expected N is.
test_that("region.N works with telemetrytype 'marking'", {
    skip_on_cran()
    fit <- secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE,
                    details = list(safeLL = TRUE, uselog = TRUE, chat = chat1, telemetryint = "GH"))
    expect_no_warning(rn <- region.N(fit))
    Dhat <- predict(fit)["D", "estimate"]
    expect_equal(rn["E.N", "estimate"], Dhat * nrow(msk) * attr(msk, "area"),
                 tolerance = 1e-6, check.attributes = FALSE)
    expect_true(is.finite(rn["E.N", "SE.estimate"]))
    expect_true(all(is.na(unlist(rn["R.N", ]))))
    expect_equal(rn["E.N", "estimate"], region.N(fit, se.N = FALSE), check.attributes = FALSE)
})

## Hidden option details$telemetryproduct (2026-11-04): the former likelihood for the pooled total
## of all sightings, with the collared animals' identified cues used twice. Retained only to show
## the effect of the independence assumption (simulations for the paper); not documented for users.
test_that("hidden option telemetryproduct: product likelihood and simulated c-hat", {
    LLp <- function (prod) {
        as.numeric(secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE, start = beta1,
                            details = list(safeLL = TRUE, uselog = TRUE, chat = chat1, telemetryint = "GH",
                                           telemetryproduct = prod, LLonly = TRUE)))
    }
    lp <- LLp(TRUE)
    expect_true(is.finite(lp))
    expect_false(isTRUE(all.equal(lp, LLp(FALSE))))
    ## simulation of c-hat is allowed only with the option (otherwise an error; see 'marking input checks')
    chatsim <- secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE, start = beta1,
                        details = list(safeLL = TRUE, uselog = TRUE, telemetryint = "GH",
                                       telemetryproduct = TRUE, nsim = 20, chatonly = TRUE))
    expect_true(is.finite(chatsim[1, 1]) && chatsim[1, 1] > 1)
})

## analytic c-hat for the same hidden option: Var(T) / E(T) of the total of all sightings
test_that("hidden option telemetryproduct: analytic c-hat agrees with simulated c-hat", {
    skip_on_cran()
    dets <- list(safeLL = TRUE, uselog = TRUE, telemetryint = "GH", telemetryproduct = TRUE, chatonly = TRUE)
    an  <- secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE, start = beta1,
                    details = c(dets, list(chatmethod = "analytic")))
    sim <- secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE, start = beta1,
                    details = c(dets, list(nsim = 2000)))
    expect_gt(unlist(an)[1], 1)
    expect_equal(unlist(an)[1], unlist(sim)[1], tolerance = 0.1, check.attributes = FALSE)
})

## 2026-11-16 second pass for analytic c-hat: secr.fit evaluates c-hat at the starting values for the
## first maximisation, then refits from the estimates so that the reported c-hat is that at the maximum
test_that("analytic c-hat is evaluated at the estimates", {
    skip_on_cran()
    det <- list(telemetryint = "GH", chatmethod = "analytic")
    bad <- beta1 + c(0, 0.4, 0)                  # deliberately poor start (lambda0 about 50% too high)
    chat.start <- secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE, start = bad,
                           details = c(det, list(chatonly = TRUE)))
    fit1 <- secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE, start = bad, details = det)
    chat.mle <- secr.fit(ch, detectfn = "HHN", mask = msk, trace = FALSE, start = fit1,
                         details = c(det, list(chatonly = TRUE)))
    ## the test discriminates: c-hat at the poor start differs from c-hat at the estimates
    expect_gt(abs(chat.start[1] - chat.mle[1]), 0.05)
    expect_equal(as.numeric(fit1$details$chat), as.numeric(chat.mle), tolerance = 1e-4)
    ## a further explicit refit leaves c-hat and the estimates unchanged
    fit2 <- secr.refit(fit1, details = list(chatmethod = "analytic"), trace = FALSE)
    expect_equal(as.numeric(fit2$details$chat), as.numeric(fit1$details$chat), tolerance = 1e-4)
    expect_equal(coef(fit2)$beta, coef(fit1)$beta, tolerance = 1e-3)
    ## the saved call is the user's call, and the internal flag has been reset
    expect_equal(as.character(fit1$call[[1]]), "secr.fit")
    expect_false(isTRUE(secr:::.localstuff$chatsecondpass))
})

test_that("analytic c-hat second pass with fixed pmix", {
    skip_on_cran()
    ## fixed pmix is held in details$fixedbeta; the second pass must keep it fixed and not warn
    expect_warning(
        fit <- secr.fit(chs, detectfn = "HHN", mask = msk, trace = FALSE, hcov = "sex",
                        model = list(sigma ~ h2, lambda0 ~ h2), fixed = list(pmix = 0.4),
                        start = c(betas, 5), details = list(telemetryint = "GH", chatmethod = "analytic")),
        NA)
    expect_equal(nrow(coef(fit)), length(betas))                 # pmix not estimated
    expect_true(is.finite(fit$details$chat[1]) && fit$details$chat[1] > 1)
    expect_true(all(is.finite(fit$beta.vcv)))
    ## secr.refit keeps the fixed pmix (previously failed: start omitted the fixed coefficient)
    fit2 <- secr.refit(fit, details = list(chatmethod = "analytic"), trace = FALSE)
    expect_equal(nrow(coef(fit2)), length(betas))
    expect_equal(as.numeric(fit2$details$chat), as.numeric(fit$details$chat), tolerance = 1e-4)
    expect_equal(coef(fit2)$beta, coef(fit)$beta, tolerance = 1e-3)
})

test_that("session with no collared animals: D ~ session shares detection parameters", {
    skip_on_cran()
    ## session A has 10 collared animals, session B has none (only unmarked cues Tu); this
    ## used to end in a segfault (no animals in the PIA of session B)
    set.seed(2024)
    tr  <- make.grid(nx = 6, ny = 6, spacing = 15, detector = "count")
    mk  <- make.mask(tr, buffer = 30, spacing = 7.5)
    dp  <- list(lambda0 = 0.5, sigma = 10)
    S   <- 4
    sim1 <- function (D, ncollar) {
        pop <- sim.popn(D = D, core = mk, buffer = 0)
        isc <- seq_len(nrow(pop)) <= ncollar
        ch  <- sim.capthist(tr, popn = pop, detectfn = "HHN", detectpar = dp, noccasions = S,
                            renumber = FALSE)
        chc <- subset(ch, rownames(ch) %in% which(isc), dropnullCH = TRUE)
        chu <- subset(ch, !(rownames(ch) %in% which(isc)), dropnullCH = TRUE)
        det <- tel <- NULL
        if (ncollar > 0) {
            w <- which(chc > 0, arr.ind = TRUE)
            w <- w[rep(seq_len(nrow(w)), chc[w]), , drop = FALSE]
            det <- data.frame(Session = "A", ID = rownames(chc)[w[,1]], Occasion = w[,2],
                              Detector = rownames(tr)[w[,3]])
            tel <- data.frame(Session = "A", ID = rep(as.character(1:ncollar), each = 20),
                              x = rnorm(20 * ncollar, rep(pop$x[1:ncollar], each = 20), dp$sigma),
                              y = rnorm(20 * ncollar, rep(pop$y[1:ncollar], each = 20), dp$sigma))
        }
        list(det = det, tel = tel, Tu = t(apply(chu, 2:3, sum)))
    }
    a <- sim1(40, 10)
    b <- sim1(60, 0)
    trs <- list(A = tr, B = tr)
    class(trs) <- c("traps", "list")
    ch <- make.marking(trs, a$det, a$tel, Tu = list(a$Tu, b$Tu), noccasions = S)
    expect_equal(sapply(ch, nrow), c(A = 10, B = 0))

    det <- list(telemetryint = "GH", chatmethod = "analytic")
    fit <- secr.fit(ch, mask = mk, detectfn = "HHN", model = D ~ session, trace = FALSE,
                    details = det)
    expect_equal(nrow(coef(fit)), 4)                 # D, D.sessionB, lambda0, sigma
    expect_true(all(is.finite(fit$beta.vcv)))
    ## session A alone gives the same estimates (session B adds only its own D)
    fitA <- secr.fit(ch[[1]], mask = mk, detectfn = "HHN", trace = FALSE, details = det)
    expect_equal(coef(fit)$beta[c(1,3,4)], coef(fitA)$beta, tolerance = 1e-3)
    ## D for session B is its total unmarked cues divided by Lambda (cues per animal per unit D)
    dv <- derived(fit)
    expect_equal(dv$B["D", "estimate"], sum(b$Tu) / dv$B["Lambda", "estimate"], tolerance = 1e-3)
    expect_equal(dv$B["D", "estimate"], predict(fit)[[2]]["D", "estimate"], tolerance = 1e-3)
    ## the analytic c-hat for session B has no collared animals and is a plain Tu dispersion
    expect_true(all(is.finite(fit$details$chat)))
    ## at least one session must have collared animals
    chB <- ch[[2]]
    expect_error(suppressWarnings(secr.fit(chB, mask = mk, detectfn = "HHN", trace = FALSE)),
                 "at least one session")
})

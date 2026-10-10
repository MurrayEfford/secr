###############################################################################
## package 'secr'
## derivedDcw.R
## 2026-10-04 Dcw for telemetry type 'marking', reported by derived();
##            sandwich variance for the density coefficients
###############################################################################

## Horvitz-Thompson estimates of density (derived, derivedDcoef, derivedDfit, derivedDsurface,
## and region.N for conditional-likelihood fits) do not apply to a CL fit of telemetry type
## 'marking': collared animals are not a sample of the population by detection, and the unmarked
## sightings that inform density are not used in CL. The full likelihood gives density.
secr_stopCLmarking <- function (object, what) {
    cap <- object$capthist
    if (isTRUE(object$CL) &&
        identical(telemetrytype(traps(if (ms(cap)) cap[[1]] else cap)), "marking"))
        stop (what, " for telemetrytype 'marking' requires the full likelihood")
    invisible(NULL)
}

## Density surface and cue-rate surfaces for a fitted model of telemetry type 'marking',
## as functions of the coefficients (beta, in the form of object$fit$par)
##   D    density per ha at each mask point
##   Lmk  list (one per latent class) of m x K matrices: expected number of cues at each
##        detector, over the sighting occasions (with usage), from an animal centred at
##        each mask point
##   pmix class proportions
##   a    cell area (ha)
secr_markingsurfaces <- function (beta, object, sessnum = 1) {
    fullbeta <- secr_fullbeta(beta, object$details$fixedbeta)
    capthist <- if (ms(object$capthist)) object$capthist[[sessnum]] else object$capthist
    mask     <- if (ms(object$mask)) object$mask[[sessnum]] else object$mask
    trps     <- traps(capthist)
    nmix     <- if (is.null(object$details$nmix)) 1 else object$details$nmix
    s        <- ncol(capthist)
    k        <- secr_getk(trps)
    K        <- if (length(k) > 1) length(k) - 1 else k
    m        <- nrow(mask)

    ## density surface at these coefficients
    object$fit$par <- beta
    D <- as.vector(secr_predictD(object, mask, group = NULL, session = sessnum, parameter = 'D'))

    ## detection parameters, hazard of detection at each detector from each mask point
    PIA0 <- object$design0$PIA[sessnum, 1, 1:s, 1:K, , drop = FALSE]   # 1 x 1 x s x K x nmix
    realparval0 <- secr_makerealparameters (object$design0, fullbeta, object$parindx,
                                            object$link, object$fixed)
    Xrealparval0 <- secr_reparameterize (realparval0, object$detectfn, object$details,
                                         mask, trps, NA, s)
    NElist  <- secr_makeNElist(object, object$mask, group = NULL, sessnum)
    ncores  <- setNumThreads(object$details$ncores)
    grain   <- if (!is.null(ncores) && ncores == 1) 0 else 1
    dettype <- secr_detectorcode(trps, MLonly = TRUE, noccasions = s)
    gkhk <- secr_makegk (dettype, object$detectfn, trps, mask, object$details, sessnum, NElist,
                         matrix(1/m, nrow = m, ncol = 1), numeric(4), Xrealparval0, grain, ncores)
    cc <- nrow(Xrealparval0)
    hk <- array(gkhk$hk, dim = c(cc, K, m))

    ## usage on sighting occasions
    markocc <- markocc(trps)
    sight   <- if (is.null(markocc)) rep(TRUE, s) else markocc < 1
    sight   <- sight & detector(trps) != 'telemetry'
    usge <- usage(trps)
    usge <- if (is.null(usge) || object$details$ignoreusage) matrix(1, nrow = K, ncol = s)
            else as.matrix(usge)[1:K, , drop = FALSE]

    ## class proportions
    pmix <- if (nmix > 1) {
        knownclass <- secr_getknownclass(capthist, nmix, object$hcov)
        attr(secr_getpmix (knownclass, PIA0, Xrealparval0), 'pmix')
    } else 1

    K1 <- K - as.integer(any(detector(trps) == "telemetry"))     # excluding the notional telemetry detector
    Lmk <- lapply(1:nmix, function(x) {
        L <- matrix(0, m, K)
        for (j in which(sight)) {
            for (kk in 1:K) {
                cx <- PIA0[1, 1, j, kk, x]
                if (cx > 0) L[, kk] <- L[, kk] + usge[kk, j] * hk[cx, kk, ]
            }
        }
        L[, 1:K1, drop = FALSE]
    })
    list(D = D, Lmk = Lmk, pmix = pmix, a = secr_getcellsize(mask), K = K1)
}

## Cue-weighted density: the average of the fitted density surface D(x) over
## the mask, weighted by the cue rate L(x) of an animal with activity centre at x
##
##     Dcw = sum_x D(x) L(x) / sum_x L(x)
##
## where L(x) is the expected number of cues, summed over detectors and sighting
## occasions (with usage), of an animal centred at x, averaged over latent classes.
## This is the density where the cameras actually sample ('local density' in the sense
## of Efford & Fletcher 2025). It equals D for a model with constant D, and, unlike the
## mean of D(x) over the whole mask, does not depend on extrapolation to places the
## cameras barely sample.
secr_Dcw <- function (beta, object, sessnum = 1) {
    sf <- secr_markingsurfaces(beta, object, sessnum)
    L  <- Reduce(`+`, Map(function(Lx, p) p * rowSums(Lx), sf$Lmk, sf$pmix))
    sum(sf$D * L) / sum(L)
}

## derived() for telemetry type 'marking': Dcw with SE by the delta method
secr_derivedDcw <- function (object, sessnum = NULL, alpha = 0.05, loginterval = TRUE) {
    if (is.null(sessnum)) sessnum <- 1
    beta <- object$fit$par
    est  <- secr_Dcw(beta, object, sessnum)
    if (is.null(object$beta.vcv) || any(is.na(object$beta.vcv))) {
        se <- NA
    }
    else {
        g  <- nlme::fdHess(beta, secr_Dcw, object, sessnum)$gradient
        se <- sqrt(as.numeric(g %*% object$beta.vcv %*% g))
    }
    z <- qnorm(1 - alpha/2)
    if (loginterval && is.finite(se)) {
        selog <- sqrt(log(1 + (se/est)^2))
        w <- exp(z * selog)
        lcl <- est / w; ucl <- est * w
    }
    else {
        lcl <- est - z * se; ucl <- est + z * se
    }
    matrix(c(est, se, lcl, ucl), nrow = 1,
           dimnames = list('Dcw', c('estimate', 'SE.estimate', 'lcl', 'ucl')))
}

###############################################################################
## Sandwich variance for the density coefficients (telemetry type 'marking')
###############################################################################

## The coefficients of D(x) other than the intercept are estimated from the shape
## of the unmarked cues across detectors, a multinomial with probabilities p_k
## proportional to the population cue rates. The model-based variance treats the
## detector counts as independent, but cues from the same animal fall at several
## detectors, so the counts are positively correlated and the variance is too small.
## The score is s = sum_k Tu_k g_k with g_k = d log p_k / d theta. Under the model
## (compound Poisson) the counts Tu_k have covariance
##    Sigma_kj = delta_kj mu_k + a sum_x D(x) L_k(x) L_j(x)
## so Var(s) = g' Sigma g and Cov(theta) = H^-1 (g' Sigma g) H^-1 with
## H = C sum_k p_k g_k g_k' (C total count). The marked animals are not among the
## unmarked animals, so their expected cue vectors (posterior means given their
## own data) are subtracted from the second-moment part of Sigma. Conditional on the
## other parameters.
## Returns data frame with the coefficient, model-based SE and sandwich SE.
secr_shapeSandwich <- function (object, sessnum = 1) {
    beta  <- object$fit$par
    fb    <- object$details$fixedbeta
    free  <- if (is.null(fb)) seq_along(beta) else which(is.na(fb))     # positions of varying coefficients
    iD    <- object$parindx$D[-1]                  # density coefficients other than the intercept
    if (length(iD) == 0) stop ("model has no density coefficients")
    ib    <- match(iD, free)                       # position in beta

    ## posterior-mean cue vectors of the marked animals at the fitted values
    .localstuff$markedpost <- NULL
    invisible(secr.refit(object, method = "none", trace = FALSE, details = list(hessian = FALSE)))
    post <- .localstuff$markedpost
    if (is.null(post)) stop ("could not evaluate marked animals")
    E <- post[, -ncol(post), drop = FALSE]         # animals x detectors

    sf <- secr_markingsurfaces(beta, object, sessnum)
    K  <- sf$K; a <- sf$a; Dx <- sf$D
    Sig <- 0; mu <- rep(0, K)
    for (x in seq_along(sf$Lmk)) {
        Lx  <- sf$Lmk[[x]]
        Sig <- Sig + sf$pmix[x] * a * crossprod(Lx, Dx * Lx)
        mu  <- mu + sf$pmix[x] * a * colSums(Dx * Lx)
    }
    lp <- function (b) {
        s <- secr_markingsurfaces(b, object, sessnum)
        m <- 0
        for (x in seq_along(s$Lmk)) m <- m + s$pmix[x] * s$a * colSums(s$D * s$Lmk[[x]])
        log(m / sum(m))
    }
    gm <- sapply(ib, function(i) {
        h <- 1e-4; bp <- beta; bm <- beta; bp[i] <- bp[i] + h; bm[i] <- bm[i] - h
        (lp(bp) - lp(bm)) / (2 * h)
    })
    gm <- matrix(gm, ncol = length(ib))
    pk <- mu / sum(mu)
    Tu <- Tu(object$capthist)
    C  <- sum(Tu)
    H  <- C * crossprod(gm, pk * gm)
    Sig_pop <- Sig + diag(mu, K)
    Sig_unm <- Sig - crossprod(E) + diag(pmax(mu - colSums(E), 1e-9 * mu), K)
    V  <- crossprod(gm, Sig_unm %*% gm)
    if (any(diag(V) <= 0)) V <- crossprod(gm, Sig_pop %*% gm)
    Hi <- solve(H)
    cv <- Hi %*% V %*% Hi
    if (any(!is.finite(sqrt(diag(cv))))) {
        V  <- crossprod(gm, Sig_pop %*% gm)
        cv <- Hi %*% V %*% Hi
    }
    nm <- rownames(coef(object))[ib]
    data.frame(coefficient = nm, beta = beta[ib],
               SE.model = sqrt(diag(object$beta.vcv))[ib],
               SE.sandwich = sqrt(diag(cv)), row.names = nm)
}

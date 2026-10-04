###############################################################################
## package 'secr'
## chatmarking.R
## 2026-10-04 analytic overdispersion (c-hat) for telemetry type 'marking'
###############################################################################

## Model-based c-hat for the sightings of unmarked animals, as the deterministic
## counterpart of the simulated c-hat of other sighting models (details$chatmethod =
## "analytic"). Called from secr_generalsecrloglikfn at the parameter values supplied.
##
## Elements returned (the same positions as details$chat):
##  [1] Tu   variance / mean of the total count of unmarked sightings. Cues of one
##           animal fall at several detectors and occasions, so the total is a
##           compound Poisson variable. For animals Poisson distributed with density D(x)
##           and L(x) the expected number of cues from an animal centred at x,
##               Var = mu + J2 (+ Vs),     J2 = a sum_x D(x) L(x)^2,
##           with J2 - J1^2 / Nexp for a fixed number of animals (binomial; J1 = a sum D L,
##           Nexp = a sum D). mu = J1 - sum_i E_i[L] is the expected number of unmarked
##           cues and Vs = sum_i Var_i[L(s_i)] is the posterior variance of the cues of
##           the marked animals, which are subtracted.
##  [2] Tm   1 (unidentified marked cues are treated as Poisson)
##  [3] shape  Pearson chi-square of the unmarked counts across detectors about their
##           expected shape, divided by its degrees of freedom (at least 1); used only
##           when density has coefficients (it scales the conditional likelihood)
secr_chatmarking <- function (hk, PIA, usge, markocc, density, cellsize, pmixpop,
                              ElamK, EL2, Tu, Tumusk, nz, n.distrib) {
    K    <- dim(PIA)[4]
    nmix <- dim(PIA)[5]
    m    <- length(density)
    cc   <- length(hk) / (K * m)
    hk3  <- array(hk, dim = c(cc, K, m))
    k1   <- ncol(ElamK)                      # detectors excluding the notional telemetry detector
    sight <- which(markocc < 1)

    L1 <- numeric(m); L2 <- numeric(m)
    for (x in 1:nmix) {
        Lx <- numeric(m)
        for (j in sight) {
            for (k in 1:k1) {
                cx <- PIA[1, 1, j, k, x]
                if (cx > 0) Lx <- Lx + usge[k, j] * hk3[cx, k, ]
            }
        }
        L1 <- L1 + pmixpop[x] * Lx
        L2 <- L2 + pmixpop[x] * Lx^2
    }
    a    <- cellsize
    J1   <- a * sum(density * L1)
    J2   <- a * sum(density * L2)
    Nexp <- a * sum(density)
    mu   <- J1 - sum(ElamK)
    Vs   <- sum(EL2 - rowSums(ElamK)^2)
    vex  <- if (n.distrib == 1) J2 - J1^2 / Nexp else J2
    chatTu <- if (mu > 0) 1 + max(0, (vex + Vs) / mu) else NA

    shape <- 1
    if (nz > 0) {
        Tuk   <- rowSums(Tu)[seq_len(k1)]
        mupop <- rowSums(Tumusk)[seq_len(k1)]
        ex <- sum(Tuk) * mupop / sum(mupop)
        ok <- ex > 0
        shape <- max(1, sum((Tuk[ok] - ex[ok])^2 / ex[ok]) / max(1, sum(ok) - 1 - nz))
    }
    c(chatTu, 1, shape)
}

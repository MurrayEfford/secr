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

###############################################################################
## Analytic overdispersion for all-sighting mark-resight models (sightmodel 5, 6)
## 
## The deterministic counterpart of the simulation in sightingchat.cpp, which it 
## reproduces in expectation: the unmarked animals are placed in mask cells 
## (multinomial given a Poisson or fixed total; the nmark marked animals are 
## removed), and the sightings in each detector x occasion cell are Poisson, or 
## Bernoulli with probability 1 - exp(-X), given the animals present. Here X_c is the
## sum over unmarked animals of the cumulative hazard h_c(location) at cell c (detector k, 
## occasion s). The estimate for the total T of unmarked sightings is Var(T) divided 
## by the naive variance (mean for counts; np pbar (1 - pbar) for binary cells).
##
## Poisson cells: with n_u the number of unmarked animals, L(x) = sum_c h_c(x),
##    E T = E n_u E L,  Var T = E T + E n_u Var L + Var n_u (E L)^2
## (Var n_u = Poisson total, 0 if fixed) where E, Var of L are over the distribution u(x) of
## unmarked animals.
## Bernoulli cells: E exp(-X_c) = G(sum_x u(x) exp(-h_c(x))) and 
## E exp(-X_c - X_c') = G(sum_x u(x) exp(-h_c(x) - h_c'(x))) where G is the probability 
## generating function of n_u. T = sum_c B_c so Var T = sum_c (M1_c - M2_cc) + sum_cc' M2_cc' 
## - (sum_c M1_c)^2 with M1_c = E exp(-X_c), M2_cc' = E exp(-X_c - X_c').
##
## Tm: marked animals are fixed at their expected distribution (as in the simulation),
## so the count has no extra-Poisson variance: 1 for counts, a ratio of Bernoulli variances 
## for binary cells. Tn: 1. 
###############################################################################
secr_chatsighting <- function (hk, PIA, usge, markocc, binomN, Nm, pimask, nmark, pmix, 
                               pID, n.distrib, sightmodel) {
    K    <- dim(PIA)[4]
    nmix <- dim(PIA)[5]
    m    <- length(Nm)
    cc   <- length(hk) / (K * m)
    hk3  <- array(hk, dim = c(cc, K, m))
    sight <- which(markocc == 0)
    if (any(markocc < 0))
        stop ("analytic c-hat is not available when there are unresolved sightings (markocc = -1)")
    if (any(markocc > 0))
        stop ("analytic c-hat requires all occasions to be sighting occasions (no marking occasions)")
    bn <- unique(binomN[sight])
    if (length(bn) != 1 || !(bn %in% c(0, -1, -2)))
        stop ("analytic c-hat requires count, proximity or multi detectors, the same on every sighting occasion")
    binary <- bn < 0
    
    ## cell hazards h[c, m] and cell indices; class-averaged hazard as in the simulation
    cells <- expand.grid(k = 1:K, s = sight)
    nc <- nrow(cells)
    h  <- matrix(0, nc, m)
    for (x in 1:nmix) {
        for (i in 1:nc) {
            cx <- PIA[1, 1, cells$s[i], cells$k[i], x]
            if (cx > 0) h[i, ] <- h[i, ] + pmix[x] * usge[cells$k[i], cells$s[i]] * hk3[cx, cells$k[i], ]
        }
    }
    ## unmarked animals
    sumNm <- sum(Nm)
    Nu    <- pmax(Nm - nmark * pimask, 0)
    nbar  <- sumNm - nmark
    if (!is.finite(nbar) || nbar <= 0) return(c(NA, NA, 1))
    u     <- Nu / sum(Nu)
    poisN <- n.distrib == 0
    
    if (!binary) {
        L   <- colSums(h)
        EL  <- sum(u * L)
        EL2 <- sum(u * L^2)
        ET  <- nbar * EL
        VT  <- ET + nbar * (EL2 - EL^2) + (if (poisN) sumNm else 0) * EL^2
        chatTu <- if (ET > 0) VT / ET else 1
        chatTm <- 1
    }
    else {
        G <- if (poisN) function(s) exp(sumNm * (s - 1)) * s^(-nmark) else function(s) s^nbar
        act <- rowSums(h) > 0                          # cells with no hazard contribute nothing
        W   <- exp(-h[act, , drop = FALSE])
        M1  <- G(as.vector(W %*% u))
        Wu  <- W * rep(sqrt(u), each = nrow(W))
        M2  <- G(tcrossprod(Wu))
        ET  <- sum(1 - M1)
        VT  <- sum(M1 - diag(M2)) + sum(M2) - sum(M1)^2
        npc <- nc                                       # all cells count towards np
        pbar <- ET / npc
        chatTu <- if (pbar > 0 && pbar < 1) VT / (npc * pbar * (1 - pbar)) else 1
        ## Tm: marked animals at their expected distribution, hazard scaled by (1 - pID)
        chatTm <- 1
        if (sightmodel == 5) {
            q1 <- 1 - sapply(cells$s, function(s) mean(pID[s, ]))
            if (any(q1 > 0)) {
                X1 <- (h %*% (nmark * pimask)) * q1
                p1 <- 1 - exp(-X1)
                pb <- mean(p1)
                if (pb > 0 && pb < 1) chatTm <- sum(p1 * (1 - p1)) / (npc * pb * (1 - pb))
            }
        }
    }
    c(chatTu, chatTm, 1)
}

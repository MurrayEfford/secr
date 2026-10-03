###############################################################################
## package 'secr'
## telemGH.R
## 2026-10-04 Gauss-Hermite integration over the activity centres of telemetered
##            animals (details$telemetryint = "GH"); alternative to the default
##            summation over mask cells (details$telemetryint = "mask")
###############################################################################

## Given an animal's telemetry fixes x_t, the likelihood of the fixes as a function
## of its activity centre s is proportional to a bivariate normal density phi(s) with
## mean = mean fix and sd = sigma / sqrt(number of fixes). The integral over s of
## g(s) * prod_t f(x_t | s) * prior(s) is computed as
##
##    (1/A) sum_j  w_j  g(s_j) * prod_t f(x_t | s_j) / phi(s_j)
##
## with nodes s_j and weights w_j from a product Gauss-Hermite rule centred on the
## mean fix and A the area of the mask. The nodes are appended to the mask as extra
## points belonging to that animal alone, and 'density' of each node is
## w_j / (A phi(s_j)), so that the usual code for histories (simplehistories) runs
## unchanged: prod_t f(x_t | s_j) is evaluated per fix and g(s_j) comes from
## detection functions evaluated at the nodes.
##
## Assumptions (differences from the mask-based method):
##   * the prior for the activity centre is uniform (no D(x), no habitat mask)
##   * nodes may lie outside the mask
##   * Euclidean distances (no userdist)

secr_ghnodes <- function (n) {
    ## probabilists' Gauss-Hermite nodes z and weights w (sum to 1)
    J <- matrix(0, n, n)
    if (n > 1) for (i in 1:(n-1)) J[i, i+1] <- J[i+1, i] <- sqrt(i)
    e <- eigen(J, symmetric = TRUE)
    list(z = e$values, w = e$vectors[1,]^2)
}

## Returns a list of replacements for the mask-based arguments of allhistsimple(),
## or NULL if there are no telemetered animals
secr_telemGH <- function (data, PIA, Xrealparval, detectfn, miscparm, gkhk,
                          pi.density, details, ngh = 7) {
    xy    <- data$xy$xy
    start <- data$xy$start
    nc    <- length(start) - 1
    tel   <- which(diff(start) > 0)      # animals with telemetry fixes
    if (length(tel) == 0) return(NULL)
    if (!is.null(details$userdist))
        stop ("details$telemetryint = 'GH' requires Euclidean distances (no userdist)")
    if (.localstuff$iter == 0 && diff(range(pi.density[,1])) > 1e-8 * mean(pi.density[,1]))
        warning ("details$telemetryint = 'GH' assumes uniform density for the activity centres ",
                 "of telemetered animals")

    mask   <- data$mask
    M      <- nrow(mask)
    K      <- dim(PIA)[4]
    telocc <- which(data$binomNcode == -3)[1]    # a telemetry occasion
    sigcol <- match('sigma', colnames(Xrealparval))
    if (is.na(sigcol)) sigcol <- 2
    Am2    <- M * secr_getcellsize(mask) * 1e4   # mask area in square metres (ha * 1e4)

    gh   <- secr_ghnodes(ngh)
    zz   <- as.matrix(expand.grid(z1 = gh$z, z2 = gh$z))
    wgh  <- as.vector(outer(gh$w, gh$w))
    G    <- nrow(zz)

    nodexy  <- matrix(0, nrow = length(tel) * G, ncol = 2)
    nodedens <- numeric(length(tel) * G)
    for (j in seq_along(tel)) {
        i    <- tel[j]
        fx   <- xy[(start[i] + 1):start[i+1], , drop = FALSE]
        xbar <- colMeans(fx)
        sig  <- Xrealparval[PIA[1, i, telocc, K, 1], sigcol]
        sdx  <- sig / sqrt(nrow(fx))
        rows <- (j-1) * G + 1:G
        nodexy[rows, ] <- cbind(xbar[1] + sdx * zz[,1], xbar[2] + sdx * zz[,2])
        phi  <- exp(-0.5 * (zz[,1]^2 + zz[,2]^2)) / (2 * pi * sdx^2)
        nodedens[rows] <- wgh / (phi * Am2)
    }
    nn <- nrow(nodexy)

    ## detection function at the nodes; node rows follow the mask rows
    dist2 <- edist2cpp(as.matrix(data$traps), nodexy)
    gknode <- makegkPointcpp (
        as.integer(detectfn),
        as.integer(details$grain),
        as.integer(details$ncores),
        as.matrix(Xrealparval),
        as.matrix(dist2),
        as.double(miscparm))
    gkhk2 <- list(gk = c(gkhk$gk, gknode$gk), hk = c(gkhk$hk, gknode$hk))

    pi.density2 <- rbind(pi.density, matrix(nodedens, nrow = nn, ncol = ncol(pi.density)))

    ## each telemetered animal gets its own row of mask indices (its nodes)
    mc      <- data$maskcond
    nrows   <- length(mc$mask_offsets) - 1
    mask_id <- rep_len(mc$mask_id, nc)
    mask_id[tel] <- nrows + seq_along(tel) - 1
    maskcond2 <- list(
        mask_id      = mask_id,
        mask_indices = c(mc$mask_indices, M + seq_len(nn) - 1),
        mask_offsets = c(mc$mask_offsets, mc$mask_offsets[nrows + 1] + G * seq_along(tel)))

    haztemp2 <- secr_gethazard (M + nn, data$binomNcode, nrow(Xrealparval), gkhk2$hk,
                                PIA, data$usge)

    ## telemetry density at the nodes only
    telemhr2 <- gethrcpp(
        as.integer(detectfn),
        as.double(start),
        as.matrix(xy),
        rbind(as.matrix(mask), nodexy),
        as.integer(M + seq_len(nn) - 1),
        as.matrix(Xrealparval))

    list(pi.density = pi.density2, gkhk = gkhk2, haztemp = haztemp2,
         maskcond = maskcond2, telemhr = telemhr2)
}

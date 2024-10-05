###############################################################################
## package 'secr'
## CV.R
## last changed 2014-01-30
###############################################################################

CV <- function (x, p, na.rm = FALSE) {
    if (missing(p)) {
        sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
    }
    else {
        p <- p/sum(p)
        xbar <- sum(x * p)
        sqrt(sum(p * (x-xbar)^2)) / xbar
    }
}

CVa0 <- function (object, ...) {
    pred <- predict(object, ...)
    nmix <- object$details$nmix
    if ((nmix == 1) | (length(pred) != nmix))
        stop ("requires 2- or 3-class mixture")
    param <- row.names(pred[[1]])
    if (!('pmix' %in% param))
        stop ("pmix not found")
    pmix <- sapply(pred, '[[', 'pmix','estimate')
    if ('a0' %in% param) {
        a0 <- sapply(pred, '[[', 'a0','estimate')
    }
    else {
        if (object$detectfn < 9) {
            if (!all(c('g0','sigma') %in% param))
                stop ("requires g0 or lambda0 and sigma")
            g0 <- sapply(pred, '[[', 'g0','estimate')
            lambda0 <- -log(1 - g0)
        }
        else {
            if (!all(c('lambda0','sigma') %in% param))
                stop ("requires lambda0 and sigma")
            lambda0 <- sapply(pred, '[[', 'lambda0','estimate')
        }
        sigma <- sapply(pred, '[[', 'sigma','estimate')
        a0 <- 2 * pi * lambda0 * sigma^2
    }
    round(CV(a0, pmix),8)
}

# can embrace fixed par?
# can embrace h3? yes - as stands
# can extend to esa?
# DOES NOT handle individual covariates
# EXCLUDES poly from 2024-09-19

CVa <- function (object, sessnum = 1, ...) {
    if (ms(object)) {
        capthists <- object$capthist[[sessnum]]
        masks <- object$mask[[sessnum]]
    }
    else {
        capthists <- object$capthist
        masks <- object$mask
    }
    pred <- predict(object, ...)
    param <- row.names(pred[[1]])
    if (!('pmix' %in% param))
        stop ("pmix not found; CVa is for mixture models with fitted pmix parameter")
    nmix <- object$details$nmix
    if ((nmix == 1) | (length(pred) != nmix))
        stop ("requires 2- or 3-class mixture")
    trps <- traps(capthists)
    dettype <- detectorcode(trps)
    if (!all(dettype %in%  c(0,1,2,5,8,13))) {
        stop ("type CVa is available only for point detectors")
    }
    binomN <- getbinomN(binomN, detector(trps))
    dpar <-  detectpar(object, ..., byclass = TRUE)
    a1 <- pdot(masks, trps, object$detectfn, detectpar = dpar[[1]],
               noccasions = ncol(capthists), binomN)
    a2 <- pdot(masks, trps, object$detectfn, detectpar = dpar[[2]],
               noccasions = ncol(capthists), binomN)
    cellsize <- attr(masks, "spacing")^2/10000
    a <- c(a1 = sum(a1), a2 = sum(a2)) * cellsize
    pmix <- sapply(pred, '[[', 'pmix','estimate')
    round(CV(a, pmix),8)
}

#msk <- make.mask(traps, buffer = buffer, ...)

CVpdot <- function (..., conditional = FALSE) {
    trps <- eval(match.call(pdot, expand.dots=T)$traps)
    dettype <- detectorcode(trps)
    if (!all(dettype %in%  c(0,1,2,5,8,13))) {
        stop ("type CVpdot is available only for point detectors")
    }
    p <- pdot (...)
    if (conditional) {
        fx <- p/sum(p)    # distribution of captured animals
        mu <- sum(p * fx)
        CV <- sqrt(sum(p^2 * fx) - mu^2)/mu
    }
    else {
        mu <- sum(p) / length(p)  # average over all X, assuming uniform
        CV <- sqrt(sum(p^2/length(p)) - mu^2)/mu
    }
    c(meanpdot = mu, CVpdot = CV)
    
}

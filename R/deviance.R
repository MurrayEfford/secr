############################################################################################
## package 'secr'
## deviance.R
## last changed 2009 09 04; 2011-11-10 N replaces D; 2024-08-20 fixed N, edited df.residual
############################################################################################

deviance.secr <- function (object, ...) {
    ## sum saturated log likelihood over sessions and groups, assuming independence
    ## object - an secr object, including capthist

    session.deviance <- function (capthist, grp, N) {
        capthist <- matrix(capthist, nrow = nrow(capthist))  # dim (n, SxK)
        groupeddata <- split.data.frame(capthist, grp)
        ngrp <- length(groupeddata)
        count <- function(x) table(make.lookup(x)$index)
        nw <- lapply(groupeddata, count)  ## list, each component a vector
        n <- sapply(nw, sum)              ## vector, each group an element
        if (object$CL)
            LLsat <- sum(sapply(n, lfactorial))
        else
            LLsat <- switch (tolower(object$details$distribution),
                poisson = sum(n * log(n) - n),   ## not 100% sure this applies to inhomog Poisson
                binomial = sum(n * log(n/N) + (N-n) * log((N-n)/N) +
                  lgamma(N+1) - lgamma(N-n+1)))
        nwfn <- function (nwg, ng) - sum(lfactorial(nwg)) + sum(nwg * log(nwg / ng ))
        LLsat + sum(sapply(1:ngrp, function(g) nwfn(nw[[g]], n[g])))
    }

    grps <- group.factor(object$capthist, object$groups)

    if (!object$CL & (object$details$distribution=='binomial')) {
        # 2024-08-20 fix code for N broken since 2011!
        N <- region.N(object)['E.N', 'estimate']
    }

    if (ms(object)) {
        nsession <- length(object$capthist)
        LLsat.allsess <- sum(sapply(1:nsession, function(r) session.deviance(
            object$capthist[[r]], grps[[r]], N[r,])))
    }
    else {
        LLsat.allsess <- session.deviance(object$capthist, grps, N)
    }

    -2*(-object$fit$value - LLsat.allsess)

}

df.residual.secr <- function (object, ...) {
    ## object - an secr object
    np <- length(object$fit$par)
    ndistinct <- function (capthist) { 
        nrow(unique(as.array(capthist)))
    }
    if (ms(object)) {
        sum.ndistinct <- sum(sapply(object$capthist, ndistinct))
    }
    else {
        sum.ndistinct <- ndistinct(object$capthist)
    }
    sum.ndistinct - np
}



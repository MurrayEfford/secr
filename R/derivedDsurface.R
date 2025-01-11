#############################################################################
## package 'secr'
## derived density from conditional (relativeD) models
## 2025-01-01
## 2025-01-03 derivedDbeta0 groups
## 2025-01-11 se.beta0
#############################################################################

# session-specific

onek <- function (beta, object, individuals, sessnum)
    ## object is a fitted secr object 
    ## individuals is vector indexing the subset of animals to be used
    # Return k-hat for given beta
    # Only 1 session
{
    out <- sum(1 / esa (object, sessnum, beta, Dweight = TRUE)[individuals])
    # cat('beta ', beta, ' out ', out, '\n')
    out
}
############################################################################################

kgradient <- function (object, individuals, sessnum, ...)
    ## object is a fitted secr object 
    ## individuals is vector indexing the subset of a to be used
{
    nlme::fdHess(object$fit$par, onek, object, individuals, sessnum, ...)$gradient
}
############################################################################################

derivedDbeta0 <- function (object, sessnum = 1, groups = NULL, se.beta0 = FALSE) {
    if (is.null(object$model$D) || is.null(object$link$D) || !object$CL) {
        warning ("not relative density model")
        return(NULL)
    }
    else {
        capthist <- object$capthist
        if (ms(capthist)) capthist <- capthist[[sessnum]]
        grp <- group.factor(capthist, groups)
        individuals <- split (1:nrow(capthist), grp)
        ngrp <- length(individuals)   ## number of groups
        
        se.derivedk <- function (selection, object, selected.a, sessnum) {
            A <-  masksize(object$mask)
            s2 <- switch (tolower(object$details$distribution),
                          poisson  = sum (1/selected.a^2),
                          binomial = sum (( 1 - selected.a / A) / selected.a^2))
            kgrad  <- kgradient (object, selection, sessnum)
            vark <- kgrad %*% object$beta.vcv %*% kgrad
            sqrt(vark + s2)
        }

        getbeta0 <- function (selection) {
            selected.a <- esa (object, sessnum, Dweight = TRUE)[selection]
            k <- sum(1 / selected.a)
            if (se.beta0) {
                se.k <- se.derivedk (selection, object, selected.a, sessnum)
            }
            else {
                se.k <- NA
            }
            c(k, se.k)
        }

        D <- predictD(object, object$mask, group = NULL, session = sessnum, parameter = 'D')
        D <- matrix(D, ncol = 1)                                # dim m x 1
        
        # NB getbeta0() gives same (estimate, SE.estimate) as derived(object, Dweight=T), 
        # and takes about the same time
        # May later switch to derived()
        
        if ( ngrp > 1) {
            # multiple groups
            kse <- sapply(individuals, getbeta0)
        }
        else {    
            # one group
            kse <- getbeta0(individuals[[1]])
        }
        
        c(
            beta0 = transform(kse[1], object$link$D),
            se.beta0 = se.transform(kse[1], kse[2], object$link$D)
        )
        
    }
}

derivedDsurface <- function (object, mask = NULL, sessnum = NULL, groups = NULL) {
    Dx <- function(object, mask, sessnum, selection) {
        D <- predictD(object, mask, group = NULL, session = sessnum, parameter = 'D')
        cellsize <- getcellsize(mask)
        px <- pxi(object, sessnum = sessnum, X = mask)   # dim N x m
        D <- matrix(D, ncol = 1)                         # dim m x 1
        intDp <- px %*% D * cellsize                     # dim N x 1
        out <- mask
        covariates(out)$D.0 <- D * sum(1/intDp[selection])
        class(out) <- c("Dsurface", "mask", "data.frame")
        out
    }
    if (ms(object) && is.null(sessnum)) {
        if (is.null(mask)) mask <- object$mask
        if (ms(mask)) mask <- mask[[sessnum]]
      out <- mapply(derivedDsurface, sessnum = 1:length(object$capthist), 
                    MoreArgs = list(object = object, mask = mask), 
                    SIMPLIFY = FALSE)   
      class (out) <- c("Dsurface", "mask", "list") 
      out
    }
    else {
        if (!object$details$relativeD) 
            stop ("derivedDsurface is for relativeD models")
        if (is.null(object$details$fixedbeta) || is.na(object$details$fixedbeta[1])) 
            stop ("derivedDsurface fixedbeta[1] expected for relativeD but not found")
        if (!(object$link$D %in% c("log", "identity")))
            warning ("derivedDsurface relativeD requires log or identity link")
        if (is.null(sessnum)) sessnum <- 1
        if (is.null(mask)) mask <- object$mask
        
        capthist <- object$capthist
        if (ms(capthist)) capthist <- capthist[[sessnum]]
        grp <- group.factor(capthist, groups)
        ind <- 1:nrow(capthist)
        if (length(ind)>0)
            individuals <- split (ind, grp)
        else
            individuals <-  split (numeric(0), grp) ## list of empty grp levels
        ngrp <- length(individuals)   ## number of groups
        if ( ngrp > 1)
            out <- mapply (
                Dx, 
                individuals, 
                MoreArgs = list(object = object, mask = mask, 
                                sessnum = sessnum), 
                SIMPLIFY = FALSE)
        else {
            if (ngrp == 1) {
                out <- Dx(object, mask = mask, sessnum = sessnum)
            }
            else {
                out <- NULL
            }
        }
        out
    }
}

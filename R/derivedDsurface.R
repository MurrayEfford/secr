#############################################################################
## package 'secr'
## derived density from conditional (relativeD) models
## 2025-01-01
## 2025-01-03 derivedDbeta0 groups
#############################################################################

# session-specific

derivedDbeta0 <- function (object, sessnum = 1, groups = NULL, se.beta0 = FALSE) {
    if (is.null(object$model$D) || is.null(object$link$D) || !object$CL) {
        warning ("not relative density model")
        return(NULL)
    }
    else {
        cellsize <- getcellsize(object$mask)
        px <- pxi(object, sessnum = sessnum, X = object$mask)   # dim n x m
        capthist <- object$capthist
        if (ms(capthist)) capthist <- capthist[[sessnum]]
                grp <- group.factor(capthist, groups)
        individuals <- split (1:nrow(capthist), grp)
        ngrp <- length(individuals)   ## number of groups
        
        # next: modify to get delta method SE
        pxk <- function (px) {
            px <- as.matrix(px)
            intDp <- px %*% D * cellsize                 # dim n x 1
            sum(1/intDp)
        }

        if (se.beta0) {
            stop("variance of beta0 not yet available")
        }  
        
        D <- predictD(object, object$mask, group = NULL, session = sessnum, parameter = 'D')
        D <- matrix(D, ncol = 1)                                # dim m x 1
        
        if ( ngrp > 1) {
            px <- split(as.data.frame(px), grp)
            k <- sapply(px, pxk)
        }
        else {    
            k <- pxk(px)
        }
        sek <- NA    # placeholder
        c(
            beta0 = transform(k[1], object$link$D),
            se.beta0 = se.transform(k[1], sek[1], object$link$D)
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

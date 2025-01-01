# session-specific
# not initially for groups

Dx <- function(object, mask, sessnum, selection) {
    D <- predictD(object, mask, group = NULL, session = sessnum, parameter = 'D')
    cellsize <- getcellsize(mask)
    px <- pxi(object, sessnum = sessnum, X = mask)              # dim N x m
    D <- matrix(D, ncol = 1)       # dim m x 1
    intDp <- px %*% D * cellsize   # dim N x 1
    out <- mask
    covariates(out)$D.0 <- D * sum(1/intDp[selection])
    class(out) <- c("Dsurface", "mask", "data.frame")
    out
}

derivedIntercept <- function (object, sessnum = 1) {
    D <- predictD(object, object$mask, group = NULL, session = sessnum, parameter = 'D')
    ch <- object$capthist
    n <- if (ms(object)) nrow(ch[[sessnum]]) else nrow(ch)
    msk <- if (ms(object)) object$mask[[sessnum]] else object$mask
    # function from regionN.R
    a <- sumDpdot (
        object   = object, 
        sessnum  = sessnum, 
        mask     = msk, 
        D        = D,
        noneuc   = NULL,
        cellsize = getcellsize(msk), 
        constant = FALSE)
    transform(n/a, object$link$D)
}

# completeDbeta <- function(object, sessnum) {
#     intercept <- derivedIntercept(object, sessnum)
#     object$details$fixedbeta[1] <- intercept
#     if (object$link$D == 'identity') {
#         object$fit$par <- object$fit$par * intercept
#     }
#     object
# }

derivedDsurface <- function (object, mask = NULL, sessnum = NULL, groups = NULL) {
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
        
        # object <- completeDbeta(object, sessnum)
        # # inefficient but reliable
        # predictDsurface(object, mask, parameter = 'D')[[sessnum]]
        
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
        browser()
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

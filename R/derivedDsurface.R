# session-specific
# not initially for groups

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

completeDbeta <- function(object, sessnum) {
    intercept <- derivedIntercept(object, sessnum)
    object$details$fixedbeta[1] <- intercept
    if (object$link$D == 'identity') {
        object$fit$par <- object$fit$par * intercept
    }
    attr(object, 'derivedIntercept') <- intercept
    object
}

derivedDsurface <- function (object, mask = NULL, sessnum = 1) {
    if (!object$details$relativeD) 
        stop ("derivedDsurface is for relativeD models")
    if (is.null(object$details$fixedbeta) || is.na(object$details$fixedbeta[1])) 
        stop ("derivedDsurface fixedbeta[1] expected for relativeD but not found")
    if (ms(object))
        warning ("derivedDsurface is not for multi-session models")
    if (!(object$link$D %in% c("log", "identity")))
        warning ("derivedDsurface relativeD requires log or identity link")
    object <- completeDbeta(object, sessnum)
    out <- predictDsurface(object, mask, parameter = 'D')
    attr(out, "derivedIntercept") <- attr(object, 'derivedIntercept') 
    out
}

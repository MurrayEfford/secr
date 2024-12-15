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
    log(n/a)
}

derivedDsurface <- function (object, mask = NULL, sessnum = 1) {
    if (object$CL || !object$details$relativeD) 
        stop ("derivedDsurface is for relativeD models")
    if (is.null(object$details$fixedbeta) || is.na(object$details$fixedbeta[1])) 
        stop ("derivedDsurface fixedbeta[1] not found")
    if (ms(object))
        warning ("derivedDsurface untested for multi-session models")
    if (!is.null(object$groups))
        stop ("derivedDsurface not enabled for groups")
    if (object$link$D != "log")
        stop ("derivedDsurface: requires log link")
    
    intercept <- derivedIntercept(object, sessnum)
    object$details$fixedbeta[1] <- intercept
    out <- predictDsurface(object, mask, parameter = 'D')
    attr(out, "derivedIntercept") <- intercept
    out
}

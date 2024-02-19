################################################################################
## package 'secr'
## list.secr.fit.R
## 2024-02-19 supercedes par.secr.fit
################################################################################

## keeping it simple
list.secr.fit <- function (..., constant = list(), names = NULL) {
    fits <- mapply(secr.fit, ..., MoreArgs = constant, SIMPLIFY = FALSE)
    nfits <- length(fits)
    defaultnames <- paste0('fit', 1:nfits)
    if (is.null(names)) names <- defaultnames
    else if (length(names) != nfits) {
        warning ("number of names does not equal number of fits")
        names <- defaultnames
    }
    secrlist(fits, names = names)
}

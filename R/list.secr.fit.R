################################################################################
## package 'secr'
## list.secr.fit.R
## 2024-02-19 supercedes par.secr.fit
## 2024-07-11 robust to partial failure
################################################################################

## keeping it simple
list.secr.fit <- function (..., constant = list(), prefix = "fit", names = NULL) {
    fits <- mapply(secr.fit, ..., MoreArgs = constant, SIMPLIFY = FALSE)
    nfits <- length(fits)
    defaultnames <- paste0(prefix, 1:nfits)
    if (is.null(names)) names <- defaultnames
    else if (length(names) != nfits) {
        warning ("number of names does not equal number of fits")
        names <- defaultnames
    }
    ok <- sapply(fits, inherits, 'secr')
    if (any(!ok)) {
        if (all(!ok)) stop ("no valid fits") 
        else warning (sum(!ok), " fits failed and were dropped")
    }
    secrlist(fits[ok], names = names[ok])
}

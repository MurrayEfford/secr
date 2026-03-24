##############################################################################
## package 'secr'
## fillNA.R
## 2026-03-25
##############################################################################

fillNA <- function (mask, cov = NULL, method = "nearest", verbose = TRUE) {
    if (ms(mask)) {
        out <- lapply(mask, fillNA, cov = cov, method = method, verbose = verbose)
        class(out) <- class(mask)
        out
    }
    else {
        method <- match.arg(method)
        if (is.null(covariates(mask))) {
            warning("covariate not found; no change")
            mask
        }
        else {
        if (is.null(cov)) cov <- names(covariates(mask))
        for (cv in cov) {
            na <- is.na(covariates(mask)[, cv])
            if (verbose) {
                if (sum(!na) == 0)
                    warning("all ", cv, " missing")
                cat (sum(na), "NA values of covariate", paste0("'",cv,"'"), 
                     "replaced with nearest non-missing value\n")
            }
        OK <- subset(mask, !na)
        pt <- nearesttrap(mask, OK)
        # use nearest non-missing value
        # this is current value when non-missing
        covariates(mask)[, cv] <- covariates(OK)[,cv][pt]
        }
        mask
        }
    }
}

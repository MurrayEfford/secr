##############################################################################
## package 'secr'
## fillNA.R
## 2026-03-25
##############################################################################

fillNA <- function (mask, cov = NULL, method = c("nearest", "value"), 
                    value = NULL, verbose = TRUE) {
    method <- match.arg(method)
    if (!is.null(value) && method != "value") {
        # warning ("value specified in fillNA; setting method = 'value'")
        method <- "value"
    }
    if (method == "value" && is.null(value)) {
        stop ("fillNA() with method = 'value' requires you to specify a replacement value")
    }
    if (ms(mask)) {
        out <- lapply(mask, fillNA, cov = cov, method = method, value = value, 
                      verbose = verbose)
        class(out) <- class(mask)
        out
    }
    else {
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
                type <- if (method == "nearest") "replaced with nearest value" else paste0("replaced with ", value)
                if (sum(na) == 0) type <- ""
                cat (sum(na), "NA values of covariate", paste0("'",cv,"'"), type, "\n")
            }
            if (method == "value") {
                covariates(mask)[is.na(covariates(mask)[, cv]), cv] <- value
            }
            else {
                OK <- subset(mask, !na)
                pt <- nearesttrap(mask, OK)
                # use nearest non-missing value
                # this is current value when non-missing
                covariates(mask)[, cv] <- covariates(OK)[,cv][pt]
            }
        }
        
        mask
        }
    }
}

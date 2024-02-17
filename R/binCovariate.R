##############################################################################
## package 'secr'
## binCovariate.R
## 2024-02-17
##############################################################################

bin <- function(x, width) {
    width * (trunc(x / width)) + width/2
}

binCovariate <- function (object, covname, width) {
    if (ms(object)) {
        # replace components while retaining attributes (e.g., names)
        object[] <- lapply(object, binCovariate, covname, width)
        object
    }
    else {
        newname <- paste0(covname, width)
        if (is.null(covariates(object)))
            stop ("object has no covariates")
        if ( !covname %in% names(covariates(object)))
            stop ("covariate ", covname, " not found")
        if (!is.numeric(covariates(object)[,covname]))
            stop ("covariate is not numeric")
        if (newname %in% names(covariates(object)))
            stop ("covariate ", newname, " already exists")
        covariates(object)[,newname] <- bin(covariates(object)[,covname], width)
        object
    }
}
##############################################################################


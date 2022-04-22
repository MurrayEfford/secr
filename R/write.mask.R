############################################################################################
## package 'secr'
## write.mask.R
## 2019-07-01, 03
############################################################################################

write.mask <- function (object, file='', header = TRUE, ndec = 2, covariates = TRUE, ...) {

    objectname <- ifelse (is.character(header),
        header, deparse(substitute(object), control=NULL))
    header <- ifelse (is.character(header), TRUE, header)

    if (!is(object, 'mask'))
        stop ("requires a 'mask' object")
    n <- nrow(object)
    object$x <- round(object$x,ndec)
    object$y <- round(object$y,ndec)
    temp <- object
    class(temp) <- "data.frame"  ## 2018-10-02

    covlist <- numeric(0)
    if (!is.null(covariates) & !is.null(covariates(object))) {
        covs <- covariates(object)
        if (is.character(covariates)) {
            covlist <- match(covariates, names(covs))
            covlist <- covlist[!is.na(covlist)]
        }
        else
            covlist <- names(covs)

        if (length(covlist)>0) {
            covnames <- paste(covlist, collapse=' ')
            covs <- covs[, covlist, drop=FALSE]
            temp <- cbind(temp, covs)
        }
    }

    if (header) {
        cat ("# Habitat mask points exported from ", objectname, "\n",
             sep = "", file = file)
        cat ('#', format(Sys.time(), "%a %b %d %X %Y"), '\n', append = TRUE,
             file = file)
    }

    suppressWarnings(write.table(temp, file = file, append = header,
        row.names = FALSE, col.names = header, quote = FALSE, ...))

}
###############################################################################

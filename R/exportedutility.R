################################################################################
## package 'secr'
## exportedutility.R
## 2024-09-16  moved from utility.R
################################################################################

## mean and SD if x numeric
getMeanSD <- function(xy) {
    MeanSD <- function (x) {
        if (is.numeric(x))
            c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
        else
            c(NA,NA)
    }
    as.data.frame (apply(xy, 2, MeanSD))
}

#-------------------------------------------------------------------------------

maskarea <- function (mask, sessnum = 1) {
    if (!ms(mask)) nrow(mask) * attr(mask,'area')
    else nrow(mask[[sessnum]]) * attr(mask[[sessnum]],'area')
}

#-------------------------------------------------------------------------------

masklength <- function (mask, sessnum = 1) {
    if (!ms(mask)) nrow(mask) * attr(mask,'spacing')/1000
    else nrow(mask[[sessnum]]) * attr(mask[[sessnum]],'spacing')/1000
}

#-------------------------------------------------------------------------------

masksize <- function (mask, sessnum = 1) {
    if (inherits(mask, 'linearmask'))
        masklength(mask, sessnum)
    else
        maskarea(mask, sessnum)
}

#-------------------------------------------------------------------------------
edist <- function (xy1, xy2) {
    nr <- nrow(xy1)
    nc <- nrow(xy2)
    x1 <- matrix(xy1[,1], nr, nc)
    x2 <- matrix(xy2[,1], nr, nc, byrow=T)
    y1 <- matrix(xy1[,2], nr, nc)
    y2 <- matrix(xy2[,2], nr, nc, byrow=T)
    sqrt((x1-x2)^2 + (y1-y2)^2)
}

#-------------------------------------------------------------------------------

## least cost paths from mask including barriers to movement
## use edist for equivalent Euclidean distances

nedist <- function (xy1, xy2, mask, inf = Inf, ...) {
    newargs <- list(...)
    if (missing(mask)) mask <- xy2
    noneuc <- covariates(mask)$noneuc
    if (is.null(noneuc)) noneuc <- rep(1, nrow(mask))
    defaultargs <- list(transitionFunction = mean, directions = 16)
    args <- replace(defaultargs, names(newargs), newargs)
    args$x <- raster(mask, values = noneuc)
    if (requireNamespace('gdistance', quietly = TRUE)) {    ## 2015-01-23
        tr <- do.call(gdistance::transition, args)
        tr <- gdistance::geoCorrection(tr, type = "c", multpl = FALSE)
        out <- gdistance::costDistance(tr, as.matrix(xy1), as.matrix(xy2))
    }
    else stop ("package gdistance is required for nedist")
    if (is.finite(inf)) out[!is.finite(out)] <- inf
    out
}

#-------------------------------------------------------------------------------

rlnormCV <- function(n, mean, cv) {
    # n simulated values from log-normal distribution with mean = mean and CV = cv
    sdlog <- log(cv^2 + 1)^0.5
    meanlog <- log(mean) - sdlog^2/2
    rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}
#-------------------------------------------------------------------------------

shareFactorLevels <- function (object, columns = NULL, stringsAsFactors = TRUE) {
    ## stringsAsFactors added 2020-05-16
    if (ms(object)) {
        if (!is.null(covariates(object))) {
            df <- do.call(rbind, covariates(object))
            if (is.null(columns)) {
                columns <- 1:ncol(df)
            }
            if (stringsAsFactors) {
                df[,columns] <- stringsAsFactors(df[,columns, drop = FALSE])
            }
            for (i in columns) {
                if (is.factor(df[,i])) {
                    levelsi <- levels(df[,i])
                    for (sess in 1:length(object)) {
                        covariates(object[[sess]])[,i] <-
                            factor(covariates(object[[sess]])[,i],
                                   levels = levelsi)
                    }
                }
            }
        }
    }
    else {
        # modified 2021-04-27 to apply to covariates, not object itself
        if (!is.null(covariates(object))) {
            if (stringsAsFactors) {
                df <- covariates(object)
                if (is.null(columns)) {
                    columns <- 1:ncol(df)
                }
                df[,columns] <- stringsAsFactors(df[,columns, drop = FALSE])
                covariates(object) <- df
            }
        }
    }
    object
}

#-------------------------------------------------------------------------------


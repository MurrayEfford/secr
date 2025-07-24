############################################################################################
## package 'secr'
## clone.R
## last changed 2016-11-29, 2018-11-16...
## 2021-05-19 sortorder for polygon data
## 2022-08-08 truncated Poisson, revised 'constant'
############################################################################################

# random truncated Poisson
rtpois <- function(n, lambda) {
    qpois(runif(n, dpois(0, lambda), 1), lambda)
}
#-------------------------------------------------------------------------------

clone <- function (object, type, ...) UseMethod("clone")

clone.default <- function (object,  type, ...)       {
    if (length(dim(object)) != 2)
        stop ("requires 2-D object")
    type <- tolower (type)
    n <- nrow(object)
    if (n == 0) {
        out <- object
    }
    else {
        if (type == 'constant') 
            freq <- rep(list(...)[[1]], length.out = n)
        else if (type == 'poisson')
            freq <- rpois(n, ...)
        else if (type == 'truncatedpoisson')
            freq <- rtpois(n, ...)
        else if (type == 'nbinom')
            freq <- rnbinom(n, ...)
        else
            stop("unrecognised type")
        freq <- round(freq)
        index <- rep(1:n, freq)
        object[index,,drop = FALSE]
    }
}

clone.popn <- function (object, type, ...) {
    if (ms(object)) {
        out <- lapply (object, clone, type, ...)
        class (out) <- c('popn','list')
        out
    }
    else {
        type <- tolower (type)
        n <- nrow(object)
        if (n == 0) {
            out <- object
        }
        else {
            if (type == 'constant')
                freq <- rep(list(...)[[1]], length.out = n)
            else if (type == 'poisson')
                freq <- rpois(n, ...)
            else if (type == 'truncatedpoisson')
                freq <- rtpois(n, ...)
            else if (type == 'nbinom')
                freq <- rnbinom(n, ...)
            else
                stop("unrecognised type")
            freq <- round(freq)
            index <- rep(1:n, freq)
            out <- object[index,,drop = FALSE]
            
            ## 2018-11-16
            seqn <- unlist(lapply(freq[freq>0], seq, from = 1))
            seqn <- leadingzero(seqn)
            rown <- paste(rep(rownames(object), freq), seqn, sep='.')
            rownames(out) <- rown
            
            if (!is.null(covariates(object))) {
                covariates(out) <- covariates(object)[index,,drop = FALSE]
                rownames(covariates(out)) <- rown
            }
            attr (out, 'freq') <- freq
        }
        out
    }
}

clone.capthist <- function (object, type, ...) {
    if (ms(object)) {
        out <- lapply (object, clone, type, ...)
        class (out) <- class(object)
        out
    }
    else {
        type <- tolower (type)
        n <- nrow(object)
        if (n == 0) {
            out <- object
        }
        else {
            if (type == 'constant')
                freq <- rep(list(...)[[1]], length.out = n)
            else if (type == 'poisson')
                freq <- rpois(n, ...)
            else if (type == 'truncatedpoisson')
                freq <- rtpois(n, ...)
            else if (type == 'nbinom')
                freq <- rnbinom(n, ...)
            else
                stop("unrecognised type")
            freq <- round(freq)
            index <- rep(1:n, freq)
            if (length(dim(object))==2)
                out <- object[index,,drop = FALSE]
            else
                out <- object[index,,,drop = FALSE]
            seqn <- unlist(lapply(freq[freq>0], seq, from = 1))
            seqn <- leadingzero(seqn)
            rown <- paste(rep(rownames(object), freq), seqn, sep='.')
            rownames(out) <- rown
            traps(out) <- traps(object)
            attr(out, 'cutval') <- attr(object, 'cutval')
            if (!is.null(covariates(object))) {
                covariates(out) <- covariates(object)[index,,drop = FALSE]
                rownames(covariates(out)) <- rown
            }
            
            if (!is.null(xy(object)) | !is.null(signalframe(object))) {
                ## these attributes are defined for each detection
                ## so we first extract the numeric animal ID for each detection
                detectionindex <- animalID(object, names = FALSE, sortorder = 'ksn')
                xy(out) <- xy(object)[index[detectionindex]]
                signalframe(out) <- signalframe(object)[index[detectionindex]]
            }
            attr (out, 'freq') <- freq
            class(out) <- class(object)
        }
        
        out
    }
}

#############################################################################
## package 'secr'
## RMarkInput.R

## 2017-03-25 RMarkInput fixed for 3D capthist only
## 2018-01-22 transferred from join.R
## 2024-02-25 unRMarkInput accepts non-numeric '.' for missing
#############################################################################

RMarkInput <- function (object, grouped = FALSE, covariates = TRUE) {
    if (!inherits(object, "capthist"))
        stop ("requires single-session capthist object")
    if (ms(object))
        stop ("requires single-session capthist object - use 'join'")
    # object <- check3D(object)
    CH <- apply(object, 1:2, function(x) as.numeric(any(abs(x)>0)))
    ntimes <- ncol(object)
    alive <- apply(object,1,function(x) all(x>=0))
    
    if (is.logical(covariates)) {
        if (covariates & !is.null(covariates(object))) {
            # if (is.null(covariates(object)))
            #     stop("no covariates in object")
            covnames <- names(covariates(object))
        }
        else
            covnames <- ""
    }
    else {
        covnames <- covariates
        if (is.character(covariates)) {
            if (is.null(covariates(object)))
                stop("no covariates in object")
        }
        found <- covnames %in% names(covariates(object))
        if (any(!found)) {
            stop(paste(covnames[!found], collapse=','), " not in covariates(object)")
        }
    }
    
    if (any(covnames != "")) {
        if (grouped) {
            warning("'grouped' is incompatible with individual covariates and will be ignored")
            grouped <- FALSE
        }
    }
    
    if (grouped)   ## bug fix 2012-07-04
        CH <- cbind(CH, alive) ## add single-digit code as last column
    CH <- data.frame(ch = apply(CH, 1, paste, collapse=''),
                     stringsAsFactors = FALSE)
    
    if (grouped) {
        temp <- table(CH$ch)
        alive <- as.numeric(substring(names(temp),ntimes+1,ntimes+1))
        CH <- data.frame(ch = substring(names(temp),1,ntimes),
                         freq = as.numeric(temp),
                         stringsAsFactors = FALSE)
        CH$freq <- ifelse(alive, CH$freq, -CH$freq)
        CH <- CH[order(CH$ch, decreasing = TRUE),]
        row.names(CH) <- 1:nrow(CH)
    }
    else {
        CH$freq <- ifelse(alive,1,-1)
        if (any(covnames != "")) {
            CH[,covnames] <- covariates(object)[,covnames]
        }
    }
    attr(CH, "intervals") <- attr(object, "intervals")
    if (is.null(attr(CH,"intervals")))
        attr(CH, "intervals") <- rep(0,ntimes-1)
    CH
}

unRMarkInput <- function(df, covariates = TRUE) {
    if (!is.data.frame(df))
        stop("requires dataframe input")
    if (!('ch' %in% names(df)))
        stop ("ch is a required field")
    if (!('freq' %in% names(df))) {
        warning ("field 'freq' not found; assuming all freq=1")
        df$freq <- rep(1, nrow(df))
    }
    nocc <- nchar(df$ch)
    if (length(unique(nocc))>1)
        stop ("ch must be a constant-length string of 0s and 1s")
    nocc <- nocc[1]
    freq <- df$freq
    alive <- sign(freq)
    freq <- abs(freq)
    freq <- rep(1:nrow(df), freq)
    alive <- alive[freq]
    ch <- df$ch[freq]
    w <- suppressWarnings(as.numeric(unlist(sapply(ch, strsplit, ''))))
    CH <- matrix(w, byrow = TRUE, ncol = nocc)
    # allow deads
    last <- function(x) {
        x[is.na(x)] <- 0
        which.max(cumsum(x))
    }
    CH[cbind(1:nrow(CH), apply(CH,1,last))] <- alive
    CH <- array(CH, dim=c(dim(CH),1))   # for version 3
    class(CH) <- 'capthist'
    if (any (is.na(CH))) warning("missing values set to NA")
    # transfer covariates if present
    if (ncol(df)>2) {
        if (is.logical(covariates)) {
            if (covariates)
                covnames <- names(df)[-match(c('ch','freq'),names(df))]
            else
                covnames <- ""
        }
        else {
            covnames <- covariates[covariates %in% names(df)]
        }
        
        if (any(covnames != ""))
            covariates(CH) <- df[freq, covnames, drop = FALSE]
    }
    CH
}

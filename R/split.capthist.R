############################################################################################
## package 'secr'
## split.capthist.R
## last changed 2009 06 11 2009 07 10 2009 10 05 2012 07 26 2012 09 04 2015-10-11
## 2021-04-24, 2023-11-23, 2024-11-08
############################################################################################

split.capthist <- function (x, f, drop = FALSE, prefix='S', bytrap = FALSE,
    byoccasion = FALSE, bysession = FALSE, ...) {
    if (!inherits(x, 'capthist'))
        stop ("argument to 'split.capthist' should have class 'capthist'")
    if (any(is.na(unlist(f)))) {
        stop ("f should not contain NA values")
    }
    if (ms(x)) {
        if (bysession) {
            if (length(f)!= length(x))
                stop ("length of f should match number of sessions")
            class(x) <- 'list'
            out <- split(x, f)
            sess <- split(session(x), f)
            for (i in 1:length(out)) {
                class(out[[i]]) <- c('capthist','list')
                session(out[[i]]) <- sess[[i]]
            }
            out
        }
        # multi-session 'capthist' added 2021-04-24
        else if (inherits(f, 'list')) {
                # make list of multisession capthists
                out <- mapply(split, x, f, MoreArgs = list(drop=drop, prefix=prefix, 
                                                           bytrap = bytrap, byoccasion=byoccasion, ...), SIMPLIFY = FALSE)
                # rearrange
                out <- unlist(out, recursive = FALSE)
                i <- unlist(sapply(f, levels))
                if (!is.matrix(i)) stop ("levels of f should be the same over sessions")
                out <- split(out, row(i))
                names(out) <- i[,1]
                reform <- function(y) {
                    y <- c(y) # combine component sessions of each split
                    class(y) <- c('capthist','list')
                    session(y) <- session(x)
                    y
                }
                lapply(out, reform)   
        }
        else {
            stop ("split of multisession capthist requires multisession (list-valued) f")
        }
    }
    else {
        if (bysession) stop ("cannot split single session capthist by session")
        oldopt <- options(warn=-1)
        
        f <- as.factor(f)  # retains unused levels
        if (any(!is.na(as.numeric(levels(f))))) {
            levels(f) <- paste (prefix, secr_leadingzero(levels(f)),sep='')
        }
        options(oldopt)
        
        if (bytrap) {
            ## 2015-10-11
            ## if (length(f)!= nrow(traps(x)))
            if (length(f)!= secr_ndetector(traps(x)))
                stop ("length of f should match number of detectors")
        }
        else if (byoccasion) {
            if (length(f)!=ncol(x))
                stop ("length of f should match number of columns in capthist")
        }
        else {
            if (length(f)!=nrow(x))
                stop ("length of f should match number of rows in capthist")
        }
        if (bytrap & byoccasion)
            stop("specify only one of bytrap and byoccasion")
        
        out <- list()
        for (i in levels(f)) {
            if (bytrap) {
                temp <- subset (x, traps = f == i, ...)
            }
            else if (byoccasion) {
                temp <- subset (x, occasions = f == i, ...)
            }
            else {
                temp <- subset (x, subset = f == i, ...)
            }
            session(temp) <- i
            if (!drop | (nrow(temp)>0))
                out[[i]] <- temp
        }
        class (out) <- c('capthist', 'list')
        out
    }
}
############################################################################################

extract.estimates <- function (x, simplify = FALSE) {
## compile a dataframe (simplify = T) or list of data.frames of session-specific real parameter estimates
## from a list of separate secr model fits
   if (!is.list(x) | !inherits(x[[1]], 'secr'))
       stop ("requires list of fitted secr models")
   temp <- lapply(x, predict)
   temp <- lapply(temp, function(x) x[,-1])  ## drop unwanted 'link' column
   temp <- lapply(temp, function(x) {x$Parameter <- row.names(x); x})
   sessions <- names(temp)
   nsessions <- length(temp)
   parnames <- row.names(temp[[1]])
   nrealpar <- nrow(temp[[1]])
   temp2 <- data.frame(abind(temp, along = 1), row.names = NULL, 
     stringsAsFactors = FALSE)
   temp2[,1:4] <- sapply(temp2[,1:4], as.numeric)
   temp2$Session <- rep(sessions, rep(nrealpar, nsessions))
   if (simplify) {
       temp3 <- temp2[order(temp2$Parameter, temp2$Session), c('Parameter','Session','estimate', 'SE.estimate', 'lcl', 'ucl')]
       row.names(temp3) <- NULL
   }
   else {
      temp3 <- split(temp2[order(temp2$Session), c('Session','estimate', 'SE.estimate', 'lcl', 'ucl')],
          temp2$Parameter)
      temp3 <- lapply(temp3, function(x) {row.names(x) <- NULL; x})
   }
   temp3
}
############################################################################################


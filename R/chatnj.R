###############################################################################
## package 'secr'
## chatnj.R
## 2025-05-20
###############################################################################

chat.nj <- function(object, bycluster = FALSE, ...) {
    ##-------------------------------------------------------------
    ## Restrict application
    if (!inherits(object, 'secr'))
        stop ("chat.nj expects fitted secr model")
    if (length(object$fixed)>0)
        stop ("chat.nj does not yet accept models with one or more fixed parameters")
    if (!is.null(object$groups))
        stop ("chat.nj does not yet accept models with grouping")
    # if (length(table(object$design0$PIA))>1 || length(table(object$design$PIA))>1)
    #     stop ("chat.nj does not yet handle varying detection probabilities")
    
    ##-------------------------------------------------------------
    
    if (bycluster) {
        if (ms(object))
            stop ("chat.nj expects single-session input when bycluster = TRUE")
        observed.nj <- cluster.counts(object$capthist)
    }
    else {
        if (!ms(object))
            stop ("chat.nj expects multi-session input when bycluster = FALSE")
        observed.nj <- sapply(object$capthist, nrow)
    }
    expected.nj <- expected.n(object, bycluster = bycluster)  # by session or cluster
    expected.nj <- unlist(expected.nj)
    np <- length(object$parindx$D)
    chat <- Fletcher.chat(
        observed    = observed.nj, 
        expected    = expected.nj, 
        np          = np, 
        type        = 'Fletcher', 
        multinomial = FALSE,
        ...)   # default verbose = TRUE
    chat
}


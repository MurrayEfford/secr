############################################################################################
## package 'secr'
## detectpar.R
## 2013-11-09 byclass option 
## 2024-08-05 pmix 
## 2024-09-07 bytrap option 
## 2024-10-03 moved from methods.R
############################################################################################

detectpar.default <- function(object, ...) {
    stop ("only for secr models")
}

detectpar.secr <- function(object, ..., byclass = FALSE, bytrap = FALSE) {
    extractpar <- function (temp, ntrap, initial = TRUE) {
        ## for one session
        if (!is.data.frame(temp))   ## assume list
        {
            if (byclass || bytrap) {
                if (initial && nclass>1 && bytrap) {
                    temp <- unlist(temp, recursive = FALSE)
                }
                out <- lapply(temp, extractpar, FALSE)
                
                if (bytrap) {
                    maketrapdf <- function(out) {
                        out <- matrix(unlist(out), nrow = ntrap, byrow = TRUE,
                                      dimnames = list(NULL, pnames))
                        as.data.frame(out)
                    }
                    if (nclass>1) {
                        ## one df for each class
                        clss <- rep(1:nclass, length.out=length(out))
                        out <- lapply(split(out, clss), maketrapdf)
                        if (byclass) out else out[[1]]
                    }
                    else {
                        maketrapdf(out)
                    }
                }
                else {
                    out
                }
            }
            else {
                extractpar(temp[[1]], ntrap = NA, initial = FALSE)
            }
        }
        else {
            if (!is.data.frame(temp) || (nrow(temp) > length(object$link)))
                stop ("unexpected input to detectpar()")
            
            temp <- temp[, 'estimate', drop = FALSE]
            temp <- split(temp[,1], rownames(temp))
            temp <- c(temp, object$fixed)
            
            if (object$details$param %in% c(4,5)) {
                Dindex <- match ('D', names(temp))
                sigmakindex <- match ('sigmak', names(temp))
                cindex <- match ('c', names(temp))
                temp[[sigmakindex]] <- temp[[sigmakindex]] / temp[[Dindex]]^0.5 + temp[[cindex]]
                names(temp)[sigmakindex] <- 'sigma'
            }
            if (object$details$param %in% c(3,5)) {
                a0index <- match ('a0', names(temp))
                sigmaindex <- match ('sigma', names(temp))
                lambda0 <- temp[[a0index]] / 2 / pi / temp[[sigmaindex]]^2 * 10000
                temp[[a0index]] <- if (object$detectfn %in% 14:19) lambda0 else 1-exp(-lambda0)
                names(temp)[a0index] <- if (object$detectfn %in% 0:8) 'g0'
                else if (object$detectfn %in% 14:19) 'lambda0'
                else stop ('invalid combination of param %in% c(3,5) and detectfn')
            }
            
            temp <- temp[pnames]
            
            # acoustic cutval argument not estimated, but add to vector
            if ((object$detectfn > 9) & (object$detectfn <14))
                temp <- c(temp, list(cutval = object$details$cutval))
            
            temp
        }
    }
    if (!inherits(object,'secr'))
        stop ("requires 'secr' object")
    pnames <- parnames(object$detectfn)
    nclass <- object$details$nmix
    if (object$details$param == 2) pnames[1] <- 'esa'
    if (nclass > 1) pnames <- c(pnames, 'pmix')
    
    tr <- traps(object$capthist)
    if (!ms(tr)) tr <- list(tr)   # force to list
    
    if (bytrap) {
        arglist <- list(...)
        if ('newdata' %in% names(arglist)) 
            newdat <- arglist$newdata
        else 
            newdat <- makeNewData (object, all.levels = TRUE)
        sessions <- levels(newdat$session)
        nsessions <- length(sessions)
        temppred <- vector('list', nsessions)
        names(temppred) <- sessions
        if (nsessions==1) {
            names(tr) <- sessions[1]
        }
        trapcovnames <- names(newdat)[names(newdat) %in% 
                                          names(covariates(tr[[1]]))]
        onepred <- function(trapno, session) {
            newdatS <- newdat[newdat$session == session, , drop = FALSE]
            if (length(trapcovnames)>0) {
                newdatS[,trapcovnames] <-
                    covariates(tr[[session]])[trapno,trapcovnames]    
            }
            predict (object, newdata = newdatS)
        }
        for (sess in sessions) {
            ntrap <- ndetector(tr[[sess]])
            temppred[[sess]] <- lapply(1:ntrap, onepred, session = sess)
        }
        temppred <- unlist(temppred, recursive = FALSE)
        names(temppred) <- apply(expand.grid(sessions, 1:ntrap),1,paste, collapse=',')
    }
    else{
        temppred <- predict (object, ...)
    }
    ntrap <-sapply(tr, ndetector)
    if (ms(object)) {
        sess <- sapply(strsplit(names(temppred), ','), '[', 1)
        temppred <- split(temppred, sess)
        temp <- mapply(extractpar, temppred, ntrap = ntrap, 
                       SIMPLIFY = FALSE)
        temp
    }
    else {
        extractpar(temppred, ntrap = ntrap)
    }
}

# new code, not finished 2024-10-04
detectparnew <- function(object, ..., byclass = FALSE, bytrap = FALSE) {
    extractpar <- function (temp, ntrap, initial = TRUE) {
        ## for one session
        if (!is.data.frame(temp))   ## assume list
        {
            if (byclass || bytrap) {
                if (initial && nclass>1 && bytrap) {
                    temp <- unlist(temp, recursive = FALSE)
                }
                out <- lapply(temp, extractpar, FALSE)
                
                if (bytrap) {
                    maketrapdf <- function(out) {
                        out <- matrix(unlist(out), nrow = ntrap, byrow = TRUE,
                                      dimnames = list(NULL, pnames))
                        as.data.frame(out)
                    }
                    if (nclass>1) {
                        ## one df for each class
                        clss <- rep(1:nclass, length.out=length(out))
                        out <- lapply(split(out, clss), maketrapdf)
                        if (byclass) out else out[[1]]
                    }
                    else {
                        maketrapdf(out)
                    }
                }
                else {
                    out
                }
            }
            else {
                extractpar(temp[[1]], ntrap = NA, initial = FALSE)
            }
        }
        else {
            if (!is.data.frame(temp) || (nrow(temp) > length(object$link)))
                stop ("unexpected input to detectpar()")
            
            temp <- temp[, 'estimate', drop = FALSE]
            temp <- split(temp[,1], rownames(temp))
            temp <- c(temp, object$fixed)
            
            if (object$details$param %in% c(4,5)) {
                Dindex <- match ('D', names(temp))
                sigmakindex <- match ('sigmak', names(temp))
                cindex <- match ('c', names(temp))
                temp[[sigmakindex]] <- temp[[sigmakindex]] / temp[[Dindex]]^0.5 + temp[[cindex]]
                names(temp)[sigmakindex] <- 'sigma'
            }
            if (object$details$param %in% c(3,5)) {
                a0index <- match ('a0', names(temp))
                sigmaindex <- match ('sigma', names(temp))
                lambda0 <- temp[[a0index]] / 2 / pi / temp[[sigmaindex]]^2 * 10000
                temp[[a0index]] <- if (object$detectfn %in% 14:19) lambda0 else 1-exp(-lambda0)
                names(temp)[a0index] <- if (object$detectfn %in% 0:8) 'g0'
                else if (object$detectfn %in% 14:19) 'lambda0'
                else stop ('invalid combination of param %in% c(3,5) and detectfn')
            }
            
            temp <- temp[pnames]
            
            # acoustic cutval argument not estimated, but add to vector
            if ((object$detectfn > 9) & (object$detectfn <14))
                temp <- c(temp, list(cutval = object$details$cutval))
            
            temp
        }
    }
    if (!inherits(object,'secr'))
        stop ("requires 'secr' object")
    pnames <- parnames(object$detectfn)
    nclass <- object$details$nmix
    if (object$details$param == 2) pnames[1] <- 'esa'
    if (nclass > 1) pnames <- c(pnames, 'pmix')
    
    tr <- traps(object$capthist)
    if (!ms(tr)) tr <- list(tr)   # force to list
    
    if (bytrap) {
        arglist <- list(...)
        if ('newdata' %in% names(arglist)) 
            newdat <- arglist$newdata
        else 
            newdat <- makeNewData (object, all.levels = TRUE, bytrap = TRUE)
        sessions <- levels(newdat$session)
        nsessions <- length(sessions)
        if (nsessions==1) {
            names(tr) <- sessions[1]
        }
        onepred <- function(session) {
            newdatS <- newdat[newdat$session == session, , drop = FALSE]
            predict (object, newdata = newdatS, se.fit = FALSE)
        }
        temppred <- lapply(sessions, onepred)
        names(temppred) <- sessions
        temppred <- unlist(temppred, recursive = FALSE)
        names(temppred) <- apply(expand.grid(sessions, 1:ntrap),1,paste, collapse=',')
    }
    else{
        temppred <- predict (object, se.fit = FALSE, ...)
    }
    ntrap <- sapply(tr, ndetector)
    if (ms(object)) {
        sess <- sapply(strsplit(names(temppred), ','), '[', 1)
        temppred <- split(temppred, sess)
        temp <- mapply(extractpar, temppred, ntrap = ntrap, 
                       SIMPLIFY = FALSE)
        temp
    }
    else {
        extractpar(temppred, ntrap = ntrap)
    }
}

############################################################################################

# return parameters for first detection
# optionally by trap and by class

detectpar0 <- function(object, ..., byclass = FALSE, bytrap = FALSE) {
    beta <- object$fit$par
    # insert fixed betas as needed
    beta <- fullbeta(beta, object$details$fixedbeta)
    # get lookup table 
    realparval0 <- makerealparameters (object$design0, beta,
                                       object$parindx, object$link, object$fixed)  # naive
    # which indices of PIA do we care about? others use only first element
    indlist <- as.list(rep(1,5))
    if (bytrap) indlist[[4]] <- TRUE 
    if (byclass) indlist[[5]] <- TRUE 
    # find indices
    arglist <- c(list(x = object$design0$PIA, drop = FALSE), indlist)
    selected <- do.call('[', arglist)
    out <- as.data.frame(realparval0[selected,, drop = FALSE])
    nmix <- dim(object$design0$PIA)[5]
    if (byclass) out$class <- rep((1:nmix)-1, each = nrow(out)/nmix)
    out
}

############################################################################################

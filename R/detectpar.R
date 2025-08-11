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
    #--------------------------------------------------------------------------
    extractpar <- function (onesessiondf) {

        if (object$details$param %in% c(4,5)) {
            Dindex <- match ('D', names(onesessiondf))
            sigmakindex <- match ('sigmak', names(onesessiondf))
            cindex <- match ('c', names(onesessiondf))
            onesessiondf[[sigmakindex]] <- onesessiondf[[sigmakindex]] / onesessiondf[[Dindex]]^0.5 + onesessiondf[[cindex]]
            names(onesessiondf)[sigmakindex] <- 'sigma'
        }
        if (object$details$param %in% c(3,5)) {
            a0index <- match ('a0', names(onesessiondf))
            sigmaindex <- match ('sigma', names(onesessiondf))
            lambda0 <- onesessiondf[[a0index]] / 2 / pi / onesessiondf[[sigmaindex]]^2 * 10000
            onesessiondf[[a0index]] <- if (object$detectfn %in% 14:19) lambda0 else 1-exp(-lambda0)
            names(onesessiondf)[a0index] <- if (object$detectfn %in% 0:8) 'g0'
            else if (object$detectfn %in% 14:19) 'lambda0'
            else stop ('invalid combination of param %in% c(3,5) and detectfn')
        }
        mixvar <- switch(nclass, character(0),'h2','h3')
        if (nclass>1) {
            if (byclass) {
            classlist <- split(onesessiondf, onesessiondf[,mixvar])
            onesessionlist <- lapply(classlist, function(x) {
                x <- x[, pnames, drop = FALSE]
                rownames(x) <- 1:nrow(x)
                as.list(x)
            })
            }
            else {
                mix <- onesessiondf[,mixvar]
                OK <- mix==levels(mix)[1]
                onesessionlist <- as.list(onesessiondf[OK,pnames, drop = FALSE])
            }
        }
        else {
            onesessionlist <- as.list(onesessiondf[,pnames, drop = FALSE])
        }    
        
        # acoustic cutval argument not estimated, but add to vector
        if ((object$detectfn > 9) & (object$detectfn <14))
            onesessionlist <- c(onesessionlist, list(cutval = object$details$cutval))
        
        onesessionlist
        
    }  # end extractpar
    #--------------------------------------------------------------------------
    
    # ----------------------------
    # preliminary
    if (!inherits(object,'secr')) stop ("requires 'secr' object")
    if (!is.null(object$fixed) && length(object$fixed)>0) warning("fixed parameters ignored") 
    
    pnames <- secr_parnames(object$detectfn)
    nclass <- object$details$nmix
    
    if (object$details$param == 2) pnames[1] <- 'esa'
    if (nclass > 1) pnames <- c(pnames, 'pmix')
    tr <- traps(object$capthist)
    if (!ms(tr)) tr <- list(tr)   # force to list
    
    # ----------------------------------------
    # make list of session-specific dataframes
    if (bytrap) {
        arglist <- list(...)
        if ('newdata' %in% names(arglist)) {
            newdat <- arglist$newdata
        }
        else {
            newdat <- makeNewData (object, all.levels = FALSE, bytrap = TRUE)
        }
        sessions <- levels(newdat$session)
        nsessions <- length(sessions)
        if (nsessions==1) {
            names(tr) <- sessions[1]
        }
        onepred <- function(session) {
            newdatS <- newdat[newdat$session == session, , drop = FALSE]
            trapID <- newdatS$trapID
            # must split by trap for latent classes
            if (nclass>1) {
                newdatS <- split(newdatS, trapID)
                out <- lapply(newdatS, predict, object = object, se.fit = FALSE)
                out <- do.call(rbind, out)
            }
            else {
                out <- predict(object, newdata = newdatS, se.fit = FALSE)
            }
            out$trapID <- trapID
            out
        }
        temppred <- lapply(sessions, onepred)
        names(temppred) <- sessions
    }
    else{
        predicted <- predict (object, se.fit = FALSE, ...)
        temppred <- split(predicted, predicted$session)
        #if (is.data.frame(temppred)) temppred <- list(temppred)
    }
    # ----------------------------------------
    
    # temppred is a list, one component per session

    #-----------------------------------------
    # 2025-08-12 apply overall tempxy
    fixxy <- function (onepred) {
        if ('sigmaxy'   %in% names(onepred)) onepred$sigma   <- onepred$sigma * onepred$sigmaxy
        if ('lambda0xy' %in% names(onepred)) onepred$lambda0 <- onepred$lambda0 * onepred$lambda0xy
        if ('a0xy'      %in% names(onepred)) onepred$a0      <- onepred$a0 * onepred$a0xy
        if ('sigmakxy'  %in% names(onepred)) onepred$sigmak  <- onepred$sigmak * onepred$sigmakxy
        onepred
    }
    if (any(.localstuff$spatialparameters %in% names(object$model))) {
        temppred <- lapply(temppred, fixxy)    
    }
    #-----------------------------------------
    
    if (ms(object)) 
        lapply(temppred, extractpar)
    else 
        extractpar(temppred[[1]])
}
############################################################################################


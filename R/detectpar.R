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
    
    pnames <- parnames(object$detectfn)
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

    if (ms(object)) 
        lapply(temppred, extractpar)
    else 
        extractpar(temppred[[1]])
}
############################################################################################

# 2024-10-04
# return parameters for first detection
# optionally by session, animal, occasion, trap or latent class

detectpar0 <- function(object,
                       bysession  = TRUE, 
                       byanimal   = FALSE, 
                       byoccasion = FALSE, 
                       bytrap     = FALSE, 
                       byclass    = FALSE) {
    beta <- object$fit$par
    # insert fixed betas as needed
    beta <- fullbeta(beta, object$details$fixedbeta)
    # get lookup table 
    realparval0 <- makerealparameters (
        object$design0, beta,
        object$parindx, object$link, object$fixed)  # naive
    param <- object$details$param
    if (param>0) {
        if (param == 3) {
            realparval0 <- reparameterize(
                realparval0, object$detectfn, object$details, object$mask, 
                traps(object$capthist), NA, NA)
        }
        else stop ("parameterisation ", param, " not available")
    }
    # which indices of PIA do we care about? others use only first element
    indlist <- as.list(rep(1,5))
    if (bysession)  indlist[[1]] <- TRUE
    if (byanimal)   indlist[[2]] <- TRUE 
    if (byoccasion) indlist[[3]] <- TRUE 
    if (bytrap)     indlist[[4]] <- TRUE 
    if (byclass)    indlist[[5]] <- TRUE 
    # find indices
    arglist <- c(list(x = object$design0$PIA, drop = FALSE), indlist)
    selected <- do.call('[', arglist)
    out <- as.data.frame(realparval0[selected,, drop = FALSE])
    
    dd <- lapply(dim(selected), function(x) 1:x)
    names(dd) <- c('session','animal','occasion','trap','class')
    id <- do.call(expand.grid, dd)[,sapply(indlist, is.logical), drop = FALSE]
    cbind(out, id)
    
}

############################################################################################

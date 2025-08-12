#############################################################################
## package 'secr'
## derived density from conditional (relativeD) models
## 2025-08-08 derivedDcoef moved from derivedDsurface and greatly modified

#-------------------------------------------------------------------------------

# exported
derivedDcoef.secrlist <- function (object, se = FALSE, ...) {
    lapply(object, se = se, ...)
}

derivedDcoef.secr <- function (object, se = FALSE, ...) {
    #-------------------------------------------------------------------------------
    onelinkk <- function (beta)
    {
        # Return k-hat for given beta
        sessfn <- function(r) sum(1 / esa (object, r, beta, Dweight = TRUE))
        kr <- sapply(1:nsessions, sessfn)  # one estimate per session
        transform(mean(kr), object$link$D)
    }
    #-------------------------------------------------------------------------------
    
    if (is.null(object$model$D) || is.null(object$link$D) || !object$CL) {
        warning ("not relative density model")
        return(NULL)
    }
    else {

        vcv.derivedk <- function (object) {
            sessfn <- function(r) 1 / esa (object, r, Dweight = TRUE)
            a <- unlist(lapply(1:nsessions, sessfn))
            # A <-  if (ms(object)) sum(sapply(object$mask, masksize)) else masksize(object)
            # s2 <- switch (tolower(object$details$distribution),
            #               poisson  = sum (1/a^2),
            #               binomial = sum (( 1 - a / A) / a^2))
            
            linkkgrad  <- nlme::fdHess(object$fit$par, onelinkk)$gradient
            vark <- linkkgrad %*% object$beta.vcv %*% linkkgrad
            np <- length(linkkgrad) + 1
            completevcv <- matrix(0, np, np)
            completevcv[2:np,2:np] <- object$beta.vcv 
            # ad hoc Poisson variance
            s2 <- se.transform(length(a),   sqrt(length(a)) ,'log')^2
            completevcv[1,1]       <- vark + s2
            completevcv[2:np,1]    <- object$beta.vcv %*% linkkgrad
            completevcv[1,2:np]    <- object$beta.vcv %*% linkkgrad
            completevcv
        }
        
        nsessions <- if (ms(object)) length(object$capthist) else 1
        beta <- object$fit$par
        linkk <- onelinkk(beta)
        if (se) {
            warning ("se(beta0) from derivedDcoef() is approximate", call. = FALSE)
            vcv <- vcv.derivedk(object) 
        }
        else {
            vcv <- matrix(NA,1,1)
        }
        
        # return on link scale
        beta0 <- linkk # transform(k, object$link$D)
        se.beta0 <- sqrt(vcv[1,1]) # se.transform(k, sqrt(vcv[1,1]), object$link$D)
        oldcoef <- coef(object)
        alpha <- attr(oldcoef, 'alpha')
        z <- abs(qnorm(1 - alpha/2))
        Dcoef <- c(beta0, se.beta0, beta0 - z * se.beta0, beta0 + z * se.beta0)
        tmp <- rbind(D = Dcoef, oldcoef)
        if (object$link$D == 'identity') {
            tmp[grepl('D.', rownames(tmp)),] <- tmp[grepl('D.', rownames(tmp)),] * tmp[1,1]
        }
        attr(tmp, 'vcv') <- vcv
        tmp
        
        # NB gives same (estimate, SE.estimate) as derived(object, Dweight=T), 
        # and takes about the same time
    }
    
}
#-------------------------------------------------------------------------------

# not exported
predictD <- function (object, ...) {
    der <- derivedDcoef(object, se = TRUE)
    object$details$relativeD <- FALSE
    object$betanames <- c('D', object$betanames)
    tmp <- attr(der, 'vcv')
    object$beta.vcv <- tmp
    object$details$fixedbeta[1] <- NA
    object$fit$par <- der$beta
    predict(object, ...)    
}
#-------------------------------------------------------------------------------

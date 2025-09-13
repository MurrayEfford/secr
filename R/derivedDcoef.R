#############################################################################
## package 'secr'
## derived density from conditional (relativeD) models
## 2025-08-08 derivedDcoef moved from derivedDsurface and greatly modified
## 2025-08-21 derivedDfit() new function
## 2025-09-09 derivedDcoef minor rearrangement
#-------------------------------------------------------------------------------

# exported
derivedDcoef.secrlist <- function (object, se = FALSE, ...) {
    lapply(object, derivedDcoef, se = se, ...)
}

derivedDcoef.secr <- function (object, se = FALSE, ...) {
    #-------------------------------------------------------------------------------
    onelinkk <- function (beta)
    {
        # Return k-hat for given beta
        sessfn <- function(r) esa (object, r, beta, Dweight = TRUE)
        ar <- lapply(1:nsessions, sessfn)  # one estimate per session
        sessionn <- sapply(ar, length)
        # weighting a_i by relative density in each session
        if (Dlam) {
            RD <- cumprod(predictDlambda(object)$estimate)
        }
        else if (nsessions>1) {
            RD <- sapply(predict(object), '[', 'D', 'estimate')
        }
        else {
            RD <- 1
        }
        Rweight <- 1 / rep(RD, sessionn)   # to match unlist(ar)
        k <- sum(1 / unlist(ar) * Rweight) / nsessions
        transform(k, object$link$D)
    }
    #-------------------------------------------------------------------------------
    nsessions <- if (ms(object)) length(object$capthist) else 1
    beta <- object$fit$par
    
    if (is.null(object$model$D) || is.null(object$link$D) || !object$CL) {
        warning ("not relative density model")
        return(NULL)
    }
    Dlam <- !is.null(object$details$Dlambda) && object$details$Dlambda
    if (nsessions > 1) {
        warning ("derivedDcoef approximate for multisession models")
        # return(NULL)
    }
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
    linkk <- onelinkk(beta)
    if (se) {
        vcv <- vcv.derivedk(object) 
    }
    else {
        vcv <- matrix(NA,1,1)
    }
    
    # return on link scale
    beta0 <- linkk 
    se.beta0 <- sqrt(vcv[1,1])
    oldcoef <- coef(object)
    alpha <- attr(oldcoef, 'alpha')
    z <- abs(qnorm(1 - alpha/2))
    
    Dcoef <- c(beta0, se.beta0, beta0 - z * se.beta0, beta0 + z * se.beta0)
    tmpcoef <- rbind(D = Dcoef, oldcoef)
    if (Dlam) rownames(tmpcoef)[1] <- 'D.D1'
    if (object$link$D == 'identity') {
        np <- nrow(tmpcoef)
        Dcoef <- grepl('D.', rownames(tmpcoef))
        tmpcoef[Dcoef,] <- tmpcoef[Dcoef,] * beta0
        if (se) {
            vcv[Dcoef, 2:np] <- vcv[Dcoef, 2:np] * beta0
            vcv[2:np, Dcoef] <- vcv[2:np, Dcoef] * beta0
        }
    }
    
    attr(tmpcoef, 'vcv') <- vcv
    tmpcoef
    
    # NB gives same (estimate, SE.estimate) as derived(object, Dweight=T), 
    # and takes about the same time
    
    
}
#-------------------------------------------------------------------------------

# function to convert fitted single-session CL relative D model to a full model
# used by region.N etc.

derivedDfit <- function(object, vcv = TRUE) {
    if (is.null(object$model$D) || is.null(object$link$D) || !object$CL) {
        return(object) # unchanged
    }
    der <- derivedDcoef(object, se = vcv)
    object$fit$par <- der$beta
    object$fit$estimate <- object$fit$par
    object$beta.vcv <- attr(der, 'vcv')
    object$details$fixedbeta[1] <- NA  # inferred, not fixed
    object$betanames <- c('D', object$betanames)
    object$CL <- FALSE
    object$details$relativeD <- FALSE
    object
}
#-------------------------------------------------------------------------------
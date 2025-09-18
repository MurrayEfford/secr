############################################################################################
## package 'secr'
## trend.R
## 2023-09-29 2025-08-12 2025-09-18 (annotation) 
############################################################################################

Dfn2 <- function (designD, beta = NULL, ...) {
    # function to infer linear predictor (lp) of log(D1), log(lambda1), log(lambda2) etc.
    # and compute the cumulative sum of these to infer session-specific log(D)
    # used for Dlambda parameterization (details$Dlambda == TRUE) by secr_getD()
    dimD <- attr(designD, 'dimD')
    designD[1:dimD[1],] <- 0
    designD <- cbind(rep(c(1,0), c(dimD[1], nrow(designD)-dimD[1])), designD[,-1])
    if (is.null(beta)) return(ncol(designD)) # number of beta parameters
    lp <- designD %*% beta
    lp <- array(lp, dim = dimD[-2])   # mask x session matrix
    t(apply(lp, 1, cumsum))           # density on link scale
}

predictDlambda <- function (object, alpha = 0.05) {
    nsessions <- length(object$capthist)   
    dimD <- attr(object$designD, 'dimD')

    beta <- secr_complete.beta(object)
    beta.vcv <- secr_complete.beta.vcv(object)
    beta.vcv[is.na(beta.vcv)] <- 0
    Dpar <- object$parindx[['D']]
    # ad hoc patch for changed beta, drop last
    if (length(Dpar) > ncol(object$designD)) {
        stop("fitted Dlambda model not compatible with 5.3.0")
    }
    beta <- beta[Dpar]
    beta.vcv <- beta.vcv[Dpar,Dpar, drop = FALSE]
    
    
    # use subset of design matrix, just 1 row per session
    mat <- object$designD[dimD[1]*(0:(nsessions-1))+1,,drop=FALSE]

    if (object$CL) warning("predictDlambda does not yet infer D1 when CL = TRUE")
    lp <- mat %*% beta
    prepost <- function(i) mat[i,, drop = FALSE] %*% beta.vcv %*% t(mat[i,, drop = FALSE])

    selp <- sapply(1:nsessions, prepost)^0.5
    z <- abs(qnorm(1-alpha/2))   
    out <- data.frame(
        estimate = exp(lp),
        SE.estimate = exp(lp) * sqrt(exp(selp^2) - 1),
        lcl = exp(lp - z*selp),
        ucl = exp(lp + z*selp)
    )
    rownames(out) <- c('D1', paste0('lambda', 1:(nrow(out)-1)))
    out
}

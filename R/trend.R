############################################################################################
## package 'secr'
## trend.R
## 2023-09-29 
############################################################################################

Dfn2 <- function (designD, beta = NULL, ...) {
    dimD <- attr(designD, 'dimD')
    designD[1:dimD[1],] <- 0
    designD <- cbind(rep(c(1,0), c(dimD[1], nrow(designD)-dimD[1])), designD)
    if (is.null(beta)) return(ncol(designD)) # number of beta parameters
    lp <- designD %*% beta
    lp <- array(lp, dim = dimD[-2])   # mask x session matrix
    t(apply(lp, 1, cumsum))           # density on link scale
}

predictDlambda <- function (object, alpha = 0.05) {
    z <- abs(qnorm(1-alpha/2))   
    Dpar <- grepl('D.', rownames(coef(object)))
    beta <- coef(object)[Dpar, 'beta']
    nsessions <- length(object$capthist)   
    dimD <- attr(object$designD, 'dimD')
    
    # use subset of design matrix, just 1 row per session
    vars <- object$designD[dimD[1]*(0:(nsessions-1))+1,,drop=FALSE]
    vars[1,1] <- 0
    mat <- cbind(rep(c(1,0), c(1, nrow(vars)-1)), vars)
    
    lp <- mat %*% beta
    vcv <- object$beta.vcv[Dpar,Dpar, drop = FALSE]
    prepost <- function(i) mat[i,, drop = FALSE] %*% vcv %*% t(mat[i,, drop = FALSE])
    selp <- sapply(1:nsessions, prepost)^0.5
    out <- data.frame(
        estimate = exp(lp),
        SE.estimate = exp(lp) * sqrt(exp(selp^2) - 1),
        lcl = exp(lp - z*selp),
        ucl = exp(lp + z*selp)
    )
    rownames(out) <- c('D1', paste0('lambda', 1:(nrow(out)-1)))
    out
}

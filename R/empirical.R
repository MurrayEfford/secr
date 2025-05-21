############################################################################################
## package 'secr'
## empirical.R
## design-based variance of derived density
## 2011 05 10 - derived.systematic quarantined in 'under development' folder
## 2012 03 17 - derived.nj extended to single n
## 2012 12 24 - clust argument required in derived.external
## 2014-10-17 revised to require spsurvey
## 2016-06-05 started derivedStratified
## 2017-11-21 renamed derived.nj to derivednj, derived.cluster to derivedCluster etc.
##            to do - treat sessions resulting from mash() clustergroups as strata
## 2018-12-16 derivednj and derivedSessions allow weighted option (now R3)
## 2018-12-17 derivedStratified moved to its own file
## 2018-12-20 derivednj and derivedSessions allow weighted option R2
## 2021-10-18 localvar suspended owing to changes in spsurvey package
## 2025-05-22 derivedCluster(), derivedSession() protected against nj = 0

################################################################################


# localvar <- function (z, xy, weights = NULL) {
#     ## z is response variable
#     ## xy is dataframe or list with x, y
#     ## choice of output verified by comparing to
#     ## sum ((nj - n/J)^2) * J / (J-1) for ovensong vector
#     ## with vartype='SRS'
#     if (requireNamespace('spsurvey', quietly = TRUE)) {
#         if (is.null(weights)) weights <- rep(1,length(z))
#         else if (length(weights) != length(z))
#             stop ("weights must be same length as z")
#         temp1 <- spsurvey::total.est(z = z, 
#                                      wgt = weights,
#                                      x = xy[,1], y = xy[,2],
#                                      vartype = "Local")
#     }
#     else {
#         stop ("package 'spsurvey' required for local variance")
#     }
# 
#     temp1[1,4]^2
# }

# LPMvar <- function (nj, xy, esa) {
#     # just playing around here
#     ## z is response variable
#     ## xy is dataframe or list with x, y
#     ## nearest neighbour variance from Robertson et al. 2018
#     
#     z <- nj/esa
#     d <- edist(xy,xy)
#     diag(d) <- Inf
#     nearest <- apply(-d, 1, nnet::which.is.max)  ## breaks ties at random
#     A <- sum(esa)
#     J <- length(nj)
#     0.5 * sum((z - z[nearest])^2) * A / (J-1)
# }

############################################################################

derivedMash <- function (object, sessnum = NULL, method = c('SRS', 'local'), 
                         alpha = 0.05, loginterval = TRUE) {
    if (ms(object) & (is.null(sessnum))) {
            ## recursive call if MS
            sessnames <- session(object$capthist)
            nsess <- length(sessnames)
            output <- vector('list', nsess )
            for (i in 1:nsess) {
                output[[i]] <- derivedMash(object, sessnum = i, method, alpha, loginterval)
            }
            names(output) <- sessnames
            output
    }
    else {
        method <- match.arg(method)
        if (is.null(sessnum)) {
            n.mash <- attr(object$capthist, 'n.mash')
            centres <- attr(object$capthist, 'centres')
        }
        else {
            if (length(sessnum)>1)
                stop ("requires single sessnum")
            n.mash <- attr(object$capthist[[sessnum]], 'n.mash')
            centres <- attr(object$capthist[[sessnum]], 'centres')
        }

        if (is.null(n.mash))
            stop ("requires capthist with 'n.mash' attribute")

        if (ms(object))
            der <-  derived(object, se.esa = TRUE, se.D = FALSE)[[sessnum]][,1:2]
        else
            der <-  derived(object, se.esa = TRUE, se.D = FALSE)[,1:2]

        derivednj(n.mash, der['esa','estimate'], der['esa','SE.estimate'],
                   method, centres, alpha, loginterval, independent.esa = FALSE)

    }
}
####################################################################################
derivedCluster <- function (object, method = c('SRS', 'R2', 'R3','local', 
                    'poisson', 'binomial'), alpha = 0.05, loginterval = TRUE) {

    ## 2018-12-17 allow different cluster geometry
    ## assume single session 
    method <- match.arg(method)
    capthist <- object$capthist
    tr <- traps(capthist)

    if (ms(object)) 
        stop("derivedCluster now does not accept multi-session models")

    if (is.null(clusterID(tr)))
        stop("derivedCluster expects traps object to have clusterID")
    
    if (!is.null(attr(capthist, 'n.mash')))
        warning ("use derivedMash for mashed data")
    
    nj <- cluster.counts(capthist)
    centres <- cluster.centres(tr)
    der <- derived(object, se.esa = TRUE, se.D = FALSE, bycluster = TRUE)    
    der <- do.call(rbind, der)
    der <- der[grepl('esa', rownames(der)), 1:2]
    
    if (any(nj==0)) {
        warning("esa imputed for nj = 0")
        der[nj==0,] <- rep(apply(der,2,mean, na.rm=TRUE), each = sum(nj==0))
    }
    
    # pass nj, esa, se.esa vectors to generic function
    # esa-hat are almost certainly non-independent
    
    derivednj(nj, esa = der[,1], se.esa = der[,2], method, centres, alpha, 
              loginterval, masksize(object$mask), independent.esa = FALSE)
}
####################################################################################

derivedExternal <- function (object, sessnum = NULL, nj, cluster,
    buffer = 100, mask = NULL, noccasions = NULL, 
    method = c('SRS','local'), xy = NULL, alpha = 0.05, loginterval = TRUE) {

    se.deriveesa <- function (selection, asess, noccasions) {
        CLesa  <- esagradient (object, selection, asess, noccasions)
        sqrt(CLesa %*% object$beta.vcv %*% CLesa)
    }

    if (ms(object) & is.null(sessnum)) {
            ## recursive call if MS
            sessnames <- session(object$capthist)
            nsess <- length(sessnames)
            output <- vector('list', nsess )
            for (i in 1:nsess) {
                output[[i]] <- derivedExternal(object, sessnum = i, nj, cluster,
                    buffer, mask, noccasions, method, xy, alpha, loginterval)
            }
            names(output) <- sessnames
            output
    }
    else {
        method <- match.arg(method)
        if (!(inherits(cluster, 'traps')) | ms(cluster))
            stop ("cluster should be a single-session traps object")

        if (is.null(sessnum)) {
            capthist <- object$capthist
            sessnum <- 1
        }
        else
            capthist <- object$capthist[[sessnum]]

        if (is.null(noccasions))
            noccasions <- ncol(capthist)

        if (is.null(mask)) {
            mask <- make.mask(cluster, buffer = buffer)
        }

        ## substitute new traps and mask into old object
        traps(object$capthist) <- cluster
        object$mask <- mask

        ## use CLmeanesa for harmonic mean in case of individual variation
        esa.local <- CLmeanesa (object$fit$par, object, 1:nrow(capthist),
            sessnum, noccasions)
        se.esa <- se.deriveesa(1, sessnum, noccasions = noccasions)

        derivednj (nj, esa.local, se.esa, 
                   method, xy, alpha, loginterval, independent.esa = FALSE)
    }
}
####################################################################################

derivedSession <- function ( object, method = c('SRS','R2', 'R3','local','poisson', 'binomial'), 
    xy = NULL, alpha = 0.05, loginterval = TRUE, area = NULL, independent.esa = FALSE ) {
    method <- match.arg(method)
    if (!inherits (object, 'secr') | !ms(object))
        stop ("requires fitted multi-session model")
    if (('session' %in% object$vars) | ('Session' %in% object$vars) |
        (!is.null(object$sessioncov)))
        stop ("derivedSession assumes detection model constant over sessions")
    nj <- sapply(object$capthist, nrow)    ## animals per session
    der <-  derived(object, se.esa = TRUE, se.D = FALSE)
    der <- do.call(rbind, der)
    der <- der[grepl('esa', rownames(der)),1:2]
    if (any(nj==0)) {
        warning("esa imputed for nj = 0")
        der[nj==0,] <- rep(apply(der,2,mean, na.rm=TRUE), each = sum(nj==0))
    }
    if (is.null(area)) 
        area <- sum(sapply(object$mask, masksize))  ## assumes non-overlap
    derivednj(nj, der[,'estimate'], der[,'SE.estimate'], 
        method, xy, alpha, loginterval, area, independent.esa)
}
############################################################################

derivednj <- function ( nj, esa, se.esa = NULL, 
                        method = c('SRS','R2', 'R3','local','poisson','binomial'), 
                        xy = NULL, alpha = 0.05, loginterval = TRUE, area = NULL, 
                        independent.esa = FALSE ) {
    method <- match.arg(method)
    esa <- unlist(esa)
    if (is.null(se.esa)) {
        esa <- matrix(esa, ncol = 2)
        se.esa <- esa[,2]
        esa <- esa[,1]
    }
    n <- sum (nj)                          ## total animals
    J <- length(nj)                        ## number of replicates
    if (length(esa) == 1) {
        esa <- rep(esa, J)
        se.esa <- rep(se.esa, J)
        independent.esa <- FALSE
    }
    if (length(esa) != J) stop ("lengths of nj and esa differ")
    A <- sum(esa)
    if (independent.esa)
        varA <- sum(se.esa^2)
    else
        varA <- J^2 * sum(esa/A * se.esa^2)

    if (J == 1) {
        method <- tolower(method)
        if (!(method %in% c('poisson','binomial'))) {
            method <- 'poisson'
            warning ("nj is of length 1; n variance component defaults to Poisson")
        }
    }
    if (method == 'binomial') {
        if (is.null(area))
            stop ("must specify 'area' for binomial variance")
        N <- area * n / A
    }
    
    if (method=='local') {
        warning("method = 'local' not available in secr 4.6 because of changes in spsurvey")
        # stop ("'local' method requires x and y coordinates")
    }
    varn <- switch (method,
        SRS = sum ((nj - n/J)^2) * J / (J-1),
        R2 = sum( esa^2 * (nj/esa - n/A)^2) * J / (J-1), ## Fewster et al 2009 R2
        R3 = A * sum( esa * (nj/esa - n/A)^2) / (J-1),   ## Fewster et al 2009 R3
        local = NA,                                      ## suppressed because total.est withdrawn spsurvey 5.0
        # local = localvar (nj, xy),                     ## from total.est in spsurvey, unweighted
        # local = localvar (nj, xy, 1/esa),              ## from total.est in spsurvey
        poisson = n,                              
        binomial = N * (A / area * (1 - A / area)),
        NA
    )
    D <- n / A 
    varD <- D^2 * (varn/n^2 + varA/A^2)
    temp <- data.frame(row.names = c('esa','D'), estimate = c(A,D), SE.estimate = c(varA,varD)^0.5)
    temp <- add.cl(temp, alpha, loginterval)
    temp$CVn <- c(NA, varn^0.5/n)
    temp$CVa <- c(NA, sqrt(varA)/A)
    temp$CVD <- c(NA, varD^0.5/D)
    attr(temp, 'nj') <- nj
    attr(temp, 'esa') <- esa
    attr(temp, 'se.esa') <- se.esa
    temp
}
############################################################################

# derivedSession(ovenbird.model.1, method='R3')


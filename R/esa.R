############################################################################################
## package 'secr'
## esa.R
## 2020-05-13 fixed bug that ignored individual covariates by providing only single-row CH0
## 2022-01-22 fixed bug that ignored details$ignoreusage
## 2024-09-22 S3 method
## 2025-01-04 experimental includes density if model$D && Dweight
############################################################################################

esa.secr <- function (object, sessnum = 1, beta = NULL, real = NULL, 
                      noccasions = NULL, ncores = NULL, Dweight = FALSE, ...)

# Return vector of 'a' for given g0, sigma, [z (if hazard fn) ] and session
# detectfn is integer code for detection function 0 = halfnormal, 1 = hazard, 2 = exponential
# 'sessnum' is integer index of session (factor level of the 'session' attribute in capthist)
# object must have at least capthist, mask, detectfn

{
    if (ms(object))
        capthists <- object$capthist[[sessnum]]
    else
        capthists <- object$capthist

    # 2023-06-01
    if (nrow(capthists) == 0) {
        return (numeric(0))    
    }
    
    if (ms(object$mask)) {
        mask <- object$mask[[sessnum]]
    }
    else {
        mask <- object$mask
    }

    if (is.null(beta) & is.null(real)) {
        beta <- object$fit$par
    }
    
    beta <- secr_fullbeta(beta, object$details$fixedbeta)
    trps <- traps(capthists)  ## need session-specific traps
    if (!all(detector(trps) %in% .localstuff$individualdetectors))
        stop ("require individual detector type for esa")
    n       <- max(nrow(capthists), 1)
    s       <- ncol(capthists)
    dettype <- secr_detectorcode(trps, noccasions = s)
    constant <- !is.null(noccasions)    ## fix 2011-04-07
    if (is.null(noccasions)) {
        noccasions <- s
    }
    markocc <- markocc(traps(capthists))
    if (is.null(markocc))
        markocc <- rep(1,s)  ## simple marking
    allsighting <- !any(markocc>0)
    if (allsighting) {
        ## drop all zero histories, consider sighting as if marking
        # capthists <- subset(capthists, apply(capthists!=0,1,sum)>0)
        warning("use of allsighting here untested")
        ## markocc[] <- 1
    }
    else {
        #not right, but getting there
        #object$capthist <- subset(object$capthist, occasions = (markocc>0))
    }
    #----------------------------------------------------------------------
    nmix    <- if (is.null(object$details$nmix)) 1 else object$details$nmix
    
    knownclass <- secr_getknownclass(capthists, nmix, object$hcov)
    k <- secr_getk(trps)
    K <- if (length(k)>1) length(k)-1 else k
    binomN <- object$details$binomN
    m      <- length(mask$x)            ## need session-specific mask...
    cellsize <- secr_getcellsize(mask)       ## length or area
    
    if (is.null(object$model$D) || !Dweight) {
        D <- rep(1,m)
    }
    else {
        D <- secr_predictD(object, object$mask, group = NULL, session = sessnum, 
                             parameter = 'D')
    }
    pi.density <- matrix(D/sum(D), ncol = 1)       # dim m x 1
    
    # Non-Euclidean distance parameters
    NElist <- secr_makeNElist(object, object$mask, group = NULL, sessnum)
    #----------------------------------------------------------------------
    if (constant) {
        ## assume constant
        if (is.null(beta))
            realparval0 <- detectpar(object)
        else {
            realparval0 <- secr_makerealparameters (object$design0, beta,
                object$parindx, object$link, object$fixed)  # naive
            # realparval0 <- as.list(realparval0)
            # names(realparval0) <- secr_parnames(object$detectfn)
            ## 2016-11-12
            realparval0 <- as.list(realparval0[1,])
            realparval0$cutval <- attr(object$capthist,'cutval')  ## 2016-05-22 may be NULL
        }
        a <- cellsize * sum(pi.density * pdot(X = mask, traps = trps, detectfn = object$detectfn,
                             detectpar = realparval0, noccasions = noccasions, 
                             ncores = ncores))
        return(rep(a,n))
    }
    else {
        if (is.null(beta)) {
            if (is.null(real))
                stop ("requires real parameter values")
            PIA0 <- array(1, dim=c(1,n,s,k,nmix))
            realparval0 <- matrix(rep(real, rep(n,length(real))), nrow = n)   ## UNTRANSFORMED
        }
        else {
            if (object$CL) {
                PIA0 <- object$design0$PIA[sessnum,1:n,1:s,1:K,,drop = FALSE]
            }
            else {
                ## use whatever columns are present as the number of groups
                ## should be constant over sessions
                PIA0 <- object$design0$PIA[sessnum,,1:s,1:K,,drop = FALSE]
                ## fill array with PI appropriate to grouping of i-th animal
                PIA0 <- PIA0[1, secr_group.factor(capthists, object$groups),,,,drop = FALSE]
            }

            realparval0 <- secr_makerealparameters (object$design0, beta,
                object$parindx, object$link, object$fixed)  # naive

        }
        ## not compatible with sigmak parameterizations
        Dtemp <- NA
        Xrealparval0 <- secr_reparameterize (realparval0, object$detectfn, object$details,
                                        mask, trps, Dtemp, s)

        usge <- usage(trps)
        # if (is.null(usge)) {
        # 2022-01-22 respect ignore usage cf Greg note
        if (is.null(usge) || object$details$ignoreusage) {
                usge <- matrix(1, nrow = K, ncol = s)
            used <- 1
        }
        else {
            ## force to binary 2012-12-17
            used <- (usge > 1e-10) * 1
        }
        if (any(used == 0))
            PIA0 <- PIA0 * rep(rep(t(used),rep(n,s*K)),nmix)

        miscparm <- numeric(4)
        miscparm[1] <- object$details$cutval

        ## not if first detector type is telemetry
        if (dettype[1] %in% c(13)) {
            return(NA)
        }
        else {
          # NE <- NULL   ## no NE covariates (yet)
          # pi.density <- matrix(1/m, nrow=m, ncol=1)
          ncores <- setNumThreads(ncores)  
          grain <- if (ncores==1) 0 else 1
          gkhk <- secr_makegk (dettype, object$detectfn, trps, mask, 
                                 object$details, sessnum, NElist, pi.density, 
                                 miscparm, Xrealparval0, grain, ncores)
          if (any(dettype==0)) {
              CH0 <- secr_nullCH (c(n,s), object$design0$individual)
          }
          else {
              CH0 <- secr_nullCH (c(n,s,K), object$design0$individual)
          }
          binomNcode <- secr_recodebinomN(dettype, binomN, secr_telemcode(trps))
          pmixn <- secr_getpmix (knownclass, PIA0, Xrealparval0)
          MRdata <- list(markocc = markocc, firstocc=rep(-1,nrow(CH0)))
          pID <- secr_getpID(PIA0, Xrealparval0, MRdata)
          pdot <- secr_integralprw1 (
              cc0 = nrow(Xrealparval0),
              haztemp = secr_gethazard(m, binomNcode, nrow(Xrealparval0), gkhk$hk, PIA0, usge),
              gkhk = gkhk,
              pi.density = pi.density,
              PIA0 = PIA0,
              CH0 = CH0,
              binomNcode = binomNcode,
              MRdata = MRdata,
              grp = rep(1,n),
              usge = usge,
              pmixn = pmixn,
              pID = pID,
              grain = grain,
              ncores = ncores)
          
          pdot * cellsize * sum(D)
        }
    }
}
############################################################################################

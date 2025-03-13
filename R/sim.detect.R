## package 'secr'
## sim.detect.R
## 2024-07-29 split from sim.secr.R and editted
## 2024-07-30 expected
## 2024-07-31 dropzeroCH
## 2024-08-02 using revised output order from Rcpp functions (i,s,k)
## 2024-08-21 edited for style
## 2025-03-14 tweaks to handling of miscparm
################################################################################

## utility function used only in sim.detect
## mixture proportions for representative animal 
## (i.e. assuming no individual or group variation)       
## assume dim(PIA)[1] == 1
getpmixall <- function(PIA, realparval)
{
    nc <- dim(PIA)[2]
    k <- dim(PIA)[4]
    nmix <- dim(PIA)[5]
    pmix <- 1   ## default single class
    if (nmix>1) {
        # index of first non-missing occasion s and detector k
        fsk <- firstsk(PIA[1,1,,,1, drop = FALSE])
        kc <- as.vector((fsk-1) %/% k + 1)[1]
        sc <- as.vector((fsk-1) %/% k + 1)[1]
        pmix <- PIA[cbind(1,1,sc,kc,1:nmix)]
    }
    pmix
}
################################################################################

sim.detect <- function (object, popnlist, maxperpoly = 100, renumber = TRUE, 
                        expected = FALSE, dropzeroCH = TRUE)
    ## popnlist is always a list of popn objects, one per session
{
    ## we use fake CH to extract parameter value dependent on prev capt
    ## BUT this is slow and clunky...
    dummycapthist<- function (capthist, pop, fillvalue = 1) {
        if (ms(capthist)) {
            output <- list()
            for (i in 1:nsession)
                output[[i]] <- dummycapthist (capthist[[i]], pop = pop[i], 
                                              fillvalue = fillvalue)
            class(output)   <- c('capthist', 'list')
            session(output) <- session(capthist)   
            output
        }
        else {
            newdim <- dim(capthist)
            newdim[1] <- nrow(pop[[1]])
            output <- array(fillvalue, dim = newdim)
            ## CAPTURE ON LAST OCCASION REGARDLESS OF FILLVALUE
            ## trick to keep array valid without misleading
            output[,,newdim[3]] <- 1
            class(output) <- 'capthist'
            traps(output) <- traps(capthist)
            session(output) <- session(capthist)
            covariates(output) <- covariates(pop[[1]])
            output
        }
    }
    ## --------------------------------------------------------------------
    ## Exclusions
    if (!is.null(object$groups) && (object$details$param>1))
        stop("simulation does not extend to groups when param>1")
    
    detectors <- unique(unlist(detector(traps(object$capthist))))
    supported <- c('single', 'multi', 'proximity', 'count', 'polygonX',
                   'transectX', 'signal', 'polygon','transect')
    unsupported <- detectors[!detectors %in% supported]
    if (length(unsupported) > 0) {
        unsupported <- paste(unsupported, collapse = ', ')
        stop("detector type ', unsupported, ' is not supported in sim.detect")
    }
    ## --------------------------------------------------------------------
    ## process behavioural responses
    Markov <- any(c('B','Bk','K') %in% object$vars)
    btype <- which (c("b", "bk", "k") %in% tolower(object$vars))
    if (length(btype) > 1)
        stop ("all behavioural responses must be of same type in sim.detect")
    if (length(btype) == 0)
        btype <- 0
  
    ## --------------------------------------------------------------------
    ## setup
    if (is.null(object$details$ignoreusage)) {
        object$details$ignoreusage <- FALSE
    }
    # 2025-03-14 set nmiscparm here
    if (is.null(object$details$miscparm)) {
        nmiscparm <- 0
        object$details$miscparm <- numeric(4)
    }
    else {
        nmiscparm <- length(object$details$miscparm)
    }
    N <- sapply(popnlist, nrow)
    nocc <- if(ms(object)) sapply(object$capthist, ncol) else ncol(object$capthist)
    ndet <- if(ms(object)) sapply(traps(object$capthist), nrow) else nrow(traps(object$capthist))
    sessionlevels <- session(object$capthist)   
    nsession <- length(sessionlevels)
    sessmask <- object$mask
    if (!ms(sessmask)) sessmask <- list(sessmask)  ## always a list
    grplevels <- group.levels(object$capthist, object$groups, sep=".")
    beta <- object$fit$par
    #2025-03-13
    # userd <- is.null(object$details$userdist)
    userd <- !is.null(object$details$userdist) && is.function(object$details$userdist)
    
    ncores <- setNumThreads()
    grain <- if (ncores==1) 0 else 1
    
    ##------------------------------------------------------------------
    ## detection parameters for each session, animal, occasion, detector
    ## realparval is lookup table of values,
    ## design$PIA is index to lookup
    
    ##----------------------------------------
    ## real parameter values for naive animals
    
    ## using existing design to save time in secr.design.MS...
    if (length(unique(object$design0$PIA)) == 1) {
        design0 <- object$design0
        dim0 <- dim(design0$PIA)
        design0$PIA <- array (1, dim = c(dim0[1], max(N), dim0[3:5]))
        dummyCH <- NULL
    }
    else {
        dummyCH <- dummycapthist(object$capthist, popnlist)
        design0 <- secr.design.MS (
            dummyCH, 
            object$model, 
            object$timecov, 
            object$sessioncov,
            object$groups, 
            object$hcov, 
            object$dframe, 
            naive       = TRUE, 
            CL          = object$CL,
            ignoreusage = object$details$ignoreusage, 
            contrasts   = object$details$contrasts)
    }
    
    ## uniform over individuals
    ## C++ code uses individual-specific PIA even in full-likelihood case
    for (i in 1:nsession) {
        if (is.null(object$groups)) {
            ## all animals the same...
            matchedgroup <- rep(1, N[i])
        }
        else {
            ## group identities for simulated animals
            ## can treat popn as capthist in group.factor if not ms(object)
            newgroup <- group.factor (popnlist[[i]], object$groups)
            ## here we assume original PIA is for a full-likelihood model with
            ## one row per group
            matchedgroup <- (1:nlevels(newgroup))[as.numeric(newgroup)]
        }
        design0$PIA[i,1:N[i],,,] <- object$design0$PIA[i,matchedgroup,,,]
    }
    realparval0 <- makerealparameters (design0, beta, object$parindx, object$link,
                                       object$fixed)
    ##----------------------------------------
    ## real parameter  values for 'caughtbefore'  animals or detectors
    ## -- this  definition of  design1 differs  from that in secr.fit()
    ## where  parameter  values  are  appropriate to  the  particular
    ## (realised) detection histories
    
    if (btype > 0) {
        if (is.null(dummyCH)) 
            dummyCH <- dummycapthist(object$capthist, popnlist)
        design1 <- secr.design.MS (
            dummyCH, 
            object$model, 
            object$timecov, 
            object$sessioncov,
            object$groups, 
            object$hcov, 
            object$dframe,
            naive       = FALSE,
            CL          = object$CL, 
            ignoreusage = object$details$ignoreusage, 
            contrasts   = object$details$contrasts)
        realparval1 <- makerealparameters (
            design1, 
            beta, 
            object$parindx, 
            object$link,
            object$fixed)  # caught before
    }
    else {   ## faster to just copy if no behavioural response
        design1 <- design0
        realparval1 <- realparval0
    }
    ##------------------------------------------
    
    # see loglikhelperfn.R for getD
    # Density
    if (!object$CL ) {
        D <- getD (object$designD, beta, sessmask, 
                   object$parindx, object$link, object$fixed,
                   grplevels, sessionlevels, 
                   parameter = 'D')
    }
    # Non-Euclidean distance parameter
    NE <- getD (object$designNE, beta, sessmask, 
                object$parindx, object$link, object$fixed,
                grplevels, sessionlevels, 
                parameter = 'noneuc')
    
    #--------------------------------------------------------------------
    
    output <- list()
    for (sessnum in 1:nsession) {
        
        ## in multi-session case must get session-specific data from lists
        if (ms(object)) {
            s <- ncol(object$capthist[[sessnum]])
            session.traps    <- traps(object$capthist)[[sessnum]]
            session.mask <- if (userd || (object$details$param %in% c(2,6)))
                object$mask[[sessnum]] else NULL
            Dtemp <- if (object$details$param %in% c(4:6))
                predict(object)[[sessnum]]['D','estimate']
            else NA
            nmash <- attr(object$capthist[[sessnum]], 'n.mash')
        }
        else {
            s <- ncol(object$capthist)
            session.traps    <- traps(object$capthist)
            session.mask <- if (userd || (object$details$param %in% c(2,6)))
                object$mask else NULL
            Dtemp <- if (object$details$param %in% c(4:6)) predict(object)['D','estimate']
                        else NA
            nmash <- attr(object$capthist, 'n.mash')
        }
     
        ## function reparameterize is in reparameterize.R
        Xrealparval0 <- reparameterize (realparval0, object$detectfn, object$details,
                                        session.mask, session.traps, Dtemp, s)
        Xrealparval1 <- reparameterize (realparval1, object$detectfn, object$details,
                                        session.mask, session.traps, Dtemp, s)
        ##------------------------------------------
        session.animals <- popnlist[[sessnum]]
        
        ## pass miscellaneous unmodelled parameter(s)
        if (nmiscparm > 0) {
            miscindx <- max(unlist(object$parindx)) + (1:nmiscparm)
            attr(session.mask, 'miscparm') <- coef(object)[miscindx, 1]
        }
        # otherwise attr(session.mask, 'miscparm') remains at pre-set value(s), if set
        
        ##------------------------------------------
        
        dettype <- detectorcode(session.traps, MLonly = FALSE, noccasions = nocc[sessnum])
        binomN <- expandbinomN(object$details$binomN, dettype)
        if (all(detector(session.traps) %in% .localstuff$polydetectors)) {
            k <- c(table(polyID(session.traps)),0)
            K <- length(k)-1
        }
        else {
            k <- nrow(session.traps)
            K <- k
        }
        
        #------------------------------------------
        ## simulate this session...
        usge <- usage(session.traps)
        
        if (is.null(usge) || object$details$ignoreusage)
            usge <- matrix(1, nrow = K, ncol = s)
        NR <- N[sessnum]
        
        if (all(detector(session.traps) %in% .localstuff$exclusivedetectors)) {
            maxdet <- NR * s
        }
        else {
            # safety margin : average detections per animal per detector per occasion
            maxdet <- NR * s * K * maxperpoly
        }
        if ((object$detectfn==12) || (object$detectfn==13)) {
            ## muN, sdN
            object$details$miscparm[2:3] <- beta[max(unlist(object$parindx))+(1:2)]
        }
        ## requires session.animals has covariate named 'hcov'
        knownclass <- getknownclass(session.animals, object$details$nmix, object$hcov)
        
        ## get fixed pmix for each class
        pmix <- getpmixall(design0$PIA, Xrealparval0)
        
        #---------------------------------------------------------------
        ## TO BE FIXED
        if (length(unique(dettype))>1 || length(unique(binomN))>1)
            stop("simulation not yet updated for varying detector type")
        #---------------------------------------------------------------
        
        ## pre-compute distances from detectors to animals
        
        if (all(dettype %in% c(-1,0,1,2,5,8))) {
            if (is.function(object$details$userdist)) {
                ## only use of NE in this fn
                noneuc <- getmaskpar(!is.null(NE), NE, nrow(session.mask), sessnum, FALSE, NULL)
                density <- getmaskpar(!object$CL, D, nrow(session.mask), sessnum, 
                                      object$details$unmash, nmash)
                distmat2 <- getuserdist(
                    session.traps, 
                    session.animals, 
                    object$details$userdist, 
                    sessnum, 
                    noneuc[,1], 
                    density[,1], 
                    object$details$miscparm)
            }
            else {
                distmat2 <- getdistmat2 (
                    session.traps, 
                    session.animals, 
                    object$details$userdist)
            }
            ## precompute gk, hk for point detectors
            gkhk0 <- makegkPointcpp (
                as.integer(object$detectfn),
                as.integer(grain),
                as.integer(ncores),
                as.matrix(Xrealparval0),
                as.matrix(distmat2),
                as.double(object$details$miscparm))
            gkhk <- makegkPointcpp (
                as.integer(object$detectfn),
                as.integer(grain),
                as.integer(ncores),
                as.matrix(Xrealparval1),
                as.matrix(distmat2),
                as.double(object$details$miscparm))
            if (any(dettype == 8)) {
                stop("sim.secr not ready for capped detectors")
            }
        }
        
        expectedarray <- NA
        
        if (all(dettype %in% c(-1,0,1,2,8))) {
            temp <- simdetectpointcpp (
                as.integer (dettype[1]),           ## detector -1 single, 0 multi, 1 proximity, 2 count,...
                as.integer (NR),
                as.integer (nrow(Xrealparval0)),   ## cc0
                as.integer (nrow(Xrealparval1)),   ## cc
                as.double  (gkhk0$gk),
                as.double  (gkhk$gk),
                as.double  (gkhk0$hk),
                as.double  (gkhk$hk),
                as.integer (design0$PIA[sessnum,1:NR,1:s,1:K,]),       ## naive
                as.integer (design1$PIA[sessnum,1:NR,1:s,1:K,]),
                as.integer (object$details$nmix),
                as.integer (knownclass),
                as.double  (pmix),
                as.matrix  (usge),
                as.integer (btype),      ## code for behavioural response  0 none etc.
                as.integer (Markov),     ## learned vs transient behavioural response 0 learned 1 Markov
                as.integer (binomN)      ## number of trials for 'count' detector modelled with binomial
            )
            if (expected) {
                expectedarray <- expdetectpointcpp (
                    as.integer (dettype[1]),           ## detector 0 multi, 1 proximity, 2 count,...
                    as.integer (NR),
                    as.integer (nrow(Xrealparval0)),   ## cc0
                    as.integer (nrow(Xrealparval1)),   ## cc
                    as.double  (gkhk0$gk),
                    as.double  (gkhk$gk),
                    as.double  (gkhk0$hk),
                    as.double  (gkhk$hk),
                    as.integer (design0$PIA[sessnum,1:NR,1:s,1:K,]),       ## naive
                    as.integer (design1$PIA[sessnum,1:NR,1:s,1:K,]),
                    as.integer (object$details$nmix),
                    as.integer (knownclass),
                    as.double  (pmix),
                    as.matrix  (usge),
                    as.integer (btype),      ## code for behavioural response  0 none etc.
                    as.integer (Markov),     ## learned vs transient behavioural response 0 learned 1 Markov
                    as.integer (binomN)      ## number of trials for 'count' detector modelled with binomial
                )
                expectedarray <- array(expectedarray, dim = c(NR, s, K))
            }
        }
        else if (all(dettype %in% c(3,4,6,7))) {
            temp <- simdetectpolycpp (
                as.integer (dettype[1]),       ## detector 3 exclusive polygon, 4 exclusive transect, 6 polygon, 7 transecr
                as.integer (object$detectfn),  ## code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 = uniform
                as.integer (object$details$nmix), ## number of classes
                as.integer (btype),            ## code for behavioural response  0 none etc.
                as.integer (Markov),           ## learned vs transient behavioural response 0 learned 1 Markov
                as.integer (k),                ## number of vertices per polygon (zero-terminated vector)
                as.matrix  (session.animals),   ## x,y points of animal range centres (first x, then y)
                as.matrix  (session.traps),     ## x,y locations of traps (first x, then y)
                as.matrix  (Xrealparval0),      ## Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal]
                as.matrix  (Xrealparval1),      ## Parameter values (matrix nr= comb of g0,sigma,b nc=3) [caught before]
                as.integer (design0$PIA[sessnum,1:NR,1:s,1:K,]), ## lookup which g0/sigma/b combination to use for given g, S, K [naive animal]
                as.integer (design1$PIA[sessnum,1:NR,1:s,1:K,]), ## lookup which g0/sigma/b combination to use for given n, S, K  [caught before]
                as.integer (knownclass),       ## known membership of 'latent' classes
                as.double  (pmix),             ## membership probabilities
                as.matrix  (usge),              ## ss x kk array of 0/1 usage codes or effort
                as.integer (binomN),            ## number of trials for 'count' detector modelled with binomial
                as.integer (maxperpoly)        ##
            )
            if ((temp$resultcode == 2) && (any(dettype %in% c(6,7)))) {
                stop (">100 detections per animal per polygon per occasion")
            }
        }
        else
            if (all(dettype %in% c(5))) {
                if (is.null(object$details$cutval))
                    stop ("sim.detect for signal model requires object$details$cutval")
                temp <- simdetectsignalcpp (
                    as.integer(dettype[1]),
                    as.integer(object$details$nmix),
                    as.integer(object$detectfn),
                    as.integer(object$details$cutval),
                    as.matrix(Xrealparval0),
                    as.integer(design0$PIA[sessnum,1:NR,1:s,1:K,]),  # index of N,S,K,x to rows in Xrealparval0
                    as.integer(pmix),
                    as.integer(knownclass),
                    as.matrix(session.animals),
                    as.matrix(session.traps),
                    as.matrix(distmat2),
                    as.matrix(usge),
                    as.double(object$details$miscparm)
                )
            }
        if (temp$resultcode != 0) {
            stop ("simulated detection failed, code ", temp$resultcode)
        }

        w <- array(temp$value, dim = c(NR, s, K), dimnames = list(NULL,1:s,NULL))
        included <- if (dropzeroCH) apply(w,1,sum)>0 else rep(TRUE, nrow(w))
        w <- w[included,,, drop = FALSE] 

        class(w)   <- 'capthist'    # NOT data.frame
        traps(w)   <- session.traps
        session(w) <- sessionlevels[sessnum]
        
        if (!is.null(covariates(popnlist))) {
            covariates(w) <- subset(covariates(popnlist[[sessnum]]), 
                                    subset = included)
        }
        
        ##---------------------
        ## signal
        if (any(dettype %in% c(5,12)) && (temp$n>0)) {
            nd <- sum(abs(w))
            signal(w) <- temp$signal[1:nd]
            if ((object$detectfn==12) || (object$detectfn==13))
                noise(w) <- temp$signal[(nd+1):(2*nd)]
            attr(w, 'cutval') <- object$details$cutval
        }
        else {
            attr(w, 'signalframe') <- NULL
            attr(w, 'cutval') <- NULL
        }
        ##---------------------
        ## polygon or transect
        if (any(dettype %in% c(3,4,6,7))) {
            nd <- sum(abs(w))
            if (nd>0)
                xymat <- temp$detectedXY[1:nd, 1:2, drop = FALSE]
            else 
                xymat <- matrix(nrow=0, ncol=2)
            detectedXY <- data.frame(xymat)
            names(detectedXY) <- c('x','y')
            attr(w, 'detectedXY') <- detectedXY
        }
        else
            attr(w, 'detectedXY') <- NULL
        ##---------------------
        
        if (renumber && (temp$n>0)) {
            rownames(w) <- 1:nrow(w)
        }
        else {
            rownames(w) <- rownames(session.animals)[included]
        }
        output[[sessnum]] <- w
        
        if (expected) {
            exparray <- expectedarray[included,,,drop = FALSE]
            rownames(exparray) <- rownames(w)
            attr(output[[sessnum]], 'expected') <- exparray
        }
        
    } ## end of session loop
    
    if (nsession == 1) output <- output[[1]]
    else {
        names(output) <- sessionlevels
        class(output) <- c('capthist', 'list')
    }
    output
}
################################################################################
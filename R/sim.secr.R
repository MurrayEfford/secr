############################################################################################
## package 'secr'
## sim.secr.R
## uses Dsurface.R, utility.R, secrloglik.R
## Copyright (C) 2019 Murray Efford
## 2019-10-05 converted to cpp

## 2019-10-06 yet to fix: groups?

## 2019-10-06 excludes: knownclass (hcov, h2)
## 2019-10-06 excludes: mixed detector types
## 2019-10-06 excludes: individual covariates
## 2019-10-06 excludes: telemetry
## 2020-11-08 bad rownames for sim.detect when renumber = FALSE
############################################################################################

disinteraction <- function (capthist, groups, sep='.') {
  ngv <- length(groups)
  grouplevels <- group.levels(capthist, groups, sep=sep)
  if (ngv>1)
    temp <- matrix(unlist(strsplit(as.character(grouplevels), sep, fixed=TRUE)),
                   byrow = T, ncol = ngv)
  else temp <- grouplevels
  temp <- data.frame(temp)
  names(temp) <- groups
  temp
}

############################################################################################

simulate.secr <- function (object, nsim = 1, seed = NULL, maxperpoly = 100, chat = 1,
                           ...)
## if CL, condition on n? what about distribution of covariates over n?
## locate according to IHP with lambda(X) controlled by f(X|covar), assuming homog Poisson
## i.e. use f(X|covar)/max(f(X|covar)) to reject until meet quota n?
## or f(X, covar | detected)?
## TOO HARD - cf MARK

## 2012-10-25
## other possible exclusions:
## mashed?

{
    ##  check input
    if (any(c("bn", "bkn", "bkc", "Bkc") %in% tolower(object$vars)))
        stop ("simulate works only with binary behavioural responses")

    if (!inherits(object,'secr'))
        stop ("requires 'secr' object")
    if (object$CL)
        stop ("not implemented for conditional likelihood")

    ## density array dim(mask points, groups, sessions)
    Darray <- getDensityArray (predictDsurface(object))

    ## setup

    ngrp <- dim(Darray)[2]
    nsession <- dim(Darray)[3]
    if (!is.null(object$groups)) {
        ## individual covariates for foundation of g
        di <- disinteraction (object$capthist, object$groups)
    }

    sesscapt <- vector('list', nsim)

    ##################
    ## set random seed
    ## copied from simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    ##################

    ## loop over replicates
    runone <- function(r) {
        ## argument r is unused dummy
        sesspopn <- list()
        for (sessnum in 1:nsession) {
            if (nsession==1) mask <- object$mask
            else mask <- object$mask[[sessnum]]
            popn <- list()
            for (g in 1:ngrp) {
                if (inherits(mask, 'linearmask')) {
                    density <- Darray[,g,sessnum]        ## vector
                    mod2D <- 'linear'                    ## linear Poisson or IHP
                }
                else {
                    if (object$model$D == ~1) {
                        density <- Darray[1,g,sessnum]   ## scalar
                        mod2D <- 'poisson'               ## homogeneous
                    }
                    else {
                        density <- Darray[,g,sessnum]    ## vector
                        mod2D <- 'IHP'                   ## inhomogeneous
                    }
                }
                if (chat > 1)
                    density <- density / chat
                ND <- switch (object$details$distribution,
                              binomial = 'fixed',
                              poisson = 'poisson',
                              'poisson')
                ##-------------------------------------------------------------------
                ## sim.popn arguments omitted:
                ## buffer          redundant when mask specified
                ## buffertype      ditto
                ## poly            ditto
                ## covariates      ignored for now...
                ## number.from = 1 fine
                ## nsession = 1    fine here within session loop as long
                ##                 as there is no turnover model
                ## details = NULL  not relevant (turnover and special model2D only)
                ## seed = NULL     default mechanism -- needs attention 2014-09-07
                popn[[g]] <- sim.popn (D = density, core = mask, model2D = mod2D, Ndist = ND)

                ## ------------------------------------------------------------------
                ## Add any needed covariates, first generating a clean dataframe
                if (is.null(object$groups) & is.null(object$hcov))
                    covariates(popn[[g]]) <- NULL
                else
                    covariates(popn[[g]]) <- data.frame(row.names = row.names(popn[[g]]))

                ## ---groups---
                if (!is.null(object$groups)) {
                    grpcov <- as.data.frame(di[rep(g, nrow(popn[[g]])),]) ## 2014-08-08
                    names(grpcov) <- object$groups
                    covariates(popn[[g]]) <- cbind(covariates(popn[[g]]),grpcov)
                }
                ## ---hcov---
                ## sample with replacement from original hcov field 2014-08-08
                if (!is.null(object$hcov)) {
                    if (ms(object))
                        oldhcov <- covariates(object$capthist[[sessnum]])[,object$hcov]
                    else
                        oldhcov <- covariates(object$capthist)[,object$hcov]
                    covariates(popn[[g]])[object$hcov] <-
                        sample(oldhcov, size = nrow(popn[[g]]), replace = TRUE)
                }
                ## ------------------------------------------------------------------
            }
            sesspopn[[sessnum]] <- do.call(rbind, popn)   ## combine groups in one popn object
        }
        sim.detect(object, sesspopn, maxperpoly)
    }
    sesscapt <- lapply(1:nsim, runone)

    ## experimental
    if (chat>1) sesscapt <- lapply(sesscapt, replicate, chat)

    attr(sesscapt,'seed') <- RNGstate   ## save random seed
    class(sesscapt) <- c('secrdata', 'list')
    sesscapt
}
############################################################################################

sim.secr <- function (object, nsim = 1, extractfn = function(x)
    c(deviance=deviance(x), df=df.residual(x)), seed = NULL,
    maxperpoly = 100, data = NULL, tracelevel = 1, hessian =
    c('none','auto','fdHess'), start = object$fit$par, ncores = NULL) {

## parametric bootstrap simulations based on a fitted secr object
## extractfn is a required function to extract values from an secr fit
## it should return a vector of named values that does not vary in length
## 'hessian' overrides value in object$details
## set hessian='auto' if extractfn requires variance-covariance matrix

    hessian <- match.arg(hessian)
    cl <- match.call(expand.dots = TRUE)
    cl <- paste(names(cl)[-1],cl[-1], sep=' = ', collapse=', ' )
    cl <- paste('sim.secr(', cl, ')')

    if (is.null(extractfn)) extractfn <- trim
    test <- extractfn(object)

    if (is.numeric(test)) {
        n.extract <- length(test)
        if (n.extract<=0)
            stop ("invalid 'extractfn'")
    }
    detectnames <- names(object$design0[[1]])   ## names of real detection parameters
    details <- replace(object$details, 'hessian', hessian)
    tracelevel <- as.integer(tracelevel)
    details$trace <- tracelevel > 1
    min.detections <- 1
    .localstuff$iter2 <- 0   

    if (any(detector(traps(object$capthist)) == 'single')) {
        memo ('multi-catch likelihood used for single-catch traps', tracelevel>0)
    }
    
    if (is.null(data)) {
        memo ('sim.secr simulating detections...', tracelevel>0)
        ## use multiple cores only for model fitting 2015-12-02
        data <- simulate(object, nsim = nsim, seed = seed, maxperpoly = maxperpoly)
    }
    else {
        if (any(class(data) != c('secrdata', 'list')))
            stop("invalid data")
    }
    fitmodel <- function (sc) {
        ## i <<- i+1
        ## use second counter so as not to interfere with secrloglik
        .localstuff$iter2 <- .localstuff$iter2 + 1   ## 2016-01-09
        if (tracelevel>0)
            memo (paste('sim.secr fitting replicate', .localstuff$iter2, '...'), TRUE)
        nc <-  sum(counts(sc)$'M(t+1)'[,'Total'])
        if (nc >= min.detections) {
            tempfit <- suppressWarnings( secr.fit(sc, model = object$model, mask = object$mask,
                CL = object$CL, detectfn = object$detectfn, binomN = details$binomN,
                start = start, link = object$link, fixed = object$fixed,
                timecov = object$timecov, sessioncov = object$sessioncov,
                groups = object$groups, dframe = object$dframe, details = details,
                method = object$fit$method, verify = FALSE, biasLimit = NA,
                ncores = ncores) )
            extractfn(tempfit)
        }
        else if (is.list(test)) list() else rep(NA, n.extract)
    }
    
    output <-  lapply (data, fitmodel)

    if (is.numeric(test)) {
        output <- do.call(rbind, output)
        output <- data.frame(output)
    }
    else {
        class(output) <- c('secrlist', 'list')
    }

    attr(output,'seed') <- attr(data,'seed')
    attr(output,'call') <- cl
    attr(output,'extractfn') <- extractfn
    output
}
############################################################################################

print.secrdata <- function(x,...) {
## suggestion of Rolf Turner 19 Jan 2009 for printing without attributes
    attributes(x) <- NULL
    print(x)
}
############################################################################################

print.secrlist <- function(x,...) {
## suggestion of Rolf Turner 19 Jan 2009 for printing without attributes
    attributes(x) <- NULL
    print(x)
}
############################################################################################

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

sim.detect <- function (object, popnlist, maxperpoly = 100, renumber = TRUE)
## popnlist is always a list of popn objects, one per session
{
    ## we use fake CH to extract parameter value dependent on prev capt
    ## BUT this is slow and clunky...
    dummycapthist<- function (capthist, pop, fillvalue=1) {
        if (ms(capthist)) {
            output <- list()
            for (i in 1:nsession)
                output[[i]] <- dummycapthist (capthist[[i]],
                    pop = pop[i], fillvalue = fillvalue)
            class(output) <- c('capthist', 'list')
            session(output) <- session(capthist)   ## 2010 03 10
            output
        }
        else {
            newdim <- dim(capthist)
            newdim[1] <- nrow(pop[[1]])
            output <- array(fillvalue, dim = newdim)
            ## CAPTURE ON LAST OCCASION
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
    ## see also below dettype %in% c(-1,0,1,2,5,6,7)
    if (!is.null(object$groups) & (object$details$param>1))
        stop("simulation does not extend to groups when param>1")

    if ('telemetry' %in% unlist(detector(traps(object$capthist))))
        stop("telemetry models are not supported in simulate.secr")
    
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
    if (is.null(object$details$ignoreusage))
        object$details$ignoreusage <- FALSE
    if (is.null(object$details$miscparm))
        object$details$miscparm <- numeric(4)
    N <- sapply(popnlist, nrow)
    nocc <- if(ms(object)) sapply(object$capthist, ncol) else ncol(object$capthist)
    ndet <- if(ms(object)) sapply(traps(object$capthist), nrow) else nrow(traps(object$capthist))
    sessionlevels <- session(object$capthist)   ## was names(capthist) 2009 08 15
    nsession <- length(sessionlevels)
    sessmask <- object$mask
    if (!ms(sessmask)) sessmask <- list(sessmask)  ## always a list
    grplevels <- group.levels(object$capthist, object$groups, sep=".")
    beta <- object$fit$par
    userd <- is.null(object$details$userdist)
    
    ncores <- setNumThreads()
    grain <- if (ncores==1) 0 else 1
    
    ##------------------------------------------
    ## detection parameters for each session, animal, occasion, detector
    ## realparval is lookup table of values,
    ## design$PIA is index to lookup

    ## real parameter values for naive animals

    ## using existing design to save time in secr.design.MS...
    if (length(unique(object$design0$PIA)) == 1) {
        design0 <- object$design0
        dim0 <- dim(design0$PIA)
        design0$PIA <- array (1, dim = c(dim0[1], max(N), dim0[3:5]))
        dummyCH <- NULL
    }
    else {
        dummyCH <- dummycapthist(object$capthist, popnlist, fillvalue = 1)
        design0 <- secr.design.MS (dummyCH, object$model, object$timecov, object$sessioncov,
                                   object$groups, object$hcov, object$dframe, 
                                   naive = TRUE, CL = object$CL,
                                   ignoreusage = object$details$ignoreusage, 
                                   contrasts = object$details$contrasts)
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

    ## real parameter  values for 'caughtbefore'  animals or detectors
    ## -- this  definition of  design1 differs  from that in secr.fit()
    ## where  parameter  values  are  approppriate to  the  particular
    ## (realised) detection histories

    if (btype > 0) {
        if (is.null(dummyCH))
            dummyCH <- dummycapthist(object$capthist, popnlist, fillvalue = 1)
        design1 <- secr.design.MS (dummyCH, object$model, object$timecov, object$sessioncov,
            object$groups, object$hcov, object$dframe, ignoreusage = object$details$ignoreusage, 
            contrasts = object$details$contrasts, CL = object$CL)
        realparval1 <- makerealparameters (design1, beta, object$parindx, object$link,
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
            session.mask <- if (userd | (object$details$param %in% c(2,6)))
                    object$mask[[sessnum]] else NULL
            Dtemp <- if (object$details$param %in% c(4:6))
                    predict(object)[[sessnum]]['D','estimate']
                else NA
            nmash <- attr(object$capthist[[sessnum]], 'n.mash')
        }
        else {
            s <- ncol(object$capthist)
            session.traps    <- traps(object$capthist)
            session.mask <- if (userd | (object$details$param %in% c(2,6)))
                object$mask else NULL
            Dtemp <- if (object$details$param %in% c(4:6)) predict(object)['D','estimate']
                else NA
            nmash <- attr(object$capthist, 'n.mash')
        }

        ## function reparameterize is in secrloglik.R
        Xrealparval0 <- reparameterize (realparval0, object$detectfn, object$details,
                                        session.mask, session.traps, Dtemp, s)
        Xrealparval1 <- reparameterize (realparval1, object$detectfn, object$details,
                                        session.mask, session.traps, Dtemp, s)
        ##------------------------------------------
        session.animals <- popnlist[[sessnum]]

        ## pre-compute distances from detectors to animals
        ## pass miscellaneous unmodelled parameter(s)
        nmiscparm <- length(object$details$miscparm)
        if (nmiscparm > 0) {
            miscindx <- max(unlist(object$parindx)) + (1:nmiscparm)
            attr(session.mask, 'miscparm') <- coef(object)[miscindx, 1]
        }
        
        ##------------------------------------------

        dettype <- detectorcode(session.traps, MLonly = FALSE, noccasions = nocc[sessnum])
        if (! all(dettype %in% c(-1,0,1,2,3,4,5,6,7)))
            stop ("detector type ",
                  paste(detector(session.traps)),
                  " not implemented")

        binomN <- expandbinomN(object$details$binomN, dettype)
        
        if (all(detector(session.traps) %in% .localstuff$polydetectors)) {
            k <- c(table(polyID(session.traps)),0)
            K <- length(k)-1
        }
        else {
            k <- nrow(session.traps)
            K <- k
        }
        ## no obvious reason to retain this; using sessnum 2019-10-20
        ## sessg <- min (sessnum, design0$R)

        #------------------------------------------
        ## simulate this session...
        usge <- usage(session.traps)

        if (is.null(usge) | object$details$ignoreusage)
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

        ## TO BE FIXED
        if (length(unique(dettype))>1 | length(unique(binomN))>1)
            stop("simulation not yet updated for varying detector type")
        
        if (all(dettype %in% c(-1,0,1,2,5,8))) {
            if (is.function(object$details$userdist)) {
                ## only use of NE in this fn
                noneuc <- getmaskpar(!is.null(NE), NE, nrow(session.mask), sessnum, FALSE, NULL)
                density <- getmaskpar(!object$CL, D, nrow(session.mask), sessnum, 
                                      object$details$unmash, nmash)
                distmat2 <- getuserdist(session.traps, session.animals, object$details$userdist, sessnum, 
                    noneuc[,1], density[,1], object$details$miscparm, object$detectfn==20)
            }
            else {
                distmat2 <- getdistmat2 (
                    session.traps, 
                    session.animals, 
                    object$details$userdist,
                    object$detectfn == 20)
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
                #     gkhk <- cappedgkhkcpp (as.integer(nrow(Xrealparval1)), as.integer(nrow(session.traps)),
                #                            as.double(getcellsize(session.mask)), as.double(density[,1]), 
                #                            as.double(gkhk$gk), as.double(gkhk$hk))  
            }
        }
  
        if (all(dettype %in% c(-1,0,1,2,8))) {
            temp <- simdetectpointcpp (dettype[1],      ## detector -1 single, 0 multi, 1 proximity, 2 count,... 
                                       as.integer(NR), 
                                       as.integer(nrow(Xrealparval1)),
                                       as.double(gkhk0$gk), 
                                       as.double(gkhk$gk), 
                                       as.double(gkhk0$hk), 
                                       as.double(gkhk$hk), 
                                       as.integer(design0$PIA[sessnum,1:NR,1:s,1:K,]),       ## naive
                                       as.integer(design1$PIA[sessnum,1:NR,1:s,1:K,]),       
                                       as.integer(object$details$nmix),
                                       as.integer(knownclass),
                                       as.double(pmix),
                                       as.matrix(usge),
                                       as.integer(btype),       ## code for behavioural response  0 none etc. 
                                       as.integer(Markov),      ## learned vs transient behavioural response 0 learned 1 Markov 
                                       as.integer(binomN)      ## number of trials for 'count' detector modelled with binomial 
            )
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
            if ((temp$resultcode == 2) & (any(dettype %in% c(6,7)))) {
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
        w <- array(temp$value, dim=c(s, K, NR), dimnames = list(1:s,NULL,NULL))
        w <- aperm(w, c(3,1,2))
        w <- w[apply(w,1,sum)>0,,, drop = FALSE] 
            
        class(w)   <- 'capthist'    # NOT data.frame
        traps(w)   <- session.traps
        session(w) <- sessionlevels[sessnum]

        if (!is.null(covariates(popnlist))) {
            covariates(w) <- subset(covariates(popnlist[[sessnum]]), subset =
                as.logical(temp$caught))
        }
        
        ##---------------------
        ## signal
        if (any(dettype %in% c(5,12)) & (temp$n>0)) {
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

        if (renumber & (temp$n>0)) {
            rownames(w) <- 1:temp$n
        }
        else {
            ## rownames(w) <- (1:N[sessnum])[as.logical(temp$caught)]
            ## changed 2020-11-08
            rownames(w) <- rownames(session.animals)[as.logical(temp$caught)]
        }
        output[[sessnum]] <- w
        
    } ## end of session loop

    if (nsession == 1) output <- output[[1]]
    else {
        names(output) <- sessionlevels
        class(output) <- c('capthist', 'list')
    }
    output
}
############################################################################################
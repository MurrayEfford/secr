## package 'secr'
## secr.fit.R

## 2019-12-03 secr.design.MS uses CL argument
## 2020-08-30 check3D restored for secrlinear arvicola example
## 2021-04-02 allcapped bug fixed  (cannot combine capped, uncapped)
## 2021-04-25 4.4.0
## 2021-06-22 global change fixedpar to fixed for consistency
###############################################################################

secr.fit <- function (capthist,  model = list(D~1, g0~1, sigma~1), mask = NULL,
                      buffer = NULL, CL = FALSE, detectfn = NULL, binomN = NULL, start = NULL,
                      link = list(), fixed = list(), timecov = NULL, sessioncov = NULL, hcov = NULL,
                      groups = NULL, dframe = NULL, details = list(), method = 'Newton-Raphson',
                      verify = TRUE, biasLimit = 0.01, trace = NULL, ncores = NULL, ...)
    
{
    # Fit spatially explicit capture recapture model
    #
    # Arguments:
    #
    #  capthist   -  capture history object (includes traps object as an attribute)
    #  model      -  formulae for real parameters in terms of effects and covariates
    #  mask       -  habitat mask object
    #  buffer     -  default buffer width should mask not be provided
    #  CL         -  logical switch : conditional likelihood (T) or full likelihood (F)
    #  detectfn   -  code for detection function 0 = halfnormal, 1 = hazard, 2 = exponential etc.
    #  start      -  start values for maximization (numeric vector link scale);
    #                if NULL then 'autoini' function is used
    #  link       -  list of parameter-specific link function names 'log', 'logit', 'identity',
    #                'sin', 'neglog'
    #  fixed      -  list of fixed values for named parameters
    #  timecov    -  data for time covariates if these are used in 'model'
    #  sessioncov -  dataframe of session-level covariates
    #  groups     -  vector of names to group fields in attr(capthist,'covariates') dataframe
    #  dframe     -  optional data frame of design data for detection model (tricky & untested)
    #  details    -  list with several additional settings, mostly of special interest
    #  method     -  optimization method (indirectly chooses
    #  verify     -  logical switch for pre-check of capthist and mask with verify()
    #  trace      -  logical; if TRUE output each likelihood as it is calculated
    #  ...        -  other arguments passed to nlm() or optim()
    
    #################################################
    ## Remember start time
    ptm  <- proc.time()
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    # gc(verbose = FALSE) ## garbage collection 2019-09-02
    
    if (is.character(capthist)) {
        capthistname <- capthist; rm(capthist)
        capthist <- get(capthistname, pos=-1)
    }
    if (is.character(mask)) {
        maskname <- mask; rm(mask)
        mask <- get(maskname, pos=-1)
    }
    if (is.character(dframe)) {
        dframename <- dframe; rm(dframe)
        dframe <- get(dframename, pos=-1)
    }
    if (is.character(details)) {
        detailsname <- details; rm(details)
        details <- get(detailsname, pos=-1)
    }
    
    if (!inherits(capthist, 'capthist'))
        stop ("requires 'capthist' object")

    ###################################################################    
    # optionally replace detector type - moved 2022-04-15
    if (!is.null(details$newdetector)) {
        warning("replacement detector type specified by user")
        if (ms(capthist)) {
            for (i in 1:length(capthist)) {
                if (inherits(details$newdetector, 'list'))
                    detector(traps(capthist[[i]])) <- details$newdetector[[i]]
                else
                    detector(traps(capthist[[i]])) <- details$newdetector
            }
        }
        else {
            detector(traps(capthist)) <- details$newdetector       
        }
    }
    ###################################################################    
    
    # restored 2020-08-30 for secrlinear arvicola example
    capthist <- check3D(capthist) 
    detectortype <- unlist(detector(traps(capthist)))
    anycount <- any(detectortype %in% .localstuff$countdetectors)
    anypoly  <- any(detectortype %in% c('polygon',  'polygonX'))
    anytrans <- any(detectortype %in% c('transect', 'transectX'))
    alltelem <- all(detectortype %in% 'telemetry')
    anytelem <- any(detectortype %in% 'telemetry')
    anysingle <- any(detectortype %in% 'single')
    anycapped <- any(detectortype %in% 'capped')
    allcapped <- all(detectortype %in% 'capped') 
    allpresence <- all(detectortype %in% 'presence')
    anysignal <- any(detectortype %in% 'signal')
    allsignal <- all(detectortype %in% 'signal')
    allsignalnoise <- all(detectortype %in% 'signalnoise')
    anysighting <- any(sighting(traps(capthist)))  # do not confuse 'anysightings' later
    anyindex <- any(detectortype %in% 'index')
    
    ##############################################
    # check valid analysis in secr
    if (anycapped) {
        if (CL) stop("CL not allowed with capped detectors")
        if (!is.null(groups)) stop("groups not implemented for capped detectors")
        if (!allcapped) stop("cannot combine capped and uncapped detector type")
    }
    details$anycapped <- anycapped
    
    #################################################
    ## Default detection function
    
    if (is.null(detectfn)) {
        if (allsignal) {
            detectfn <- 10
            warning ("detectfn not specified; using signal strength (10)")
        }
        else if (allsignalnoise) {
            detectfn <- 12
            warning ("detectfn not specified; using signal-noise (12)")
        }
        else if (anytelem) {
            detectfn <- 14
            warning ("detectfn not specified; using hazard half-normal (14)")
        }
        else if (anypoly | anytrans) {
            detectfn <- 14
        }
        else {
            detectfn <- 0
        }
    }
    else {
        if (anytelem)
            detectfn <- valid.detectfn(detectfn, c(14,16))
        else if (allpresence)
            detectfn <- valid.detectfn(detectfn, 0:8)
        else  if (anypoly | anytrans)
            detectfn <- valid.detectfn(detectfn, 14:20)  ## 2017-04-04, 2021-03-30
        else
            detectfn <- valid.detectfn(detectfn)
    }
    #################################################
    if (anysingle) warning ("multi-catch likelihood used for single-catch traps")
    if (anycapped) warning ("capped likelihood is an approximation")
    
    #################################################
    ## Use input 'details' to override various defaults
    
    defaultdetails <- list(
        distribution = 'poisson',
        hessian = 'auto',
        trace = TRUE,
        LLonly = FALSE,
        centred = FALSE,
        binomN = 0,                      ## Poisson
        cutval = 0,
        minprob = 1e-200,                ## before 2015-05-24 minprob = 1e-50
        tx = 'identity',
        param = 0,
        unmash = FALSE,
        telemetryscale = 1,
        ignoreusage = FALSE,
        debug = 0,
        intwidth2 = 0.8,
        usecov = NULL,
        userdist = NULL,
        autoini = 1,
        knownmarks = TRUE,
        nsim = 0,
        chatonly = FALSE,
        chat = NULL,
        savecall = TRUE,
        newdetector = NULL,
        contrasts = NULL,
        grain = 1,
        maxdistance = NULL,
        stackSize = "auto",   ## ignored on Windows
        fastproximity = TRUE,
        f = NULL              ## optional function f(x)
    )
    if (!is.null(attr(capthist,'cutval'))) {
        defaultdetails$cutval <- attr(capthist,'cutval')
    }
    else if (ms(capthist)) {
        if (!is.null(attr(capthist[[1]],'cutval')))   ## 2012-09-04
            defaultdetails$cutval <- attr(capthist[[1]],'cutval')
    }
    if (is.logical(details$hessian))
        details$hessian <- ifelse(details$hessian, 'auto', 'none')
    details <- replace (defaultdetails, names(details), details)
    if (!is.null(trace)) details$trace <- trace
    
    if (!is.null(binomN)) {
        if (tolower(binomN) == 'usage')
            binomN <- 1   ## code for 'binomial size from usage' 2012-12-22
        if (tolower(binomN) == 'poissonhazard')
            binomN <- -1  ## interpret g() as true detection function 2013-01-07
        details$binomN <- binomN   ## 2011 01 28
    }
    
    details$ncores <- setNumThreads(ncores, stackSize = details$stackSize)
    if (details$LLonly)  {
        details$trace <- FALSE
    }
    
    if (details$param == 1) {
        warning ("Gardner & Royle parameterisation discontinued in secr 3.0, using param = 0")
        details$param <- 0
    }
    
    if (details$savecall) {
        cl   <- match.call(expand.dots=TRUE)
        cl[[1]] <- quote(secr.fit)
    }
    else cl <- NULL
    
    # moved up from after LLonly 2022-03-06
    if (is.null(details$chat))
        details$chat <- matrix(1, nrow = length(session(capthist)), ncol = 3)
    else
        details$chat <- rep(details$chat,3)[1:3]  ## duplicate scalar
    
    #################################################
    ## MS - indicator TRUE if multi-session (logical)

    ## sessionlevels - names of sessions (character)
    MS <- ms(capthist)
    sessionlevels <- session(capthist)
    if (is.null(sessionlevels)) sessionlevels <- '1'
    if (alltelem) CL <- TRUE

    #################################################
    ## Optional data check
    if (verify) {
        memo ('Checking data', details$trace)
        test <- verify(capthist, report = 1)
        if (test$errors)
            stop ("'verify' found errors in 'capthist' argument")
        
        if (!is.null(mask)) {
            if (MS & ms(mask)) {
                ## list of masks
                test <- lapply(mask, verify, report = 1)
                notOK <- any(unlist(test))
            }
            else notOK <- verify(mask, report = 1)$errors
            if (notOK)
                stop ("'verify' found errors in 'mask' argument")
        }
    }
    
    #################################################
    ## Ensure valid mask
    ## assume traps(capthist) will extract a list of trap layouts
    ## if multi-session (MS == TRUE)
    
    usebuffer <- is.null(mask)    ## flag for later check
    if (usebuffer) {
        if (is.null(buffer)) {
            buffer <- 100
            if (!allpresence)
                warning ("using default buffer width 100 m")
        }
        makemaskCH <- function (CH, ...) {
            tr <- traps(CH)
            # specific to a session, don't use global anytelem
            if (any(detector(tr)=='telemetry')) {
                centroids <- t(sapply(telemetryxy(CH), apply, 2, mean))
                tmpxy <- rbind(centroids, data.frame(tr))
                make.mask(tmpxy, ...)
            }
            else make.mask(tr, ...)
        }
        # default is to buffer around both traps and centroids if telemetry
        if (MS) mask <- lapply (capthist, makemaskCH, buffer = buffer, type = "trapbuffer")
        else    mask <- makemaskCH(capthist, buffer = buffer, type = "trapbuffer")
    }
    else {
        if (MS & !ms(mask)) {
            if (inherits(mask, 'linearmask'))
                newclass <- c('linearmask', 'mask', 'list')
            else
                newclass <- c('mask', 'list')
            ## inefficiently replicate mask for each session!
            mask <- lapply(sessionlevels, function(x) mask)
            class (mask) <- newclass
            names(mask) <- sessionlevels
        }
    }
    #################################################
    
    nc <- ifelse (MS, sum(sapply(capthist, nrow)), nrow(capthist))
    
    if (nc < 1) {
        warning (nc, " detection histories")
    }
    if (is.null(details$userdist) & inherits(mask, 'linearmask')) {
        warning ("using Euclidean distances with linear mask")
    }
    
    #################################################
    ## mark-resight
    ## the checks here do not take account of details$markresight 2016-12-05
    if (MS) {
        sighting <- sighting(traps(capthist[[1]]))
        Tu <- Tu(capthist[[1]])
        Tm <- Tm(capthist[[1]])
    }
    else {
        sighting <- sighting(traps(capthist))
        Tu <- Tu(capthist)
        Tm <- Tm(capthist)
    }
    if (('pID' %in% names(fixed)) & !is.null(Tm)){
        ## if ((fixed$pID == 1) & (sum(Tm)>0) & any(markocc(traps(capthist))==0))
        ## 2019-12-16
        if ((fixed$pID == 1) & (sum(Tm)>0) & any(unlist(markocc(traps(capthist)))==0))
            warning ("mark-resight nonID sightings ignored when fixed$pID = 1")
    }
    if (sighting & CL & !is.null(Tu)) {
        warning ("mark-resight unmarked (but not nonID) sightings ignored when CL = TRUE")
    }
    
    #################################################
    ## optional centring of traps and mask 2010 04 27
    if (details$centred) {
        centre <- function (xy, dxy) {
            xy[,] <- sweep(xy, MARGIN = 2, FUN='-', STATS = dxy)
            xy
        }
        if (MS) {
            nsess <- length(traps(capthist))
            offsetxy <- lapply(traps(capthist), function(xy) apply(xy, 2, mean))
            for (i in 1:nsess) {
                temptraps <- centre(traps(capthist[[i]]), offsetxy[[i]])
                traps(capthist[[i]]) <- temptraps
                mask[[i]] <- centre(mask[[i]], offsetxy[[i]])
                attr(mask[[i]], 'meanSD')[1,1:2] <- attr(mask[[i]], 'meanSD')[1,1:2] -
                    offsetxy[[i]]
                attr(mask[[i]], 'boundingbox') <- centre(attr(mask[[i]], 'boundingbox'),
                                                         offsetxy[[i]])
            }
        }
        else {
            offsetxy <- apply(traps(capthist), 2, mean)
            traps(capthist) <- shift(traps(capthist), -offsetxy)
            mask <- shift.traps(mask, -offsetxy)
            attr(mask, 'meanSD')[1,1:2] <- attr(mask, 'meanSD')[1,1:2] - offsetxy
            attr(mask, 'boundingbox') <- centre(attr(mask, 'boundingbox'), offsetxy)
        }
    }
    
    #################################################
    ## standardize user model and parameterisation
    #################################################
    
    if ('formula' %in% class(model)) model <- list(model)
    model <- stdform (model)  ## named, no LHS
    if (CL) model$D <- NULL
    if (all(detectortype %in% c('telemetry'))) {
        model$g0 <- NULL
        model$lambda0 <- NULL
    }
    details$param <- new.param(details, model, CL)
    ## intercept and fix certain models with bad defaults
    model <- updatemodel(model, detectfn, 9, c('g0', 'sigma'), c('b0', 'b1'))
    model <- updatemodel(model, detectfn, 10:13, c('g0', 'sigma'), c('beta0','beta1'))
    model <- updatemodel(model, detectfn, 14:20, 'g0', 'lambda0')
    
    allvars <- unlist(lapply(model, all.vars))
    learnedresponse <- any(.localstuff$learnedresponses %in% allvars) ## || !is.null(dframe)
    timevarying <- any(c('t', 'T', 'tcov', names(timecov)) %in% allvars) ## || !is.null(dframe)
    
    ## 2022-01-05
    timevarying <- timevarying || any(names(timevaryingcov(traps(capthist))) %in% allvars)
    
    ## 2023-03-07 detect multi-session list of time covariate dataframes
    if (inherits(timecov, 'list')) {
        timevarying <- timevarying || any(names(timecov[[1]]) %in% allvars)
    }
    
    ##############################################
    ## 2019-09-01 fast proximity option
    details$fastproximity <- !is.null(details$fastproximity) &&
        details$fastproximity  && 
        ## all(detectortype %in% c('proximity', 'count', 'capped')) && 
        all(detectortype %in% c('proximity', 'count')) && 
        !learnedresponse && !timevarying && !anysighting &&
        is.null(groups)
    
    if (details$fastproximity) {
      if (any(unlist(detector(traps(capthist))) == 'count') && 
          (is.null(details$binomN) || any(details$binomN == 0)))
        count.distrib <- 'poisson'
      else
        count.distrib <- 'binomial'
      ## ensure single occasion 
      ## do it here so that design etc. use reduced capthist
      ## 2022-01-04 added verify = FALSE, dropunused = FALSE
      ## 2022-01-23 enforce all-ones usage if ignoreusage
      if (details$ignoreusage) {
        capthist <- uniformusage(capthist)
      }
      capthist <- reduce(capthist, by = 'all', outputdetector = 'count', 
        verify = FALSE, dropunused = FALSE)
      
      if (count.distrib == 'binomial') {
        if (!is.null(details$binomN) && (details$binomN[1]>1)) {
          ## use explicit binomN
          if (ms(capthist)) {
            for (r in 1:length(capthist)) {
              usage(traps(capthist[[r]])) <- usage(traps(capthist[[r]])) * details$binomN[1]
            }
          }
          else {
            usage(traps(capthist)) <- usage(traps(capthist)) * details$binomN[1] 
          }
        }
        details$ignoreusage <- FALSE  ## 2022-01-22
        details$binomN <- 1   ## binomial size from usage
      }
      loglikefn <- fastsecrloglikfn
    }
    else {
      loglikefn <- generalsecrloglikfn
    }
    
    #################################################
    ## which real parameters are fixed?
    #################################################
    
    ## c fixed by default in sigmak parameterisation
    if (details$param %in% 4:6) {
      if (! ("c" %in% names(model))) {
        ## default to fixed c = 0
        if (!("c" %in% names(fixed)))
          fixed$c <- 0
      }
      if (! ("d" %in% names(model))) {
        ## default to fixed d = 0
        if (!("d" %in% names(fixed)))
          fixed$d <- 0
      }
    }
    
    if (alltelem & !("lambda0" %in% names(fixed))) {
        ## default to fixed lambda0 = 1
        fixed$lambda0 <- 1.0
    }
    
    fnames <- names(fixed)
    
    #################################################
    ## build default model and update with user input
    #################################################
    
    defaultmodel <- list(D=~1, g0=~1, lambda0=~1,  esa=~1, a0=~1,
                         sigma=~1, sigmak=~1, z=~1, w=~1, c=~1, d=~1,
                         noneuc=~1, beta0=~1, beta1=~1,
                         sdS=~1, b0=~1, b1=~1, pID=~1, pmix=~1)
    defaultmodel <- replace (defaultmodel, names(model), model)
    
    #################################################
    # finite mixtures 
    #################################################
    
    nmix <- get.nmix(model, capthist, hcov)
    if (nmix > 3)
        stop ("number of latent classes exceeds 3")
    if ((nmix>1) & !is.null(hcov) & !is.null(groups))
        stop ("hcov mixture model incompatible with groups")
    if ((nmix == 1) & ('pmix' %in% c(fnames,names(model))))
        stop ("pmix specified for invariant detection model")
    
    if ((nmix>1) & !('pmix' %in% fnames)) {
        if (is.null(model$pmix)) model$pmix <- ~1
        pmixvars <- all.vars(model$pmix)
        if (!any (pmixvars %in% c('h2','h3'))) ## add mixing h2 or h3
        {
            defaultmodel$pmix <- if (nmix == 2)
                update(model$pmix, ~. + h2)
            else
                update(model$pmix, ~. + h3)
        }
        else {
            defaultmodel$pmix <- model$pmix   ## use as-is
            # 2021-06-17
            # badvar <- !(pmixvars %in% c('session','Session',sessioncov,'h2','h3'))
            badvar <- !(pmixvars %in% c('session','Session', names(sessioncov),'h2','h3'))
            if (any(badvar))
                stop ("formula for pmix may not include ", pmixvars[badvar])
        }
    }
    details$nmix <- nmix
    
    #################################################
    ## parameter names
    #################################################
    
    pnames <- valid.pnames (details, CL, detectfn, alltelem, sighting, nmix)
    
    #################################################
    ## test for irrelevant parameters in user's model
    #################################################
    
    OK <- names(model) %in% pnames
    if (any(!OK))
        stop ("parameters in model not consistent with detectfn etc. : ",
              paste(names(model)[!OK], collapse = ', '))
    OK <- fnames %in% pnames
    if (any(!OK))
        stop ("attempt to fix parameters not in model : ",
              paste(fnames[!OK], collapse = ', '))
    
    #################################################
    ## finalise model
    #################################################
    
    pnames <- pnames[!(pnames %in% fnames)]   ## drop fixed real parameters
    model <- defaultmodel[pnames]             ## select real parameters
    valid.model(model, CL, detectfn, hcov, details$userdist, names(sessioncov))
    vars <-  unlist(lapply(model, all.vars))
    
    #################################################
    ## Specialisations
    #################################################
    if (CL & !is.null(groups)) {
        groups <- NULL
        warning ("groups not valid with CL; groups ignored")
    }
    if (CL & var.in.model('g', model))
        stop ("'g' is not a valid effect when 'CL = TRUE'")
    if ((length(model) == 0) & !is.null(fixed))
        stop ("all parameters fixed")     ## assume want only LL
    
    ## mark-resight
    if ('pID' %in% names(model)) {
        pIDvars <- all.vars(model[['pID']])
        if (!all(pIDvars %in% c('session','Session',
                                names(sessioncov),
                                names(covariates(traps)))
        ))
            warning ("predictors in model for pID may be invalid")
    }
    #################################################
    # Link functions (model-specific)
    #################################################
    
    defaultlink <- list(D='log', g0='logit', lambda0='log', esa='log',
                        a0='log', sigma='log', sigmak='log', z='log',
                        w='log', c='identity', d='log', noneuc='log',
                        beta0='identity', beta1='neglog', sdS='log',
                        b0='log', b1='neglog',  pID='logit',
                        pmix='logit', cut='identity')
    
    # if (anycount) defaultlink$g0 <- 'log'
    link <- replace (defaultlink, names(link), link)
    link[!(names(link) %in% c(fnames,pnames))] <- NULL

    ##############################################
    # Prepare detection design matrices and lookup
    ##############################################
    memo ('Preparing detection design matrices', details$trace)
    design <- secr.design.MS (capthist, model, timecov, sessioncov, groups, hcov,
                              dframe, ignoreusage = details$ignoreusage, naive = FALSE,
                              CL = CL, contrasts = details$contrasts)
    design0 <- secr.design.MS (capthist, model, timecov, sessioncov, groups, hcov,
                               dframe, ignoreusage = details$ignoreusage, naive = TRUE,
                               CL = CL, contrasts = details$contrasts)
    ############################################
    # Prepare density design matrix
    ############################################
    D.modelled <- !CL & is.null(fixed$D)
    smoothsetup <- list(D = NULL, noneuc = NULL)
    if (!D.modelled) {
        designD <- matrix(nrow = 0, ncol = 0)
        grouplevels <- 1    ## was NULL
        attr(designD, 'dimD') <- NA
        nDensityParameters <- integer(0)
    }
    else {
        grouplevels  <- group.levels(capthist,groups)
        if (!is.null(details$userDfn)) {
            ## may provide a function used by getD in functions.R
            ## userDfn(mask, beta[parindx$D], ngrp, nsession)
            designD <- details$userDfn
            if (!is.function(designD))
                stop ("details$userDfn should be a function")
            ## this form of call returns only coefficient names
            Dnames <- designD('parameters', mask)
        }
        else {
            memo ('Preparing density design matrix', details$trace)
            ## 2021-06-17 tentative inclusion of session covariates 
            ## if (!all (all.vars(model$D) %in%
            ## c('session', 'Session','g')) & details$param %in% c(4,5)) {
            if (!all (all.vars(model$D) %in%
                    c('session', 'Session','g', names(sessioncov))) & details$param %in% c(4,5)) {
                if (is.null(details$userdist))
                    stop ("only session and group models allowed for density when details$param = ",
                        details$param)
            }
            temp <- D.designdata( mask, model$D, grouplevels, session(capthist), sessioncov)
            if (any(smooths(model$D))) {
                smoothsetup$D <- gamsetup(model$D, temp)
            }
            ## otherwise, smoothsetup$D remains NULL
            envD <- attr(model$D, '.Environment')
            if (!is.null(envD)) {
              assign('Dfn', identity, envir = envD)
            }
            designD <- general.model.matrix(model$D, data = temp, gamsmth = smoothsetup$D, 
                                            contrasts = details$contrasts)
            attr(designD, 'dimD') <- attr(temp, 'dimD')
            attr(designD, 'Dfn') <- details[['Dfn']]
            if(!is.null(details[['Dfn']])) {
                nDbeta <- details[['Dfn']](designD)
                Dnames <- paste0('D', 1:nDbeta)
            }
            else {
                Dnames <- colnames(designD)
            }
        }
        nDensityParameters <- length(Dnames)
    }
    ############################################
    # Prepare non-Euclidean design matrix
    ############################################
    NE.modelled <- ('noneuc' %in% getuserdistnames(details$userdist)) &
        is.null(fixed$noneuc)
    if (!NE.modelled) {
        designNE <- matrix(nrow = 0, ncol = 0)
        grouplevelsNE <- 1    ## was NULL
        attr(designNE, 'dimD') <- NA
        nNEParameters <- integer(0)
    }
    else {
        grouplevelsNE  <- group.levels(capthist,groups)
        memo ('Preparing non-Euclidean parameter design matrix', details$trace)
        temp <- D.designdata( mask, model$noneuc, grouplevelsNE, session(capthist), sessioncov)
        if (any(smooths(model$noneuc)))
            smoothsetup$noneuc <- gamsetup(model$noneuc, temp)
        ## otherwise, smoothsetup$NE remains NULL
        designNE <- general.model.matrix(model$noneuc, data = temp, gamsmth = smoothsetup$noneuc, 
                                         contrasts = details$contrasts)
        attr(designNE, 'dimD') <- attr(temp, 'dimD')
        NEnames <- colnames(designNE)
        nNEParameters <- length(NEnames)
    }

    ############################################
    # Parameter mapping (general)
    ############################################
    if (!is.null(design$designMatrices)) {
        np <- sapply(design$designMatrices, ncol)
    }
    else {
        np <- c(detectpar = 0)
    }
    np <- c(D = nDensityParameters, np, noneuc = nNEParameters)
    NP <- sum(np)
    parindx <- split(1:NP, rep(1:length(np), np))
    names(parindx) <- names(np)[np>0]
    if (!D.modelled) parindx$D <- NULL
    if (!NE.modelled) parindx$noneuc <- NULL
    
    data <- prepareSessionData(capthist, mask, details$maskusage, design, design0, detectfn, 
                               groups, fixed, hcov, details)

    ############################################
    # code for start vector shifted to separate function from 4.5.4
    ############################################
    start <- makeStart (start, parindx, capthist, mask, detectfn, link, 
        details, fixed, CL, anypoly, anytrans, alltelem, sighting) 
    if (is.null(start)) return(list(call = cl, fit = NULL))
    
    ############################################
    # Single evaluation option
    ############################################
    .localstuff$iter <- 0
    if (details$LLonly) {
        if (is.null(start) || any(is.na(start)))   ## is.list(start))
            stop ("must provide full vector of beta parameter values in 'start'")
        prep <- proc.time()
        LL <- - loglikefn (
            beta       = start,
            parindx    = parindx,
            link       = link,
            fixed      = fixed,
            designD    = designD,
            designNE   = designNE,
            design     = design,
            design0    = design0,
            detectfn   = detectfn,
            learnedresponse = learnedresponse,
            sessionlevels = sessionlevels,
            CL         = CL,
            data       = data,
            details    = details
        )
        out <- c(logLik=LL)
        attr(out, "npar") <- length(unlist(parindx))
        attr(out, "preptime") <- (prep-ptm)[3]
        attr(out, "LLtime") <- (proc.time() - prep)[3]
        return(out)
    }
    ############################################
    
    ############################################
    ## estimate overdispersion by simulation
    ############################################
    
    ## if (nsim > 0) and details$chatonly then exit, returning only chat
    if (details$nsim > 0) {
        TuTm <- function(x) !(is.null(Tu(x)) & is.null(Tm(x)))
        anysightings <- if (MS)
            any (sapply(capthist, TuTm))
        else TuTm(capthist)
        if (anysightings) {
            memo('Simulating sightings to estimate overdispersion...', details$trace)
            chat <- loglikefn (
                beta       = start,
                parindx    = parindx,
                link       = link,
                fixed      = fixed,
                designD    = designD,
                designNE   = designNE,
                design     = design,
                design0    = design0,
                CL         = CL,
                detectfn   = detectfn,
                learnedresponse = learnedresponse,
                sessionlevels = sessionlevels,
                data = data,
                details = details
            )
            if (details$chatonly)
                return(chat)   ## and no more!
            else {
                details$chat <- apply(chat,2,mean)  ## collapse multisession matrix
                details$nsim <- 0
                ## and proceed to call loglikefn again
            }
        }
    }
  
    ############################################
    ## ad hoc fix for experimental parameters
    ############################################
    if (is.null(details$miscparm))
        nmiscparm <- 0
    else {
        nmiscparm <- length(details$miscparm)
        if (length(start) < max(unlist(parindx)) + nmiscparm )
            start <- c(start, details$miscparm)
    }
    if (allsignalnoise)
        nmiscparm <- 2
    
    NP <- NP + nmiscparm
    stopifnot (length(start) == NP)
    
    ############################################
    # Fixed beta parameters
    ############################################
    fb <- details$fixedbeta
    if (!is.null(fb)) {
        if (!(length(fb)== NP))
            stop ("invalid fixed beta - require NP-vector")
        if (sum(is.na(fb))==0)
            stop ("cannot fix all beta parameters")
        ## drop unwanted betas; remember later to adjust parameter count
        start <- start[is.na(fb)]
        NP <- length(start)
    }
    ############################################
    # Variable names (general)
    ############################################
    betanames <- unlist(sapply(design$designMatrices, colnames))
    names(betanames) <- NULL
    realnames <- names(model)
    ## coefficients for D precede all others
    if (D.modelled) betanames <- c(paste('D', Dnames, sep='.'), betanames)
    ## coefficients for noneuc follow all others (except model-specific in para below)
    if (NE.modelled) betanames <- c(betanames, paste('noneuc', NEnames, sep='.'))
    betanames <- sub('..(Intercept))','',betanames)
    
    ############################################
    # Variable names (model-specific)
    ############################################
    
    if (detectfn %in% c(12,13)) {
        betanames <- c(betanames, 'muN', 'sdN')
        realnames <- c(realnames, 'muN', 'sdN')
    }
    else if (nmiscparm>0) {
        miscnames <- names(details$miscparm)
        if (is.null(miscnames))
            miscnames <- paste('miscparm', 1:nmiscparm, sep='')
        betanames <- c(betanames, miscnames)
    }
    
    ## retain betanames only for non-fixed beta (i.e. NA fixedbeta)
    if (!is.null(details$fixedbeta))
        betanames <- betanames[is.na(details$fixedbeta)]
    betaw <- max(max(nchar(betanames)),8)   # for 'trace' formatting
    names(start) <- betanames
    
    ############################################
    # Maximize likelihood
    ############################################
    
    lcmethod <- tolower(method)
    if (lcmethod != 'none') {
        memo('Maximizing likelihood...', details$trace)
        if (details$trace)
            cat('Eval     Loglik', formatC(betanames, format='f', width=betaw), '\n')
    }
    ## arguments always passed to loglikefn
    secrargs <- list(
        parindx    = parindx,
        link       = link,
        fixed      = fixed,
        designD    = designD,
        designNE   = designNE,
        design     = design,
        design0    = design0,
        CL         = CL,
        detectfn   = detectfn,
        learnedresponse = learnedresponse,
        sessionlevels = session(capthist),
        data       = data,
        details    = details,
        betaw      = betaw    # for trace format
    )   
    ############################################
    ## calls for specific maximisation methods
    ## 2013-04-21
    if (NP == 1) {
        lcmethod <- "optimise"
        signs <- c(-1,1) * sign(start)
        args <- list (f         = loglikefn,
                      interval  = start * (1 + details$intwidth2 * signs))
        args <- c(args, secrargs)
        args <- replace (args, names(list(...)), list(...))
        this.fit <- try(do.call (optimise, args))
        if (inherits(this.fit, 'try-error'))
            warning ("univariate search for minimum failed")
        this.fit$par <- this.fit$minimum
        this.fit$value <- this.fit$objective
        if (details$hessian != "none")
            details$hessian <- "fdHess"
    }
    else
        if (lcmethod %in% c('newton-raphson')) {
            args <- list (p         = start,
                          f         = loglikefn,
                          hessian   = tolower(details$hessian)=='auto',
                          stepmax   = 10)
            args <- c(args, secrargs)
            args <- replace (args, names(list(...)), list(...))
            
            this.fit <- do.call (nlm, args)

            this.fit$par <- this.fit$estimate     # copy for uniformity
            this.fit$value <- this.fit$minimum    # copy for uniformity
            if (this.fit$code > 2)
                warning ("possible maximization error: nlm returned code ",
                         this.fit$code, ". See ?nlm")
        }
    #-----------------------------------------------------------------
    else if (method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B",
                           "SANN", "Brent")) {
        args <- list(par     = start,
                     fn      = loglikefn,
                     hessian = tolower(details$hessian)=='auto',
                     method  = method)
        args <- c(args, secrargs)
        args <- replace (args, names(list(...)), list(...))
        this.fit <- do.call (optim, args)
        if (this.fit$convergence != 0)
            warning ("probable maximization error: optim returned convergence ",
                     this.fit$convergence, ". See ?optim")
    }
    #-----------------------------------------------------------------
    # Hessian-only 2013-02-23
    else if (lcmethod %in% 'none') {
        memo ('Computing Hessian with fdHess in nlme', details$trace)
        loglikfn <- function (beta) {
            do.call(loglikefn, c(list(beta=beta), secrargs))
        }
        grad.Hess <- nlme::fdHess(start, fun = loglikfn, .relStep = 0.001, minAbsPar=0.1)
        this.fit <- list (value = loglikfn(start), par = start,
                          gradient = grad.Hess$gradient,
                          hessian = grad.Hess$Hessian)
        biasLimit <- NA   ## no bias assessment
    }
    else stop ("maximization method", method, "not recognised")
    ############################################################################
    
    this.fit$method <- method         ## remember which method we used...
    covar <- NULL
    N <- NULL
    if (this.fit$value > 1e9) {     ## failed
        # this.fit$beta[] <- NA
        this.fit$par[] <- NA
    }
    else {
        
        ############################################
        ## Variance-covariance matrix
        ############################################
        
        if (tolower(details$hessian)=='fdhess') {
            memo ('Computing Hessian with fdHess in nlme', details$trace)
            loglikfn <- function (beta) {
                do.call(loglikefn, c(list(beta=beta), secrargs))
            }
            grad.Hess <- nlme::fdHess(this.fit$par, fun = loglikfn,
                                      .relStep = 0.001, minAbsPar = 0.1)
            this.fit$hessian <- grad.Hess$Hessian
        }
        if (!is.null(this.fit$hessian)) {
            covar <- try(MASS::ginv(this.fit$hessian))
            if (inherits(covar, "try-error")) {
                warning ("could not invert Hessian to compute ",
                         "variance-covariance matrix")
                covar <- matrix(nrow = NP, ncol = NP)
            }
            else if (any(diag(covar)<0)) {
                warning ("at least one variance calculation failed ")
            }
            dimnames(covar) <- list(betanames, betanames)
        }
    }
    
    ############################################
    ## form output list
    ############################################
    desc <- packageDescription("secr")  ## for version number

    ## if density modelled then smoothsetup already contains D,noneuc
    ## otherwise an empty list; now add detection parameters from
    ## secr.design.MS
    smoothsetup <- c(smoothsetup, design$smoothsetup)
    
    ## 2019-08-21
    ## experimentally do not save environment of model formulae
    
    model <- lapply(model, function(x) {environment(x) <- NULL; x})

    output <- list (call = cl,
                    capthist = capthist,
                    mask = mask,
                    detectfn = detectfn,
                    CL = CL,
                    timecov = timecov,
                    sessioncov = sessioncov,
                    hcov = hcov,
                    groups = groups,
                    dframe = dframe,
                    designD = designD,        # new 2019-10-12
                    designNE = designNE,      # new 2019-10-12
                    design = design,
                    design0 = design0,
                    start = start,
                    link = link,
                    fixed = fixed,
                    parindx = parindx,
                    model = model,
                    details = details,
                    vars = vars,
                    betanames = betanames,
                    realnames = realnames,
                    fit = this.fit,
                    beta.vcv = covar,
                    smoothsetup = smoothsetup,
                    learnedresponse = learnedresponse,
                    version = desc$Version,
                    starttime = starttime,
                    proctime = (proc.time() - ptm)[3]
    )
    
    class(output) <- 'secr'
    
    if (usebuffer & !is.na(biasLimit)) {
        test <- try(bufferbiascheck(output, buffer, biasLimit), silent = TRUE)
        if (inherits(test, 'try-error'))
            warning("test for mask truncation bias could not be performed")
    }
    
    memo(paste('Completed in ', round(output$proctime,2), ' seconds at ',
               format(Sys.time(), "%H:%M:%S %d %b %Y"),
               sep=''), details$trace)

    output
}
################################################################################

############################################################################################
## package 'secr'
## scoretest.R
## score test for secr models
############################################################################################

prepare <- function (secr, newmodel) {

    ## Prepare detection design matrices and lookup
    ## secr provides raw data etc.

    capthist <- secr$capthist
    mask     <- secr$mask
    timecov  <- secr$timecov
    sessioncov <- secr$sessioncov
    groups   <- secr$groups
    hcov     <- secr$hcov
    dframe   <- secr$dframe  ## never used?
    CL       <- secr$CL
    detectfn <- secr$detectfn
    link     <- secr$link
    details  <- secr$details
    fixed    <- secr$fixed
    smoothsetup <- secr$smoothsetup
    
    sessionlevels <- session(capthist)
    if (is.null(sessionlevels)) sessionlevels <- '1'
    grouplevels  <- group.levels(capthist, groups)

    design <- secr.design.MS (capthist, newmodel, timecov, sessioncov, groups, hcov, dframe,
                              ignoreusage = details$ignoreusage, CL = CL || details$relativeD, contrasts = details$contrasts)
    design0 <- secr.design.MS (capthist, newmodel, timecov, sessioncov, groups, hcov, dframe,
                               ignoreusage = details$ignoreusage, CL = CL || details$relativeD, contrasts = details$contrasts,
                               naive = T)
    D.modelled <- !CL & is.null(fixed$D)
    NE.modelled <- is.function(details$userdist) & is.null(fixed$noneuc)           
    sessionlevels <- session(capthist)
    grouplevels  <- group.levels(capthist, groups)
    D.designmatrix <- designmatrix (D.modelled, mask, newmodel$D,
                            grouplevels, sessionlevels, sessioncov, smoothsetup$D, details$contrasts)
    NE.designmatrix <- designmatrix (NE.modelled, mask, newmodel$noneuc,
                            grouplevels, sessionlevels, sessioncov, smoothsetup$noneuc, details$contrasts)        
            
    #############################
    # Parameter mapping (general)
    #############################

    np <- sapply(design$designMatrices, ncol)
    if (D.modelled) np <-  c(D = ncol(D.designmatrix), np)
    if (NE.modelled) np <-  c(np, noneuc = ncol(NE.designmatrix))
    
    NP <- sum(np)
    parindx <- split(1:NP, rep(1:length(np), np))
    names(parindx) <- names(np)[np>0]
    if (!D.modelled) parindx$D <- NULL
    if (!NE.modelled) parindx$noneuc <- NULL
    
    ## DOES THIS DEAL WITH SESSION COVARIATES OF DENSITY? 2009 08 16

    ################
    # Variable names
    ################

    betanames <- unlist(sapply(design$designMatrices, colnames))
    names(betanames) <- NULL
    if (D.modelled) betanames <- c(paste('D', colnames(D.designmatrix), sep='.'), betanames)
    if (NE.modelled) betanames <- c(betanames, paste('noneuc', colnames(NE.designmatrix), sep='.'))
    betanames <- sub('..(Intercept))','',betanames)
    # if (detector(traps(capthist)) %in% .localstuff$simpledetectors)
    # 2023-12-12
    if (all(unlist(detector(traps(capthist))) %in% .localstuff$simpledetectors))
        savedlogmultinomial <- logmultinom(capthist, group.factor(capthist, groups))
    else
        savedlogmultinomial <- 0

    ################
    list (
       link      = link,
       fixed     = fixed,       ## 2015-04-15
       parindx   = parindx,
       capthist  = capthist,
       mask      = mask,
       CL        = CL,
       detectfn  = detectfn,
       D.design  = D.designmatrix,
       NE.design = NE.designmatrix,
       design    = design,
       design0   = design0,
       hcov      = hcov,
       groups    = groups,
       betanames = betanames,
       details   = details,
       logmult   = savedlogmultinomial
     )

 }   # end of 'prepare'
############################################################################################

score.test <- function (secr, ..., betaindex = NULL, trace = FALSE, ncores = NULL,
                        .relStep = 0.001, minAbsPar = 0.1) {
    if (!inherits(secr, 'secr'))
        stop ("only for 'secr' objects")
    models <- list(...)
    setNumThreads(ncores)
    if (length(models)>1) { ##  ) {
        # apply to each component of 'model' list

        score.list <- lapply (models, score.test, secr=secr,
                            betaindex=betaindex)
        class(score.list) <- c('score.list', class(score.list))
        score.list
    }
    else {
        if (is.list(models) & !('formula' %in% class(models[[1]])) )
            model <- models[[1]]
        else
            model <- models

        if (inherits (model, 'secr')) model <- model$model
        # model may be incompletely specified
        model <- stdform (model)

        model <- replace (secr$model, names(model), model)

        ## OBSOLETE: 2014-10-25, BUT DOES IT NEED REPLACING?
        if (secr$CL) model$D <- NULL
        if (secr$detectfn!=1) model$z <- NULL

        #########################################################
        # check components on which model differs from secr$model
        
        if (length(secr$model) != length(model))
            stop ("models differ in number of real parameters")
        if (any(names(secr$model) != names(model)))
            stop ("real parameters do not match")
        differ <- logical(length(model))
        for (j in 1:length(model)) differ[j] <- secr$model[[j]] != model[[j]]
        if (sum(differ)==0)
            stop ("models are identical")

        #########################################################
        # construct essential parts of new secr model
        newsecr <- prepare (secr, model)
        newsecr$details <- replace (secr$details, 'trace', trace)  ## override
        newsecr$details$ncores <- setNumThreads(ncores)  ## 2021-10-17
        beta0 <- secr$fit$par
        names(beta0) <- secr$betanames
        beta1 <- mapbeta(secr$parindx, newsecr$parindx, beta0, betaindex)
        names(beta1) <- newsecr$betanames

        allvars <- unlist(lapply(model, all.vars))
        learnedresponse <- any(.localstuff$learnedresponses %in% allvars) 
        
        if (trace) {
            ## .localstuff is environment for temp variables in namespace
            .localstuff$iter <- 0
            cat('Eval     logLik', formatC(newsecr$betanames, format='f', width=11), '\n')
        }
        loglikfn <- function (beta, design) {
          # Return the negative log likelihood for spatial capture-recapture model
            -generalsecrloglikfn(
               beta     = beta,
               parindx  = design$parindx,
               link     = design$link,
               fixed    = design$fixed,
               designD  = design$D.design,
               designNE = design$NE.design, 
               design   = design$design,
               design0  = design$design0,
               CL       = design$CL,
               detectfn = design$detectfn,
               learnedresponse = learnedresponse,
               sessionlevels = session(design$capthist),
               data = data,
               details  = design$details,
               dig      = 4,
               betaw    = 11)
        }
        data <- prepareSessionData(newsecr$capthist, newsecr$mask, newsecr$details$maskusage, 
                                   newsecr$design, newsecr$design0, newsecr$detectfn, newsecr$groups, 
                                   newsecr$fixed, newsecr$hcov, newsecr$details)
        
        # cat('logLik',loglikfn(beta1, design=newsecr), '\n')   ## testing
        # cat('logmult', newsecr$logmult, '\n')
        
        ## fdHess is from package nlme
        grad.Hess <- nlme::fdHess(beta1, fun = loglikfn, design = newsecr,
                            .relStep = .relStep, minAbsPar = minAbsPar)
        u.star <- grad.Hess$gradient
        i.star <- -grad.Hess$Hessian

        ## score <- try (solve(i.star), silent = TRUE)
        ## use MASS function 2016-03-10
        score <- try (MASS::ginv(i.star), silent = TRUE)
        if (inherits(score, "try-error")) {
            warning ("could not invert information matrix")
            score <- NA
        }
        else  {
            score <- u.star %*% score %*% u.star
        }

        ## score statistic cf -2(logLik1 - logLik0)
        statistic <- c('X-squared' = as.numeric(score))

        ## number caught, for AICc
        n <- ifelse (is.list(secr$capthist),
            sum(sapply (secr$capthist, nrow)),
            nrow(secr$capthist) )

        np0 <- length(beta0)
        np1 <- length(beta1)
        np  <- c(np0, np1)
        parameter <- c(df = np1-np0)

        H0    <- model.string(secr$model, secr$details$userDfn)
        H1    <- model.string(model, secr$details$userDfn)
        call0 <-  paste('[', format(secr$call), ']', sep = '')
        AIC   <- c(0, -statistic + 2 * parameter)
        AICc  <- ifelse ((n-np-1)>0,
            AIC + 2 * np * (np+1) / (n-np-1),
            NA)
        names(AIC) <- c('AIC.0','AIC.1')
        names(AICc) <- c('AICc.0','AICc.1')

        temp <- list(
             statistic = statistic,
             parameter = parameter,
             p.value = 1 - pchisq(score, parameter),
             method = 'Score test for two SECR models',
             data.name = paste(H0, 'vs', H1, call0),
             np0 = np0,
             H0 = H0,
             H1 = H1,
             H1.beta = beta1,
             AIC = AIC,
             AICc = AICc
        )
        class(temp) <- c('score.test', 'htest')
        temp

    }
}
############################################################################################

score.table <- function (object, ..., sort = TRUE, dmax = 10) {

    if ('score.test' %in% class(object)) object <- list(object)
    if (length(match.call(expand.dots=F)$...)==0) score.list <- object
    else score.list <- c(object, list(...))
    if ( any( sapply(score.list, function(x) !('score.test' %in% class(x)) )) )
        stop ("arguments must be 'score.test' objects")

    output <- data.frame(model= sapply (score.list, function(x) x$H1))
    output$X.squared <- sapply(score.list, function(x) x$statistic)
    output$df <- sapply(score.list, function(x) x$parameter)
    output$p.value <- round(sapply(score.list, function(x) x$p.value),6)
    output$npar <- score.list[[1]]$np0 + sapply(score.list, function(x) x$parameter)
    output$AIC <- sapply(score.list, function(x) x$AIC[2])
    output$AICc <- sapply(score.list, function(x) x$AICc[2])

    output <- rbind(data.frame(model=score.list[[1]]$H0,
        X.squared = NA,
        df = NA,
        p.value = NA,
        npar = score.list[[1]]$np0,
        AIC = score.list[[1]]$AIC[1],
        AICc = score.list[[1]]$AICc[1]), output)

    # if (nrow(output)>1) {
    #     output$dAICc <- output$AICc - min(output$AICc)
    #     OK <- abs(output$dAICc) < abs(dmax)
    #     sumdAICc <- sum(exp(-output$dAICc[OK]/2))
    #     output$AICwt <- ifelse ( OK, exp(-output$dAICc/2) / sumdAICc, 0)
    # }
    if (nrow(output)>1) {
        output$dAIC <- output$AIC - min(output$AIC)
        OK <- abs(output$dAIC) < abs(dmax)
        sumdAIC <- sum(exp(-output$dAIC[OK]/2))
        output$AICwt <- ifelse ( OK, exp(-output$dAIC/2) / sumdAIC, 0)
    }
    
    output[,c('AIC','AICc','dAIC')] <- round(output[,c('AIC','AICc','dAIC')],3)
    if (!is.null(output$AICwt)) output[,'AICwt'] <- round(output[,'AICwt'],3)

    if (sort) output <- output [order(output$AICc),]
    else      output <- output [,]
    row.names(output) <- NULL

    output
}
############################################################################################

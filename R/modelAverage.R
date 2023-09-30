############################################################################################
## package 'secr'
## modelAverage.R
## replacing model.average
############################################################################################

MATA <- function (wt, est, se, alpha) {
    ## Turek and Fletcher 2012 "model-averaged tail area" interval
    fn <- function (theta, lcl = TRUE) {
        sum(wt * pnorm((est-theta) / se, lower.tail = lcl )) - alpha/2
    }
    lcl <- uniroot (fn, interval = c(min(est - 10*se), max(est)), lcl = FALSE)$root
    ucl <- uniroot (fn, interval = c(min(est), upper = max(est + 10*se)), lcl = TRUE)$root
    cbind(lcl, ucl)
}

modelAverage <- function (object, ..., 
    realnames = NULL, 
    betanames = NULL,
    newdata = NULL, 
    alpha = 0.05, 
    dmax = 10, 
    covar = FALSE, 
    average = c('link', 'real'), 
    criterion = c('AICc','AIC'), 
    CImethod = c('Wald', 'MATA'),
    chat = NULL) 
{
    UseMethod("modelAverage") 
} 

modelAverage.default <- function (object, ..., realnames = NULL, betanames = NULL,
    newdata = NULL, alpha = 0.05, dmax = 10, covar = FALSE, average =
        c('link', 'real'), criterion = c('AICc','AIC'), CImethod =
        c('Wald', 'MATA'), chat = NULL) 
{
    cat ('no modelAverage method for objects of class', class(object), '\n')
}

modelAverage.secr <- function (object, ..., realnames = NULL, betanames = NULL,
    newdata = NULL, alpha = 0.05, dmax = 10, covar = FALSE, average =
        c('link', 'real'), criterion = c('AICc','AIC'), CImethod =
        c('Wald', 'MATA'), chat = NULL) 
{
    allargs <- list(...)
    modelnames <- (c ( as.character(match.call(expand.dots=FALSE)$object),
        as.character(match.call(expand.dots=FALSE)$...) ))
    allargs <- secrlist(object, allargs)
    names(allargs) <- modelnames
    modelAverage(allargs,  realnames = realnames, betanames = betanames,
        newdata = newdata, alpha = alpha, dmax = dmax, covar = covar, 
        average = average, criterion = criterion, CImethod = CImethod,
        chat = NULL)
}

 modelAverage.secrlist <- function (object, ..., realnames = NULL, betanames = NULL,
    newdata = NULL, alpha = 0.05, dmax = 10, covar = FALSE, average =
        c('link', 'real'), criterion = c('AICc','AIC'), CImethod =
        c('Wald', 'MATA'), chat = NULL) 
{
     if (length(list(...)) > 0) {
         warning ("... argument ignored in 'modelAverage.secrlist'")
     }
     
    #########
    ## SETUP

    ## match character arguments
    CImethod <- match.arg(CImethod)
    criterion <- match.arg(criterion)
    quasi <- !is.null(chat)
    Q <- if (quasi) "Q" else ""
    average <- match.arg(average)

    if (!is.null(names(object))) {
        modelnames <- names(object)
    }
    else {
        modelnames <- as.character(match.call(expand.dots=FALSE)$...)
    }
    
    ## checks
    if ( any (!sapply(object, function (x) inherits(x, 'secr'))) )
        stop ("require fitted 'secr' objects")
    if ( length(object) < 2 )
        warning ("only one model")
    if (!is.list(object) | !inherits(object[[1]], 'secr'))
        stop ("object must be secr or list of secr")

    ## preliminaries
    type <- 'real'                     ## default
    parnames <- object[[1]]$realnames  ## default
    links <- object[[1]]$link
    if (!is.null(realnames))
        parnames <- realnames
    else if (!is.null(betanames)) {
        type <- 'beta'
        average <- 'beta'   ## override
        parnames <- betanames
    }
    np <- length(parnames)
    nsecr <- length(object)
    
    ## rudimentary checks for compatible models
    if (nsecr > 1) {
        objnames <- function(i) switch (type, real=object[[i]]$realnames, beta=object[[i]]$betanames)
        test <- sapply (2:nsecr, function(i) sum(match(parnames, objnames(i), nomatch=0)>0) == np)
        if (!all(test))
            stop ("requested parameters not found in all models, ",
                  "or models incompatible")
    }

    ## NEWDATA construct if not provided
    if (is.null(newdata)) {
        #############################################################
        ## form unified 'newdata' containing all necessary predictors
        # tempnewdata <- lapply (object, secr.make.newdata)
        tempnewdata <- lapply (object, makeNewData)
        ## extract list of columns from all components of 'newdata'
        ## modified 2010 02 14
        column.list <- list(0)
        for (i in 1:nsecr) column.list[[i]] <- as.list(tempnewdata[[i]])
        column.list <- unlist(column.list, recursive = FALSE)

        column.list <- column.list[unique(names(column.list))]
        column.list <- lapply(column.list, unique)
        common <- names(column.list)[names(column.list) %in% names(newdata)]
        column.list[common] <- newdata[common]   ## user input

        ## modified 2010 02 23 to match levels of session covariates to session
        ## not tested with different session covariates in diff sessions
        sessioncovs <- lapply(object, function(x)
            if(!is.null(x$sessioncov)) data.frame(session=session(x$capthist), x$sessioncov)
            else NULL)
        sessioncovs <- sessioncovs[!sapply(sessioncovs, is.null)]
        scn <- as.vector(sapply(sessioncovs, names))
        scn <- match(unique(scn),scn)
        sessioncovs <- as.data.frame(sessioncovs)[,scn]
        sessioncovs <- as.data.frame(sessioncovs[!sapply(sessioncovs, is.null)])
        sessioncovnames <- unlist(lapply(object, function(x) names(x$sessioncov)))
        sessioncovariate <- names(column.list) %in% sessioncovnames
        newdata <- expand.grid (column.list[!sessioncovariate])
        ## replaced 2010 03 09 to retain sort order
        #   if (nrow(sessioncovs)>0) newdata <- merge(newdata, sessioncovs)
        if (nrow(sessioncovs)>0) {
            for (i in names(sessioncovs)) {
                if (i != 'session') newdata[,i] <- sessioncovs[newdata$session,i]
            }
        }
    }
    
    nr <- nrow(newdata)
    rownames <- apply(newdata, 1, function(x) paste(names(newdata), '=', x, sep='', collapse=','))

    ## VECTOR OF MODEL WEIGHTS

    AICdf <- AIC(object, sort = FALSE, dmax = dmax, criterion = criterion, chat = chat)
    wt <- unlist(AICdf[,paste0(Q, criterion, 'wt')])

    ## AVERAGE MODELS

    if (type == 'beta') {
        nr   <- 1
        ests <- lapply (object, function(x) coef(x)[parnames, 'beta'])
        ests <- array(unlist(ests), dim=c(nr, np, nsecr))
        estse <- lapply (object, function(x) coef(x)[parnames, 'SE.beta'])   ## for MATA
        estse <- array(unlist(estse), dim=c(nr, np, nsecr))
        ests <- array(unlist(ests), dim=c(nr, np, nsecr))
        estse <- array(unlist(estse), dim=c(nr, np, nsecr))
        wtd <- sweep(ests, MARGIN = 3, STATS = wt, FUN = '*')
        wtd <- apply(wtd, 1:2, sum)
        ## fails if no dimnames on x$beta.vcv:
        vcvs <- lapply (object, function(x) { vcov(x)[parnames, parnames] })
        vcv1 <- array(unlist(vcvs), dim=c(np, np, nsecr))  ## between beta parameters
        vcv <- sweep (vcv1, MARGIN = 3, STATS = wt, FUN = '*')
        vcv <- apply(vcv, 1:2, sum)
        sewtd  <- diag(vcv)^0.5
    }
    else { ## type == 'real'
        getLP <- function (object1) {  ## predicted values of real parameters
            getfield <- function (x) {
                secr.lpredictor (
                    formula = object1$model[[x]], 
                    newdata = newdata,
                    indx = object1$parindx[[x]], 
                    beta = object1$fit$par,
                    beta.vcv = object1$beta.vcv, field = x,
                    smoothsetup = object1$smoothsetup[[x]], 
                    contrasts = object1$details$contrasts,
                    Dfn = attr(object1$designD, 'Dfn')
                )
            }
            ## check added 2016-06-16, fixed 2016-11-08
            if (any(unlist(nclusters(object1$capthist))>1))
                warning("model.average is ignoring n.mashed", call. = FALSE)
            sapply (names(object1$model), getfield, simplify = FALSE)
        }
        predicted <- lapply (object, getLP)

        ests <- lapply (predicted, function(x) lapply(x[parnames], function(y) y[, 'estimate'] ))
        estse <- lapply (predicted, function(x) lapply(x[parnames], function(y) y[, 'se'] ))
        if (average == 'real') {
            ests <- lapply (ests, function (x) Xuntransform(unlist(x),
                varnames=rep(parnames, rep(nr, np)), linkfn=links))
            estse <- mapply (function (x,se.x) se.Xuntransform(unlist(x), unlist(se.x),
                varnames=rep(parnames, rep(nr, np)), linkfn=links), ests, estse)
            vcvs <- lapply (object, vcov, realnames = parnames, newdata=newdata)
        }
        else {
            vcvs <- lapply (predicted, function(x) lapply(x[parnames], function(y) attr(y, 'vcv')))
        }
        ests <- array(unlist(ests), dim=c(nr, np, nsecr))
        estse <- array(unlist(estse), dim=c(nr, np, nsecr))
        wtd <- sweep(ests, MARGIN = 3, STATS = wt, FUN = '*')
        wtd <- apply(wtd, 1:2, sum)
        ## variance from Burnham & Anderson 2004 eq (4)
        vcv1 <- array(unlist(vcvs), dim=c(nr, nr, np, nsecr))  ## between rows of newdata
        betak <- sweep(ests, MARGIN = 1:2, STATS = wtd, FUN='-')
        vcv2 <- array(apply(betak, 2:3, function(x) outer(x,x)), dim=c(nr,nr,np,nsecr))
        vcv <- sweep (vcv1 + vcv2, MARGIN = 4, STATS = wt, FUN = '*')
        vcv <- apply(vcv, 1:3, sum)
        sewtd  <- apply(vcv, 3, function (x) diag(x)^0.5)

    }

    ## OUTPUT
    sewtd <- array(sewtd, dim = c(nr,np))
    z <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard z!
    output <- array (dim=c(nr, 4, np))
    if (type=='beta') {
        # extra dimension just to share later code
        output[1,1,] <- wtd
        output[1,2,] <- sewtd
        output[1,3,] <- wtd - z*sewtd
        output[1,4,] <- wtd + z*sewtd
        ## overwrite confidence intervals if MATA requested
        if (CImethod == 'MATA') {
            for (p in 1:np)
                output[1,3:4,p]  <- MATA(wt, ests[1,p,], estse[1,p,], alpha)
        }
        dimnames(output) <- list(NULL,
            c('beta', 'SE.beta', 'lcl', 'ucl'),
            parnames)
    }
    else { ## type=='real'
        ## 2010 02 15, 2010 02 23, 2015-09-30
        for (m in 1: length(object))
        if (object[[m]]$details$param > 1) {
            stop ("'model.average' cannot handle nonstandard parameterisations")
        }
        for (i in 1:nr) {
            if (average == 'real') {
                output[i,1,] <- wtd[i,]
                output[i,2,] <- sewtd[i,]
                lpwtd <- Xtransform (wtd[i,], links, parnames)
                selpwtd <- se.Xtransform (wtd[i,], sewtd[i,], links, parnames)
                output[i,3,] <- Xuntransform(lpwtd-z*selpwtd, links, parnames)
                output[i,4,] <- Xuntransform(lpwtd+z*selpwtd, links, parnames)
            }
            else {
                output[i,1,] <- Xuntransform(wtd[i,], links, parnames)
                output[i,2,] <- se.Xuntransform(wtd[i,], sewtd[i,], links, parnames)
                output[i,3,] <- Xuntransform(wtd[i,]-z*sewtd[i,], links, parnames)
                output[i,4,] <- Xuntransform(wtd[i,]+z*sewtd[i,], links, parnames)
            }
            ## overwrite confidence intervals if MATA requested
            if (CImethod == 'MATA') {
                for (p in 1:np)
                    output[i,3:4,p] <- untransform(MATA(wt,ests[i,p,], estse[i,p,],alpha), links[[p]])
            }
        }
        dimnames(output) <- list(rownames,
            c('estimate', 'SE.estimate', 'lcl', 'ucl'),
            parnames)
    }

    ## collapse if possible
    if (dim(output)[1] == 1) {
        if (np==1) output <- output[1,,]
        else output <- t(output[1,,])
    }
    else if (np==1) output <- output[,,1]

    if (covar) {
        ## 2011-11-05; 2012-02-16 bug fix
        if (type == 'real') {
            stop("sorry - covar not available at present",
                 " for model-averaged real parameters")
            dimnames(vcv) <- list (rownames, rownames, parnames)
            output <- list (ma = output, rowvcv = vcv)
        }
        else {
            dimnames(vcv) <- list (parnames, parnames)
            output <- list (ma = output, linkvcv = vcv)
        }
    }
    return(output)
}


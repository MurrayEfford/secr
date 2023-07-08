############################################################################################
## package 'secr'
## collate.R
## Moved to own file 2023-07-08
## 2023-07-08 bug fixed: betanames failed with multiple sessions
############################################################################################

collate <- function (..., realnames = NULL, betanames = NULL, newdata = NULL,
                     alpha = 0.05, perm = 1:4, fields = 1:4) {
    ## 2013-04-20
    object <- secrlist(...)
    if (!is.null(names(object)))
        modelnames <- names(object)
    else
        modelnames <- as.character(match.call(expand.dots=FALSE)$...)
    
    if ( any (!sapply(object, function (x) inherits(x, c('secr','ipsecr'))) ))
        stop ("require fitted secr or ipsecr objects")
    if ( length(object) < 2 )
        warning ("only one model")
    if (!is.list(object) | !inherits(object[[1]], c('secr','ipsecr')))
        stop("object must be secr or list of secr")
    
    type <- 'real'                     ## default
    
    parnames <- unique(as.vector(unlist(sapply(object,
                                               function(x) x$realnames))))  ## default
    
    if (!is.null(realnames))
        parnames <- realnames
    else if (!is.null(betanames)) {
        type <- 'beta'
        parnames <- betanames
    }
    
    np <- length(parnames)
    nsecr <- length(object)
    
    ## rudimentary checks for compatible models
    if (nsecr > 1) {
        objnames <- function(i) switch (type,
                                        real = object[[i]]$realnames, beta = object[[i]]$betanames)
        test <- sapply (2:nsecr, function(i)
            sum(match(parnames, objnames(i), nomatch=0)>0) == np)
        if (!all(test))
            stop ("parameters not found in all models, or incompatible models")
    }
    
    getLP <- function (object1) {  ## for predicted values of real parameters
        getfield <- function (x) {
            if (is.null(object1$beta)) object1$beta <-  object1$fit$par         
            secr.lpredictor (
                formula = object1$model[[x]], 
                newdata = newdata,
                indx = object1$parindx[[x]], 
                beta = object1$beta,
                beta.vcv = object1$beta.vcv, 
                field = x,
                smoothsetup = object1$smoothsetup[[x]],
                contrasts = object1$details$contrasts,
                f = object1$details[['f']]
            )
        }
        ## check added 2016-06-16, fixed 2016-11-08
        if (any(unlist(nclusters(object1$capthist))>1))
            warning("collate is ignoring n.mashed", call. = FALSE)
        sapply (names(object1$model), getfield, simplify = FALSE)
    }
    
    if (is.null(newdata)) {
        
        ## form unified 'newdata' containing all necessary predictors
        
        ## start with list of model-specific newdata --
        ## each component of tempnewdata is a data.frame of newdata
        ## for the corresponding model
        # tempnewdata <- lapply (object, secr.make.newdata)
        tempnewdata <- lapply (object, makeNewData)
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
    
    z <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard z!
    if (type == 'real') {
        nr <- nrow(newdata)
        rownames <- apply(newdata, 1, function(x) paste(names(newdata), '=', x, sep='', collapse=','))
        predict <- lapply (object, getLP)
        stripped <- lapply(predict, function(x) lapply(x[parnames], function(y) y[, c('estimate','se')] ))
        stripped <- array(unlist(stripped), dim=c(nr, 2, np, nsecr))
    }
    else {
        nr <- 1
        coefs <- lapply (object, coef)
        stripped <- lapply(coefs, function(x) x[parnames, c('beta','SE.beta')] )
        stripped <- array(unlist(stripped), dim=c(nr, np, 2, nsecr))
        stripped <- aperm(stripped, c(1,3,2,4))
    }
    output <- array (dim=c(nr, 4, np, nsecr))
    if (type=='real') {
        output[,1:2,,] <- stripped
        for (i in 1:nr) {
            for (m in 1:nsecr) {
                ## changed object[[1]]$link to object[[m]]$link following lines 2010 02 14 - seems right
                output[i,1,,m] <- Xuntransform(stripped[i,1,,m, drop=FALSE], object[[m]]$link, parnames)
                output[i,2,,m] <- se.Xuntransform(stripped[i,1,,m, drop=FALSE], stripped[i,2,,m, drop=FALSE], object[[m]]$link, parnames)
                output[i,3,,m] <- Xuntransform(stripped[i,1,,m, drop=FALSE]-z*stripped[i,2,,m, drop=FALSE], object[[m]]$link, parnames)
                output[i,4,,m] <- Xuntransform(stripped[i,1,,m, drop=FALSE]+z*stripped[i,2,,m, drop=FALSE], object[[m]]$link, parnames)
            }
        }
    }
    else { ## type=='beta'
        output[1,1:2,,] <- stripped
        output[1,3,,] <- stripped[1,1,,,drop=FALSE]-z*stripped[1,2,,,drop=FALSE]
        output[1,4,,] <- stripped[1,1,,,drop=FALSE]+z*stripped[1,2,,,drop=FALSE]
    }
    
    if (type=='real') {
        dimnames(output) <- list(rownames,
                                 c('estimate', 'SE.estimate', 'lcl', 'ucl'),
                                 parnames,
                                 modelnames)
    }
    else {
        dimnames(output) <- list(NULL,
                                 c('beta', 'SE.beta', 'lcl', 'ucl'),
                                 parnames,
                                 modelnames)
    }
    ## default dimensions:
    ## row, model, statistic, parameter
    output <- aperm(output, c(1,4,2,3))
    
    return(aperm(output[,,fields,,drop=FALSE], perm))
    
}
############################################################################################


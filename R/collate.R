############################################################################################
## package 'secr'
## collate.R
## Moved to own file 2023-07-08
## 2023-07-08 bug fixed: betanames failed with multiple sessions
## 2025-07-17 bug fixed: allows for fixedbeta via complete.beta()
## 2025-07-20 improved check for compatibility
############################################################################################

collate <- function (object, ..., 
                     realnames = NULL, betanames = NULL, newdata = NULL,
                     alpha = 0.05, perm = 1:4, fields = 1:4) 
{
    UseMethod("collate") 
} 

collate.default <- function (object, ..., 
                             realnames = NULL, betanames = NULL, newdata = NULL,
                             alpha = 0.05, perm = 1:4, fields = 1:4) 
{
    cat ('no collate method for objects of class', class(object), '\n')
} 

collate.secr <- function (object, ..., realnames = NULL, betanames = NULL, newdata = NULL,
                            alpha = 0.05, perm = 1:4, fields = 1:4) {
    allargs <- list(...)
    modelnames <- (c ( as.character(match.call(expand.dots=FALSE)$object),
                       as.character(match.call(expand.dots=FALSE)$...) ))
    allargs <- secrlist(object, allargs)
    names(allargs) <- modelnames
    collate(allargs,  realnames = realnames, betanames = betanames,
            newdata = newdata, alpha = alpha, perm = perm, fields = fields)
}

collate.ipsecr <- function (object, ..., 
                            realnames = NULL, betanames = NULL, newdata = NULL,
                            alpha = 0.05, perm = 1:4, fields = 1:4) 
{
    allargs <- list(...)
    modelnames <- (c ( as.character(match.call(expand.dots=FALSE)$object),
                       as.character(match.call(expand.dots=FALSE)$...) ))
    class(object) <- 'secr'
    allargs <- secrlist(object, allargs)
    names(allargs) <- modelnames
    collate(allargs, realnames = realnames, betanames = betanames, 
            newdata = newdata, alpha = alpha, perm = perm, fields = fields)
}

collate.secrlist <- function (object, ..., realnames = NULL, betanames = NULL, newdata = NULL,
                                alpha = 0.05, perm = 1:4, fields = 1:4) {
    if (length(list(...)) > 0) {
        warning ("... argument ignored in 'collate.secrlist'")
    }
    if (!is.null(names(object)))
        modelnames <- names(object)
    else
        modelnames <- as.character(match.call(expand.dots=FALSE)$...)
    
    if ( length(object) < 2 )
        warning ("only one model")

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
                                        real = object[[i]]$realnames, 
                                        beta = object[[i]]$betanames)
        test <- sapply (1:nsecr, function(i)
            sum(match(parnames, objnames(i), nomatch=0)>0) == np)
        if (!all(test))
            stop (call. = FALSE, 
                  "Parameter(s) not found in all models, or incompatible models; \n",
                  "specify common parameter(s) with realnames or betanames?")
    }
    
    getLP <- function (object1) {  ## for predicted values of real parameters
        getfield <- function (x) {
            secr_lpredictor (
                formula = object1$model[[x]], 
                newdata = newdata,
                indx = object1$parindx[[x]], 
                beta = beta, 
                beta.vcv = beta.vcv, 
                field = x,
                smoothsetup = object1$smoothsetup[[x]],
                contrasts = object1$details$contrasts,
                Dfn = attr(object1$designD, 'Dfn')
            )
        }
        if (any(unlist(secr_nclusters(object1$capthist))>1))
            warning("collate is ignoring n.mashed", call. = FALSE)
        # 2025-07-17 allow for fixedbeta
        beta <- secr_complete.beta(object1)
        beta.vcv <- secr_complete.beta.vcv(object1)
        sapply (names(object1$model), getfield, simplify = FALSE)
    }
    
    if (is.null(newdata)) {
        
        ## form unified 'newdata' containing all necessary predictors
        
        ## start with list of model-specific newdata --
        ## each component of tempnewdata is a data.frame of newdata
        ## for the corresponding model
        tempnewdata <- lapply (object, makeNewData)
        column.list <- list(0)
        for (i in 1:nsecr) column.list[[i]] <- as.list(tempnewdata[[i]])
        column.list <- unlist(column.list, recursive = FALSE)
        
        column.list <- column.list[unique(names(column.list))]
        column.list <- lapply(column.list, unique)
        common <- names(column.list)[names(column.list) %in% names(newdata)]
        column.list[common] <- newdata[common]   ## user input
        
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


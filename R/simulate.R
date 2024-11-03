################################################################################
## package 'secr'
## simulate.R

## 2024-07-29 simulate.R split from sim.secr.R
################################################################################

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

################################################################################

sim.onepopn <- function (object, Darray, chat = 1) {
    ngrp <- dim(Darray)[2]
    nsession <- dim(Darray)[3]
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
                # 2024-10-30 suppress to enable masking
                # if (object$model$D == ~1) {
                #     density <- Darray[1,g,sessnum]   ## scalar
                #     mod2D <- 'poisson'               ## homogeneous
                # }
                # else {
                    density <- Darray[,g,sessnum]    ## vector
                    mod2D <- 'IHP'                   ## inhomogeneous
                # }
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
            if (is.null(object$groups) && is.null(object$hcov))
                covariates(popn[[g]]) <- NULL
            else
                covariates(popn[[g]]) <- data.frame(row.names = row.names(popn[[g]]))
            
            ## ---groups---
            if (!is.null(object$groups)) {
                grpcov <- as.data.frame(object$di[rep(g, nrow(popn[[g]])),]) ## 2014-08-08
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
    sesspopn   # list
}

simulate.secr <- function (object, nsim = 1, seed = NULL, maxperpoly = 100, chat = 1,
                           poponly = FALSE, ...)
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
    if (!is.null(object$groups)) {
        ## individual covariates for foundation of g
        ## pass with object 2024-09-09
        object$di <- disinteraction (object$capthist, object$groups)
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
    runone <- function() {
        sesspopn <- sim.onepopn(object, Darray, chat)
        if (poponly) 
            sesspopn
        else
            sim.detect(object, sesspopn, maxperpoly)
    }
    sesscapt <- replicate(nsim, runone(), simplify = FALSE)
    if (poponly) {
        attr(sesscapt,'seed') <- RNGstate   ## save random seed
        class(sesscapt) <- c('popn', 'list')
        sesscapt
    }
    else {
        ## experimental
        if (chat>1) sesscapt <- lapply(sesscapt, replicate, chat)
        
        attr(sesscapt,'seed') <- RNGstate   ## save random seed
        class(sesscapt) <- c('secrdata', 'list')
        sesscapt
    }
}
################################################################################

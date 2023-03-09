###############################################################################
## package 'secr' 4.5
## secr.design.MS.R

## 2019-12-03 replaced bygroup with CL
## 2023-03-10 individualcovariates() moved from utility.R

################################################################################

## Based on Charles C. Berry on R-help 2008-01-13
## drop = FALSE 2018-11-22
n.unique.rows <- function(x) {
    order.x <- do.call(order, as.data.frame(x))
    equal.to.previous <- rowSums(x[tail(order.x,-1),,drop = FALSE] != 
            x[head(order.x,-1),,drop = FALSE])==0 
    1 + sum(!equal.to.previous)
}

individualcovariates <- function (PIA) {
    pia <- matrix(aperm(PIA, c(2:5,1)), nrow = dim(PIA)[2])
    n.unique.rows(pia) > 1
}
#-------------------------------------------------------------------------------

secr.design.MS <- function (capthist, models, timecov = NULL, sessioncov = NULL,
                            groups = NULL, hcov = NULL, dframe = NULL, naive = FALSE,
                            CL = FALSE, keep.dframe = FALSE, full.dframe = FALSE,
                            ignoreusage = FALSE, contrasts = NULL, ...) {

## Generate design matrix, reduced parameter array, and parameter index array (PIA)
## for detection function parameters
## 'capthist' must be of class 'capthist' or 'list'
## uses pad1 to pad session-specific covar to constant length with first value,
## pad1 defined in 'utility.R'
##
## groups is a vector of factor names whose intersection defines group
## groups only defined for CL = FALSE  
## use of 'g' requires valid groups definition
## grouping variables are also added individually to dframe

    #--------------------------------------------------------------------------------
    findvars.MS <- function (cov, vars, dimcov) {
        ## function to add covariates to a design data frame 'dframe'
        ## cov may be a dataframe or list of dataframes, one per session (R > 1),
        ## if list, then require predictors to appear in all sessions
        ## uses insertdim from utility.R
        ## NOT to be used to add group variables
        ## Does NOT standardize numeric covariates:

        if (is.null(cov) | (length(cov)==0) | (length(vars)==0)) return()
        else {
            found <- ''
            if (!is.data.frame(cov)) {
                ## therefore assume cov is a list, one per session
                ## care required as sessions may differ in n, S, K
                ## solution is to pad to max [n,S,K] over sessions
                if (!is.list(cov) | (R==1))
                    stop ("irregular covariates; check multisession structure")
                
                #############################################
                ## 2020-05-15 convert all character to factor
                cov <- lapply(cov, stringsAsFactors)
                #############################################
                
                covnames <- lapply(cov, names)
                varincov <- sapply(covnames, function(nam) vars %in% nam)
                if (length(vars)>1) found <- vars[apply(varincov,1,all)]
                else found <- vars[all(varincov)]
                for (variable in found) {
                    ## if factor, check all levels the same 2013-03-11
                    ## 2020-05-15 NOTE no check for character covariates
                    if (any(sapply(cov, function(x) is.factor(x[,variable])))) {
                        baselevels <- levels(cov[[1]][,variable])
                        lev <- lapply(cov, function(x) levels(x[,variable]))
                        if (!all(sapply(lev[-1], function(x) identical(x,baselevels)))) {
                            print(lev)
                            stop("covariate factor levels differ between sessions")
                        }
                    }
                    ## pad on first dimension
                    vals <- lapply(cov, function(x) pad1(x[,variable], dims[dimcov[1]]))
                    vals <- unlist(vals)
                    dframe[,variable] <<- insertdim (vals, dimcov, dims)
                }
            }
            else
            {
                #############################################
                ## 2020-05-15 convert all character to factor
                cov <- stringsAsFactors(cov)
                #############################################
                
                found <- names(cov) %in% vars
                if (is.data.frame(cov) & any(found)) {
                    found <- names(cov)[found]
                    values <- as.data.frame(cov[,found])
                    names(values) <- found
                    if (length(values)>0) {
                        for (variable in found) {
                            vals <- values[,variable]
                            dframe[,variable] <<- insertdim (vals, dimcov, dims)
                        }
                    }
                }
            }
            vars <<- vars[!(vars %in% found)]
        }
    }
    #--------------------------------------------------------------------------------

    findvars.traptime <- function (covindices, vars) {
        ## function to add time-specific trap covariates to a design data frame 'dframe'
        ## covindices should be a list or list of lists, one per session (R > 1),
        ## if list, then require predictors to appear in all sessions
        ## uses pad1 and insertdim from utility.R

        found <- ''
        dimcov <- c(1,4,3) ## session, trap, time
        ## covindices is list of numeric or character index vectors, one component per session
        if (length(covindices) != R)
            stop ("require one set of indices per session")
        if (is.data.frame(trapcov))   ## single-session
            trapcov <- list(trapcov)
        covnames <- unique(sapply(covindices,names))
        found <- vars[vars %in% covnames]
        vars <<- vars[!(vars %in% found)]

        for (variable in found) {
            firstcol <- trapcov[[1]][,covindices[[1]][[1]][1]]
            factorlevels <- NULL
            if (is.factor(firstcol)) {
                ## all must have same levels!!
                factorlevels <- levels(firstcol)
            }
            vals <- vector(mode = 'list', length = R)
            for (i in 1:R) {  ## over sessions
                getvals <- function (indices, trcov) {
                    notOK <- is.na(trcov[,indices])
                    if (any(notOK)) {
                        warning ("detector covariate(s) missing values set to -1")
                        trcov[,indices][notOK] <- -1
                    }
                    mat <- as.matrix(trcov[,indices]) ## detectors x occasions
                    padarray(mat, dims[c(4,3)])
                }
                covs <- covindices[[i]][[variable]]  ## indices this cov, this session
                vals[[i]] <- getvals(covs, trapcov[[i]])
            }
            vals <- unlist(vals)  ## concatenate sessions
            if (!is.null(factorlevels))
                vals <- factor(vals, factorlevels)

            dframe[,variable] <<- insertdim (vals, dimcov, dims)
        }
    }
    #--------------------------------------------------------------------------------

    models$D <- NULL                          # drop density model
    models$noneuc <- NULL                     # drop non-Euclidean parameter model
    npar     <- length(models)                # real parameters
    grouplevels  <- group.levels(capthist,groups)
    ngrp    <- max(1,length(grouplevels))

    ## 'session-specific' list if MS
    MS   <- ms(capthist) # logical for multi-session
    sessionlevels <- session(capthist)
    if (is.null(sessionlevels)) sessionlevels <- '1'

    if (MS) {
        R <- length(capthist)
        n <- max(sapply(capthist, nrow))                         # max over sessions
        S <- max(sapply(capthist, ncol))                         # max over sessions
        K <- max(sapply(traps(capthist), ndetector))             # max over sessions
    }
    else {
        R <- 1
        n <- nrow(capthist)
        S <- ncol(capthist)
        K <- ndetector(traps(capthist))
    }
    ## cover unmarked case
    if (n == 0) n <- 1

    if (npar == 0) {
        ## 2014-01-25
        ## no detection parameters estimated
        ## return null design object (list with no contents)
        constantPIA <- array(1, dim=c(R,n,S,K,1))
        return(list(designMatrices = NULL, parameterTable = NULL, PIA = constantPIA, R = R))
    }

    parnames <- names(models)                 # typically c('g0', 'sigma', 'z')
    vars     <- unique (unlist(sapply (models, all.vars)))
    vars     <- vars[!(vars %in% groups)]     # groups treated separately
    nmix     <- get.nmix(models, capthist, hcov)
    trps    <- traps(capthist)                 # session-specific trap array
    used    <- usage(trps)                     # session-specific usage
    if (ignoreusage) used <- NULL
    zcov    <- covariates(capthist)            # session-specific individual covariates
    trapcov <- covariates(trps)                # session-specific trap covariates

    if (('g' %in% vars) & is.null(groups))
        stop ("requires valid 'groups' covariate")

    #--------------------------------------------------------------------------
    # timecov may be a vector or a dataframe or a list of vectors or a list of data frames
    # conform to either dataframe or list of dataframes (1 per session)
    if (!is.data.frame(timecov)) {
        if (is.list(timecov)) {
            if (length(timecov) != R)
                stop ("wrong number of sessions in 'timecov' list")
            timecov <- lapply(timecov, as.data.frame)
        }
        else timecov <- as.data.frame(timecov)
    }

    #--------------------------------------------------------------------------
    # session covariates
    if (!is.null(sessioncov)) {
        sessioncov <- as.data.frame(sessioncov)
        if (nrow(sessioncov) != R)
            stop("number of rows in 'sessioncov' should equal ",
                 "number of sessions")
    }

    #--------------------------------------------------------------------------
    dims <- c(R,n,S,K,nmix)    # 'virtual' dimensions
    dframenrow <- prod(dims)   # number of rows
    autovars <- c('session','Session','g','t','T', 'ts', 'tt',
                  'b','bn','B','bk','bkn','Bk', 'k', 'K', 'bkc', 'Bkc',
                  'kcov','tcov', 'h2', 'h3')
    #--------------------------------------------------------------------------
    # user-specified dframe
    # 2011-11-27
    if (is.null(dframe)) {
        ## data frame with dframenrow rows and no columns 2013-03-11
        dframe <- data.frame(row.names = 1:dframenrow)
        dframevars <- ""
    }
    else {
        tempn <- n
        if (nrow(dframe) !=  dframenrow )
            stop ("dframe should have ", R*tempn*S*K*nmix, " rows ( R*n*S*K*nmix )")
        dframevars <- names(dframe)
    }
    #--------------------------------------------------------------------------
    ## session, Session

    dframe$session <- factor( insertdim (sessionlevels, 1, dims),
                             levels = sessionlevels)

    if ('Session' %in% vars) {
        dframe$Session <- insertdim (0:(R-1), 1, dims)
    }
    
    #--------------------------------------------------------------------------
    ## t, T, tcov, ts

    if ('t' %in% vars) {
        dframe$t <- factor( insertdim (1:S, 3, dims) )
    }
    if ('T' %in% vars) {
        dframe$T <- insertdim (0:(S-1), c(3,1), dims)
    }
    if ('ts' %in% vars) {
        markocc <- markocc(traps(capthist))
        ## reworked 2017-09-01
        if (setequal(markocc, 0:1) | setequal(markocc, -1:1)) {
            ## treat -1 (unresolved) and 0 (sighting) occasions the same
            markocc <- pmax(markocc,0)  
            ts <- factor(c('sighting','marking')[markocc+1])   ## markocc 0 vs markocc 1
        }
        else if (setequal(markocc, -1:0)) {
            ts <- factor(c('unresolved','sighting')[markocc+2])   ## markocc -1 vs markocc 0
        }
        else
            stop ("ts is used only for mark-resight data")
        dframe$ts <- insertdim (ts, c(3,1), dims)
    }
    if ('tt' %in% vars) {
        detect <- detector(traps(capthist))
        telem <- detect == 'telemetry'
        if (!(any(telem) & !all(telem)))
            stop ("tt is appropriate for telemetry data only when mixed with another detector type")
        tt <- factor(c('nontelem','telem'))[telem+1] 
        dframe$tt <- insertdim (tt, c(3,1), dims)
    }
    if ('tcov' %in% vars) {
        if (is.null(timecov))
            stop ("requires valid time covariate 'timecov'")
        if (length(unlist(timecov)) > S)                                             ### CHECK
            warning ("length of 'timecov' exceeds number of occasions")
        if (is.data.frame(timecov)) {
            if (nrow(timecov) < S)
                stop ("requires valid time covariate 'timecov'")
            timecov <- timecov[,1,drop=F]  ## retain only first
        }
        else {
            if (any(sapply(timecov, nrow) < S))
                stop ("requires valid time covariate 'timecov'")
            timecov <- lapply (timecov, function(x) pad1(x[,1], S))
        }
        dframe$tcov <- insertdim (unlist(timecov), c(3,1), dims)
    }

    #--------------------------------------------------------------------------
    if (!is.null(groups)) {

        ################
        # add g factor

        gvar <- group.factor(capthist, groups)          # list if MS
        if (MS) gvar <- lapply(gvar, pad1, n)           # constant length
        # by animal within session
        dframe$g <- insertdim ( unlist(gvar), c(2,1), dims)  ## unlist works on factors, too

        #################################
        # also add separate group factors

        # Get group membership from covariates(capthist) for each session,
        # and pad to max(n) if needed

        for (i in groups) {
            if (MS) {
                grouping <- lapply(zcov, function(x) x[,i])
                grouping <- unlist(lapply(grouping, pad1, n))
            }
            else grouping <- zcov[,i]
            ## 2011-11-28 these insertdim seem to do nothing - should assign to dframe column
            insertdim(grouping, c(2,1), dims)
        }
    }
    #--------------------------------------------------------------------------

    ## behavioural response fields

    if (sum(c('b','bn','B','bk','bkn','Bk','bkc','Bkc', 'k', 'K') %in% vars) > 1)
    stop ("model should not use more than one type of behavioural response")

    ## assume sessions are same type
    ## 2016-10-12
    makeb <- function (caphist) {      ## global response
        # condition added 2016-10-01
        if (nrow(caphist)==0) 
            array(dim = c(0,S))
        else {
            temp0 <- apply(abs(caphist), 1:2, sum)
            t(apply(temp0, 1, prevcapt))
        }
    }
    makek <- function (caphist) {      ## trap responds to capture of any animal
        temp <- apply(abs(caphist), c(2,3), sum) # occasion x trap
        apply(temp, 2, prevcapt)
    }
    makebk <- function (caphist) {     ## individual trap-specific response
        # condition added 2016-10-01
        if (nrow(caphist)==0) 
            array(dim = c(0,S,K))
        else {
            temp <- apply(abs(caphist), c(1,3), prevcapt)
            aperm(temp, c(2,1,3))
        }
    }
    makebkc <- function (caphist) {     ## trap-specific multi-level response
        # condition added 2016-10-01
        if (nrow(caphist)==0) 
            array(dim = c(0,S,K))
        else {
            ## same animal
            caphist[,,] <- abs(caphist)
            b1 <- apply(caphist, c(1,3), prevcapt)
            b1 <- aperm(b1, c(2,1,3))
            
            ## other animal
            b2 <- array(dim = dim(caphist))
            for (i in 1:nrow(caphist)) {
                temp <- apply(caphist[-i,,], 2:3, sum)
                b2[i,,] <- apply(temp, 2, prevcapt)
            }
            ## output array dim (n,S,K)
            ## 1 none, 2 this animal, 3 other animal, 4 both
            b1 + 2 * b2 + 1
        }
    }

    #--------------------------------------------------------------------------

    if ('b' %in% vars) {
        if (naive) dframe$b <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, cumsum(x[-S])>0)
            if (MS) {
                temp <- lapply(capthist, makeb)
                temp <- lapply(temp, padarray, c(n,S))
                temp <- unlist(temp)
            }
            else temp <- makeb(capthist)  # one session
            dframe$b <- insertdim (as.vector(temp), c(2,3,1), dims)
        }
    }

    #------------------------------------------------

    if ('bn' %in% vars) {
        if (naive) dframe$bn <- rep(0, dframenrow)
        else {
            prevcapt <- function(x) c(0, cumsum(x[-S]))
            if (MS) {
                temp <- lapply(capthist, makeb)
                temp <- lapply(temp, padarray, c(n,S))
                temp <- unlist(temp)
            }
            else temp <- makeb(capthist)  # one session
            dframe$bn <- insertdim (as.vector(temp), c(2,3,1), dims)
        }
    }

    #------------------------------------------------

    if ('B' %in% vars) {
        if (naive) dframe$B <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, x[-S]>0)
            if (MS) {
                temp <- lapply(capthist, makeb)
                temp <- lapply(temp, padarray, c(n,S))
            }
            else temp <- makeb(capthist)  # one session
            dframe$B <- insertdim (as.vector(unlist(temp)), c(2,3,1), dims)
       }
    }

    #------------------------------------------------

    if ('bk' %in% vars) {
        if (naive) dframe$bk <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, cumsum(x[-S])>0)
            if (MS) {
                temp <- lapply(capthist, makebk)
                temp <- lapply(temp, padarray, c(n,S,K))
            }
            else temp <- makebk(capthist)  # one session
            dframe$bk <- insertdim(as.vector(unlist(temp)), c(2,3,4,1), dims)
        }
    }

    #------------------------------------------------
    if ('bkn' %in% vars) {
        if (naive) dframe$bkn <- rep(0, dframenrow)
        else {
            prevcapt <- function(x) c(0, cumsum(x[-S]))
            if (MS) {
                temp <- lapply(capthist, makebk)
                temp <- lapply(temp, padarray, c(n,S,K))
            }
            else temp <- makebk(capthist)  # one session
            dframe$bkn <- insertdim(as.vector(unlist(temp)), c(2,3,4,1), dims)
        }
    }

    #------------------------------------------------
    if ('Bk' %in% vars) {
        if (naive) dframe$Bk <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, x[-S]>0)
            if (MS) {
                temp <- lapply(capthist, makebk)
                temp <- lapply(temp, padarray, c(n,S,K))
            }
            else temp <- makebk(capthist)  # one session
            dframe$Bk <- insertdim(as.vector(unlist(temp)), c(2,3,4,1), dims)
       }
    }
    #--------------------------------------------------------------------------

    if ('bkc' %in% vars) {
        if (naive)
            dframe$bkc <- rep(1, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, cumsum(x[-S])>0)
            if (MS) {
                temp <- lapply(capthist, makebkc)
                temp <- lapply(temp, padarray, c(n,S,K))
            }
            else temp <- makebkc(capthist)  # one session
            temp <- unlist(temp)
            dframe$bkc <- insertdim(temp, c(2,3,4,1), dims)
        }
        dframe$bkc <- factor(dframe$bkc)
        levels(dframe$bkc) <- c('None','Self','Other','Both')
    }

    #------------------------------------------------
    if ('Bkc' %in% vars) {
        if (naive)
            dframe$Bkc <- rep(1, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, x[-S]>0)
            if (MS) {
                temp <- lapply(capthist, makebkc)
                temp <- lapply(temp, padarray, c(n,S,K))
            }
            else temp <- makebkc(capthist)  # one session
            temp <- unlist(temp)
            dframe$Bkc <- insertdim(temp, c(2,3,4,1), dims)
        }
        dframe$Bkc <- factor(dframe$Bkc)
        levels(dframe$Bkc) <- c('None','Self','Other','Both')
    }
    #--------------------------------------------------------------------------

    if ('k' %in% vars) {
        if (naive) dframe$k <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, cumsum(x[-S])>0)
            if (MS) {
                temp <- lapply(capthist, makek)
                temp <- lapply(temp, padarray, c(S,K))
            }
            else temp <- makek(capthist)  # one session
            dframe$k <- insertdim(as.vector(unlist(temp)), c(3,4,1), dims)
       }
    }
    #--------------------------------------------------------------------------

    if ('K' %in% vars) {
        if (naive) dframe$K <- rep(FALSE, dframenrow)
        else {
            prevcapt <- function(x) c(FALSE, x[-S]>0)
            if (MS) {
                temp <- lapply(capthist, makek)
                # temp <- lapply(temp, padarray, c(n,S,K))
                # bug fixed 2012-09-10
                temp <- lapply(temp, padarray, c(S,K))
            }
            else temp <- makek(capthist)  # one session
            dframe$K <- insertdim(as.vector(unlist(temp)), c(3,4,1), dims)
       }
    }
    #--------------------------------------------------------------------------
    if ('kcov' %in% vars) {
        if (is.null(trapcov))
            stop ("model uses trap covariates, but valid covariate ",
                  "data not found")
        if (is.data.frame(trapcov)) trapcov <- trapcov[,1,drop=F]  ## retain only first
        else trapcov <- lapply (trapcov, function(x) pad1(x[,1], K))
        dframe$kcov <- insertdim(unlist(trapcov), c(4,1), dims)
    }

    #--------------------------------------------------------------------------

    ## h2 or h3
    if (nmix > 1) {
        mixture <- paste('h',nmix,sep='')
        classnames <- h.levels(capthist, hcov, nmix)
        tempclass <- insertdim(classnames, 5, dims)
        dframe[,mixture] <- factor(tempclass)
    }

    #--------------------------------------------------------------------------

    ## all autovars should have now been dealt with
    vars <- vars[!(vars %in% c(autovars, dframevars))]

    #--------------------------------------------------------------------------
    # add zcov, sessioncov, timecov, trapcov

    if (CL) findvars.MS (zcov, vars, c(2,1)) ## CL only
    findvars.MS (sessioncov, vars, 1)
    findvars.MS (timecov, vars, c(3,1))      ## session-specific list
    findvars.MS (trapcov, vars, c(4,1))      ## session-specific list

    #--------------------------------------------------------------------------
    # time-varying trap covariates
    if (MS)
        tvc <- timevaryingcov(trps)
    else
        tvc <- list(timevaryingcov(trps))

    if (!is.null(tvc) & (length(vars)>0)) {
        
        findvars.traptime (tvc, vars)
    }
    #--------------------------------------------------------------------------

    if (length(vars)>0) {
        if (!is.null(zcov)) {
            if (is.data.frame(zcov))
                znames <- names(zcov)
            else
                znames <- unlist(lapply(zcov, names))
            if (any (vars %in% znames))
                stop ("seems you are trying to use individual covariates ",
                      " in a full-likelihood model")
        }
        stop ("covariate(s) ", paste(vars,collapse=","), " not found")
    }

    #=========================================================================
    # OK, we now have a complete dframe
    # Next construct the design matrices, PIA, and deal with usage
    #=========================================================================

    #------------------------------------------------------------
    make.designmatrix <- function (formula, prefix, ...) {
     # combine formula and dframe to generate design matrix
        if (is.null(formula)) {
            list (model = NULL, index = rep(1,dframenrow))
        }
        else {
            ## see utility.R for general.model.matrix
            ## allows regression splines (mgcv)
            tempmat <- general.model.matrix(formula, data = dframe, contrasts = contrasts, ...)  ## secr.design.MS contrasts
            ## drop pmix beta0 column from design matrix
            if (prefix=='pmix') {
                tempmat <- tempmat[,-1,drop=FALSE]
            }
            temp <- make.lookup (tempmat)   # retain unique rows
            list (model=temp$lookup, index=temp$index)
        }
    }
    #------------------------------------------------------------

    ## replace NA with dummy value '0' to stop model.matrix dropping rows added as padding
    ## might be more efficient to allow it to drop and then match later
    dframe[is.na(dframe)] <- 0
    # list with one component per real parameter
    # each of these is a list with components 'model' and 'index'
    designMatrices <- sapply (1:length(models), simplify=FALSE,
        function (x) make.designmatrix(models[[x]], names(models[x])), ...)
    names(designMatrices) <- names(models)

    ## dim(indices) = c(R*n*S*K*nmix, npar)
    indices <- sapply (designMatrices, function(x) x$index)
    indices <- matrix(unlist(indices), ncol = npar)

    # retain just the 'model' components of 'designMatrices'
    designMatrices <- lapply (designMatrices, function(x)x$model )

    # prefix column names in 'designMatrices' with parameter name
    for (i in 1:npar)
        colnames(designMatrices[[i]]) <- paste (parnames[i], '.',
            colnames(designMatrices[[i]]), sep='')

    # repackage indices to define unique combinations of parameters
    indices2 <- make.lookup(indices)

    #--------------------------------------------------------------------
    # PIA = Parameter Index Array
    #       index to row of parameterTable for a given R,n,S,K,nmix
    # dim(parameterTable) = c(uniqueparcomb, npar)
    #       index to row of designMatrix for each real parameter
    #--------------------------------------------------------------------

    PIA <- array(indices2$index, dim = dims)

    parameterTable <- indices2$lookup

    colnames(parameterTable) <- parnames

    #--------------------------------------------------------------------
    # Zero the index of trap+time pairs that were 'not set'
    # the external C code checks for this and sets p(detection) to zero
    #--------------------------------------------------------------------

    # 'used' is list if MS

    if ((!is.null(used)) & (length(used)>0)) {
        if (!is.null(unlist(used)))
        if (any(!unlist(used))) {
            if (!MS) {
                if (all(dim(used) == c(K,S)))
                    used <- rep(t(used),rep(n,S*K))
                else {
                    ## allowance for animal-specific 3-D usage!
                    ## otherwise assume K,S,n
                    used <- aperm(used, c(3,1,2))
                }
                # 2012-12-17
                PIA[1, , , ,] <- PIA[1, , , ,] * rep(used>0,nmix)
            }
            else for (r in 1:R) {
                if (!is.null(used[[r]])) {
                    use <- array (0, dim=c(S,K))
                    temp <- t(used[[r]])  # dim (S',K')
                    if (length(dim(temp)) == 3)
                        stop("3-D usage not available with multiple sessions")
                    use[1:nrow(temp), 1:ncol(temp)] <- temp  # padding
                    PIA[r, , , ,] <- PIA[r, , , ,] * rep(rep(use>0,rep(n,S*K)),nmix)
                }
            }
        }
    }

    #--------------------------------------------------------------------
    ## 2014-08-21

    smoothsetup <- vector(length(parnames), mode = 'list')
    names(smoothsetup) <- parnames
    for (i in parnames) {
        if (any(smooths(models[[i]]))) {
            ## collapse data to unique rows
            temp <- make.lookup (dframe)
            smoothsetup[[i]] <- gamsetup(models[[i]], temp$lookup)
        }
    }
    
    #--------------------------------------------------------------------
        individual <- individualcovariates(PIA)
    #--------------------------------------------------------------------
    if (keep.dframe) {
        ## 2013-03-11
        purge.dframe <- function (validdim) {
            ## drop padding rows in dframe
            tmps <- split(dframe, dframe$session)
            ia <- array(1:prod(dims[-1]), dim = dims[-1])
            getsess  <- function (x,valid) {
                inddf <- do.call(expand.grid, lapply(valid, seq_len))
                names(inddf) <- c('animal','occasion','detector','mixture')
                ind <- as.matrix(inddf)
                x <- x[as.vector(ia[ind]),]
                x <- cbind(inddf, x)
                row.names(x) <- 1:nrow(x)
                x
            }
            dfs <- mapply(getsess, tmps, validdim, SIMPLIFY = FALSE)
            do.call(rbind, dfs)
        }

        if (MS) {
            validdim <- as.list(data.frame(t(matrix(c(sapply(capthist, nrow),
                          sapply(capthist, ncol),
                          sapply(trps, ndetector),
                          rep(nmix,R)), nrow = R))))
            names(validdim) <- sessionlevels
        }
        else
            validdim <- list(dim(PIA)[2:5])

        ## optionally purge rows that pad each session to constant nrows
        if (!full.dframe)
            dframe <- purge.dframe(validdim)
        else {
            inddf <- do.call(expand.grid, lapply(dims, seq_len))[,-1]
            names(inddf) <- c('animal','occasion','detector','mixture')
            dframe <- cbind(inddf,dframe)
        }
        ## 2014-08-22
        ## dframe <- dframe[,c(5,1:4,6:ncol(dframe))]
        dframe[,1:5] <- dframe[,c(5,1:4)]
        list(designMatrices = designMatrices, parameterTable = parameterTable, PIA = PIA, R = R,
             dframe = dframe, validdim = validdim, smoothsetup = smoothsetup,
             individual = individual)
    }
    else
        list(designMatrices = designMatrices, parameterTable = parameterTable, PIA = PIA, R = R,
             smoothsetup = smoothsetup, individual = individual)
}
############################################################################################



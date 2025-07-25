############################################################################################
## package 'secr'
## makeNewData.R
## last changed
## 2009 12 13 (mixtures)
## 2010 03 10 'T'
## 2010 06 17 'Session'
## 2010 06 21 'x2', 'y2', 'xy'
## 2010 08 28 fix bug with T
## 2011 11 28 user dframe factors now covered
## 2015-10-08 'ts'
## 2017-12-18 all.levels argument
## 2021-03-24 fix all.levels = FALSE bug
## 2024-10-04 renamed file
## Create (neutral) design data suitable for 'predict'
## generic method makeNewData
############################################################################################

makeNewData <- function (object, all.levels = FALSE, ...) UseMethod("makeNewData")

makeNewData.default <- function (object, all.levels = FALSE, ...) {
    cat ('no makeNewData method for objects of class', class(object), '\n')
}

makeNewData.secr <- function (object, all.levels = FALSE, bytrap = FALSE, ...) {
    # Session treated separately later
    autovars <- c('g','x','y','x2','y2','xy','session',
                  't','T','ts','b','B','bk','Bk','Br','bkc','Bkc','k','K','tcov','kcov','h2','h3')
    capthist <- object$capthist
    trps <- traps(capthist)
    mask <- object$mask
    vars <- object$vars
    groups <- object$groups
    timecov <- object$timecov
    sessioncov <- object$sessioncov
    nmix <- object$details$nmix
    hcov <- object$hcov
    
    if(is.null(nmix)) nmix <- 1
    mixvar <- switch(nmix, character(0),'h2','h3')
    
    nocc <- if (ms(capthist)) sapply(capthist, ncol) else ncol(capthist)
    nocc <- max(nocc)
    grouplevels <- secr_group.levels(capthist, groups)
    ngrp <- max(1, length(grouplevels))
    sessions <- session(capthist)
    R <- length(sessions)
    dims <- c(R, ngrp, nmix)
    
    onesession <- function(session) {
        findvars <- function (basevars, cov, source) {
            ## function to add covariates to a list
            ## cov should be dataframe or list of dataframes, one per session (R > 1),
            if (!is.data.frame(cov)) cov <- cov[[session]] ## assume multisession list
            if (is.null(cov) | (length(cov)==0) | (length(sessvars)==0)) return(basevars)
            else {
                found <- ''
                for (v in sessvars) {
                    if (v %in% names(cov)) {
                        vals <- cov[,v]
                        if (is.character(vals)) vals <- factor(vals)
                        basevars[[v]] <- if (is.factor(vals))
                            factor(levels(vals), levels = levels(vals))
                        else
                            unique(vals)
                        if (source == 'mask') {
                            if (is.factor(vals))
                                basevars[[v]] <- basevars[[v]][1]
                            else
                                basevars[[v]] <- mean(basevars[[v]])
                        }
                        
                        found <- c(found, v)
                    }
                }
                sessvars <<- sessvars[!(sessvars %in% found)]
                return(basevars)
            }
        }
        
        ## missing timevarying...
        ## 2021-07-30
        
        sessvars <- vars
        if (ms(capthist)) 
            trps <- traps(capthist[[session]])
        else 
            trps <- traps(capthist)
        
        basevars <- list(session = factor(sessions[session], levels=sessions))
        if (ngrp>1) basevars$g <- factor(grouplevels)
        if (nmix>1) basevars[mixvar] <- list(secr_h.levels(capthist, hcov, nmix))
        
        for (v in sessvars) {
            if (v=='x')  basevars$x <- 0     # mean attr(mask,'meanSD')[1,'x']
            if (v=='y')  basevars$y <- 0     # mean attr(mask,'meanSD')[1,'y']
            if (v=='x2') basevars$x2 <- 0   # mean attr(mask,'meanSD')[1,'x']
            if (v=='y2') basevars$y2 <- 0   # mean attr(mask,'meanSD')[1,'y']
            if (v=='xy') basevars$xy <- 0   # mean attr(mask,'meanSD')[1,'x']
            if (v=='T')  basevars$T <- 0   
            
            if (v=='t')  basevars$t <- factor(1:nocc)
            if (v=='ts') basevars$ts <- factor(c('marking','sighting'))
            if (v=='b')  basevars$b <- factor(0:1)
            if (v=='B')  basevars$B <- factor(0:1)
            if (v=='bk') basevars$bk <- factor(0:1)
            if (v=='Bk') basevars$Bk <- factor(0:1) 
            if (v=='Br') basevars$Br <- 0:1
            if (v=='k')  basevars$k <- factor(0:1)
            if (v=='K')  basevars$K <- factor(0:1)
            NSOB <- c('None','Self','Other','Both')
            if (v=='bkc') basevars$bkc <- factor(NSOB, levels = NSOB)
            if (v=='Bkc') basevars$Bkc <- factor(NSOB, levels = NSOB)
            
            if (v=='tcov') {
                timecov <- object$timecov
                if (is.factor(timecov)) {
                    basevars$tcov <- unique(timecov)
                }
                else
                    basevars$tcov <- 0        # ideally use mean or standardize?
            }
            if (v=='kcov') {
                kcov <- covariates(trps)[,1]
                if (is.factor(kcov)) {
                    basevars$kcov <- unique(kcov)
                }
                else {
                    basevars$kcov <- 0   
                }
            }
        }
        ## all autovars except Session should now have been dealt with
        sessvars <- sessvars[!sessvars %in% autovars]
        if (ngrp==1) 
            basevars <- findvars (basevars, covariates(capthist), 'CH') ## individual covariates
        
        basevars <- findvars (basevars, timecov, 'time')
        if (bytrap) {
            basevars <- c(basevars, list(trapID = 1:ndet[session]))
            trapvars <- sessvars[sessvars %in% trapcovnames]
            sessvars <- sessvars[!(sessvars %in% trapvars)]
        }
        else {
            trapvars <- character(0)
            basevars <- findvars (basevars, covariates(trps), 'traps')
        }
        basevars <- findvars (basevars, covariates(mask), 'mask')
        
        ## revert to first level (original default)
        ## 2021-03-24 repaired in 4.3.4
        for (v in names(basevars)) {
            # if (length(v)>0) { 
            if (!all.levels & !(v %in% c('session', 'g', 'h2','h3','trapID'))) {
                basevars[[v]] <- basevars[[v]][1] 
            }
        }
        sessnewdat <- expand.grid(basevars)
        if (bytrap) {
            if (length(trapvars)>0) {
            # need to match by session
            sessnewdat <- cbind(sessnewdat, covariates(trps)[sessnewdat$trapID, trapvars])
            }
        }
        sessnewdat
    }   # end of onesession

    # some setup
    if (bytrap) {
        if (ms(capthist)) {
            ndet <- sapply(traps(capthist), secr_ndetector)
            trapcovnames <- names(covariates(traps(capthist)[[1]]))
        }
        else {
            ndet <- secr_ndetector(trps)
            trapcovnames <- names(covariates(traps(capthist)))
        }
    }
    
    newdata <- lapply(1:length(sessions), onesession)
    newdata <- do.call(rbind, newdata)

    ## 2020-08-09
    if (!is.null(sessioncov)) {
        for (i in names(sessioncov)) {
            if ((i %in% vars) & !(i %in% names(newdata)))
            newdata[,i] <- sessioncov[newdata$session,i]
        }
    }
    
    if ('Session' %in% vars) 
        newdata$Session <- as.numeric(newdata$session) - 1   
    newdata
    
}
############################################################################################


################################################################################
## package 'secr'
## utility.R
################################################################################

## 2019-07-27 secr 4.0.0
## 2019-07-27 makelookupcpp replaces makelookup
## 2019-08-12 individualcovariates
## 2019-08-14 setNumThreads
## 2019-11-29 secr 4.1.0
## 2020-01-08 distancetotrap revised for polygon detectors
## 2020-01-26 selectCHsession
## 2020-02-21 secr 4.2.0
## 2020-05-15 stringsAsFactors function
## 2020-07-14 secr 4.3.0 distmat
## 2020-09-05 getknownclass factor bug fixed
## 2021-04-25 4.4.0
## 2021-12-16 tidy up transformations, allow arbitrary link X(), invX(), se.invX()
## 2022-01-04 4.5.0
## 2022-01-18 4.5.1
## 2022-01-23 uniformusage() function to assign all-ones usage matrix
## 2022-02-10 4.5.2 shift to sf where possible
## 2022-04-22 4.5.4
## 2022-05-31 4.5.5
## 2022-10-08 4.5.7
## 2022-11-29 4.5.8
## 2023-03-09 4.5.9
## 2023-03-10 4.5.11
## 2023-03-10 distancetotrap and nearesttrap moved to separate file
## 2023-03-10 setNumThreads moved to separate file
## 2023-05-21 4.6.0
## 2024-03-19 rlnormCV
## 2024-07-31 secr_addzeroCH tweaked to allow zero-row covariate df + drop = FALSE
## 2024-09-25 purged a couple of unused fn, moved xy2CH to xy2CH.R
## 2024-10-09 span() (now in plot.capthist.R)
## 2025-01-07 secr_allzero bug fixed
## 2025-03-18 secr_saveprogress()
## 2025-06-17 filterw(), captinhood()
## 2025-06-25 detectfn 20 OU
## 2025-07-20 5.3.0
## 2025-07-24 secr_ prefix attached to most functions used in other .R files
################################################################################

# Global variables in namespace
#
## define a local environment for temporary variables e.g. iter
## e.g. Roger Peng https://stat.ethz.ch/pipermail/r-devel/2009-March/052883.html

.localstuff <- new.env()

.localstuff$packageType <- ' pre-release'
##.localstuff$packageType <- ''

.localstuff$validdetectors <- c('single','multi','proximity','count', 
    'polygonX', 'transectX', 'signal', 'polygon', 'transect', 
    'capped', 'null','null','null','null', 'telemetry', 'signalnoise')
.localstuff$simpledetectors <- c('single','multi','proximity','count', 'capped')
.localstuff$individualdetectors <- c('single','multi','proximity','count',
    'polygonX', 'transectX', 'signal', 'signalnoise', 'polygon', 'transect',
                                     'telemetry', 'capped')
.localstuff$pointdetectors <- c('single','multi','proximity','count',
    'signal', 'signalnoise', 'unmarked','presence','capped')
.localstuff$polydetectors <- c('polygon','transect','polygonX','transectX')
.localstuff$exclusivedetectors <- c('single','multi','polygonX','transectX')
.localstuff$countdetectors <- c('count','polygon','transect','unmarked','telemetry')
.localstuff$iter <- 0   ## counter 1
.localstuff$iter2 <- 0  ## counter 2
.localstuff$detectionfunctions <-
  c('halfnormal',
    'hazard rate',
    'exponential',
    'compound halfnormal',
    'uniform',
    'w exponential',
    'annular normal',
    'cumulative lognormal',
    'cumulative gamma',
    'binary signal strength',
    'signal strength',
    'signal strength spherical',
    'signal-noise',
    'signal-noise spherical',
    'hazard halfnormal',
    'hazard hazard rate',
    'hazard exponential',
    'hazard annular normal',
    'hazard cumulative gamma',
    'hazard variable power',
    'Ornstein-Uhlenbeck')

.localstuff$DFN <- c('HN', 'HR', 'EX', 'CHN', 'UN', 'WEX', 'ANN', 'CLN', 'CG',
  'BSS', 'SS', 'SSS', 'SN', 'SNS',
  'HHN', 'HHR', 'HEX', 'HAN', 'HCG', 'HVP', 'OU')

## Bk added 2020-02-26; Br added 2025-06-17
.localstuff$learnedresponses <- c('b', 'bk', 'B', 'k', 'Bk', 'Br') 

## 2025-07-22

.localstuff$spatialparameters <- c('noneuc','sigmaxy','lambda0xy','a0xy','sigmakxy')
.localstuff$spatialparametersD <- c('D','noneuc','sigmaxy','lambda0xy','a0xy','sigmakxy')
#-------------------------------------------------------------------------------

secr_detectionfunctionname <- function (fn) {
    .localstuff$detectionfunctions[fn+1]
}

#-------------------------------------------------------------------------------

secr_detectionfunctionnumber <- function (detname) {
    dfn <- match (toupper(detname), .localstuff$DFN)
    if (is.na(dfn))
        dfn <- match (tolower(detname), .localstuff$detectionfunctions)
    if (is.na(dfn))
        stop ("unrecognised detection function ", detname)
    dfn-1
}

#-------------------------------------------------------------------------------

secr_parnames <- function (detectfn) {
    switch (detectfn+1,
        c('g0','sigma'),   ## 0
        c('g0','sigma','z'),
        c('g0','sigma'),
        c('g0','sigma','z'),
        c('g0','sigma'),
        c('g0','sigma','w'),
        c('g0','sigma','w'),
        c('g0','sigma','z'),
        c('g0','sigma','z'),
        c('b0','b1'),
        c('beta0','beta1', 'sdS'),    ## include cutval?
        c('beta0','beta1', 'sdS'),    ## include cutval?
        c('beta0','beta1', 'sdS','muN','sdN'),
        c('beta0','beta1', 'sdS','muN','sdN'),
        c('lambda0','sigma'),
        c('lambda0','sigma','z'),
        c('lambda0','sigma'),
        c('lambda0','sigma','w'),
        c('lambda0','sigma','z'),
        c('lambda0','sigma','z'),
        c('epsilon','sigma','tau')    ## 20   2025-06-25
    )
}

#-------------------------------------------------------------------------------

secr_getdfn <- function (detectfn) {
    switch (detectfn+1, HN, HR, EX, CHN, UN, WEX, ANN, CLN, CG, BSS, SS, SSS,
                       SN, SNS, HHN, HHR, HEX, HAN, HCG, HVP, HHN)
}

#-------------------------------------------------------------------------------

secr_valid.detectfn <- function (detectfn, valid = c(0:3,5:19)) {
# exclude 4 uniform: too numerically flakey
    if (is.null(detectfn))
        stop ("requires 'detectfn'")
    if (is.character(detectfn))
        detectfn <- secr_detectionfunctionnumber(detectfn)
    if (any(!(detectfn %in% valid)))    # allow vector of detectfn 2024-02-12
        stop ("invalid detection function")
    detectfn
}

#-------------------------------------------------------------------------------

secr_valid.detectpar <- function (detectpar, detectfn) {
    if (is.null(detectpar) | is.null(detectfn))
        stop ("requires valid 'detectpar' and 'detectfn'")

    ## 2013-07-19, 2013-10-22
    ## replace a0 with g0 or lambda0 as appropriate to detectfn
    if ('a0' %in% names(detectpar)) {
        aname <- if (detectfn %in% 0:8) 'g0' else 'lambda0'
        lambda0 <- detectpar[['a0']] / 2 / pi / detectpar[[2]]^2 * 10000
        detectpar[[aname]] <- if (detectfn %in% 0:8) 1-exp(-lambda0) else lambda0
    }

    if (!all(secr_parnames(detectfn) %in% names(detectpar)))
        stop ("requires 'detectpar' ", paste(secr_parnames(detectfn), collapse=','),
            " for ", secr_detectionfunctionname(detectfn), " detectfn")
    detectpar[secr_parnames(detectfn)]
}

#-------------------------------------------------------------------------------

secr_valid.model <- function(model, CL, detectfn, hcov, userdist, sessioncovnames) {
    badsmooths <- function (formula) {
        ## does smooth specification conform to secr requirements?
        ## returns TRUE/FALSE
        labels <- attr(terms(formula), "term.labels")
        if (length(labels) > 0) {
            smoothterms <- sapply(labels, function (x)
                any(sapply(c("s\\(", "te\\("), grepl, x)))
            labels <- labels[smoothterms]
            any(sapply(labels, function(x)
                grepl("s\\(", x) & !grepl("k =", x))) |
                any(sapply(labels, function(x)
                    grepl("te\\(", x) & (!grepl("fx = TRUE", x) | !grepl("k =", x))))
        }
        else
            FALSE
    }
    
    if (any(sapply(model, badsmooths))) {
        warning ("smooth term may be unsuitable for secr: ",
                 "does not specify k or fx where required")
    }
}

#-------------------------------------------------------------------------------

secr_getuserdistnames <- function (userdist) {
    ## return the names of any supplementary arguments of user-provided function
    ## for non-euclidean distance computations
    if (is.function(userdist)) {
        distnames <- try(userdist(), silent = TRUE)
        if (!is.character(distnames))
            stop("invalid userdist function - ",
                 "should return parameter names when called with no arguments")
        distnames
    }
    else
        character(0)
}

#-------------------------------------------------------------------------------

secr_getuserdist <- function (traps, mask, userdist, sessnum, NElist, density, ...) {
    ## Apply user-provided distance function or basic distance function secr_getdistmat2()
    if (is.null(userdist)) {
        secr_getdistmat2(traps, mask, NULL)
    }
    else {
        userdistnames <- secr_getuserdistnames(userdist)
        m <- nrow(mask)
        if (is.null(covariates(mask)))
            covariates(mask) <- data.frame(row.names = 1:m)
        
        if (length(NElist)>0) {
            noneucpar <- lapply(NElist, secr_getmaskpar, OK = TRUE, m, sessnum, FALSE, NULL)
            noneucpar <- as.data.frame(noneucpar)
            covariates(mask) <- cbind(covariates(mask), noneucpar)
        }
        
        if (('D' %in% userdistnames) && !is.null(density))
            covariates(mask)$D <- density
        
        ## pass miscellaneous unmodelled parameter(s)
        extra <- list(...)
        for (n in names(extra)) attr(mask, n) <- extra[n]
        
        distmat2 <- secr_valid.userdist (userdist,
                                         detector(traps),
                                         xy1 = traps,
                                         xy2 = mask,
                                         mask = mask,
                                         sessnum = sessnum)^2
        baddist <- (!is.finite(distmat2)) | (distmat2<0) | is.na(distmat2)
        if (any(baddist)) {
            warning ("replacing infinite, negative and NA userdist values with 1e10")
            distmat2[baddist] <- 1e10
        }
        distmat2
    }
}
#--------------------------------------------------------------------------------

secr_valid.pnames <- function (details, CL, detectfn, alltelem, sighting, nmix) {
    ## modelled parameters
    pnames <- switch (detectfn+1,
        c('g0','sigma'),           # 0 halfnormal
        c('g0','sigma','z'),       # 1 hazard rate
        c('g0','sigma'),           # 2 exponential
        c('g0','sigma','z'),       # 3
        c('g0','sigma'),           # 4
        c('g0','sigma','w'),       # 5
        c('g0','sigma','w'),       # 6
        c('g0','sigma','z'),       # 7
        c('g0','sigma','z'),       # 8
        c('b0','b1'),              # 9
        c('beta0','beta1','sdS'),  # 10
        c('beta0','beta1','sdS'),  # 11
        c('beta0','beta1','sdS'),  # 12  cf secr_parnames() in utility.R: muN, sdN?
        c('beta0','beta1','sdS'),  # 13  cf secr_parnames() in utility.R: muN, sdN?
        c('lambda0','sigma'),      # 14 hazard halfnormal
        c('lambda0','sigma','z'),  # 15 hazard hazard rate
        c('lambda0','sigma'),      # 16 hazard exponential
        c('lambda0','sigma','w'),  # 17
        c('lambda0','sigma','z'),  # 18
        c('lambda0','sigma','z'),  # 19
        c('epsilon','sigma','tau'))  # 20 OU 2025-06-25

    if (details$param %in% c(2,6))
        pnames[1] <- 'esa'
    if (details$param %in% c(3,5))
        pnames[1] <- 'a0'
    if (details$param %in% 4:6) {
        pnames[2] <- 'sigmak'
        pnames <- c(pnames, 'c')
        pnames <- c(pnames, 'd')
    }
    if (!CL || details$relativeD) {
        # include density D if needed
        pnames <- c('D', pnames)
    }
    ## 'noneuc', 'sigmaxy', 'lambda0xy', 'a0xy' etc.
    for (parm in secr_getuserdistnames(details$userdist)) {
        parm <- parm[parm != 'D']   # drop unwanted
        pnames <- c(pnames, parm)
    }
    if (sighting)
      pnames <- c(pnames, 'pID')
    # if (alltelem) {
    #     rnum <- match(c('D','lambda0','a0','esa','g0'), pnames)
    #     rnum[is.na(rnum)] <- 0
    #     pnames <- pnames[-rnum]
    # }
    if (nmix>1)
        pnames <- c(pnames, 'pmix')
    pnames
}
#-------------------------------------------------------------------------------

secr_valid.userdist <- function (userdist, detector, xy1, xy2, mask, sessnum) {
    if (is.null(userdist)) {
        ## default to Euclidean distance
        result <- edist(xy1, xy2)
    }
    else {
        if (any(detector %in% .localstuff$polydetectors)) {
            stop ("userdist cannot be used with polygon detector types;")
        }
        if (is.function(userdist))
        {
            OK <- secr_getuserdistnames(userdist) %in% names(covariates(mask))
            if ((length(OK)>0) & !all(OK))
                stop ("covariates required by userdist function not in mask : ",
                      paste(secr_getuserdistnames(userdist)[!OK], collapse=','))
            # 2023-02-06 selected columns 1:2 only (mask passes miscparm)
            result <- do.call(userdist, c(list(xy1[,1:2], xy2[,1:2], mask)))
        }
        else {
            if (is.character(userdist)) {
                userdist <- get(userdist, pos=-1)
            }

            if (is.list(userdist) & !is.data.frame(userdist)) {
                if (missing(sessnum) || is.na(sessnum))
                    stop("This use does not yet allow for session-specific userdist")
                else
                    result <- userdist[[sessnum]]
            }
            else {
                result <- userdist
            }
        }
        if (!all(dim(result) == c(nrow(xy1), nrow(xy2))))
            stop ("invalid distance matrix dim = ", dim(result)[1], ',', dim(result)[2])
        baddist <- (!is.finite(result)) | (result<0) | is.na(result)
        if (any(baddist)) {
            warning ("replacing infinite, negative and NA userdist values with 1e10", 
                     call. = FALSE)
            result[baddist] <- 1e10
        }
    }
    result
}
#-------------------------------------------------------------------------------

secr_new.param <- function (details, model, CL) {
    esa <- 'esa' %in% names(model)
    a0 <- 'a0' %in% names(model)
    sigmak <- 'sigmak' %in% names(model)
    newparam <- details$param
    if (esa & !sigmak) {
        newparam <- 2
    }
    if (a0 & !sigmak) {
        newparam <- 3
    }
    if (sigmak) {
        if (esa) {
            newparam <- 6
        }
        else {
            if (CL)
                stop ("sigmak parameterization requires full model, not CL, unless also 'esa'")
            newparam <- ifelse(a0, 5, 4)
        }
    }
    if (newparam  != details$param)
        warning ("Using parameterization details$param = ", newparam, call. = FALSE)
    newparam
}

#-------------------------------------------------------------------------------
## MULTI-SESSION FORM?

secr_detectorcode <- function (object, MLonly = TRUE, noccasions = NULL) {
    ## numeric detector code from a traps object
    detcode <- sapply(detector(object), switch,
        single      = -1,
        multi       = 0,
        proximity   = 1,
        count       = 2,
        polygonX    = 3,
        transectX   = 4,
        signal      = 5,
        polygon     = 6,
        transect    = 7,
        capped      = 8,
        unmarked    = 10,
        presence    = 11,
        signalnoise = 12,
        telemetry   = 13,
        -2)
    
    if (MLonly) {
        detcode <- ifelse (detcode==-1, rep(0,length(detcode)), detcode)
        if (any(detcode<0))
            stop ("Unrecognised detector type")
    }

    if (!is.null(noccasions) & (length(detcode)==1))
        detcode <- rep(detcode, noccasions)
    detcode
}
#-------------------------------------------------------------------------------

secr_expanddet <- function(CH) {
    trps <- traps(CH)
    if (is.null(trps))
        return ('nonspatial')
    else {
        det <- detector(trps)
        if (length(det)<ncol(CH))
            rep(det[1], ncol(CH))
        else det
    }
}

#-------------------------------------------------------------------------------

secr_replacedefaults <- function (default, user) replace(default, names(user), user)

#-------------------------------------------------------------------------------

secr_memo <- function (text, trace) {
    ## could use message(text), but does not immediately flush console
    if (trace) { cat (text, '\n')
    flush.console() }
}

#-------------------------------------------------------------------------------

## regularize a list of formulae
secr_stdform <- function (flist) {
    LHS <- function (form) {
        trms <- as.character (form)
        if (length(trms)==2) '' else trms[2]
    }
    RHS <- function (form) {
        trms <- as.character (form)
        ## 2020-05-14 for compatibility with R 4.0
        if (length(trms)==3) as.formula(paste(trms[c(1,3)], collapse = " ")) else form
    }
    lhs <- sapply(flist, LHS)
    temp <- lapply(flist, RHS)
    if (is.null(names(flist))) names(temp) <- lhs
    else names(temp) <- ifelse(names(flist) == '', lhs, names(flist))
    temp
}

#-------------------------------------------------------------------------------

secr_lnbinomial <- function (x,size,prob) {
    # dbinom allowing non-integer x, forcing log = TRUE
    if (x <= size) {
        lgamma (size+1) - lgamma (size-x+1) - lgamma (x+1) +
            x * log(prob) + (size-x) * log (1-prob)
    }
    else {
        -Inf
    }
}

#-------------------------------------------------------------------------------

secr_model.string <- function (model, userDfn) {
    # 2023-04-16 Note: model should be a list
    if (!is.null(userDfn)) {
        if (!is.null(model$D))
            model$D <- paste('~userD', userDfn('name'), sep='.')
    }
    temp <- paste (names(model), as.character(model), collapse=' ', sep='')
    temp
}

#-------------------------------------------------------------------------------

secr_fixed.string <- function (fixed) {
    if (is.null(fixed) | length(fixed)==0) 'none'
    else paste (names(fixed), as.character(fixed), collapse=', ', sep=' = ')
}

#-------------------------------------------------------------------------------

secr_var.in.model <- function(v,m) v %in% unlist(lapply(m, all.vars))

#-------------------------------------------------------------------------------

secr_get.nmix <- function (model, capthist, hcov) {
    model$D <- NULL  ## ignore density model
    model$pmix <- NULL ## pmix alone cannot make this a mixture model
    nmix <- 1
    if (any(secr_var.in.model('h2', model))) {
        nmix <- 2
        if (any(secr_var.in.model('h3', model)))
            stop ("do not combine h2 and h3")
    }
    if (any(secr_var.in.model('h3', model))) {
        nmix <- 3
    }
    if ((nmix == 1) & (!is.null(hcov))) {
        if (ms(capthist))
            capthist <- capthist[[1]]
        if (is.factor(covariates(capthist)[,hcov]))
            lev <- levels(covariates(capthist)[,hcov])
        else
            lev <- levels(factor(covariates(capthist)[,hcov]))
        if (all(is.na(covariates(capthist)[,hcov])))
            stop ("hcov missing for all individuals, but detection model invariant")
        if (length(lev) < 2)
            stop ("hcov covariate not found or has fewer than 2 levels")
        if (length(lev) > 2)
            warning ("hcov covariate has more than 2 levels; using only first two", call. = FALSE)
        nmix <- 2
    }
    nmix
}

#-------------------------------------------------------------------------------

secr_add.cl <- function (df, alpha, loginterval, lowerbound = 0) {

## add lognormal or standard Wald intervals to dataframe with columns
## 'estimate' and 'SE.estimate'
## lowerbound added 2011-07-15
    z <- abs(qnorm(1-alpha/2))
    if (loginterval) {
        delta <- df$estimate - lowerbound
        df$lcl <- delta / exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
        df$ucl <- delta * exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
    }
    else {
        df$lcl <- pmax(lowerbound, df$estimate - z * df$SE.estimate)
        df$ucl <- df$estimate + z * df$SE.estimate
    }
    df
}

#-------------------------------------------------------------------------------

secr_spatialscale <- function (object, detectfn, sessnum = 1) {
    if (inherits(object, 'secr')) {
        if (ms(object))
            detpar <- detectpar(object)[[sessnum]]
        else
            detpar <- detectpar(object)
        cutval <- object$details$cutval
    }
    else {
        detpar <- object
        cutval <- object$cutval
    }
    
    if (!is.null(detpar$sigma)) detpar$sigma
    else if (detectfn == 10) {
        (cutval - detpar$beta0) / detpar$beta1
    }
    else if (detectfn == 11) {
        d11 <- function(d, beta0, beta1, c) beta0 +
            beta1 * (d-1) - 10 * log10(d^2) - c
        interval <- c(0,10 * (cutval - detpar$beta0) / detpar$beta1)
        uniroot (d11, interval, detpar$beta0, detpar$beta1, cutval)$root
    }
    else if (detectfn == 9) {
        - 1 / detpar$b1   
    }
    else stop ("unrecognised detectfn")
}

#-------------------------------------------------------------------------------

## logical for whether object specifies userDfn
secr_userD <- function (object) {
  if (!inherits(object, c('secr','ipsecr')))
    stop ("requires fitted model")
  !is.null(object$details$userDfn)
}

#-------------------------------------------------------------------------------

secr_nclusters <- function (capthist) {
    if (ms(capthist)) {
        lapply(capthist, secr_nclusters)
    }
    else 	{
        nmash <- attr(capthist, 'n.mash')
        ifelse (is.null(nmash), 1, length(nmash))
    }
}

#-------------------------------------------------------------------------------

## clunky but effective re-write 2012-09-04, improved 2016-02-20, 2016-05-10
secr_leadingzero <- function (x) {
    xc <- as.character(x)
    w <- max(nchar(xc))
    n0 <- function(n) paste(rep('0',n), collapse='')
    paste(sapply(w-nchar(xc), n0), x, sep='')

    ## or, 2016-01-15, 2016-02-20 BUT DOESN'T HANDLE NON-INTEGER 2016-05-10
    #     if (is.character(x)) x <- as.numeric(x)
    #     sprintf(paste("%0", w, "d", sep = ""), x)
}

#-------------------------------------------------------------------------------

secr_group.levels <- function (capthist, groups, sep='.') {
    # 2016-06-05 use also for trap strata
    if (inherits(capthist, 'list')) {
        temp <- lapply(capthist, secr_group.levels, groups, sep)
        unique(unlist(temp))  ## vector of global levels
    }
    else {
        if (is.null(groups)) 0
        else {
            if (!all(groups %in% names(covariates(capthist))))
                stop ("one or more grouping variables is missing ",
                      "from covariates")
            temp <- as.data.frame(covariates(capthist)[,groups])
            levels(interaction(temp, drop = TRUE, sep = sep, lex.order = FALSE))
        }
    }
}
#-------------------------------------------------------------------------------

secr_group.factor <- function (capthist, groups)
    ## convert a set of grouping factors to a single factor (g)
    ## levels common to all sessions
{
    if (inherits(capthist, 'list')) {
        temp <- lapply(capthist, secr_group.factor, groups)  ## recursive call
        grouplevels <- secr_group.levels(capthist, groups)
        if (length(grouplevels)<2)
            temp
        else
            # list; force shared factor levels on each component
            lapply (temp, factor, levels = grouplevels)
    }
    else {
        if (is.null(groups) | (length(groups)==0) ) {
            return (factor(rep(1, nrow(capthist)), levels = 1))  
        }
        temp <- as.data.frame(covariates(capthist)[,groups])
        if (ncol(temp) != length(groups)) {
            stop ("one or more grouping variables is missing from ",
                  "covariates(capthist)")
        }
        temp <- interaction(temp, drop = TRUE, sep = '.', lex.order = FALSE) 
        temp
    }
}

#-------------------------------------------------------------------------------

secr_getgrpnum <- function (capthist, groups) {
    # vector of group factor values for individuals in single-session capthist
    if (is.null(groups))
        rep(1, nrow(capthist))
    else
        match(secr_group.factor(capthist, groups), secr_group.levels(capthist, groups))
}

#-------------------------------------------------------------------------------

secr_h.levels <- function (capthist, hcov, nmix) {
    ## determine the first nmix levels of a factor individual covariate
    if (is.null(hcov))
        as.character(1:nmix)
    else {
        if (ms(capthist)) {
            ## take first session as we can assume factor covariates have same levels in
            ## all sessions
            capthist <- capthist[[1]]
        }
        hcov <- covariates(capthist)[,hcov]
        if (!is.factor(hcov)) {
            warning ("hcov was coerced to a factor", call. = FALSE)
            hcov <- factor(hcov)
        }
        levels(hcov)[1:nmix]
    }
}

#-------------------------------------------------------------------------------

## Return an integer vector of class membership defined by a categorical
## individual covariate in a capthist object. Individuals of unknown
## class (including those with class exceeding nmix) are coded 1,
## others as (class number + 1). When no mixture is specified (nmix == 1)
## all are coded as unknown.

## knownclass 1 'unknown' 
## knownclass 2 'latent class 1' 
## knownclass 3 'latent class 2' 

secr_getknownclass <- function(capthist, nmix, hcov) {
    if (ms(capthist)) {
        lapply(capthist, secr_getknownclass, nmix = nmix, hcov = hcov)
    }
    else {
        if ((nmix>1) & (!is.null(hcov))) {
          ## 2020-09-05 use as.factor() instead of factor() to coerce 
          ## (if already factor, coercing with factor() loses old levels)
          var <- as.factor(covariates(capthist)[,hcov])
          tmp <- as.numeric(var) + 1
          tmp[is.na(tmp) | (tmp>(nmix+1))] <- 1
          attr(tmp,'levels') <- levels(factor(covariates(capthist)
            [,hcov]))[1:nmix]
          tmp
        }
        else
            rep(1,nrow(capthist))
    }
}

#-------------------------------------------------------------------------------

## inflate a convex outline along all radii by linear factor 'rmult'
secr_inflate <- function (xy, rmult = 1) {
    xy <- as.matrix(xy)
    centre <- apply(xy, 2, mean)
    xy <- sweep(xy, MARGIN = 2, STATS = centre, FUN = '-')
    r <- apply(xy, 1, function(z) sqrt(sum(z^2)))
    theta <- atan2 (xy[,2], xy[,1])
    r <- r * rmult
    xy <- cbind(r * cos(theta), r * sin(theta))
    sweep(xy, MARGIN = 2, STATS = centre, FUN = '+')
}

#-------------------------------------------------------------------------------

## moved from pdot.R 2013-11-09
## scalar 2016-10-14
secr_getbinomN <- function (binomN, detectr) {
    if (any(detectr %in% .localstuff$countdetectors)) {
        if (is.null(binomN))
            return(0)
        else if (binomN == 'usage')
            return(1)
        else
            return(binomN)
    }
    else
        return(1)
}

#-------------------------------------------------------------------------------

## expand beta parameter vector using template of 'fixed beta'
## fixed beta fb input is missing (NA) for estimated beta parameters
secr_fullbeta <- function (beta, fb) {
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta  ## partial beta (varying only)
        beta <- fb             ## complete beta
    }
    beta
}
#-------------------------------------------------------------------------------

secr_NEmodelled <- function (details, fixed, NEnames) {
    userdistnames <- secr_getuserdistnames(details$userdist)
    sapply(NEnames, function(NEname)
        (NEname %in% userdistnames) && is.null(fixed[[NEname]])
    )
}
#-------------------------------------------------------------------------------

secr_fullbetanames <- function (object) {
    # 2024-12-23
    betanames <- unlist(sapply(object$design$designMatrices, colnames))
    names(betanames) <- NULL
    if(!is.null(attr(object$designD, 'Dfn'))) {
        nDbeta <- attr(object$designD, 'Dfn')(object$designD)
        Dnames <- paste0('D', 1:nDbeta)
    }
    else {
        Dnames <- colnames(object$designD)
    }
    ## coefficients for D precede all others
    D.modelled <- (!object$CL || object$details$relativeD) && is.null(object$fixed$D)
    # NULL happens when no density beta (relativeD with D~1)
    if (D.modelled && !is.null(Dnames)) {
        betanames <- c(paste('D', Dnames, sep='.'), betanames)
    }
    NE.modelled <- secr_NEmodelled(object$details, object$fixed, names(object$designNE))
    if (any(NE.modelled)) {
        exnames <- mapply(paste, names(object$designNE), sapply(object$designNE, colnames), sep='.')
        exnames <- unname(unlist(exnames))
        betanames <- c(betanames, exnames)
    }
    
    betanames <- sub('..(Intercept))','',betanames)
    betanames
}

#-------------------------------------------------------------------------------

secr_complete.beta <- function (object) {
    fb <- object$details$fixedbeta
    if (inherits(object, 'secr')) 
        beta <- setNames(object$fit$par, object$betanames) 
    else 
        beta <- object$beta
    secr_fullbeta(beta, fb)
}
#-------------------------------------------------------------------------------

secr_complete.beta.vcv <- function (object) {
    fb <- object$details$fixedbeta
    if (!is.null(fb) && !is.null(object$beta.vcv)) {
        nbeta <- length(fb)
        beta.vcv <- matrix(NA, nrow = nbeta, ncol = nbeta)
        beta.vcv[is.na(fb[row(beta.vcv)]) & is.na(fb[col(beta.vcv)])] <- object$beta.vcv
    }
    else {
        beta.vcv <- object$beta.vcv
    }
    beta.vcv
}
#-------------------------------------------------------------------------------

secr_smooths <- function (formula) {
    ## which terms in formula are smooths?
    ## returns logical vector
    labels <- attr(terms(formula), "term.labels")
    if (length(labels) > 0)
        sapply(labels, function (x) any(sapply(c("s\\(", "te\\(", "poly\\("), grepl, x)))
    else
        logical(0)
}

#-------------------------------------------------------------------------------

secr_gamsetup <- function(formula, data, ...) {
    ## use 'session' column as dummy LHS so gam does not gag
    ## (cf secrgam:::make.density.design.matrix)
    ## session is always present in detection data, must be added for D
    if (is.null(data$session)) data$session <- rep(1,nrow(data))
    formula <- update.formula(formula, session ~ .)
    setup <- gam(formula, data = data, fit = FALSE, ...)
    colnames(setup$X) <- setup$term.names
    setup
}
#-------------------------------------------------------------------------------

secr_general.model.matrix <- function (formula, data, gamsmth = NULL, 
    contrasts = NULL, ...) {

    ## A function to compute the design matrix for the model in
    ## 'formula' given the data in 'data'. This is merely the result
    ## of model.matrix() unless 'formula' includes smooth terms -- s()
    ## or te() as described in mgcv ?formula.gam.

    ## If smooth terms are present then the matrix may be based on a
    ## previous gam setup (provided in the argument 'gamsmth') or
    ## computed de novo with gam(..., fit = FALSE)

    ## note 2014-08-24
    ## orthogonal polynomials e.g. poly(x,2) are handled by model.matrix,
    ## but the information needed for prediction at new data is not
    ## saved by secr.fit, so predict.secr generally fails with message
    ## "'degree' must be less than number of unique points"

    ##  head(eval(parse(text = attr(terms(~ poly(x,y, degree=2)),
    ##  'term.labels')[1]), env=possummask))

    ## 2014-08-24, 2014-09-09, 2017-11-30
    ## 2019-10-12 drop row names
    ## 2021-12-09 f optional argument
    
    #--------------------------------------------------------
    polys <- function (formula) {
        ## which terms in formula are orthogonal polynomials?
        ## returns logical vector
        labels <- attr(terms(formula), "term.labels")
        if (length(labels) > 0)
            sapply(labels, grepl, pattern = "poly\\(")
        else
            logical(0)
    }
    #--------------------------------------------------------
    
    dots <- list(...)

    if (any(polys(formula)))
        stop ("orthogonal polynomials are temporarily blocked")  ## 2014-09-12
    if (any(secr_smooths(formula))) {
        if (is.null(gamsmth)) {
            ## setup knots etc from scratch
            mat <- secr_gamsetup(formula, data, ...)$X
        }
        else {
            ## fool predict.gam into generating the necessary
            ## predictor matrix from previous setup
            class (gamsmth) <- 'gam'
            gamsmth$coefficients <- rep(NA, ncol(gamsmth$X))
            mat <- mgcv::predict.gam(gamsmth, newdata = data, type = 'lpmatrix')
            colnames(mat) <- colnames(gamsmth$X)
        }
    }
    else {
        ## model.matrix(formula, data, ...)
        mat <- model.matrix(formula, data = data, contrasts.arg = contrasts)
    }
    rownames (mat) <- NULL
    mat
}

#-------------------------------------------------------------------------------

secr_makerealparameters <- function (design, beta, parindx, link, fixed) {
    modelfn <- function(i) {
        ## linear predictor for real parameter i
        Yp <- design$designMatrices[[i]] %*% beta[parindx[[i]]]
        if (names(link)[i] == 'pmix') {
            ## 2013-04-14 index of class groups (pmix sum to 1.0 within latentmodel)
            cols <- dimnames(design$designMatrices[[i]])[[2]]
            h2 <- grep('.h2', cols, fixed=T)
            h3 <- grep('.h3', cols, fixed=T)
            h2c <- grep(':h2', cols, fixed=T)
            h3c <- grep(':h3', cols, fixed=T)
            h.cols <- c(h2,h3,h2c,h3c)
            tmp <- design$designMatrices[[i]][,-h.cols, drop = FALSE]
            tmph <- design$designMatrices[[i]][,h.cols, drop = FALSE]
            ## 2018-02-23 why as.numeric()? 
            latentmodel <- as.numeric(factor(apply(tmp,1,paste, collapse='')))
            refclass <- apply(tmph,1,sum) == 0
            Yp[refclass] <- NA
            Yp <- mlogit.untransform(Yp, latentmodel)
            Yp[design$parameterTable[,i]]
        }
        else {
            Yp <- untransform(Yp, link[[i]])
            Yp[design$parameterTable[,i]]   ## replicate as required
        }
    }
    ## construct matrix of detection parameters
    nrealpar  <- length(design$designMatrices)
    nondetect <- .localstuff$spatialparametersD
    for (i in nondetect) {
        parindx[[i]] <- NULL ## detection parameters only
        link[[i]]    <- NULL ## detection parameters only
    }
    detectionparameters <- names(link)
    fixed.dp <- fixed[detectionparameters[detectionparameters %in% names(fixed)]]
    
    if (length(fixed.dp)>0)
        for (a in names(fixed.dp))  ## bug fixed by adding this line 2011-09-28
            link[[a]] <- NULL
    if (length(link) != nrealpar)
        stop ("number of links does not match design matrices")
    
    if (nrealpar == 0) {
        return(matrix(unlist(fixed.dp),nrow = 1))
    }
    
    temp <- sapply (1:nrealpar, modelfn)
    if (nrow(design$parameterTable)==1) temp <- t(temp)
    nrw <- nrow(temp)
    ## make new matrix and insert columns in right place
    temp2 <- as.data.frame(matrix(nrow = nrw, ncol = length(detectionparameters)))
    names(temp2) <- detectionparameters
    temp2[ , names(design$designMatrices)] <- temp          ## modelled
    if (!is.null(fixed.dp) & length(fixed.dp)>0)
        temp2[ , names(fixed.dp)] <- sapply(fixed.dp, rep, nrw)    ## fixed
    as.matrix(temp2)
    
}

#-------------------------------------------------------------------------------

secr_lpredictor <- function (formula, newdata, indx, beta, field, beta.vcv=NULL,
    smoothsetup = NULL, contrasts = NULL, Dfn = NULL) {
    ## form linear predictor for a single 'real' parameter
    ## smoothsetup should be provided whenever newdata differs from
    ## data used to fit model and the model includes smooths from gam
    vars <- all.vars(formula)
    OK <- vars %in% names(newdata)
    if (any(!OK)) {
        missingvars <- paste(vars[!OK], collapse = ', ')
        if (sum(!OK) == 1)
            stop ("model covariate ", missingvars, " not found in 'newdata'")
        else
            stop ("model covariates ", missingvars, " not found in 'newdata'")
    }
    newdata <- as.data.frame(newdata)
    lpred <- matrix(ncol = 2, nrow = nrow(newdata), dimnames = list(NULL,c('estimate','se')))

    if (!is.null(Dfn) && field == 'D') {
        warning("secr_lpredictor is not ready for D as function -  do not use estimates")
       nsess <- length(unique(newdata$session))
       Yp <- Dfn(newdata[,vars[1]], beta = beta[indx], dimD = c(nrow(newdata)/nsess,1,nsess)) 
       mat <- as.matrix(newdata[,vars[1], drop = FALSE])
    }
    else {
        
        mat <- secr_general.model.matrix(formula, data = newdata, gamsmth = smoothsetup, 
            contrasts = contrasts)
        if (nrow(mat) < nrow(newdata))
            warning ("missing values in predictors?", call. = FALSE)
        
        nmix <- 1
        if (field=='pmix') {
            ## drop pmix beta0 column from design matrix (always zero)
            mat <- mat[,-1,drop=FALSE]
            if ('h2' %in% names(newdata)) nmix <- 2
            if ('h3' %in% names(newdata)) nmix <- 3
            mixfield <- c('h2','h3')[nmix-1]
        }
        
        ###############################
        Yp <- mat %*% beta[indx]
        ###############################
        
        ## A latent model comprises one row for each latent class.
        ## Back transformation of pmix in mlogit.untransform() requires all rows of 
        ## each latent model. That function splits vector Yp by latent model.
        
        if (field == 'pmix') {
            nonh <- newdata[, names(newdata) != mixfield, drop = FALSE]
            latentmodel <- factor(apply(nonh, 1, paste, collapse = ''))
            refclass <- as.numeric(newdata[, mixfield]) == 1
            Yp[refclass] <- NA
            Yp <- mlogit.untransform(Yp, latentmodel)
            Yp <- logit(Yp)  # return to logit scale for later untransform!
            if (nmix==2) {
                h2.1 <- as.numeric(newdata$h2)==1
                h2.2 <- as.numeric(newdata$h2)==2
            }
        }
    }

    lpred[,1] <- Yp
    if (is.null(beta.vcv) || (any(is.na(beta[indx])))) return ( cbind(newdata,lpred) )
    else {
        if (is.null(Dfn) || field != 'D') {
            vcv <- beta.vcv[indx,indx, drop = FALSE]
            vcv[is.na(vcv)] <- 0
            nrw <- nrow(mat)
            vcv <- apply(expand.grid(1:nrw, 1:nrw), 1, function(ij)
                mat[ij[1],, drop=F] %*% vcv %*% t(mat[ij[2],, drop=F])) 
            
            vcv <- matrix (vcv, nrow = nrw)
            if (field=='pmix') {
                if (nmix==2)
                    vcv[h2.1,h2.1] <- vcv[h2.2,h2.2]
                else
                    vcv[,] <- NA
            }
            lpred[,2] <- diag(vcv)^0.5
        }
        else {
            vcv <- NULL
        }
        
        temp <- cbind(newdata,lpred)
        attr(temp, 'vcv') <- vcv
        return(temp)
    }
}

#-------------------------------------------------------------------------------

secr_getcellsize <- function (mask) {
    if (inherits(mask, 'linearmask'))
        cell <- attr(mask, 'spacing') / 1000  ## per km
    else
        cell <- attr(mask, 'area')            ## per ha
    if (is.null(cell))
        stop ("mask lacks valid cell size (area or spacing)")
    cell
}

#-------------------------------------------------------------------------------

## intercept and fix certain models with bad defaults
secr_updatemodel <- function (model, detectfn, detectfns, oldvar, newvar, warn = FALSE) {
    if (detectfn %in% detectfns) {
        for (i in 1:length(oldvar)) {
            if (oldvar[i] %in% names(model)) {
                names(model)[names(model) == oldvar[i]] <- newvar[i]
                if (warn)
                    warning ("replacing ", oldvar[i], " by ", newvar[i],
                             " in model for detectfn ", detectfn)
            }
        }
    }
    model
}

#-------------------------------------------------------------------------------


secr_nparameters <- function (object) {
    Npar <- max(unlist(object$parindx))
    Npar <- Npar + length(object$details$miscparm)
    ## allow for fixed beta parameters
    if (!is.null(object$details$fixedbeta))
        Npar <- Npar - sum(!is.na(object$details$fixedbeta))
    Npar
}

#-------------------------------------------------------------------------------

secr_mapbeta <- function (parindx0, parindx1, beta0, betaindex, default = 0)

    ## Extend beta vector from simple model (beta0) to a more complex (i.e. general)
    ## model, inserting neutral values (zero) as required.
    ## For each real parameter, a 1:1 match is assumed between
    ## beta values until all beta values from the simpler model are
    ## used up. THIS ASSUMPTION MAY NOT BE JUSTIFIED.
    ## betaindex is a user-controlled alternative.
    ## 2025-07-19 explicit default 0

{
    ## list of zeroed vectors, one per real parameter
    beta1 <- lapply(parindx1, function (x) {x[]<-default; x})
    if (is.null(beta0)) {
        unlist(beta1)
    }
    else {
        if (!is.null(betaindex)) {
            beta1 <- unlist(beta1)
            if (sum(betaindex>0) != length(beta0))
                stop ("invalid 'betaindex'")
            beta1[betaindex] <- beta0
            beta1
        }
        else {
            ## indx is within-parameter rather than absolute index
            ## for each _original_ real parameter
            indx <- lapply(parindx0, function(x) x-x[1]+1)
            was <- function (parname) {
                parname %in% names(beta0) && !parname %in% names(beta1) 
            }
            wasnt <- function (parname) {
                !parname %in% names(beta0)
            }
            for (j in names(beta1)) {
                
                # new = old
                if (j %in% names(beta0)) {
                    beta1[[j]][indx[[j]]] <- beta0[parindx0[[j]]]
                }
                
                # transfers between overlapping parameters 2025-07-27
                
                # xy > base
                if (j == 'sigma' && was('sigmaxy')) {
                    beta1[[j]][1] <- beta0[parindx0[['sigmaxy']][1]]
                }
                if (j == 'lambda0' && was ('lambda0xy')) {
                    beta1[[j]][1] <- beta0[parindx0[['lambda0xy']][1]]
                }
                if (j == 'sigma' && was('sigmakxy')) {
                    beta1[[j]][1] <- beta0[parindx0[['sigmakxy']][1]] +
                        log(100) - beta0[parindx0[['D']][1]] / 2
                }
                if (j == 'lambda0' && was ('a0xy')) {
                    beta1[[j]][1] <- beta0[parindx0[['a0xy']][1]] - 
                        log(2*pi) - 2 * beta0[parindx0[['sigma']][1]]
                }
                
                # base > xy
                if (j == 'sigmaxy' && wasnt('sigmaxy')) {
                    beta1[[j]][1] <- beta0[parindx0[['sigma']][1]]
                    beta1[['sigma']] <- 0
                }
                if (j == 'lambda0xy' && wasnt ('lambda0xy')) {
                    beta1[[j]][1] <- beta0[parindx0[['lambda0']][1]]
                    beta1[['lambda0']] <- 0
                }
                if (j == 'sigmakxy' && wasnt('sigmakxy')) {
                    beta1[[j]][1] <- beta0[parindx0[['sigma']][1]] -
                        log(100) + beta0[parindx0[['D']][1]] / 2
                    beta1[['sigma']] <- 0
                }
                if (j == 'a0xy' && wasnt('a0xy')) {
                    beta1[[j]][1] <- beta0[parindx0[['lambda0']][1]] +
                        log(2*pi) + 2 * beta0[parindx0[['sigma']][1]]
                    beta1[['lambda0']] <- 0
                }
            }
            unlist(beta1)
        }
    }
}

#-------------------------------------------------------------------------------

secr_xyinpoly <- function (xy, trps) {
    ptinside <- function (i,k) {
        ## is point i inside poly k?
        polyxy <- as.matrix(lxy[[k]])
        polyxy <- rbind(polyxy, polyxy[1,])   ## close 2014-08-28
        nr <- nrow(polyxy)
        temp <- insidecpp(unlist(xy[i,]), 0, nr-1, as.matrix(polyxy))
    }
    lxy <- split (trps, polyID(trps))
    firstinside <- function (i) {
        frstk <- 0
        for (k in 1:length(lxy)) {
            if (ptinside(i,k)) {
                frstk <- k
                break
            }
        }
        frstk
    }
    sapply(1:nrow(xy), firstinside)
}

#-------------------------------------------------------------------------------

## including pre-marked animals never sighted
## cov is optional dataframe of covariates
secr_addzeroCH <- function (CH, nzero, cov = NULL, prefix = 'Z') {
    if (nzero == 0)
        return(CH)
    else {
        nc <- nrow(CH)
        chdim <- dim(CH)
        chdim[1] <- nzero
        extra <- array(0, dim=chdim)
        dimnames(extra) <- c(list(paste(prefix, 1:nzero, sep='')), dimnames(CH)[2:3])
        CH2 <- abind(CH, extra, along = 1)
        class(CH2) <- 'capthist'
        traps(CH2) <- traps(CH)
        xy(CH2) <- xy(CH)  ## order is not affected by adding zero histories
        # added ncol>0 check 2024-09-04
        if (!is.null(covariates(CH)) && nrow(covariates(CH))>0 && 
            ncol(covariates(CH))>0 && (nrow(CH)>0)) {
            if (is.null(cov)) {
                cov <- covariates(CH)[rep(1,nzero),,drop = FALSE]
                cov[,] <- NA   ## covariates are unknown
            }
            covariates(CH2) <- rbind(covariates(CH), cov[1:nzero,,drop = FALSE])
        }
        ## ... and other essential attributes?
        CH2
    }
}

#-------------------------------------------------------------------------------

secr_expandbinomN <- function (binomN, detectorcodes) {
    # assumes detectorcodes is a vector of length = noccasions
    binomN <- ifelse (detectorcodes %in% c(2,6,7), binomN, 1)
    if (any(is.na(binomN))) stop ("NA value in binomN")
    binomN
}

#-------------------------------------------------------------------------------

secr_check3D <- function (object) {
    
    if (ms(object)) {
        out <- lapply(object, secr_check3D)
        class(out) <- class(object)
        out
    }
    else {
        if (is.matrix(object)) {
            warning("secr >= 3.0 requires 3-D capthist; ",
                    "using updateCH() to convert", call. = FALSE)
            updateCH(object)
        }
        else {
            object
        }
    }
}

#-------------------------------------------------------------------------------

secr_allzero <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, secr_allzero)
    }
    else {
        telemocc <- detector(traps(object))=='telemetry'
        # abs() applied 2025-01-07
        apply(abs(object[,!telemocc,,drop=FALSE]),1,sum)==0
    }
}

#-------------------------------------------------------------------------------

secr_selectCHsession <- function(capthist, sessnum) {
    if (ms(capthist)) 
        capthist[[sessnum]]
    else 
        capthist
}

#-------------------------------------------------------------------------------

secr_stringsAsFactors <- function (DF) {
    # convert any character columns of a data.frame (or list) to factor
    if (is.list(DF) && length(DF)>0) {    ## bug fix 2020-08-14
        chr <- sapply(DF, is.character)
        DF[chr] <- lapply(DF[chr], as.factor)
    }
    DF
}

#-------------------------------------------------------------------------------

secr_getdistmat2 <- function (traps, mask, userdist) {
    ## Static distance matrix
    if (is.function(userdist)) {
        NULL   ## compute dynamically later
    }
    else {
        if (is.matrix(userdist)) {
            if (nrow(userdist) != nrow(traps) || 
                ncol(userdist) != nrow(mask))
                stop("dimensions of userdist matrix should be c(nrow(traps), nrow(mask))")
            userdist
            
        }
        else if (any(detector(traps) %in% .localstuff$polydetectors)) {
            ## do not use result if detector is one of
            ## polygonX, polygon, transectX, transect, OR telemetry?
            matrix(0, nrow = nrow(traps), ncol = nrow(mask))
        }
        else {
            # Euclidean distance
            edist2cpp(as.matrix(traps), as.matrix(mask))
        }
    }
}

#-------------------------------------------------------------------------------

## function to assign all-ones usage matrix
secr_uniformusage <- function(object, noccasions) {
  if (inherits(object, 'capthist')) {
    if (ms(object)) {
      for (r in 1:length(object)) {
        ndet <- dim(object[[r]])[3]
        noccasions <- dim(object[[r]])[2]
        usage(traps(object[[r]])) <- matrix(1, ndet, noccasions)
      }
    }
    else {
      ndet <- dim(object)[3]
      noccasions <- dim(object)[2]
      usage(traps(object)) <- matrix(1, ndet, noccasions)
    }
  }
  else if (inherits(object, 'traps')) {
    if (missing(noccasions)) {
      stop ('noccasions should be specified for traps input')
    }
    if (ms(object)) {
      for (r in 1:length(object)) {
        ndet <- secr_ndetector(object[[r]])
        usage(object[[r]]) <- matrix(1, ndet, noccasions)
      }
    }
    else {
      ndet <- secr_ndetector(object)
      usage(object) <- matrix(1, ndet, noccasions)
    }
  }
  object
}

#-------------------------------------------------------------------------------

secr_saveprogress <- function (beta, loglik, filename) {
    log <- data.frame(
        eval = .localstuff$iter,
        loglik = loglik,
        time = format(Sys.time(), "%H:%M:%S %d %b %Y"))
    names(beta) <- names(.localstuff$savedinputs$start)
    log <- cbind(log, as.list(beta))
    attr(.localstuff$savedinputs, 'log') <- rbind(attr(.localstuff$savedinputs, 'log'), log)
    saveRDS(.localstuff$savedinputs, file = filename)
}
#-------------------------------------------------------------------------------

secr_sigmaxydistfn <- function (xy1, xy2, mask) {
    if (missing(xy1)) return("sigmaxy")
    sig <- covariates(mask)$sigmaxy   # sigma(x,y) at mask points
    sig <- matrix(sig, byrow = TRUE, nrow = nrow(xy1), ncol = nrow(xy2))
    euc <- edist(xy1, xy2) 
    euc / sig
}
#-------------------------------------------------------------------------------

secr_Dsigmakxydistfn <- function (xy1, xy2, mask) {
    if (missing(xy1)) return(c("D", "sigmakxy"))
    D   <- covariates(mask)$D   # D(x,y) at mask points
    sigk <- covariates(mask)$sigmakxy   # sigma(x,y) at mask points
    sig <- 100 * sigk / sqrt(D)
    sig <- matrix(sig, byrow = TRUE, nrow = nrow(xy1), ncol = nrow(xy2))
    euc <- edist(xy1, xy2) 
    euc / sig
}
#-------------------------------------------------------------------------------

secr_Dsigmakxya0xydistfn <- function (xy1, xy2, mask) {
    if (missing(xy1)) return(c("D", "sigmakxy","a0xy"))
    D   <- covariates(mask)$D   # D(x,y) at mask points
    sigk <- covariates(mask)$sigmakxy   # sigma(x,y) at mask points
    a0 <- covariates(mask)$a0xy * 10000 # a0(x,y) at mask points, sq. m
    sig <- 100 * sigk / sqrt(D)
    sig <- matrix(sig, byrow = TRUE, nrow = nrow(xy1), ncol = nrow(xy2))
    a0 <- matrix(a0, byrow = TRUE, nrow = nrow(xy1), ncol = nrow(xy2))
    euc <- edist(xy1, xy2) 
    
    detectfn <- attr(mask, 'detectfn')
    if (detectfn == 14)
        sqrt((euc / sig)^2 - 2 * log(a0/(2 * pi * sig^2)))   # HHN
    else if (detectfn == 16)
        euc / sig - log(a0/(2 * pi * sig^2))   # HEX
    else 
        stop ("detectfn not implemented for secr_sigmakxya0xydistfn")
}
#-------------------------------------------------------------------------------

# speculative alternative for full spatial model of lambda0 
secr_siglamxydistfn <- function (xy1, xy2, mask) {
    if (missing(xy1)) return(c("sigmaxy","lambda0xy"))
    sig <- covariates(mask)$sigmaxy   # sigma(x,y) at mask points
    sig <- matrix(sig, byrow = TRUE, nrow = nrow(xy1), ncol = nrow(xy2))
    lam <- covariates(mask)$lambda0xy   # lambda0(x,y) at mask points
    lam <- matrix(lam, byrow = TRUE, nrow = nrow(xy1), ncol = nrow(xy2))
    euc <- edist(xy1, xy2) 
    detectfn <- attr(mask, 'detectfn')
    if (detectfn == 14)
        sqrt((euc / sig)^2 - 2 * log(lam))   # HHN
    else if (detectfn == 16)
        euc / sig - log(lam)   # HEX
    else 
        stop ("detectfn not implemented for secr_siglamxydistfn")
}
#-------------------------------------------------------------------------------

secr_siga0xydistfn <- function (xy1, xy2, mask) {
    if (missing(xy1)) return(c("sigmaxy","a0xy"))
    sig <- covariates(mask)$sigmaxy   # sigma(x,y) at mask points
    sig <- matrix(sig, byrow = TRUE, nrow = nrow(xy1), ncol = nrow(xy2))
    a0 <- covariates(mask)$a0xy * 10000      # a0(x,y) at mask points, sq. m
    a0 <- matrix(a0, byrow = TRUE, nrow = nrow(xy1), ncol = nrow(xy2))
    euc <- edist(xy1, xy2) 
    detectfn <- attr(mask, 'detectfn')
    if (detectfn == 14)
        sqrt((euc / sig)^2 - 2 * log(a0/(2 * pi * sig^2)))   # HHN
    else if (detectfn == 16)
        euc / sig - log(a0/(2 * pi * sig^2))   # HEX
    else 
        stop ("detectfn not implemented for secr_siga0xydistfn")
}
#-------------------------------------------------------------------------------

secr_lambda0xydistfn <- function (xy1, xy2, mask, scale = 1) {
    if (missing(xy1)) return("lambda0xy")
    lam <- covariates(mask)$lambda0xy   # lambda0(x,y) at mask points
    lam <- matrix(lam, byrow = TRUE, nrow = nrow(xy1), ncol = nrow(xy2))
    euc <- edist(xy1, xy2)
    euc - scale * log(lam)       # HEX
}
#-------------------------------------------------------------------------------

secr_setfixedbeta <- function (fb, parindx, link, CL) {
    if (is.null(fb)) 
        fb <- rep(NA, max(unlist(parindx)))
    if (CL && !is.null(parindx$D)) {     # details$relativeD
        if (!(link$D %in% c('log','identity')))
            warning ("density link ", link$D, " not implemented for relativeD")
        if (!is.na(fb[parindx$D[1]]))
            warning ("overriding provided fixedbeta[1] for D")
        fb[parindx$D[1]] <- if (link$D == 'log') 0 else 1
    }
    if (!is.null(parindx$sigmaxy) || !is.null(parindx$sigmakxy)) {
        if (!(link$sigma %in% c('log')))
            warning ("sigma link should be log")
        fb[parindx$sigma[1]] <- 0
    }
    if (!is.null(parindx$lambda0xy)) {
        if (!(link$lambda0xy %in% c('log')))
            warning ("lambda0 link should be log")
        fb[parindx$lambda0[1]] <- 0
    }
    if (!is.null(parindx$a0xy)) {
        if (!(link$a0xy %in% c('log')))
            warning ("a0 link should be log")
        fb[parindx$lambda0[1]] <- 0   # because a0 is not in model
    }
    fb
}
#-------------------------------------------------------------------------------

# 2025-07-23 moved from preparedata.R
secr_getk <- function(traps) {
    if (!all(detector(traps) %in% .localstuff$polydetectors)) {
        nrow(traps)
    }
    else {
        if (all(detector(traps) %in% c('polygon','polygonX'))) {        
            k <- table(polyID(traps))       
        }
        else if (all(detector(traps) %in% c('transect','transectX'))) {        
            k <- table(transectID(traps))   # transectX, transect
        }
        else stop ("unrecognised poly detector type")
        c(k,0) ## zero terminate
    }
}
#--------------------------------------------------------------------------------

# 2025-07-23 moved from preparedata.R
secr_getxy <- function(dettype, capthist) {
    if (all(detector(traps(capthist)) %in% .localstuff$polydetectors)) {
        xy <- xy(capthist)
        ## start[z] indexes the first row in xy (or element in signal)
        ## for each possible count z (including zeros), where z is w-order (isk) 
        start <- abs(capthist)
        start <- head(cumsum(c(0,start)),length(start))
    }
    else {
        if (any(dettype == 13)) {
            ## ensure order matches
            ## should have null histories in capthist
            telem <- telemetryxy(capthist)
            ord <- match(names(telem), rownames(capthist), nomatch = 0)
            newtelem <- vector('list',nrow(capthist))
            newtelem[ord] <- telem
            xy <- do.call(rbind, newtelem)
            tmp <- sapply(newtelem, nrow)
            tmp[sapply(tmp, is.null)] <- 0
            start <- cumsum(c(0,tmp))
        }
        else {
            xy <- 0
            start <- 0
        }
    }
    list(xy = xy, start = start)
}
#--------------------------------------------------------------------------------

secr_makeNElist <- function (object, mask, group, sessnum) {
## 2025-07-24 construct list of design matrices for non-euclidean parameters
## object is a previously fitted secr model
## used in esa() etc.
    
    if (length(object$designNE) == 0) return(NULL)
    else if (length(object$designNE) > 3) {
        # assume legacy fit, designNE not a list
        designNE <- list(noneuc = designNE)
    }
    NElist <- secr_predictD(object = object, regionmask = mask, 
                            group = NULL, session = sessnum,
                            parameter = names(object$designNE), aslist = TRUE)
    ## convert vectors to 3-D arrays for historical reasons
    sapply(NElist, function(x) array(x, dim = c(length(x),1,1)), 
           simplify = FALSE, USE.NAMES = TRUE)
    
}
#-------------------------------------------------------------------------------

secr_telemcode <- function(object, ...) {
    if (inherits(object, 'traps') && !ms(object))
        switch (telemetrytype(object), none = 0, 
                independent = 1, dependent = 2, concurrent = 3, 0)
    else 
        NA
}

#-------------------------------------------------------------------------------

secr_ndetector <- function (traps) {
    if (is.null(traps))
        return(1)
    else if (all(detector(traps) %in% .localstuff$polydetectors))
        length(levels(polyID(traps)))
    else
        nrow(traps)
}

#-------------------------------------------------------------------------------

secr_pad1 <- function (x, n) {
    ## pad x to length n with dummy (first value)
    if (is.factor(x)) {
        xc <- as.character(x)
        xNA <- c(xc, rep(xc[1], n-length(xc)))
        out <- factor(xNA, levels=levels(x))
    }
    else out <- c(x, rep(x[1], n-length(x)))
    out
}

#-------------------------------------------------------------------------------

secr_primarysessions <- function(intervals) {
    primarysession <- cumsum(c(0,intervals))
    match(primarysession, unique(primarysession))
}

#-------------------------------------------------------------------------------

secr_secondarysessions <- function(intervals) {
    primary <- secr_primarysessions(intervals)
    unname(unlist(sapply(table(primary), seq_len)))  
}

#-------------------------------------------------------------------------------

## return indices of first occasion and detector for which PIAx is non-zero 
secr_firstsk <- function (PIAx) {
    ## PIAx dim n,s,k
    wh <- function(d2) {
        match(TRUE, d2>0)
    }
    apply(PIAx,1,wh)
}

#-------------------------------------------------------------------------------

secr_maskboolean <- function (ch, mask, threshold) {
    if (ms(ch)) {
        if (!ms(mask)) stop ("masklookup: multisession ch requires multisession mask")
        outlist <- mapply(secr_maskboolean, ch, mask, MoreArgs = list(threshold = threshold), SIMPLIFY = FALSE)
        outlist
    }
    else {
        id <- animalID(ch, names = FALSE, sortorder = 'snk')
        tr <- trap(ch, names = FALSE, sortorder = 'snk')
        trps <- traps(ch)
        m <- nrow(mask)
        if (!is.null(threshold) && all(detector(trps) %in% .localstuff$pointdetectors)) {
            df <- data.frame(id = id, x = trps$x[tr], y = trps$y[tr])
            x <- tapply(df$x, df$id, mean, na.rm=T)
            y <- tapply(df$y, df$id, mean, na.rm=T)
            xy <- data.frame(x=x,y=y)
            d2 <- edist2cpp(as.matrix(xy), as.matrix(mask))
            out <- (d2 <= threshold^2)
        }
        else {
            ## NULL option
            out <- matrix(TRUE, nrow = nrow(ch), ncol = m)
        }
        out
    }
}

#-------------------------------------------------------------------------------

# captinhood and filterw are used for novel Br behavioural response 2025-06-18

captinhood <- function (CH, maxd = NULL) {
    oneisk  <- function (isk) {
        # which detectors are in neighbourhood?
        nhood <- adj[[isk[3]]]  
        # did any detect animal i on occasions s?
        out <- CH[cbind(isk[1], isk[2], nhood)]
        # if (inherits(out, 'try-error')) browser() else
        any(out>0)
    }
    if (ms(CH)) {
        stop (" captinhood not ready for multi-session capthist")
    }
    else {
        tr <- traps(CH)
        dmat <- as.matrix(dist(tr))
        if (is.null(maxd)) maxd <- 1.5 * spacing(tr)  # kings's move on square grid
        dmat[dmat>maxd] <- 0
        if (!requireNamespace("igraph", quietly = TRUE)) {
            stop ("captinhood requires package igraph")
        }
        g <- igraph::graph_from_adjacency_matrix(dmat, weighted = TRUE, mode = "undirected")
        adj <- igraph::adjacent_vertices(g, 1:nrow(tr))
        # include focal detector in each list of adjacent indices
        adj <- mapply(c, 1:nrow(tr), adj)
        dimCH <- dim(CH)
        linearIndices <- 1:length(CH)
        # each row of isk has indices of one cell in CH
        isk <- arrayInd(linearIndices, dimCH)
        ch2 <- array(apply(isk,1,oneisk), dim = dimCH)
        dimnames(ch2)[[1]] <- dimnames(CH)[[1]]
        # artificially define as a capthist object for plotting
        class(ch2) <- 'capthist'
        traps(ch2) <- tr
        ch2
    }
}
#-------------------------------------------------------------------------------

# exponentially weight past, with padding
filterw <- function (x, w = 5, lambda = 0.6) {
    if (lambda==1)
        weights <- rep(1,w)
    else
        weights <- (1 - lambda) * lambda^(0:(w - 1))
    weights <- weights / sum(weights)
    xp <- c(rep(0,w), x)
    filter(xp, filter = weights, sides = 1)[-(1:w)]
}

#-------------------------------------------------------------------------------


## miscellaneous functions

invlogit <- function (y) 1/(1+exp(-y))   # plogis(y)
logit    <- function (x) log(x/(1-x))    # qlogis(x), except for invalid argument
sine     <- function (x) asin (x*2-1)
invsine  <- function (y) (sin(y)+1) / 2
odds     <- function (x) x / (1-x)
invodds  <- function (y) y / (1+y)

#-------------------------------------------------------------------------------

## Detection functions

HN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]
    g0 * exp (-r^2 / 2 / sigma^2)
}
HR <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    g0 * (1 - exp (-(r / sigma)^-z))
}
EX <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]
    g0 * exp (-r / sigma)
}
UN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]
    ifelse (r<=sigma, g0, 0)
}
CHN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    g0 * ( 1 - (1 - exp (-r^2 / 2 / sigma^2)) ^ z )
}
WEX <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
    ifelse(r<=w, g0, g0*exp (-(r-w) / sigma))
}
ANN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
    g0 * exp (-(r-w)^2 / 2 / sigma^2)
}
CLN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    CV2 <- (z/sigma)^2
    sdlog <- log(1 + CV2)^0.5
    meanlog <- log(sigma) - sdlog^2/2
    g0 * plnorm(r, meanlog, sdlog, lower.tail = FALSE)
}
CG <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    g0 * pgamma(r, shape=z, scale=sigma/z, lower.tail = FALSE)
}
CN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    x <- z * (r - sigma)
    g0 * (1 + (1 - exp(x)) / (1 + exp(x)))/2
}
BSS <- function (r, pars, cutval) {
    b0 <- pars[1]; b1 <- pars[2]
    gam <- -(b0 + b1 * r);
    pnorm (gam, mean=0, sd=1, lower.tail=FALSE)
}
SS <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
    if (is.null(cutval))
        stop ("require 'details$cutval' for signal strength plot")
    mu <- beta0 + beta1 * r
    1 - pnorm (q=cutval, mean=mu, sd=sdS)
}
SSS <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
    if (is.null(cutval))
        stop ("require 'details$cutval' for signal strength plot")
    ## spherical so assume distance r measured from 1 m
    mu <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
    mu[r<1] <- beta0
    1 - pnorm (q=cutval, mean=mu, sd=sdS)
}
SN <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3];
    muN <- pars[4]; sdN <- pars[5]
    muS <- beta0 + beta1 * r
    1 - pnorm (q=cutval, mean=muS-muN, sd=sqrt(sdS^2+sdN^2))
}
SNS <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3];
    muN <- pars[4]; sdN <- pars[5]
    ## spherical so assume distance r measured from 1 m
    muS <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
    muS[r<1] <- beta0
    1 - pnorm (q=cutval, mean=muS-muN, sd=sqrt(sdS^2+sdN^2))
}
HHN <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]
    1 - exp(-lambda0 * exp (-r^2 / 2 / sigma^2))
}
HHR <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    1 - exp(-lambda0 * ( 1 - exp (-(r / sigma)^-z)))
}
HEX <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]
    1 - exp(-lambda0 * exp (-r / sigma))
}
HAN <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
    1 - exp(-lambda0 * exp (-(r-w)^2 / 2 / sigma^2))
}
HCG <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    lambda0 * pgamma(r, shape=z, scale=sigma/z, lower.tail = FALSE)
}
HVP <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    1 - exp(-lambda0 * exp(-(r/sigma)^z))
}

#-------------------------------------------------------------------------------

# transformation tidy up 2021-12-16
# arbitrary link function specified with functions X, invX, se.invX

transform <- function (x, link) {
    switch (link,
            identity = x,
            i1000 = x * 1000,
            log = log(x),
            neglog = log(-x),
            logit = logit(x),
            odds = odds(x),
            sin = sine(x),
            do.call(link, list(x))
    )
}
#-------------------------------------------------------------------------------

# used only in model.average, modelAverage
se.transform <- function (real, sereal, link) {
    switch (link,
            identity = sereal,
            i1000 = sereal / 1000,
            log = log((sereal/real)^2 + 1)^0.5,
            neglog = log((sereal/-real)^2 + 1)^0.5,
            logit = sereal / real / (1 - real),
            sin = NA,
            do.call(paste0('se.',link), list(real, sereal) )
    )
}
#-------------------------------------------------------------------------------

untransform <- function (beta, link) {
    switch (link,
            identity = beta,
            i1000 = beta / 1000,
            log = exp(beta),
            neglog = -exp(beta),
            logit = invlogit(beta),
            odds = invodds(beta),
            sin = invsine(beta),
            do.call(paste0('inv',link), list(beta))
    )
}
#-------------------------------------------------------------------------------

se.untransform <- function (beta, sebeta, link) {
    # Approximate translation of SE to untransformed scale
    # Delta method cf Lebreton et al 1992 p 77
    switch (link,
            identity = sebeta,
            i1000 = sebeta / 1000,
            log = exp(beta) * sqrt(exp(sebeta^2)-1),
            neglog = exp(beta) * sqrt(exp(sebeta^2)-1),
            logit = invlogit(beta) * (1-invlogit(beta)) * sebeta,
            sin = NA,                ####!!!!
            do.call(paste0('se.inv', link), list(beta=beta, sebeta=sebeta))
    )
}
#-------------------------------------------------------------------------------

# vectorized transformations
Xtransform <- function (real, linkfn, varnames) {
    mapply(transform, real, linkfn[varnames])
}
se.Xtransform <- function (real, sereal, linkfn, varnames) {
    mapply(se.transform, real, sereal, linkfn[varnames])
}
Xuntransform <- function (beta, linkfn, varnames) {
    mapply(untransform, beta, linkfn[varnames])
}
se.Xuntransform <- function (beta, sebeta, linkfn, varnames)
{
    if (length(beta)!=length(sebeta))
        stop ("'beta' and 'sebeta' do not match")
    if (!all(varnames %in% names(linkfn)))
        stop ("'linkfn' component missing for at least one real variable")
    mapply(se.untransform, beta, sebeta, linkfn[varnames])
}
#-------------------------------------------------------------------------------

mlogit.untransform <- function (beta, latentmodel) {
    if (!missing(latentmodel)) {
        for (i in unique(latentmodel))
            beta[latentmodel==i] <- mlogit.untransform(beta[latentmodel==i])
        beta
    }
    else {
        ## beta should include values for all classes (mixture components)
        nmix <- length(beta)
        if (sum(is.na(beta)) != 1) {
            ## require NA for a single reference class
            rep(NA, length(beta))
        }
        else {
            nonreference <- !is.na(beta)   # not reference class
            b <- beta[nonreference]
            pmix <- numeric(nmix)
            pmix[nonreference] <- exp(b) / (1+sum(exp(b)))
            pmix[!nonreference] <- 1 - sum(pmix[nonreference])
            pmix
        }
    }
}

#-------------------------------------------------------------------------------

clean.mlogit <- function(x) {
    ## 2014-08-19 for robustness...
    if (is.na(x[2])) x[2] <- 1-x[1]
    x[1] <- NA   ## assumed reference class
    logit(mlogit.untransform(x))
}

#-------------------------------------------------------------------------------

mlogit <- function (x) {
    ## return the mlogit of an unscaled vector of positive values
    ## 2013-04-14
    logit(x/sum(x))
}

## End of miscellaneous functions

#-------------------------------------------------------------------------------


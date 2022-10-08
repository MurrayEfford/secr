#######################################################################################
## utility.R
#######################################################################################

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
#######################################################################################

# Global variables in namespace
#
## define a local environment for temporary variables e.g. iter
## e.g. Roger Peng https://stat.ethz.ch/pipermail/r-devel/2009-March/052883.html

.localstuff <- new.env()

.localstuff$packageType <- ' pre-release'
#.localstuff$packageType <- ''

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
    'hazard pixelar')

.localstuff$DFN <- c('HN', 'HR', 'EX', 'CHN', 'UN', 'WEX', 'ANN', 'CLN', 'CG',
  'BSS', 'SS', 'SSS', 'SN', 'SNS',
  'HHN', 'HHR', 'HEX', 'HAN', 'HCG', 'HVP','HPX')

.localstuff$learnedresponses <- c('b', 'bk', 'B', 'k', 'Bk')   ## Bk added 2020-02-26

detectionfunctionname <- function (fn) {
    .localstuff$detectionfunctions[fn+1]
}

detectionfunctionnumber <- function (detname) {
    dfn <- match (toupper(detname), .localstuff$DFN)
    if (is.na(dfn))
        dfn <- match (tolower(detname), .localstuff$detectionfunctions)
    if (is.na(dfn))
        stop ("unrecognised detection function ", detname)
    dfn-1
}
parnames <- function (detectfn) {
    parnames <- switch (detectfn+1,
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
        c('lambda0')    ## 20
    )
}
getdfn <- function (detectfn) {
    switch (detectfn+1, HN, HR, EX, CHN, UN, WEX, ANN, CLN, CG, BSS, SS, SSS,
                       SN, SNS, HHN, HHR, HEX, HAN, HCG, HVP, HPX)
}

valid.detectfn <- function (detectfn, valid = c(0:3,5:19, 20)) {
# exclude 4 uniform: too numerically flakey
    if (is.null(detectfn))
        stop ("requires 'detectfn'")
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)
    if (!(detectfn %in% valid))
        stop ("invalid detection function")
    detectfn
}

valid.detectpar <- function (detectpar, detectfn) {
    if (is.null(detectpar) | is.null(detectfn))
        stop ("requires valid 'detectpar' and 'detectfn'")

    ## 2013-07-19, 2013-10-22
    ## replace a0 with g0 or lambda0 as appropriate to detectfn
    if ('a0' %in% names(detectpar)) {
        aname <- if (detectfn %in% 0:8) 'g0' else 'lambda0'
        lambda0 <- detectpar[['a0']] / 2 / pi / detectpar[[2]]^2 * 10000
        detectpar[[aname]] <- if (detectfn %in% 0:8) 1-exp(-lambda0) else lambda0
    }

    if (!all(parnames(detectfn) %in% names(detectpar)))
        stop ("requires 'detectpar' ", paste(parnames(detectfn), collapse=','),
            " for ", detectionfunctionname(detectfn), " detectfn")
    detectpar[parnames(detectfn)]
}

valid.model <- function(model, CL, detectfn, hcov, userdist, sessioncovnames) {
    ## 2014-08-25
    if (any(sapply(model, badsmooths)))
        warning("smooth term may be unsuitable for secr: does not specify k or fx where required")
}

getuserdistnames <- function (userdist) {
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

valid.pnames <- function (details, CL, detectfn, alltelem, sighting, nmix) {
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
        c('beta0','beta1','sdS'),  # 12  cf parnames() in utility.R: muN, sdN?
        c('beta0','beta1','sdS'),  # 13  cf parnames() in utility.R: muN, sdN?
        c('lambda0','sigma'),      # 14 hazard halfnormal
        c('lambda0','sigma','z'),  # 15 hazard hazard rate
        c('lambda0','sigma'),      # 16 hazard exponential
        c('lambda0','sigma','w'),  # 17
        c('lambda0','sigma','z'),  # 18
        c('lambda0','sigma','z'),  # 19
        c('lambda0','sigma'))      # 20 hazard pixelar 2021-03-25    

    if (details$param %in% c(2,6))
        pnames[1] <- 'esa'
    if (details$param %in% c(3,5))
        pnames[1] <- 'a0'
    if (details$param %in% 4:6) {
        pnames[2] <- 'sigmak'
        pnames <- c(pnames, 'c')
        pnames <- c(pnames, 'd')
    }
    if (!CL)
      pnames <- c('D', pnames)
    if ('noneuc' %in% getuserdistnames(details$userdist))
      pnames <- c(pnames, 'noneuc')
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

valid.userdist <- function (userdist, detector, xy1, xy2, mask, sessnum) {
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
            OK <- getuserdistnames(userdist) %in% names(covariates(mask))
            if ((length(OK)>0) & !all(OK))
                stop ("covariates required by userdist function not in mask : ",
                      paste(getuserdistnames(userdist)[!OK], collapse=','))
            result <- do.call(userdist, c(list(xy1, xy2, mask)))
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
            warning ("replacing infinite, negative and NA userdist values with 1e10")
            result[baddist] <- 1e10
        }
    }
    result
}
#-------------------------------------------------------------------------------

new.param <- function (details, model, CL) {
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
        warning ("Using parameterization details$param = ", newparam)
    newparam
}

# secr 3.0 2016-10-05
## NEEDS MULTI-SESSION FORM
detectorcode <- function (object, MLonly = TRUE, noccasions = NULL) {
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

expanddet <- function(CH) {
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

## 2013-06-16
ndetectpar <- function (detectfn) {
    length(parnames(detectfn))
}

replacedefaults <- function (default, user) replace(default, names(user), user)

discreteN <- function (n, N) {
    tN <- trunc(N)
    if (N != tN) tN + sample (x = c(1,0), prob = c(N-tN, 1-(N-tN)),
        replace = T, size = n)
    else rep(tN,n)
}

ndetector <- function (traps) {
    if (is.null(traps))
        return(1)
    else if (all(detector(traps) %in% .localstuff$polydetectors))
        length(levels(polyID(traps)))
    else
        nrow(traps)
}

memo <- function (text, trace) {
    ## could use message(text), but does not immediately flush console
    if (trace) { cat (text, '\n')
    flush.console() }
}

insertdim <- function (x, dimx, dims) {
  ## make vector of values
  ## using x repeated so as to fill array
  ## with dim = dims and the x values occupying dimension(s) dimx
  olddim <- 1:length(dims)
  olddim <- c(olddim[dimx], olddim[-dimx])
  temp <- array (dim=c(dims[dimx], dims[-dimx]))
  tempval <- array(dim=dims[dimx])
  if (length(x) > length(tempval))
      tempval[] <- x[1:length(tempval)]
  else
      tempval[] <- x     ## repeat as needed
  temp[] <- tempval  ## repeat as needed
  if (is.factor(x))
    factor(levels(x), levels=levels(x))[aperm(temp, order(olddim))]   ## 2010 02 25
  else
    as.vector(aperm(temp, order(olddim)))
}

pad1 <- function (x, n) {
## pad x to length n with dummy (first value)
    if (is.factor(x)) {
        xc <- as.character(x)
        xNA <- c(xc, rep(xc[1], n-length(xc)))
        out <- factor(xNA, levels=levels(x))
    }
    else out <- c(x, rep(x[1], n-length(x)))
    out
}

padarray <- function (x, dims) {
    temp <- array(dim=dims)
    dimx <- dim(x)
    if (all(dimx>0)) {
        if (length(dimx)<2 | length(dimx)>3)
            stop ("invalid array")
        if (length(dimx)>2) temp[1:dimx[1], 1:dimx[2], 1:dimx[3]] <- x
        else temp[1:dimx[1], 1:dimx[2]] <- x
    }
    temp
}

## regularize a list of formulae
stdform <- function (flist) {
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

## Start of miscellaneous functions

invlogit <- function (y) 1/(1+exp(-y))   # plogis(y)
logit    <- function (x) log(x/(1-x))    # qlogis(x), except for invalid argument
sine     <- function (x) asin (x*2-1)
invsine  <- function (y) (sin(y)+1) / 2
odds     <- function (x) x / (1-x)
invodds  <- function (y) y / (1+y)

lnbinomial <- function (x,size,prob) {
  lgamma (size+1) - lgamma (size-x+1) - lgamma (x+1) +
      x * log(prob) + (size-x) * log (1-prob)
}
############################################################################################
## moved from methods.r 2012-10-28

model.string <- function (model, userDfn) {
    if (!is.null(userDfn)) {
        if (!is.null(model$D))
            model$D <- paste('~userD', userDfn('name'), sep='.')
    }
    temp <- paste (names(model), as.character(model), collapse=' ', sep='')
    temp
}
fixed.string <- function (fixed) {
    if (is.null(fixed) | length(fixed)==0) 'none'
    else paste (names(fixed), as.character(fixed), collapse=', ', sep=' = ')
}
############################################################################################

var.in.model <- function(v,m) v %in% unlist(lapply(m, all.vars))

############################################################################################

get.nmix <- function (model, capthist, hcov) {
    model$D <- NULL  ## ignore density model
    model$pmix <- NULL ## pmix alone cannot make this a mixture model
    nmix <- 1
    if (any(var.in.model('h2', model))) {
        nmix <- 2
        if (any(var.in.model('h3', model)))
            stop ("do not combine h2 and h3")
    }
    if (any(var.in.model('h3', model))) {
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
            warning ("hcov covariate has more than 2 levels; using only first two")
        nmix <- 2
    }
    nmix
}
############################################################################################

fixpmix <- function(x, nmix) {

    ## x is a list with component pmix that is a matrix (dataframe)
    ## with columns 'estimate' and 'se' (and possibly others)
    ## referring to the linear predictor of pmix (i.e. on mlogit
    ## scale) and rows corresponding to rows in newdata
    ## (i.e. arbitrary combinations of predictors, including mixture
    ## class h2 or h3)

    ####################################################
    ## It is necessary that newdata include all levels
    ## of the mixture class.
    ####################################################

    ## 2013-10-29
    ## assuming mixture is always last dimension...

    ## previously used in collate, model.average and predict.secr
    ## 2015-09-30 incorporated in secr.lpredictor

    temp <- matrix(x$pmix[,'estimate'], ncol = nmix)
    if (nmix==2) temp[,x$pmix[,'h2']] <- x$pmix[,'estimate']
    if (nmix==3) temp[,x$pmix[,'h3']] <- x$pmix[,'estimate']
    temp2 <- apply(temp, 1, clean.mlogit)
    x$pmix[,'estimate'] <- as.numeric(t(temp2))
    if (nmix==2)
        x$pmix[as.numeric(x$pmix$h2)==1,'se'] <- x$pmix[as.numeric(x$pmix$h2)==2,'se']
    else
        x$pmix[,'se'] <- rep(NA, nrow(x$pmix))   ## don't know how
    x
}

############################################################################################

add.cl <- function (df, alpha, loginterval, lowerbound = 0) {

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

###############################################################################

spatialscale <- function (object, detectfn, session = '') {
    if (inherits(object, 'secr')) {
        if (ms(object))
            detpar <- detectpar(object)[[session]]
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
#        (0.5 - detpar$b0) / detpar$b1
        - 1 / detpar$b1   ## 2010-11-01
    }
    else stop ("unrecognised detectfn")
}

###############################################################################

###############################################################################
## logical for whether object specifies userDfn

userD <- function (object) {
  if (!inherits(object, c('secr','ipsecr')))
    stop ("requires fitted model")
  !is.null(object$details$userDfn)
}

###############################################################################

## mean and SD if x numeric
## was statfn 2011-11-08
getMeanSD <- function(xy) {
    MeanSD <- function (x) {
        if (is.numeric(x))
            c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
        else
            c(NA,NA)
    }
    as.data.frame (apply(xy, 2, MeanSD))
}
###############################################################################

nclusters <- function (capthist) {
    if (ms(capthist)) {
	lapply(capthist, nclusters)
    }
    else 	{
        nmash <- attr(capthist, 'n.mash')
        ifelse (is.null(nmash), 1, length(nmash))
    }
}
###############################################################################
## moved here from make.grid 2012-09-04

# leadingzero <- function (x) {
#     formatC(x, width=max(nchar(x)), flag='0')  ## returns character value
# }

## clunky but effective re-write 2012-09-04, improved 2016-02-20, 2016-05-10
leadingzero <- function (x) {
    xc <- as.character(x)
    w <- max(nchar(xc))
    n0 <- function(n) paste(rep('0',n), collapse='')
    paste(sapply(w-nchar(xc), n0), x, sep='')

    ## or, 2016-01-15, 2016-02-20 BUT DOESN'T HANDLE NON-INTEGER 2016-05-10
    #     if (is.character(x)) x <- as.numeric(x)
    #     sprintf(paste("%0", w, "d", sep = ""), x)
}

###############################################################################

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
HPX <- function (r, pars, cutval) {
    g0 <- 1-exp(-pars[1])
    radius <- pars[2]
    ifelse (r<=radius, g0, 0)  # circular, not square! crude approx
}

############################################################################################

gradient <- function (pars, fun, eps=0.001, ...)
## quick & dirty 2009 09 14
## used by plot.secr for delta method limits
{
  est <- pars
  g   <- pars
  for (i in 1:length(est))
  {
      temp     <- est[i]
      if (temp != 0.0) delta <- eps * abs(temp)
      else             delta <- eps
      est[i]  <- temp - delta
      fminus  <- fun (est, ...)
      est[i]  <- temp + delta
      fplus   <- fun (est, ...)
      g[i]    <- (fplus - fminus) / (2.0 * delta)
      est[i]  <- temp;
  }
  g
}
############################################################################################

distancetopoly <- function (X, traps) {
  ## X should be 2-column dataframe, mask, matrix or similar
  ## with x coord in col 1 and y coord in col 2
  
  if (is.null(X)) return (NULL) 
  
  detecttype <- detector(traps)
  detecttype <- ifelse (is.null(detecttype), "", detecttype)
  if (!all(detecttype %in% c('polygon', 'polygonX')))
    stop("distancetopoly is for polygon detectors only")
  
  xy <- st_as_sf(data.frame(X), coords=1:2)
  trps <- split(traps, polyID(traps))
  polys <- lapply(trps, boundarytoSF)

  dlist <- lapply(polys, st_distance, x = xy)
  matrix(unlist(dlist), ncol = length(dlist))
}

distancetotrap <- function (X, traps) {
    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coor in col 2
    
    if (is.null(X)) return (NULL)  ## 2022-01-04
    
    # X <- matrix(unlist(X), ncol = 2)
    # 2022-02-18
    X <- as.data.frame(X)
    xy <- st_as_sf(X, coords = 1:2)     # POINTS
    
    nxy <- nrow(X)
    detecttype <- detector(traps)
    detecttype <- ifelse (is.null(detecttype), "", detecttype)
    
    if (all(detecttype %in% c('polygon', 'polygonX'))) {
        trps <- split(traps, polyID(traps))
        polys <- lapply(trps, boundarytoSF)
        dlist <- lapply(polys, st_distance, x = xy)
        dmat <- matrix(unlist(dlist), ncol = length(dlist))
        d <- apply(dmat,1,min)
        return (d)
    }
    
    if (inherits(traps, 'SpatialPolygons')) {
        traps <- st_as_sf(traps)
        d <- st_distance(xy, traps)
        return (d)
    }
    else if (all(detecttype %in% .localstuff$polydetectors)) {
        ## approximate only
        
        traps <- split(traps, polyID(traps))
        trpi <- function (i, n = 100) {
            intrp <- function (j) {
                ## 2020-01-08 dodge issue with polyID in as.data.frame
                ## tmp <- as.data.frame(traps[[i]][j:(j+1),])[,-1]   
                tmp <- data.frame(x = traps[[i]]$x[j:(j+1)], y = traps[[i]]$y[j:(j+1)])
                if (tmp$x[1] == tmp$x[2])
                    data.frame(x=rep(tmp$x[1], n),
                        y=seq(tmp$y[1], tmp$y[2], length=n))
                else {
                    ## 2019-11-30 suppress warnings such as :
                    ## In regularize.values(x, y, ties, missing(ties)) :
                    ## collapsing to unique 'x' values
                    suppressWarnings(data.frame(approx(tmp, n = n)))
                }
            }
            tmp <- lapply(1:(nrow(traps[[i]])-1),intrp)
            do.call(rbind, tmp)
        }
        trps <- do.call(rbind, lapply(1:length(traps), trpi))
        trps <- matrix(unlist(trps), ncol = 2)
    }
    else {
        ## 2015-10-18 added protection
        trps <- matrix(unlist(traps), ncol = 2)
    }
    
    temp <- nearestcpp(as.matrix(X), as.matrix(trps))
    if (all(detecttype %in% c('polygon', 'polygonX'))) {
        inside <- lapply(traps, pointsInPolygon, xy = X)
        inside <- do.call(rbind, inside)
        temp$distance [apply(inside,2,any)] <- 0
    }
    temp$distance
}

nearesttrap <- function (X, traps) {
  ## X should be 2-column dataframe, mask, matrix or similar
  ## with x coord in col 1 and y coord in col 2
  
  if (is.null(X)) return (NULL)  ## 2022-01-04
  
  X <- matrix(unlist(X), ncol = 2)
  nxy <- nrow(X)
  if (inherits(traps, 'SpatialPolygons')) {
    stop ("nearesttrap currently does not accept SpatialPolygons (from 4.5.3)")
    # traps <- sp::coordinates(traps@polygons[[1]]@Polygons[[1]])
    # warning("using only first polygon of SpatialPolygons")
  }
  temp <- nearestcpp(as.matrix(X), as.matrix(traps))
  temp$index
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

clean.mlogit <- function(x) {
    ## 2014-08-19 for robustness...
    if (is.na(x[2])) x[2] <- 1-x[1]
    x[1] <- NA   ## assumed reference class
    logit(mlogit.untransform(x))
}

mlogit <- function (x) {
    ## return the mlogit of an unscaled vector of positive values
    ## 2013-04-14
    logit(x/sum(x))
}

## End of miscellaneous functions
############################################################################################

group.levels <- function (capthist, groups, sep='.') {
    # 2016-06-05 use also for trap strata
    if (inherits(capthist, 'list')) {
        temp <- lapply(capthist, group.levels, groups, sep)
        sort(unique(unlist(temp)))  ## vector of global levels
    }
    else {
        if (is.null(groups)) 0
        else {
            if (!all(groups %in% names(covariates(capthist))))
                stop ("one or more grouping variables is missing ",
                      "from covariates")
            temp <- as.data.frame(covariates(capthist)[,groups])
            # omit null combinations, sort as with default of factor levels
            sort(levels(interaction(temp, drop=T, sep=sep)))
        }
    }
}
############################################################################################

h.levels <- function (capthist, hcov, nmix) {
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
            warning("hcov was coerced to a factor")
            hcov <- factor(hcov)
        }
        levels(hcov)[1:nmix]
    }
}
############################################################################################

n.occasion <- function (capthist) {
## return the number of sampling occasions for each session in capthist
    if (inherits(capthist, 'list')) {
        sapply(capthist, n.occasion)
    }
    else {
        ncol(capthist)
    }
}

############################################################################################

group.factor <- function (capthist, groups, sep='.')
    ## convert a set of grouping factors to a single factor (g)
    ## levels common to all sessions
{
    if (inherits(capthist, 'list')) {
        temp <- lapply(capthist, group.factor, groups)  ## recursive call
        grouplevels <- group.levels(capthist, groups)
        if (length(grouplevels)<2)
            temp
        else
            # list; force shared factor levels on each component
            lapply (temp, factor, levels=grouplevels)
    }
    else {
        if (is.null(groups) | (length(groups)==0) )
            return (factor(rep(1, nrow(capthist)), levels = 1))  # added levels 2017-04-18
        temp <- as.data.frame(covariates(capthist)[,groups])
        if (ncol(temp) != length(groups))
            stop ("one or more grouping variables is missing from ",
                  "covariates(capthist)")
        temp <- interaction(temp, drop=T, sep=sep)  # omit null combinations
        temp
    }
}
############################################################################################

getgrpnum <- function (capthist, groups) {
    if (is.null(groups))
        rep(1, nrow(capthist))
    else
        match(group.factor(capthist, groups), group.levels(capthist, groups))
}
############################################################################################

## adapted for cpp 2017-07-27
make.lookup <- function (tempmat) {

    ## should add something to protect make.lookup from bad data...
    nrw <- nrow(tempmat)
    ncl <- ncol(tempmat)
    nam <- colnames(tempmat)

    df <- is.data.frame(tempmat)
    if (df) {
       lev <- lapply(tempmat, levels)
       tempmat[] <- sapply(tempmat, as.numeric)
       tempmat <- as.matrix(tempmat)
    }
    dimnames(tempmat) <- NULL

    temp <- makelookupcpp(tempmat)
    
    lookup <- temp$lookup
    colnames(lookup) <- nam
    if (df) {
        lookup <- as.data.frame(lookup)
        ## restore factors
        for (i in 1: length(lev))
            if (!is.null(lev[[i]]))
                lookup[,i] <- factor(lev[[i]][lookup[,i]], levels = lev[[i]])
    }
    list (lookup=lookup, index=temp$index)
}
###############################################################################

## Return an integer vector of class membership defined by a categorical
## individual covariate in a capthist object. Individuals of unknown
## class (including those with class exceeding nmix) are coded 1,
## others as (class number + 1). When no mixture is specified (nmix == 1)
## all are coded as unknown.

## knownclass 1 'unknown' 
## knownclass 2 'latent class 1' 
## knownclass 3 'latent class 2' 

getknownclass <- function(capthist, nmix, hcov) {
    if (ms(capthist)) {
        lapply(capthist, getknownclass, nmix = nmix, hcov = hcov)
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

###############################################################################

getnmix <- function (details) {
    if (is.null(details$nmix))
       1
    else
       details$nmix
}

###############################################################################

## expand beta parameter vector using template of 'fixed beta'
## fixed beta fb input is missing (NA) for estimated beta parameters
fullbeta <- function (beta, fb) {
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta  ## partial beta (varying only)
        beta <- fb             ## complete beta
    }
    beta
}
###############################################################################

## inflate a convex outline along all radii by linear factor 'rmult'
## 2013-06-15
inflate <- function (xy, rmult = 1) {
    xy <- as.matrix(xy)
    centre <- apply(xy, 2, mean)
    xy <- sweep(xy, MARGIN = 2, STATS = centre, FUN = '-')
    r <- apply(xy, 1, function(z) sqrt(sum(z^2)))
    theta <- atan2 (xy[,2], xy[,1])
    r <- r * rmult
    xy <- cbind(r * cos(theta), r * sin(theta))
    sweep(xy, MARGIN = 2, STATS = centre, FUN = '+')
}
###############################################################################
## moved from pdot.R 2013-11-09
## scalar 2016-10-14
getbinomN <- function (binomN, detectr) {
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
###############################################################################

## convert telemetryxy attribute of a combined dataset into a standalone capthist

## TO DO: option of telemetry or polygon output

xy2CH <- function (CH, inflation = 1e-8) {
    if (ms(CH)) {
        out <- lapply(CH, xy2CH, inflation)
        class(out) <- c('capthist', 'list')
        out
    }
    else {
        xylist <- telemetryxy(CH)
        if (is.null(xylist))
            stop ("requires 'telemetryxy' attribute")
        n <- length(xylist)
        neach <- sapply(xylist, nrow)
        allxy <- do.call(rbind, xylist)

        trps <-  allxy[chull(allxy),]
        trps <- rbind(trps, trps[1,,drop=F])
        trps <- inflate(trps, 1 + inflation)  ## see also telemetry.R

        trps <- as.data.frame(trps)
        dimnames(trps) <- list(1:nrow(trps), c('x','y'))
        class(trps) <- c("traps","data.frame")
        detector(trps) <- "polygon"
        polyID(trps) <- factor(rep(1,nrow(trps)))

        rown <- rep(names(xylist), neach)
        newCH <- array(neach, dim = c(n, 1, 1))
        attr(newCH, "detectedXY") <- allxy
        if (!is.null(covariates(CH))) {
            rowlookup <- match(names(xylist), rownames(CH))
            covariates(newCH) <- covariates(CH)[rowlookup,, drop=FALSE]
        }
        class(newCH) <- "capthist"
        traps(newCH) <- trps
        newCH
    }
}

###############################################################################
## moved from mask.check.r 2014-08-28

inflatechull <- function (poly, r, ntheta = 60) {
    theta <- (2*pi) * (1:ntheta) / ntheta
    ## add supernumerary vertices
    temp  <- data.frame(x = apply(expand.grid(poly$x, r * cos(theta)),1,sum),
                   y = apply(expand.grid(poly$y, r * sin(theta)),1,sum))
    hull <- chull(temp)
    temp[c(hull,hull[1]), ]
}

###############################################################################
## used by sim.capthist to update telemetry boundary polygon 2013-11-21
## OBSOLETE 2017-01-11
refreshMCP <- function (CH, tol) {
    if (all(detector(traps(CH)) %in% c('polygon','polygonX')))
        allxy <- xy(CH)
    else
        stop ("requires polygon detector type")
    trps <-  allxy[chull(allxy),]
    trps <- inflatechull(trps, r = tol)   ## 2014-08-28
    class(trps) <- c("traps","data.frame")
    names(trps) <- c('x','y')
    detector(trps) <- detector(traps(CH))
    polyID(trps) <- rep(1,nrow(trps))
    traps(CH) <- trps
    CH
}
###############################################################################

maskarea <- function (mask, sessnum = 1) {
    if (!ms(mask)) nrow(mask) * attr(mask,'area')
    else nrow(mask[[sessnum]]) * attr(mask[[sessnum]],'area')
}
###############################################################################

masklength <- function (mask, sessnum = 1) {
    if (!ms(mask)) nrow(mask) * attr(mask,'spacing')/1000
    else nrow(mask[[sessnum]]) * attr(mask[[sessnum]],'spacing')/1000
}
###############################################################################

masksize <- function (mask, sessnum = 1) {
    if (inherits(mask, 'linearmask'))
        masklength(mask, sessnum)
    else
        maskarea(mask, sessnum)
}
###############################################################################

complete.beta <- function (object) {
    fb <- object$details$fixedbeta
    # modified 2022-04-02 for consistency with ipsecr
    beta <- if (inherits(object, 'secr')) object$fit$par else object$beta
    if (!is.null(fb)) {
        nbeta <- length(fb)
        fb[is.na(fb)] <- beta
        beta <- fb
    }
    beta
}
###############################################################################

complete.beta.vcv <- function (object) {
    fb <- object$details$fixedbeta
    if (!is.null(fb)) {
        nbeta <- length(fb)
        beta.vcv <- matrix(NA, nrow = nbeta, ncol = nbeta)
        beta.vcv[is.na(fb[row(beta.vcv)]) & is.na(fb[col(beta.vcv)])] <- object$beta.vcv
    }
    else {
        beta.vcv <- object$beta.vcv
    }
    beta.vcv
}
###############################################################################

smooths <- function (formula) {
    ## which terms in formula are smooths?
    ## returns logical vector
    labels <- attr(terms(formula), "term.labels")
    if (length(labels) > 0)
        sapply(labels, function (x) any(sapply(c("s\\(", "te\\(", "poly\\("), grepl, x)))
    else
        logical(0)
}
############################################################################################

polys <- function (formula) {
    ## which terms in formula are orthogonal polynomials?
    ## returns logical vector
    labels <- attr(terms(formula), "term.labels")
    if (length(labels) > 0)
        sapply(labels, grepl, pattern = "poly\\(")
    else
        logical(0)
}
############################################################################################

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
############################################################################################

gamsetup <- function(formula, data, ...) {
    ## use 'session' column as dummy LHS so gam does not gag
    ## (cf secrgam:::make.density.design.matrix)
    ## session is always present in detection data, must be added for D
    if (is.null(data$session)) data$session <- rep(1,nrow(data))
    formula <- update.formula(formula, session ~ .)
    setup <- gam(formula, data = data, fit = FALSE, ...)
    colnames(setup$X) <- setup$term.names
    setup
}
############################################################################################

general.model.matrix <- function (formula, data, gamsmth = NULL, 
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
    
    dots <- list(...)

    if (any(polys(formula)))
        stop ("orthogonal polynomials are temporarily blocked")  ## 2014-09-12
    if (any(smooths(formula))) {
        if (is.null(gamsmth)) {
            ## setup knots etc from scratch
            mat <- gamsetup(formula, data, ...)$X
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
###############################################################################

## shifted from secrloglik 2016-10-16
makerealparameters <- function (design, beta, parindx, link, fixed) {
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
    parindx$D <- NULL ## detection parameters only
    link$D    <- NULL ## detection parameters only
    parindx$noneuc <- NULL ## detection parameters only
    link$noneuc    <- NULL ## detection parameters only
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
############################################################################################

secr.lpredictor <- function (formula, newdata, indx, beta, field, beta.vcv=NULL,
    smoothsetup = NULL, contrasts = NULL, f = NULL) {
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

    if (!is.null(f) && field == 'D') {
       Yp <- f(newdata[,vars[1]], beta = beta[indx]) 
       mat <- as.matrix(newdata[,vars[1], drop = FALSE])
    }
    else {
        
        mat <- general.model.matrix(formula, data = newdata, gamsmth = smoothsetup, 
            contrasts = contrasts)
        if (nrow(mat) < nrow(newdata))
            warning ("missing values in predictors?")
        
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
        if (is.null(f) || field != 'D') {
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
############################################################################################

edist <- function (xy1, xy2) {
    nr <- nrow(xy1)
    nc <- nrow(xy2)
    x1 <- matrix(xy1[,1], nr, nc)
    x2 <- matrix(xy2[,1], nr, nc, byrow=T)
    y1 <- matrix(xy1[,2], nr, nc)
    y2 <- matrix(xy2[,2], nr, nc, byrow=T)
    sqrt((x1-x2)^2 + (y1-y2)^2)
}

############################################################################################

## least cost paths from mask including barriers to movement
## use edist for equivalent Euclidean distances

nedist <- function (xy1, xy2, mask, inf = Inf, ...) {
    newargs <- list(...)
    if (missing(mask)) mask <- xy2
    noneuc <- covariates(mask)$noneuc
    if (is.null(noneuc)) noneuc <- rep(1, nrow(mask))
    defaultargs <- list(transitionFunction = mean, directions = 16)
    args <- replace(defaultargs, names(newargs), newargs)
    args$x <- raster(mask, values = noneuc)
    if (requireNamespace('gdistance', quietly = TRUE)) {    ## 2015-01-23
        tr <- do.call(gdistance::transition, args)
        tr <- gdistance::geoCorrection(tr, type = "c", multpl = FALSE)
        out <- gdistance::costDistance(tr, as.matrix(xy1), as.matrix(xy2))
    }
    else stop ("package gdistance is required for nedist")
    if (is.finite(inf)) out[!is.finite(out)] <- inf
    out
}

############################################################################################

getcellsize <- function (mask) {
    if (inherits(mask, 'linearmask'))
        cell <- attr(mask, 'spacing') / 1000  ## per km
    else
        cell <- attr(mask, 'area')            ## per ha
    if (is.null(cell))
        stop ("mask lacks valid cell size (area or spacing)")
    cell
}
############################################################################################

## 2014-10-25, 2017-01-24
## intercept and fix certain ,models with bad defaults
updatemodel <- function (model, detectfn, detectfns, oldvar, newvar, warn = FALSE) {
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
############################################################################################

## Manually remove some mask points
# simplified 2022-02-03

deleteMaskPoints <- function (mask, onebyone = TRUE, add = FALSE, poly = NULL,
                              poly.habitat = FALSE, ...) {
    ## interface does not work properly in RStudio

    if (ms(mask)) {         ## a list of mask objects
        if (inherits(poly, 'list') & (!is.data.frame(poly)))
            stop ("lists of polygons not implemented in 'make.mask'")
        temp <- lapply (mask, deleteMaskPoints, onebyone = onebyone, add = add,
                        poly = poly, poly.habitat = poly.habitat, ...)
        class (temp) <- c('mask', 'list')
        temp
    }
    else {
        plot(mask, add = add, ...)
        if (!is.null(poly)) {
            if (poly.habitat)
                pointstodrop <- (1:nrow(mask))[!pointsInPolygon(mask, poly)]
            else
                pointstodrop <- (1:nrow(mask))[pointsInPolygon(mask, poly)]
        }
        else if (onebyone) {
            cat ('Click to select points; right-click to stop\n')
            flush.console()
            xy <- locator(type = 'p', pch=1, col='red')
            pointstodrop <- if (length(xy$x)==0)
                numeric(0)
            else
                nearesttrap(xy, mask)
        }
        else {
            cat ('Click to select polygon vertices; right-click to stop\n')
            flush.console()
            xy <- locator(type = 'l', col='red')
            xy <- as.data.frame(xy)
            xy <- rbind(xy, xy[1,])
            if (poly.habitat)
                pointstodrop <- (1:nrow(mask))[!pointsInPolygon(mask, xy)]
            else
                pointstodrop <- (1:nrow(mask))[pointsInPolygon(mask, xy)]
        }
        npts <- length(pointstodrop)
        if (npts>0) {
            points(mask[pointstodrop,], pch = 16, col = 'red')
            if(.Platform$OS.type == "windows") {
                pl <- if (npts>1) 's' else ''
                msg <- paste ('Delete ', npts, ' red point',pl, '?', sep='')
                response <-  utils::winDialog(type = "okcancel", msg)
            } else {
                response <- 'OK'
            }
            if (response == 'OK') {
                mask <- subset(mask, -pointstodrop)
            if (npts==1)
                message("1 point deleted")
            else
                message(npts, " points deleted")
            }
        else
            message ("point(s) not deleted")
        }
        else
            message ("no points to delete")
        plot(mask, col='green')
        mask
    }
}
############################################################################################

nparameters <- function (object) {
    Npar <- max(unlist(object$parindx))
    Npar <- Npar + length(object$details$miscparm)
    ## allow for fixed beta parameters
    if (!is.null(object$details$fixedbeta))
        Npar <- Npar - sum(!is.na(object$details$fixedbeta))
    Npar
}
############################################################################################

mapbeta <- function (parindx0, parindx1, beta0, betaindex)

    ## Extend beta vector from simple model (beta0) to a more complex (i.e. general)
    ## model, inserting neutral values (zero) as required.
    ## For each real parameter, a 1:1 match is assumed between
    ## beta values until all beta values from the simpler model are
    ## used up. THIS ASSUMPTION MAY NOT BE JUSTIFIED.
    ## betaindex is a user-controlled alternative.

{
    ## list of zeroed vectors, one per real parameter
    beta1 <- lapply(parindx1, function (x) {x[]<-0; x})

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
        ## for (j in 1:length(beta1))
        ## improved replace by name2015-11-17
        for (j in names(beta1)) {
            if (j %in% names(beta0))
                beta1[[j]][indx[[j]]] <- beta0[parindx0[[j]]]
        }
        unlist(beta1)
    }
}
############################################################################################

xyinpoly <- function (xy, trps) {
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
############################################################################################
##
addzerodf <- function (df, oldCH, sess) {
    ## add dummy detection records to dataframe for 'all-zero' case
    ## that arises in sighting-only mark-resight with known marks
    allzero <- apply(oldCH,1,sum)==0
    naz <- sum(allzero)
    if (naz > 0) {
        df0 <- expand.grid(
          newID = rownames(oldCH)[allzero], 
          newocc = NA,
          newtrap = trap(oldCH)[1], 
          alive = TRUE, 
          sess = sess,
          stringsAsFactors = FALSE)
        df$x <- NULL; df$y <- NULL  ## 2021-04-08
        df <- rbind(df,df0)
        if (!is.null(xy(oldCH))) {
            df$x <- c(xy(oldCH)$x, rep(NA, naz))
            df$y <- c(xy(oldCH)$y, rep(NA, naz))
        }
        if (!is.null(signal(oldCH)))  {
            df$signal <- c(signal(oldCH), rep(NA, naz))
        }
    }
    df
}
############################################################################################

## including pre-marked animals never sighted
## cov is optional dataframe of covariates
addzeroCH <- function (CH, nzero, cov = NULL) {
    if (nzero == 0)
        return(CH)
    else {
        nc <- nrow(CH)
        chdim <- dim(CH)
        chdim[1] <- nzero
        extra <- array(0, dim=chdim)
        dimnames(extra) <- c(list(paste('Z', 1:nzero, sep='')), dimnames(CH)[2:3])
        CH2 <- abind(CH, extra, along = 1)
        class(CH2) <- 'capthist'
        traps(CH2) <- traps(CH)
        xy(CH2) <- xy(CH)  ## order is not affected by adding zero histories
        if (!is.null(covariates(CH)) & (nrow(CH)>0)) {
            if (is.null(cov)) {
                cov <- covariates(CH)[rep(1,nzero),]
                cov[,] <- NA   ## covariates are unknown
            }
            covariates(CH2) <- rbind(covariates(CH), cov[1:nzero,])
        }
        ## ... and other essential attributes?
        CH2
    }
}
############################################################################################

expandbinomN <- function (binomN, detectorcodes) {
    # assumes detectorcodes is a vector of length = noccasions
    binomN <- ifelse (detectorcodes %in% c(2,6,7), binomN, 1)
    if (any(is.na(binomN))) stop ("NA value in binomN")
    binomN
}
############################################################################################

check3D <- function (object) {
    if (ms(object)) {
        out <- lapply(object, check3D)
        class(out) <- class(object)
        out
    }
    else {
        if (is.matrix(object)) {
            warning("secr >= 3.0 requires 3-D capthist; using updateCH() to convert")
            updateCH(object)
        }
        else {
            object
        }
    }
}
############################################################################################

updateCH <- function(object) {
    if (!inherits(object, 'capthist'))
        stop ("requires capthist object")
    # following replaces this old code 2020-08-29
    # reduce(object, dropunused = FALSE)
    if (ms(object)) {
        out <- lapply(object, updateCH)
        class (out) <- c("capthist", "list")
        out
    }
    else {
        if (length(dim(object)) == 3) {
            return(object)
        }
        else {
            K <- ndetector(traps(object))
            ch <- array(0, dim = c(dim(object), K), dimnames = 
                    list(rownames(object), colnames(object), 1:K))
            OK <- as.logical(object!=0)
            animal <- row(object)[OK]
            occ <- col(object)[OK] 
            detn <- object[OK]
            ch[cbind(animal, occ, detn)] <- 1
            traps(ch) <- traps(object)
            class (ch) <- "capthist"
            session(ch) <- session(object)
            ch
        }
    }
}
############################################################################################

newstr <-function (strings) {
    ## compress a character vector
    ## use run length encoding function
    rl <- rle(strings)
    st <- rl$values
    le <- paste0(' (',as.character(rl$lengths), ')')
    le[le==' (1)'] <- ''
    paste(paste0(st, le), collapse = ', ')
}
# newstr(c("single", rep("proximity",4)))
############################################################################################

## moved from Telemetry 2016-12-28
outsidemask <- function(CH, mask, threshold = spacing(mask) / sqrt(2)) {
    xylist <- telemetryxy(CH)
    dfun <- function(xy) {
        centres <- matrix(apply(xy, 2, mean), ncol = 2)
        distancetotrap(centres, mask)
    }
    sapply(xylist, dfun) > threshold
}
############################################################################################

shareFactorLevels <- function (object, columns = NULL, stringsAsFactors = TRUE) {
    ## stringsAsFactors added 2020-05-16
    if (ms(object)) {
        if (!is.null(covariates(object))) {
            df <- do.call(rbind, covariates(object))
            if (is.null(columns)) {
                columns <- 1:ncol(df)
            }
            if (stringsAsFactors) {
                df[,columns] <- stringsAsFactors(df[,columns, drop = FALSE])
            }
            for (i in columns) {
                if (is.factor(df[,i])) {
                    levelsi <- levels(df[,i])
                    for (sess in 1:length(object)) {
                        covariates(object[[sess]])[,i] <-
                            factor(covariates(object[[sess]])[,i],
                                   levels = levelsi)
                    }
                }
            }
        }
    }
    else {
        # modified 2021-04-27 to apply to covariates, not object itself
        if (!is.null(covariates(object))) {
            if (stringsAsFactors) {
                df <- covariates(object)
                if (is.null(columns)) {
                    columns <- 1:ncol(df)
                }
                df[,columns] <- stringsAsFactors(df[,columns, drop = FALSE])
                covariates(object) <- df
            }
        }
    }
    object
}
############################################################################################

allzero <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, allzero)
    }
    else {
        telemocc <- detector(traps(object))=='telemetry'
        apply(object[,!telemocc,,drop=FALSE],1,sum)==0
    }
}
############################################################################################

primarysessions <- function(intervals) {
    primarysession <- cumsum(c(0,intervals))
    match(primarysession, unique(primarysession))
}
############################################################################################

secondarysessions <- function(intervals) {
    primary <- primarysessions(intervals)
    unname(unlist(sapply(table(primary), seq_len)))  
}
############################################################################################

boundarytoSF <- function (poly) {
  if (is.null(poly)) {
    NULL
  }
  else if(inherits(poly, c('sf','sfc'))) {
    poly <- st_geometry(poly) # extract sfc if not already sfc
    geomtype <- st_geometry_type(poly, by_geometry = FALSE)
    if (!geomtype %in% c("POLYGON", "MULTIPOLYGON")) {
      stop ("poly sfc should be of type POLYGON or MULTIPOLYGON")
    }
    poly
  }
  else if (inherits(poly, 'SpatialPolygons')) {   # also SPDF?
    st_as_sfc(poly)
  }
  else if (inherits(poly, 'SpatVector')) {
    st_as_sfc(as(poly,"Spatial"))
  }
  else if (inherits(poly, c('matrix', 'data.frame'))) {
    ## input is 2-column matrix for a single polygon
    poly <- matrix(unlist(poly), ncol = 2)
    poly <- rbind (poly, poly[1,])  ## force closure of polygon
    st_sfc(st_polygon(list(poly)))
  }
  else stop (class(poly), " not valid input to boundarytoSF")
}
###############################################################################

pointsInPolygon <- function (xy, poly, logical = TRUE) {
  # xy is 2-column matrix or data.frame of coordinates
  if (inherits(poly, 'mask')) { 
    if (ms(poly))
      stop ("multi-session masks not supported")
    sp <- spacing(poly)
    minx <- min(poly$x, na.rm = TRUE)
    miny <- min(poly$y, na.rm = TRUE)
    mask <- sweep(poly, MARGIN = 2, FUN = '+', STATS = c(-minx, -miny))
    mask <- round(mask/sp) + 1
    xy <- matrix(unlist(xy), ncol = 2)  ## in case dataframe
    xy <- sweep(xy, MARGIN = 2, FUN = '+', STATS = c(-minx, -miny))
    xy <- round(xy/sp) + 1
    xy[xy<=0] <- NA
    xy[,1][xy[,1]>max(mask$x, na.rm = TRUE)] <- NA
    xy[,2][xy[,2]>max(mask$y, na.rm = TRUE)] <- NA
    
    maskmatrix <- matrix(0, ncol = max(mask$y, na.rm = TRUE), nrow = max(mask$x, na.rm = TRUE))
    maskmatrix[as.matrix(mask)] <- 1:nrow(mask)
    inside <- maskmatrix[as.matrix(xy)]
    inside[is.na(inside)] <- 0
    if (logical)
      inside <- inside > 0
    inside
  }
  else {
    poly <- boundarytoSF(poly)
    if (inherits(poly, c('sf','sfc'))) {
      xy <- st_as_sf(data.frame(xy), coords = 1:2)
      st_crs(xy) <- st_crs(poly)
      apply(st_within(xy, poly, sparse = FALSE), 1, any)
    }
    else {
      stop ("unknown input to pointsInPolygon")
    }
  }
}
###############################################################################

## return indices of first occasion and detector for which PIAx is non-zero 
firstsk <- function (PIAx) {
  ## PIAx dim n,s,k
  wh <- function(d2) {
    match(TRUE, d2>0)
  }
  apply(PIAx,1,wh)
}

###############################################################################

maskboolean <- function (ch, mask, threshold) {
  if (ms(ch)) {
    if (!ms(mask)) stop ("masklookup: multisession ch requires multisession mask")
    outlist <- mapply(maskboolean, ch, mask, MoreArgs = list(threshold = threshold), SIMPLIFY = FALSE)
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

# mskl <- maskboolean(captdata, msk, 100)
# par(mfrow=c(4,4), mar = c(1,1,1,1))
# for (i in 1:16) {
#   plot(traps(captdata))
#   plot(subset(msk,mskl[i,]), add=T)
#   plot(subset(captdata,i), add=T)
# }
##############################################################################

## Based on Charles C. Berry on R-help 2008-01-13
## drop = FALSE 2018-11-22

## used in secr.design.MS

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
##############################################################################

setNumThreads <- function (ncores, ...) {
    ## environment variable RCPP_PARALLEL_NUM_THREADS is set by RcppParallel::setThreadOptions
    ## current is NA if variable not previously set
    current <- as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", ""))
    if (missing(ncores)) ncores <- NULL
    if (is.na(current)) {
        if (is.null(ncores)) {
            ncores <-  min(RcppParallel::defaultNumThreads(), 2)
        }
    }
    else if (is.null(ncores) || (ncores == current)) {
        return(current)   ## no need to change 
    }
    if (ncores > RcppParallel::defaultNumThreads()) 
        stop("requested ncores exceeds number available")
    if (ncores<1)
        stop ("specified ncores < 1")
    ncores <- min(ncores, RcppParallel::defaultNumThreads())
    RcppParallel::setThreadOptions(ncores, ...) 
    return(as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "")))
}
##############################################################################

telemcode <- function(object, ...) {
    if (inherits(object, 'traps') && !ms(object))
        switch (telemetrytype(object), none = 0, 
                independent = 1, dependent = 2, concurrent = 3, 0)
    else 
        NA
}
##############################################################################

uniquerownames <- function (capthist) {
    if (!ms(capthist)) {
        return(capthist)
    }
    else {
        last <- 0
        for (i in 1:length(capthist)) {
            nr <- nrow(capthist[[i]])
            if (nr > 0) {
            rownames(capthist[[i]]) <- last + (1:nr)
            last <- last+nr
            }
        }
        capthist
    }
}
##############################################################################

selectCHsession <- function(capthist, sessnum) {
    if (ms(capthist)) 
        capthist[[sessnum]]
    else 
        capthist
}
##############################################################################

stringsAsFactors <- function (DF) {
    # convert any character columns of a data.frame (or list) to factor
    if (is.list(DF) && length(DF)>0) {    ## bug fix 2020-08-14
        chr <- sapply(DF, is.character)
        DF[chr] <- lapply(DF[chr], as.factor)
    }
    DF
}
##############################################################################

# see also getuserdist in loglikhelperfn.R
# 2021-03-30 moved from preparedata.R 
# 2021-03-30 integrate HPX 

getdistmat2 <- function (traps, mask, userdist, HPX = FALSE) {
    ## Static distance matrix
    if (is.function(userdist)) {
        NULL   ## compute dynamically later
    }
    else {
        if (HPX) {
            if (any(detector(traps) %in% .localstuff$polydetectors)) {
                trps <- split(traps, polyID(traps))
                inside <- t(sapply(trps, pointsInPolygon, xy = mask))
                d2 <- 1-inside      # 0 inside, 1 outside
                d2[d2>0] <- 1e10    # 0 inside, 1e10 outside
                d2
            }
            else {
                # maximum of squared distance in x- or y- directions
                xydist2cpp(as.matrix(traps), as.matrix(mask))
            }
        }
        else {
            if (any(detector(traps) %in% .localstuff$polydetectors)) {
                ## do not use result if detector is one of
                ## polygonX, polygon, transectX, transect, OR telemetry?
                matrix(0, nrow = nrow(traps), ncol = nrow(mask))
            }
            else {
                # Euclidean distance
              # edist(as.matrix(traps), as.matrix(mask))^2
              edist2cpp(as.matrix(traps), as.matrix(mask))
            }
        }
    }
}
#--------------------------------------------------------------------------------

## 2022-01-04 for subset.capthist
rownum <- function (x) {
  if (length(dim(x)) < 1 || dim(x)[1] == 0) NULL
  else 1: (dim(x)[1])
}
colnum <- function (x) {
  if (length(dim(x)) < 2 || dim(x)[2] == 0) NULL
  else 1: (dim(x)[2])
}

## 2022-01-23
## function to assign all-ones usage matrix

uniformusage <- function(object, noccasions) {
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
        ndet <- ndetector(object[[r]])
        usage(object[[r]]) <- matrix(1, ndet, noccasions)
      }
    }
    else {
      ndet <- ndetector(object)
      usage(object) <- matrix(1, ndet, noccasions)
    }
  }
  object
}

sfrotate <- function (x, degrees, centrexy = NULL, usecentroid = FALSE) {
    rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
    gx <- st_geometry(x)
    if (is.null(centrexy)) {
        if (usecentroid) {
            centrexy <- st_centroid(gx)[1,]   # unique centre
        }
        else {
            centrexy <- st_centroid(st_as_sfc(st_bbox(x)))
        }
    } 
    else {
        centrexy <- st_sfc(st_point(centrexy) )
    }
    (gx - centrexy) * rot(degrees/360*2*pi) + centrexy
}

# Based on Tim Salabim stackoverflow Jul 12 2018
# https://stackoverflow.com/questions/51292952/snap-a-point-to-the-closest-point-on-a-line-segment-using-sf

snap_points <- function(x, y, max_dist = 1000) {
    
    if (inherits(x, "sf")) n = nrow(x)
    if (inherits(x, "sfc")) n = length(x)
    
    out = do.call(c,
        lapply(seq(n), function(i) {
            nrst = st_nearest_points(st_geometry(x)[i], y)
            nrst_len = st_length(nrst)
            nrst_mn = which.min(nrst_len)
            if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
            return(st_cast(nrst[nrst_mn], "POINT")[2])
        })
    )
    return(out)
}

rtpois <- function(n, lambda) {
    qpois(runif(n, dpois(0, lambda), 1), lambda)
}
